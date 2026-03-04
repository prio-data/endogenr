
#' Spatial Lag Model
#'
#' Computes a spatially-lagged variable (W·y) as a weighted average of a variable
#' across neighboring geographic units at each time step. No model fitting is
#' required — it is a pure cross-sectional transformation applied at every \code{t}.
#'
#' Use \code{lag()} on the RHS to reference t-1 values (always available during
#' simulation). Without \code{lag()}, the source variable is read at the current
#' time \code{t} (valid when the source variable is already computed at \code{t}
#' by an earlier model in the topological execution order).
#'
#' @param formula Formula of the form \code{sl_y ~ lag(y)} or \code{sl_y ~ y}.
#'   The LHS names the new spatial-lag column; the RHS names the source variable.
#' @param nb Neighbor list (e.g., from \code{sfdep::st_contiguity()}).
#' @param wt Weights list (e.g., from \code{sfdep::st_weights()}).
#' @param unit_ids Ordered vector of geographic unit IDs matching the positional
#'   indexing of \code{nb} and \code{wt}.
#' @param island_default Value assigned to units with no neighbors. Defaults to
#'   \code{NA_real_}. Supply \code{0} or a global mean if preferred.
#'
#' @return An endogenmodel of class \code{c("spatial_lag", "endogenmodel")}.
#' @export
spatial_lag_model <- function(formula, nb, wt, unit_ids, island_default = NA_real_) {
  if (!requireNamespace("sfdep", quietly = TRUE)) {
    stop("Package 'sfdep' is required for spatial lag models. Install it with install.packages('sfdep').")
  }

  model <- new_endogenmodel(formula)
  model$independent <- FALSE

  # Extract source variable name from RHS
  model$source_var <- base::all.vars(rlang::f_rhs(formula))[1]

  # Determine if the RHS uses lag() — func_in_term returns c(lhs_bool, rhs_bool...)
  model$use_lag <- any(func_in_term(formula, func = "lag")[-1])

  # Store spatial weights
  model$nb <- nb
  model$wt <- wt
  model$unit_ids <- unit_ids
  model$island_default <- if (is.null(island_default)) NA_real_ else island_default

  # Pre-compute island mask (units with no neighbors)
  model$is_island <- lengths(nb) == 0

  # Set outcome name
  model$outcome <- base::all.vars(rlang::f_lhs(formula))

  class(model) <- c("spatial_lag", class(model))
  return(model)
}


#' Predict method for spatial lag models
#'
#' @param model A \code{spatial_lag} model object.
#' @param t The current time step.
#' @param data A tsibble with the full simulation data.
#' @param ... Not used.
#'
#' @return A tsibble keyed by the group variables with the spatial lag values at
#'   time \code{t}.
#' @export
predict.spatial_lag <- function(model, t, data, ...) {
  grp <- tsibble::key_vars(data)
  idx <- tsibble::index_var(data)
  outcome <- model$outcome

  # Separate the geographic unit key from the simulation index key
  geo_var <- setdiff(grp, "sim")
  sim_var <- intersect(grp, "sim")

  t_source <- if (model$use_lag) t - 1L else t

  source_data <- dplyr::as_tibble(data) |>
    dplyr::filter(!!rlang::sym(idx) == t_source)

  # Compute the spatial lag for one data.frame (single sim, all geo units at t_source)
  compute_sl <- function(df) {
    # Reorder values to match the positional indexing in nb/wt
    reorder_idx <- match(model$unit_ids, df[[geo_var]])
    ordered_vals <- df[[model$source_var]][reorder_idx]

    if (is.null(ordered_vals)) {
      stop(sprintf(
        "spatial_lag: source variable '%s' not found in simulation data at t=%s. ",
        model$source_var, unique(df[[idx]])
      ))
    }
    if (!is.numeric(ordered_vals)) {
      stop(sprintf(
        "spatial_lag: source variable '%s' is %s, not numeric (t=%s). ",
        model$source_var, class(ordered_vals)[1], unique(df[[idx]])
      ))
    }
    if (length(ordered_vals) != length(model$unit_ids)) {
      stop(sprintf(
        "spatial_lag: ordered_vals has length %d but unit_ids has length %d (t=%s). ",
        length(ordered_vals), length(model$unit_ids), unique(df[[idx]])
      ))
    }

    ordered_vals <- as.double(ordered_vals)

    # Compute weighted spatial lag
    sl_vals <- sfdep::st_lag(ordered_vals, model$nb, model$wt, na_ok = FALSE, allow_zero = TRUE)

    # Override island NAs with the user-specified default
    sl_vals[model$is_island] <- model$island_default

    # Map sl_vals (indexed by unit_ids position) back to df row order
    sl_lookup <- stats::setNames(sl_vals, as.character(model$unit_ids))
    sl_mapped <- sl_lookup[as.character(df[[geo_var]])]
    names(sl_mapped) <- NULL
    sl_mapped
  }

  if (length(sim_var) > 0) {
    sims <- unique(source_data[[sim_var]])
    result_list <- lapply(sims, function(s) {
      df <- source_data[source_data[[sim_var]] == s, , drop = FALSE]
      df[[outcome]] <- compute_sl(df)
      df[[idx]] <- t  # assign result at time t (matters for the lagged variant)
      df[, c(grp, idx, outcome), drop = FALSE]
    })
    result <- do.call(rbind, result_list)
  } else {
    source_data[[outcome]] <- compute_sl(source_data)
    source_data[[idx]] <- t
    result <- source_data[, c(grp, idx, outcome), drop = FALSE]
  }

  tsibble::as_tsibble(result, key = grp, index = idx)
}


#' Build spatial weights from an sf object
#'
#' A convenience helper that computes contiguity-based neighbor and weight lists
#' from an sf data.frame. The returned list can be passed directly to
#' \code{build_model("spatial_lag", ...)}.
#'
#' @param sf_data An sf data.frame with one row per geographic unit.
#' @param id_col Character. Name of the column in \code{sf_data} containing unit IDs.
#' @param contiguity_args A named list of additional arguments passed to
#'   \code{\link[sfdep]{st_contiguity}} (e.g., \code{list(queen = FALSE)}).
#' @param weights_args A named list of additional arguments passed to
#'   \code{\link[sfdep]{st_weights}} (e.g., \code{list(style = "W")}).
#'
#' @return A list with elements \code{nb}, \code{wt}, and \code{unit_ids}.
#' @export
st_weights_from_sf <- function(sf_data, id_col, contiguity_args = list(), weights_args = list()) {
  if (!requireNamespace("sfdep", quietly = TRUE)) {
    stop("Package 'sfdep' is required. Install it with install.packages('sfdep').")
  }
  nb <- do.call(sfdep::st_contiguity, c(list(sf_data), contiguity_args))
  wt <- do.call(sfdep::st_weights, c(list(nb), weights_args))
  list(nb = nb, wt = wt, unit_ids = sf_data[[id_col]])
}
