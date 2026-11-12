
#' Validate a cross_section formula at spec build time
#'
#' Requires `outcome ~ I(expr)` with a single `I()` term and rejects
#' time-series functions (the canonical `.pt_ts_fns` registry), which have no
#' meaning when the expression is evaluated across units at one time step.
#'
#' @param formula The formula passed to [build_model()].
#' @return `NULL`, invisibly; errors on an invalid formula.
#' @keywords internal
#' @noRd
.check_cross_section_formula <- function(formula) {
  rhs <- if (inherits(formula, "formula") && length(formula) == 3L) {
    rlang::f_rhs(formula)
  }
  if (!(is.call(rhs) && identical(rhs[[1L]], as.name("I")) && length(rhs) == 2L)) {
    stop("A cross_section formula must have the form `outcome ~ I(expr)` with a single I() term.",
         call. = FALSE)
  }
  ts_used <- intersect(all.vars(rhs, functions = TRUE), .pt_ts_fns)
  if (length(ts_used) > 0L) {
    stop("cross_section formulas are evaluated across units at a single time step, ",
         "so time-series functions (", paste(ts_used, collapse = ", "),
         ") are not allowed. Reference an already-lagged column instead, e.g. ",
         "produce `x_l1` with build_model(\"deterministic\", x_l1 ~ I(lag(x))).",
         call. = FALSE)
  }
  invisible(NULL)
}

#' @exportS3Method
fit_model.cross_section_spec <- function(spec, ctx = NULL, ...) {
  cross_section_model(formula = spec$formula)
}

#' Cross-section model
#'
#' Computes an outcome from an expression evaluated across all units at a
#' single simulated time step (e.g. a cross-sectional quantile). The formula
#' has already been validated by [build_model()].
#'
#' @param formula A two-sided formula of the form `outcome ~ I(expr)`.
#'
#' @return An endogenmodel of class `cross_section`.
#' @keywords internal
cross_section_model <- function(formula) {
  model <- new_endogenmodel(formula)
  class(model) <- c("cross_section", class(model))
  model$independent <- FALSE
  model$outcome <- parse_formula(model)$outcome
  # Argument of I(); shape validated by build_model().
  model$expr <- rlang::f_rhs(formula)[[2L]]
  model
}

#' Predict function for a cross-section model
#'
#' Evaluates the model expression over all unit rows at time `t`, separately
#' within each simulation draw. A length-1 result is copied to every unit; a
#' result with one value per unit is kept as is.
#'
#' @param model A `cross_section` endogenmodel.
#' @param t Time step to predict.
#' @param data A data.table.
#' @param ctx A panel_context object.
#' @param ... Ignored, accepted for S3 generic consistency.
#'
#' @return A data.table with key + index + outcome columns.
#' @family simulation
#' @export
predict.cross_section <- function(model, t, data, ctx, ...) {
  idx <- ctx_time(ctx)
  sim_var <- ctx_sim(ctx)
  all_keys <- ctx_keys(ctx)
  y <- model$outcome

  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(as.data.frame(data))
  }

  # .dt_rows(): `data[data[[idx]] == t]` breaks when the time column is named `t`.
  rows <- .dt_rows(data, data[[idx]] == t)
  env <- environment(model$formula)
  result_cols <- c(all_keys, idx)

  eval_slice <- function(df) {
    v <- eval(model$expr, envir = df, enclos = env)
    n <- nrow(df)
    if (length(v) == 1L) {
      v <- rep(v, n)
    } else if (length(v) != n) {
      stop(sprintf(
        "cross_section model '%s': expression returned length %d; expected 1 or %d (one per unit at %s = %s).",
        y, length(v), n, idx, t), call. = FALSE)
    }
    out <- df[, ..result_cols]
    # as.vector() drops names such as "90%" from quantile().
    data.table::set(out, j = y, value = as.vector(v))
    out
  }

  if (is.null(sim_var)) return(eval_slice(rows))
  # One pass: split by sim instead of filtering per sim id.
  data.table::rbindlist(lapply(split(rows, by = sim_var, sorted = TRUE), eval_slice))
}
