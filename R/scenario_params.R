# Scenario parameter framework -----------------------------------------------
#
# Time fixed-effect draws and coefficient overrides for linear models.
# Extension seams (scenario_terms S3 generic) allow other model types to gain
# support by implementing their own methods.

# --- Internal: detect factor(timevar) in a fitted linear model ---------------

# Scan fit_formula for a term of the form factor(<timevar>). Extract the
# estimated year effects from fitted_lm and return metadata, or NULL.
.detect_time_fe <- function(fit_formula, fitted_lm, timevar) {
  term_labels <- attr(stats::terms(fit_formula), "term.labels")

  # Find the label whose parse is factor(<timevar>)
  fe_label <- NULL
  for (lbl in term_labels) {
    expr <- tryCatch(str2lang(lbl), error = function(e) NULL)
    if (is.null(expr) || !is.call(expr)) next
    fn <- as.character(expr[[1L]])
    if (!fn %in% c("factor", "stats::factor")) next
    if (length(expr) < 2L) next
    arg1 <- expr[[2L]]
    if (!is.symbol(arg1)) next
    if (as.character(arg1) != timevar) next
    fe_label <- lbl
    break
  }

  if (is.null(fe_label)) return(NULL)

  # Guard: timevar must not appear in any other RHS term
  other_labels <- setdiff(term_labels, fe_label)
  for (lbl in other_labels) {
    vars_in <- tryCatch(all.vars(str2lang(lbl)), error = function(e) character(0L))
    if (timevar %in% vars_in) {
      stop(
        "time fixed-effect scenario handling requires the time variable to enter ",
        "only through `factor(", timevar, ")` ",
        "(found it used elsewhere in the formula).",
        call. = FALSE
      )
    }
  }

  # xlevels for this term — present only when the model has been fitted
  xl <- fitted_lm$xlevels[[fe_label]]
  if (is.null(xl)) return(NULL)

  # Build effects vector: baseline (xl[1]) -> 0, remaining levels -> coefficient
  effects <- stats::setNames(numeric(length(xl)), xl)
  for (lv in xl[-1L]) {
    nm  <- paste0(fe_label, lv)
    val <- stats::coef(fitted_lm)[[nm]]
    if (!is.null(val)) effects[[lv]] <- val
  }
  ref_value <- as.numeric(xl[1L])

  list(
    term      = fe_label,
    timevar   = timevar,
    ref_value = ref_value,
    effects   = unname(effects)
  )
}

# --- Internal: policy resolution ---------------------------------------------

# Return a fully-specified time-FE policy list.  Missing fields are filled from
# the estimated effects so the draw kernel always sees concrete numbers.
.resolve_policy <- function(scenario, outcome, time_fe) {
  default_pol <- list(type = "resample")
  if (is.null(scenario)) return(default_pol)
  oc_entry <- scenario[[outcome]]
  if (is.null(oc_entry)) return(default_pol)
  pol <- oc_entry$policy
  if (is.null(pol)) return(default_pol)

  if (identical(pol$type, "distribution")) {
    if (is.null(pol$mean)) pol$mean <- mean(time_fe$effects)
    if (is.null(pol$sd))   pol$sd   <- stats::sd(time_fe$effects)
  } else if (identical(pol$type, "fixed")) {
    if (is.null(pol$value)) pol$value <- mean(time_fe$effects)
  }
  pol
}

# --- Internal: per-row time-FE draw -----------------------------------------

# One draw per unique sim id (year effect is shared within a trajectory);
# result is mapped back to rows.
.draw_time_fe_offset <- function(effects, sim_ids, policy) {
  sims <- unique(sim_ids)
  d <- switch(
    policy$type,
    resample     = sample(effects, length(sims), replace = TRUE),
    distribution = stats::rnorm(length(sims), policy$mean, policy$sd),
    fixed        = rep(policy$value, length(sims)),
    stop("unknown time-FE policy type: '", policy$type, "'", call. = FALSE)
  )
  d[match(sim_ids, sims)]
}

# --- Internal: scalar for the expectation path ------------------------------

.policy_mean <- function(policy, effects) {
  switch(
    policy$type,
    resample     = mean(effects),
    distribution = policy$mean,
    fixed        = policy$value,
    stop("unknown time-FE policy type: '", policy$type, "'", call. = FALSE)
  )
}

# --- S3 generic: scenario_terms ---------------------------------------------

#' Retrieve time fixed-effect metadata for scenario draws
#'
#' Returns the time-FE metadata list stored on a fitted model, or `NULL` when
#' the model carries no time fixed effect. This is the per-model-type extension
#' seam consumed by [setup_param()]. Currently only `linear` models support
#' time-FE scenario draws; all other model types return `NULL` via the default
#' method. Add support for additional model types by implementing
#' `scenario_terms.<class>`.
#'
#' @param model A fitted endogenr model object.
#' @return A list with elements `term`, `timevar`, `ref_value`, and `effects`,
#'   or `NULL` when no time fixed effect is present.
#' @family simulation
#' @export
scenario_terms <- function(model) UseMethod("scenario_terms")

#' @rdname scenario_terms
#' @exportS3Method scenario_terms default
scenario_terms.default <- function(model) NULL

#' @rdname scenario_terms
#' @exportS3Method scenario_terms linear
scenario_terms.linear <- function(model) model$time_fe

# --- Internal: resolve one coefficient override at forecast step h -----------

.resolve_coef_override <- function(spec, h, beta_hat) {
  if (is.function(spec)) {
    val <- as.numeric(spec(h, beta_hat))
  } else if (length(spec) == 1L) {
    val <- as.numeric(spec)
  } else if (length(spec) > 1L) {
    val <- as.numeric(spec[[h]])
  } else {
    stop(
      "coefficient override must be a scalar, a length-horizon numeric vector, ",
      "or a function(h, beta_hat)",
      call. = FALSE
    )
  }
  if (!is.finite(val)) {
    stop("coefficient override resolved to a non-finite value at step h = ", h, ".",
         call. = FALSE)
  }
  val
}

# --- setup_param: exported --------------------------------------------------

#' Set up scenario parameters for a fitted endogenr system
#'
#' Builds an editable list of scenario parameters, one entry per model outcome
#' that carries adjustable parameters. Inspect the returned object before
#' [simulate_system()] to see which coefficients and time-FE policies can be
#' overridden.
#'
#' Each entry (keyed by outcome name) contains:
#'
#' \describe{
#'   \item{`type`}{Model type string (e.g. `"linear"`).}
#'   \item{`time_fe`}{`NULL`, or a list describing the time fixed effect:
#'     `term`, `timevar`, `n_levels`, `effects` (numeric vector of estimated
#'     year effects with baseline fixed at 0), and a `summary` vector of
#'     mean/sd/min/max.}
#'   \item{`policy`}{Time-FE draw policy. Set `$type` to `"resample"` (the
#'     default; resample uniformly from estimated year effects),
#'     `"distribution"` (Normal fit; optionally set `$mean` and `$sd`), or
#'     `"fixed"` (constant; set `$value`). Fields left `NULL` are filled from
#'     the estimated effects at simulation time. Consumed only when `time_fe`
#'     is non-`NULL`.}
#'   \item{`coef_table`}{A `data.frame` of non-FE coefficient estimates and
#'     standard errors, or `NULL`. Time-FE level rows are excluded (managed
#'     by `policy`).}
#'   \item{`estimates`}{Named numeric vector of overridable coefficient
#'     estimates. Names match those used as keys in `$coefficients`.}
#'   \item{`coefficients`}{An empty `list()` for the user to fill with
#'     coefficient overrides. Each key is a coefficient name from `estimates`;
#'     each value is one of:
#'     \itemize{
#'       \item a scalar (constant across all forecast steps),
#'       \item a length-`horizon` numeric vector (per-step trajectory), or
#'       \item a `function(h, beta_hat)` returning the override β* for forecast
#'             step `h` given the per-draw fitted estimate.
#'     }
#'     Coefficient overrides are only supported for `linear` models;
#'     [simulate_system()] errors if they are set on any other model type.}
#' }
#'
#' @param fitted_system An `endogenr_fitted_system` from [fit_system()].
#' @return An `endogenr_scenario_params` object (named list keyed by outcome).
#'   An empty list is valid — simulation then finds nothing to configure.
#' @seealso [simulate_system()], [scenario_terms()]
#' @family simulation
#' @export
setup_param <- function(fitted_system) {
  if (!inherits(fitted_system, "endogenr_fitted_system")) {
    stop("`fitted_system` must be the output of fit_system().", call. = FALSE)
  }

  entries <- list()

  for (model in fitted_system$fitted_models) {
    has_coefs   <- !is.null(model$coefs)
    tfe         <- scenario_terms(model)
    has_time_fe <- !is.null(tfe)

    if (!has_coefs && !has_time_fe) next

    oc <- model$outcome

    time_fe_info <- if (has_time_fe) {
      list(
        term     = tfe$term,
        timevar  = tfe$timevar,
        n_levels = length(tfe$effects),
        effects  = tfe$effects,
        summary  = c(
          mean = mean(tfe$effects),
          sd   = stats::sd(tfe$effects),
          min  = min(tfe$effects),
          max  = max(tfe$effects)
        )
      )
    } else NULL

    # coef_table: strip FE-level rows and NA-estimate rows
    coef_tbl <- if (has_coefs) {
      tbl <- as.data.frame(model$coefs)[, c("term", "estimate", "std.error"),
                                        drop = FALSE]
      if (has_time_fe) {
        tbl <- tbl[!startsWith(tbl$term, tfe$term), , drop = FALSE]
      }
      tbl <- tbl[!is.na(tbl$estimate), , drop = FALSE]
      rownames(tbl) <- NULL
      if (nrow(tbl) == 0L) NULL else tbl
    } else NULL

    estimates <- if (!is.null(coef_tbl)) {
      stats::setNames(coef_tbl$estimate, coef_tbl$term)
    } else numeric(0L)

    entries[[oc]] <- list(
      type         = class(model)[1L],
      time_fe      = time_fe_info,
      policy       = list(type = "resample"),
      coef_table   = coef_tbl,
      estimates    = estimates,
      coefficients = list()
    )
  }

  structure(entries, class = "endogenr_scenario_params")
}

# --- print.endogenr_scenario_params -----------------------------------------

#' @exportS3Method print endogenr_scenario_params
print.endogenr_scenario_params <- function(x, ...) {
  cat("<endogenr_scenario_params>\n")
  if (length(x) == 0L) {
    cat("  (no parameterised outcomes)\n")
    return(invisible(x))
  }

  for (oc in names(x)) {
    e <- x[[oc]]
    cat("\n-- outcome:", oc, " [", e$type, "] --\n")

    if (!is.null(e$time_fe)) {
      smry <- e$time_fe$summary
      cat(sprintf(
        "  time FE: term='%s', n_levels=%d, policy='%s'\n",
        e$time_fe$term, e$time_fe$n_levels, e$policy$type
      ))
      cat(sprintf(
        "    effects: mean=%.3f, sd=%.3f, min=%.3f, max=%.3f\n",
        smry["mean"], smry["sd"], smry["min"], smry["max"]
      ))
    }

    if (!is.null(e$coef_table)) {
      tbl <- e$coef_table
      ovr_col <- vapply(tbl$term, function(nm) {
        spec <- e$coefficients[[nm]]
        if (is.null(spec)) return("\u2014")  # em dash
        if (is.function(spec)) return("function")
        if (length(spec) == 1L) return(format(as.numeric(spec)))
        paste0("[", length(spec), "-step vector]")
      }, character(1))
      tbl$override <- ovr_col
      print(tbl, row.names = FALSE)
    }
  }

  # Hint lines for the first outcome
  oc1 <- names(x)[1L]
  cat("\nHints:\n")
  cat("  dp$", oc1, "$policy$type <- \"distribution\"  # or 'fixed'\n", sep = "")
  if (!is.null(x[[oc1]]$coef_table) && nrow(x[[oc1]]$coef_table) > 0L) {
    term1 <- x[[oc1]]$coef_table$term[1L]
    cat("  dp$", oc1, "$coefficients$`", term1,
        "` <- 0  # scalar, length-horizon vector, or function(h, beta_hat)\n",
        sep = "")
  }

  invisible(x)
}
