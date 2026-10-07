# Unified predictive-draw interface ------------------------------------------
#
# One generic replaces the per-family draw logic (the old getpi, getpi_glm,
# predict.heterolm's inline rnorm, .sample_from_fitdist point draws, ...) with
# an explicit separation of PARAMETER uncertainty (n_param) from INNOVATION
# uncertainty (n_innov). Family methods live next to their predict methods.

#' Draw from a model's predictive distribution
#'
#' Unified predictive-draw interface for endogenr models. Explicitly separates
#' parameter (coefficient / distribution-parameter) uncertainty from
#' innovation (residual / response-scale) uncertainty.
#'
#' @details
#' ## Contract
#'
#' `newdata` contains rows on which to draw, **already materialized**: the
#' columns the fitted stage-2 object needs (what the `predict.*` methods hold
#' after `.apply_ts_map()` + `.pt_apply_aliases()`). Stage-1 panel
#' materialization stays in `predict.*`; `draw_predictive()` is the seam below
#' it.
#'
#' The return value is a numeric **matrix** with `nrow(newdata)` rows and
#' `P * K` columns, where `P = max(n_param, 1)` and `K = max(n_innov, 1)`.
#' Column `(j - 1) * K + k` holds parameter draw `j`, innovation draw `k`.
#' Always a matrix, even at 1 x 1 (wrap with `as.vector()` if needed).
#'
#' - `n_param = 0`: condition on the point estimates (no parameter
#'   uncertainty), `P = 1`. `n_param >= 1`: that many independent parameter
#'   draws.
#' - `n_innov = 0`: return the conditional mean per parameter draw (this
#'   subsumes `predict(..., what = "expectation")`), `K = 1`. `n_innov >= 1`:
#'   independent innovations per row within each parameter block.
#' - `param_scope` (linear, glm, and glmmTMB methods only): `"draw"` (default)
#'   uses one parameter deviate per column block, shared across rows —
#'   coherent replicates for analysis. `"row"` uses an independent parameter
#'   deviate per row — the simulation engine's row-expansion convention where
#'   each grid row is its own sim world; the engine's `predict.*` calls use
#'   this scope so ensemble statistics are unchanged.
#'
#' ## Per-family support
#'
#' | family | parameter draw | innovation draw |
#' |---|---|---|
#' | `linear` | scaled-t (se.fit-based; `param_scope`) | Gaussian at the drawn residual scale |
#' | `glm` | t on the link scale (se.fit; `param_scope`) | family response draw at dispersion |
#' | `glmmTMB` | normal on the link scale (se.fit; `param_scope`) | family response draw + zero-inflation mask |
#' | `heterolm` | none (warns once; `n_param` ignored) | `N(mu_i, sigma_i)` |
#' | `gamlss` | none (warns once; `n_param` ignored) | family `r<FAM>` at predicted parameters |
#' | `parametric_distribution` | MVN(estimate, vcov) per block | `r<dist>` at the drawn parameters; rejects `n_innov = 0` |
#' | `univariate_fable` | none (warns once; `n_param` ignored) | coherent `fabletools::generate()` paths; rejects `n_innov = 0` |
#' | `deterministic`, `cross_section`, `exogen`, `spatial_lag` | — | no stochastic predictive distribution (error) |
#'
#' @param model A fitted endogenmodel (from [fit_model()] / [fit_system()]).
#' @param newdata A data.frame or data.table of prediction rows (see Details).
#' @param n_param Integer >= 0. Number of parameter draws (0 = point
#'   estimates).
#' @param n_innov Integer >= 0. Number of innovation draws per parameter draw
#'   (0 = conditional mean).
#' @param param_scope Parameter-deviate granularity for the linear, glm, and
#'   glmmTMB methods: `"draw"` (default) shares one deviate per column block
#'   (coherent replicates); `"row"` draws an independent deviate per row (the
#'   engine's row-expansion convention).
#' @param ctx A [panel_context()] object (univariate_fable method only) naming
#'   the unit/time columns of `newdata`.
#' @param horizon Integer number of forecast steps to generate
#'   (univariate_fable method only); defaults to the number of distinct times
#'   in `newdata`.
#' @param ... Further method-specific arguments.
#'
#' @return A numeric matrix, `nrow(newdata)` x `max(n_param, 1) *
#'   max(n_innov, 1)`.
#' @family simulation
#' @export
draw_predictive <- function(model, newdata, n_param = 1L, n_innov = 1L, ...) {
  UseMethod("draw_predictive")
}

#' @rdname draw_predictive
#' @export
draw_predictive.default <- function(model, newdata, n_param = 1L, n_innov = 1L, ...) {
  stop("No draw_predictive method for class '", class(model)[1],
       "': deterministic, cross_section, exogen, and spatial_lag models have ",
       "no stochastic predictive distribution.", call. = FALSE)
}

# Validate the draw-count arguments shared by every method.
.check_draw_counts <- function(n_param, n_innov) {
  ok <- function(v) is.numeric(v) && length(v) == 1L && is.finite(v) && v >= 0
  if (!ok(n_param) || !ok(n_innov)) {
    stop("`n_param` and `n_innov` must each be a single non-negative finite number.",
         call. = FALSE)
  }
  invisible(NULL)
}

# One-time-per-session warning keyed by an options() flag (same pattern as
# .glm_warn_unsupported).
.warn_once <- function(key, msg) {
  if (!isTRUE(getOption(key))) {
    warning(msg, call. = FALSE)
    do.call(options, stats::setNames(list(TRUE), key))
  }
  invisible(NULL)
}

# One-time-per-session warning for families without a parameter-uncertainty
# draw path.
.warn_no_param_draw <- function(label) {
  .warn_once(paste0(".endogenr_no_param_draw_", label),
             paste0(label, " models carry no parameter-uncertainty draw; ",
                    "n_param is ignored (point estimates used)."))
}
