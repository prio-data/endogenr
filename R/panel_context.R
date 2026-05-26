# Required for data.table to work correctly inside packages
.datatable.aware <- TRUE

#' Create a panel context object
#'
#' A lightweight metadata carrier that stores panel dimension names (unit, time,
#' sim) separately from the data. This decouples panel metadata from the data
#' carrier, preventing silent metadata loss on data operations.
#'
#' @param unit Character. Column name identifying panel units.
#' @param time Character. Column name identifying time periods.
#' @param sim Character or NULL. Column name identifying simulation draws.
#'
#' @return A `panel_context` object.
#' @export
panel_context <- function(unit, time, sim = NULL) {
  stopifnot(is.character(unit), length(unit) == 1L)
  stopifnot(is.character(time), length(time) == 1L)
  structure(
    list(unit = unit, time = time, sim = sim),
    class = "panel_context"
  )
}

#' @rdname panel_context
#' @param ctx A `panel_context` object.
#' @export
ctx_unit <- function(ctx) ctx$unit

#' @rdname panel_context
#' @export
ctx_time <- function(ctx) ctx$time

#' @rdname panel_context
#' @export
ctx_sim <- function(ctx) ctx$sim

#' @rdname panel_context
#' @export
ctx_keys <- function(ctx) c(ctx$unit, ctx$sim)

#' Infer a panel_context from objects
#'
#' Walks `...` looking for either an existing `panel_context` or an object
#' carrying `panel_unit` / `panel_time` attributes (the convention used by
#' `paneltools::as_panel()`). Returns the first match, or `NULL` if none
#' of the inputs supply panel metadata.
#'
#' @param ... Objects to inspect (typically simulation results and/or truth data).
#' @param sim Character or NULL. Optional simulation column name to attach.
#'
#' @return A `panel_context` object or `NULL`.
#' @keywords internal
.infer_ctx <- function(..., sim = NULL) {
  args <- list(...)
  for (x in args) {
    if (is.null(x)) next
    if (inherits(x, "panel_context")) return(x)
    unit <- attr(x, "panel_unit", exact = TRUE)
    time <- attr(x, "panel_time", exact = TRUE)
    if (is.character(unit) && length(unit) == 1L &&
        is.character(time) && length(time) == 1L) {
      return(panel_context(unit = unit, time = time, sim = sim))
    }
  }
  NULL
}
