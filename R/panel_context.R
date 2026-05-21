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
