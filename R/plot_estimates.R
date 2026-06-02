#' Plot regression estimates from dynamic simulation models
#'
#' Creates a coefficient plot with pointranges for all fitted linear models
#' in a cross-validated simulation setup, optionally appended with a table
#' of goodness-of-fit statistics (R-squared and N).
#'
#' @param models A named list (keyed by test start year) of [fit_system()]
#'   outputs (each carrying a representative `$fitted_models`). Typically
#'   loaded from a saved `models.rds` file.
#' @param outcome_labels An optional named character vector mapping outcome
#'   variable names to display labels (e.g.,
#'   `c(dem = "Democracy", yjbest = "Conflict")`). If `NULL`, raw outcome
#'   names are used.
#' @param show_gof Logical. If `TRUE` (default), a goodness-of-fit table
#'   (R-squared and N) is appended to the right of the coefficient plot.
#' @param gof_stats Character vector of columns from [broom::glance()] to
#'   display. Defaults to `c("r.squared", "nobs")`.
#' @param gof_labels Named character vector of display labels for `gof_stats`.
#'   Defaults to `c(r.squared = "R\u00b2", nobs = "N")`.
#' @param gof_formats Named character vector of [sprintf()] format strings for
#'   `gof_stats`. Defaults to `c(r.squared = "%.3f", nobs = "%.0f")`.
#' @param widths Numeric vector of length 2 controlling the relative widths of
#'   the coefficient plot and the goodness-of-fit table. Default `c(7, 1)`.
#' @param base_size Numeric. Base font size passed to [ggplot2::theme_bw()].
#'   Default `9`.
#'
#' @return A `patchwork` object (if `show_gof = TRUE`) or a `ggplot` object.
#' @family postprocess
#' @export
plot_estimates <- function(
    models,
    outcome_labels = NULL,
    show_gof       = TRUE,
    gof_stats      = c("r.squared", "nobs"),
    gof_labels     = c(r.squared = "R\u00b2", nobs = "N"),
    gof_formats    = c(r.squared = "%.3f", nobs = "%.0f"),
    widths         = c(7, 1),
    base_size      = 9
) {
  # ── Extract coefficient estimates ──────────────────────────────────────
  estimates <- do.call(rbind, lapply(names(models), function(ts) {
    do.call(rbind, lapply(models[[ts]]$fitted_models, function(m) {
      if (!is.null(m$coefs)) {
        m$coefs |>
          dplyr::mutate(
            test_start = as.integer(ts),
            outcome    = m$outcome,
            .before    = 1
          )
      }
    }))
  }))

  estimates$est_low  <- estimates$estimate - 1.96 * estimates$std.error
  estimates$est_high <- estimates$estimate + 1.96 * estimates$std.error

  # ── Apply outcome labels ───────────────────────────────────────────────
  if (!is.null(outcome_labels)) {
    estimates$outcome <- dplyr::recode(estimates$outcome, !!!outcome_labels)
    lvls <- unname(outcome_labels)
  } else {
    lvls <- unique(estimates$outcome)
  }
  estimates$outcome <- factor(estimates$outcome, levels = lvls)

  est_no_int <- estimates[estimates$term != "(Intercept)", ]
  max_terms <- max(tapply(est_no_int$term, est_no_int$outcome,
                          function(x) length(unique(x))))

  # ── Coefficient plot ───────────────────────────────────────────────────
  p_coef <- ggplot2::ggplot(
    est_no_int,
    ggplot2::aes(
      y        = factor(.data$test_start),
      x        = .data$estimate,
      xmin     = .data$est_low,
      xmax     = .data$est_high,
      linetype = .data$p.value >= 0.05,
      color    = factor(.data$test_start)
    )
  ) +
    ggplot2::geom_pointrange() +
    ggplot2::geom_vline(xintercept = 0, linetype = 5) +
    ggh4x::facet_nested_wrap(
      outcome ~ term,
      scales     = "free_x",
      ncol       = max_terms,
      nest_line  = FALSE,
      trim_blank = FALSE,
      strip      = ggh4x::strip_nested(
        background_x = ggh4x::elem_list_rect(fill = c("grey85", "white")),
        text_x       = ggh4x::elem_list_text(face = c("bold", "plain")),
        by_layer_x   = TRUE
      )
    ) +
    ggplot2::labs(y = "Test start", x = "Estimate") +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(legend.position = "none")

  if (!show_gof) return(p_coef)

  # ── Goodness-of-fit table ──────────────────────────────────────────────
  fit_stats <- do.call(rbind, lapply(names(models), function(ts) {
    do.call(rbind, lapply(models[[ts]]$fitted_models, function(m) {
      if (!is.null(m$gof)) {
        m$gof |>
          dplyr::select(dplyr::any_of(gof_stats)) |>
          dplyr::mutate(
            test_start = as.integer(ts),
            outcome    = m$outcome,
            .before    = 1
          )
      }
    }))
  }))

  if (!is.null(outcome_labels)) {
    fit_stats$outcome <- dplyr::recode(fit_stats$outcome, !!!outcome_labels)
  }
  fit_stats$outcome <- factor(fit_stats$outcome, levels = lvls)

  fit_long <- tidyr::pivot_longer(
    fit_stats,
    dplyr::any_of(gof_stats),
    names_to  = "stat",
    values_to = "value"
  )

  fit_long$label <- mapply(
    function(s, v) sprintf(gof_formats[s], v),
    fit_long$stat, fit_long$value,
    USE.NAMES = FALSE
  )

  stat_labels <- gof_labels[gof_stats]
  fit_long$stat <- factor(fit_long$stat, levels = gof_stats, labels = stat_labels)

  p_gof <- ggplot2::ggplot(
    fit_long,
    ggplot2::aes(
      x     = .data$stat,
      y     = factor(.data$test_start),
      label = .data$label
    )
  ) +
    ggplot2::geom_text(size = 2.5) +
    ggh4x::facet_nested_wrap(
      outcome ~ "",
      ncol      = 1,
      nest_line = TRUE,
      axes      = "x",
      strip     = ggh4x::strip_nested(
        background_x = ggh4x::elem_list_rect(fill = c("grey85", "white")),
        text_x       = ggh4x::elem_list_text(face = c("bold", "plain")),
        by_layer_x   = TRUE
      )
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )

  pw <- patchwork::wrap_plots(p_coef, p_gof, widths = widths) &
    ggplot2::theme(axis.title.x = ggplot2::element_blank())
  pw
}

#' Extract across-draw coefficients from a fitted system
#'
#' Walks the stored coefficient draws of a [fit_system()] output and returns a
#' tidy table of the regression coefficients that were actually used in the
#' simulation, tagged by draw and outcome. Only `linear` and `glm` models
#' carry coefficients (via [broom::tidy()]); other model types are skipped.
#'
#' @param fitted_system An `endogenr_fitted_system` from [fit_system()].
#'
#' @return A `data.table` with columns `.draw` (integer, `1..nsim`),
#'   `.window_start`/`.window_end` (the training window that produced the
#'   draw), `outcome`, `term`, `estimate`, `std.error`, `statistic`, and
#'   `p.value`. For shared/non-refit fits the window is the full
#'   `[train_start, test_start - 1]` range. Zero rows if no model carries
#'   coefficients.
#' @seealso [fit_system()], [plot_coefficients()], [plot_estimates()]
#' @family postprocess
#' @export
get_coefficients <- function(fitted_system) {
  if (!inherits(fitted_system, "endogenr_fitted_system")) {
    stop("`fitted_system` must be the output of fit_system().", call. = FALSE)
  }

  # Window for fits that did not draw a random one (subset = NULL): the full
  # training range [train_start, test_start - 1].
  full_start <- fitted_system$train_start
  full_end   <- fitted_system$test_start - 1L

  draws <- fitted_system$fitted_draws
  rows <- vector("list", length(draws))

  for (d in seq_along(draws)) {
    draw_rows <- list()
    for (model in draws[[d]]) {
      if (is.null(model) || is.null(model$coefs)) next
      win <- model$subset
      ct <- data.table::as.data.table(model$coefs)
      ct[, `:=`(
        .draw = d,
        .window_start = if (is.null(win$start)) full_start else win$start,
        .window_end   = if (is.null(win$end))   full_end   else win$end,
        outcome = model$outcome
      )]
      draw_rows[[length(draw_rows) + 1L]] <- ct
    }
    if (length(draw_rows) > 0L) {
      rows[[d]] <- data.table::rbindlist(draw_rows, use.names = TRUE, fill = TRUE)
    }
  }

  out <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (nrow(out) == 0L) return(out)

  lead <- c(".draw", ".window_start", ".window_end", "outcome", "term",
            "estimate", "std.error", "statistic", "p.value")
  data.table::setcolorder(out, c(intersect(lead, names(out)),
                                 setdiff(names(out), lead)))
  out[]
}

#' Plot the across-draw distribution of system coefficients
#'
#' Visualises the spread of each regression coefficient across the `nsim`
#' draws stored by [fit_system()]. With refit specs (`min_window` set), the
#' draws differ and the boxplots show genuine sampling spread; without refits
#' every draw is identical and each box collapses to a point. The intercept is
#' dropped. Built on [get_coefficients()], so it covers `linear`/`glm` models.
#'
#' @param fitted_system An `endogenr_fitted_system` from [fit_system()].
#' @param outcome_labels An optional named character vector mapping outcome
#'   variable names to display labels. If `NULL`, raw outcome names are used.
#' @param base_size Numeric. Base font size passed to [ggplot2::theme_bw()].
#'   Default `9`.
#'
#' @return A `ggplot` object.
#' @seealso [fit_system()], [get_coefficients()], [plot_estimates()]
#' @family postprocess
#' @export
plot_coefficients <- function(fitted_system, outcome_labels = NULL, base_size = 9) {
  coefs <- get_coefficients(fitted_system)
  if (nrow(coefs) == 0L) {
    stop("No coefficients to plot: only `linear`/`glm` models carry coefficients.",
         call. = FALSE)
  }

  coefs <- coefs[coefs$term != "(Intercept)", ]
  if (nrow(coefs) == 0L) {
    stop("No non-intercept coefficients to plot.", call. = FALSE)
  }

  oc <- as.character(coefs$outcome)
  if (!is.null(outcome_labels)) {
    hit <- oc %in% names(outcome_labels)
    oc[hit] <- outcome_labels[oc[hit]]
    lvls <- unname(outcome_labels)
  } else {
    lvls <- unique(oc)
  }
  coefs[, outcome := factor(oc, levels = lvls)]

  ggplot2::ggplot(coefs, ggplot2::aes(x = .data$estimate, y = "")) +
    ggplot2::geom_vline(xintercept = 0, linetype = 5) +
    ggplot2::geom_boxplot() +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$outcome, .data$term),
      scales = "free_x"
    ) +
    ggplot2::labs(x = "Estimate (across draws)", y = NULL) +
    ggplot2::theme_bw(base_size = base_size)
}
