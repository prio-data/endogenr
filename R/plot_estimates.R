#' Plot regression estimates from dynamic simulation models
#'
#' Creates a coefficient plot with pointranges for all fitted linear models
#' in a cross-validated simulation setup, optionally appended with a table
#' of goodness-of-fit statistics (R-squared and N).
#'
#' @param models A named list (keyed by test start year) of
#'   [setup_simulator()] outputs. Typically loaded from a saved
#'   `models.rds` file.
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
  estimates <- purrr::map_dfr(names(models), function(ts) {
    purrr::map_dfr(models[[ts]]$fitted_models, function(m) {
      if (!is.null(m$coefs)) {
        m$coefs |>
          dplyr::mutate(
            test_start = as.integer(ts),
            outcome    = m$outcome,
            .before    = 1
          )
      }
    })
  })

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
  fit_stats <- purrr::map_dfr(names(models), function(ts) {
    purrr::map_dfr(models[[ts]]$fitted_models, function(m) {
      if (!is.null(m$gof)) {
        m$gof |>
          dplyr::select(dplyr::any_of(gof_stats)) |>
          dplyr::mutate(
            test_start = as.integer(ts),
            outcome    = m$outcome,
            .before    = 1
          )
      }
    })
  })

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

  fit_long$label <- purrr::map2_chr(
    fit_long$stat, fit_long$value,
    function(s, v) sprintf(gof_formats[s], v)
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
