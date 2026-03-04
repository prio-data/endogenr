# ==============================================================================
# comparison.R
#
# Compares two approaches to long-horizon GDP forecasting over H = 12 years:
#
#   1. Dynamic simulation (endogenr): a system of three endogenous equations
#      for growth (ΔY), conflict (X), and democracy (D) is simulated year-by-
#      year, compounding uncertainty over the horizon.
#
#   2. Long-horizon prediction (longhorizon): for each horizon h, a direct OLS
#      model predicts GDP at t+h from baseline covariates at t — avoiding the
#      need to forecast inherently uncertain events like conflict.
#
# The central question: does conditioning on observed baseline conflict and
# democracy (LH) outperform propagating them forward through the simulation?
#
# Evaluation: CRPS in asinh(GDP per capita) space, cross-validated over three
# forecast origins (2000, 2005, 2010), horizons h = 1 … 12.
# ==============================================================================

library(endogenr)
library(dplyr)
library(tsibble)
library(zoo)

# ── 1. Parameters ─────────────────────────────────────────────────────────────

TRAIN_START <- 1970          # earliest training observation
HORIZONS    <- 1:12          # forecast horizons (years ahead)
TEST_STARTS <- c(2000, 2005, 2010)  # CV fold origins; also LH train_end
H           <- max(HORIZONS)

NSIM        <- 25L           # outer simulations  ← increase for final runs
INNER_SIMS  <- 10L           # inner draws per simulation
MIN_WINDOW  <- 10L           # minimum bootstrap training window (years)

GROUPVAR <- "gwcode"
TIMEVAR  <- "year"

# ── 2. Data ───────────────────────────────────────────────────────────────────

df <- endogenr::example_data |>
  tsibble::as_tsibble(key = GROUPVAR, index = TIMEVAR) |>
  # Pre-compute asinh(gdppc) — used as the evaluation scale for both approaches
  dplyr::mutate(asinh_gdppc = asinh(gdppc))

# ── 3. Simulation model system ────────────────────────────────────────────────
#
# Equations 1–3 of the paper. Each endogenous variable (ΔY, X, D) is predicted
# one step ahead from the previous period's values, a 5-year rolling mean
# (medium-run dynamics), and the other two endogenous variables.
# Population and education are treated as exogenous (future values supplied).
# GDP per capita is updated deterministically from the simulated growth rate.
#
#   ΔY_{t+1} = β₁ ΔY_t  + β₂ Ȳ_t  + γ X_t  + δ D_t  + θ Z_t + ε^Y
#   X_{t+1}  = β₁ X_t   + β₂ X̄_t  + γ ΔY_t + δ D_t  + θ Z_t + ε^X
#   D_{t+1}  = β₁ D_t   + β₂ D̄_t  + γ X_t  + δ ΔY_t + θ Z_t + ε^D
#
# where Ȳ_t = rollmean(Y, k = 5) and Z_t = {log(population), psecprop}.
# All three equations use residual bootstrap to propagate estimation uncertainty.

f_growth <- gdppc_grwt ~
  lag(gdppc_grwt) +
  lag(zoo::rollmeanr(gdppc_grwt, k = 5, fill = NA)) +
  lag(yjbest) +
  lag(dem) +
  lag(log(population)) +
  lag(psecprop)

f_conflict <- yjbest ~
  lag(yjbest) +
  lag(zoo::rollmeanr(yjbest, k = 5, fill = NA)) +
  lag(gdppc_grwt) +
  lag(dem) +
  lag(log(population)) +
  lag(psecprop)

f_democracy <- dem ~
  lag(dem) +
  lag(zoo::rollmeanr(dem, k = 5, fill = NA)) +
  lag(yjbest) +
  lag(gdppc_grwt) +
  lag(log(population)) +
  lag(psecprop)

model_system <- list(
  build_model("linear", formula = f_growth,    boot = "resid"),
  build_model("linear", formula = f_conflict,  boot = "resid"),
  build_model("linear", formula = f_democracy, boot = "resid"),
  # Deterministic GDP update: Y_t = |Y_{t-1} * (1 + ΔY_t)|
  build_model("deterministic",
              formula = gdppc ~ I(abs(lag(gdppc) * (1 + gdppc_grwt)))),
  build_model("exogen", formula = ~population),
  build_model("exogen", formula = ~psecprop)
)

# ── 4. Long-horizon formulas ──────────────────────────────────────────────────
#
# Equation 4 of the paper. For each h, a direct OLS model estimates:
#
#   asinh(Y_{i,t+h}) = α^h + β^h_Y asinh(Y_{it})
#                           + β^h_X X_{it}          [full only]
#                           + β^h_D D_{it}          [full only]
#                           + γ^h Z_{it}            [full only]
#
# baseline: GDP-only benchmark (captures convergence dynamics)
# full:     adds observed conflict intensity and democracy at the baseline year
#
# Key asymmetry vs. simulation: X_{it} and D_{it} are OBSERVED at t, so the
# model avoids forecasting inherently unpredictable future conflict.

lh_formulas <- list(
  baseline = asinh(gdppc) ~ asinh(gdppc),
  full     = asinh(gdppc) ~ asinh(gdppc) + yjbest + dem + psecprop
)

# ── 5. Cross-validation: simulation ──────────────────────────────────────────

message("=== Running simulation cross-validation ===")

sim_results <- list()

for (ts in TEST_STARTS) {
  message("  Simulation fold: test_start = ", ts)

  sim_setup <- setup_simulator(
    models      = model_system,
    data        = df,
    train_start = TRAIN_START,
    test_start  = ts,
    horizon     = H,
    groupvar    = GROUPVAR,
    timevar     = TIMEVAR,
    inner_sims  = INNER_SIMS,
    min_window  = MIN_WINDOW
  )

  set.seed(42)
  res <- simulate_endogenr(
    nsim             = NSIM,
    simulator_setup  = sim_setup,
    parallel         = TRUE
  )

  sim_results[[as.character(ts)]] <- res |>
    dplyr::mutate(asinh_gdppc = asinh(gdppc))
}

# Compute accuracy per horizon per fold.
# get_accuracy() expects a single-period tsibble; we iterate over each
# forecast year and map year → horizon.

message("  Computing simulation accuracy by horizon...")

sim_acc_rows <- list()

for (ts in TEST_STARTS) {
  res <- sim_results[[as.character(ts)]]

  for (h in HORIZONS) {
    year_h <- ts + h

    res_h <- res |>
      tsibble::as_tsibble(key = c(GROUPVAR, ".sim"), index = TIMEVAR) |>
      dplyr::filter(!!rlang::sym(TIMEVAR) == year_h)

    if (nrow(res_h) == 0L) next

    acc_h <- get_accuracy(res_h, "asinh_gdppc", df) |>
      dplyr::summarise(crps = mean(crps, na.rm = TRUE),
                       mae  = mean(mae,  na.rm = TRUE))

    sim_acc_rows[[paste(ts, h, sep = "_")]] <- dplyr::mutate(
      acc_h, test_start = ts, horizon = h
    )
  }
}

sim_acc <- dplyr::bind_rows(sim_acc_rows) |>
  dplyr::group_by(horizon) |>
  dplyr::summarise(crps = mean(crps, na.rm = TRUE),
                   mae  = mean(mae,  na.rm = TRUE),
                   .groups = "drop")

# ── 6. Cross-validation: long-horizon ────────────────────────────────────────

message("=== Running long-horizon cross-validation ===")

# Truth in asinh space — must match the scale of the LH draws
truth_asinh <- df |>
  dplyr::as_tibble() |>
  dplyr::select(dplyr::all_of(c(GROUPVAR, TIMEVAR, "asinh_gdppc"))) |>
  tsibble::as_tsibble(key = GROUPVAR, index = TIMEVAR)

lh_all_forecasts <- list()

for (ts in TEST_STARTS) {
  message("  Long-horizon fold: test_start = ", ts)

  lh_setup <- setup_long_horizon(
    data      = df,
    formulas  = lh_formulas,
    horizons  = HORIZONS,
    groupvar  = GROUPVAR,
    timevar   = TIMEVAR,
    train_end = ts,
    boot      = "resid"
  )

  set.seed(42)
  lh_fcst <- forecast_long_horizon(
    lh_setup   = lh_setup,
    data       = df,
    test_start = ts,
    nsim       = NSIM,
    inner_sims = INNER_SIMS
  )

  lh_all_forecasts[[as.character(ts)]] <- lh_fcst
}

lh_combined <- dplyr::bind_rows(lh_all_forecasts)

lh_acc <- get_lh_accuracy(
  lh_forecasts = lh_combined,
  truth        = truth_asinh,
  outcome      = "asinh_gdppc",
  groupvar     = GROUPVAR,
  timevar      = TIMEVAR
)

# ── 7. Compare ────────────────────────────────────────────────────────────────

comparison <- compare_approaches(lh_acc, sim_acc)

cat("\n=== Accuracy: CRPS in asinh(GDP per capita) space ===\n\n")
print(
  comparison |>
    dplyr::select(approach, variant, horizon, crps, mae) |>
    dplyr::arrange(horizon, approach, variant)
)

# ── 8. Plot ───────────────────────────────────────────────────────────────────

if (requireNamespace("ggplot2", quietly = TRUE)) {
  p <- ggplot2::ggplot(
    comparison |>
      dplyr::mutate(
        label = dplyr::case_when(
          approach == "simulation"   ~ "Dynamic simulation",
          variant  == "baseline"     ~ "LH: baseline (GDP only)",
          variant  == "full"         ~ "LH: full (+ conflict + democracy)",
          TRUE                       ~ paste(approach, variant, sep = " / ")
        )
      ),
    ggplot2::aes(x = horizon, y = crps, colour = label, linetype = approach)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_x_continuous(breaks = HORIZONS) +
    ggplot2::scale_linetype_manual(
      values = c("long_horizon" = "solid", "simulation" = "dashed"),
      guide  = "none"
    ) +
    ggplot2::labs(
      title   = "Dynamic simulation vs. long-horizon prediction",
      subtitle = paste0(
        "CRPS by forecast horizon, CV over test_starts = ",
        paste(TEST_STARTS, collapse = ", ")
      ),
      x       = "Horizon h (years ahead)",
      y       = "CRPS  [asinh(GDP per capita)]",
      colour  = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")

  print(p)

  if (!is.null(getOption("endogenr.save_plots")) &&
        getOption("endogenr.save_plots")) {
    ggplot2::ggsave("paper/comparison_crps.pdf", p,
                    width = 7, height = 4.5)
  }
}
