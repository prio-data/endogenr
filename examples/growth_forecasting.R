library(endogenr)
library(paneltools)
library(countrycode)
df <- endogenr::example_data
df <- as_panel(df, "gwcode", "year")
df <- paneltools::p_balance(df)
df <- df[year >= 1970 & year <= 2024]
df <- paneltools::p_complete(df) # only keep complete data

df$region <- countrycode::countrycode(
  df$gwcode,
  origin = "gwn",
  destination = "region",
  custom_match = c("6" = "North America",
                   "699" = "Middle East & North Africa",
                   "816" = "South Asia")
  )


yj1 <- scales::yj_trans(p = 0.1)
df$yjbest <- yj1$transform(df$best)
df$dem <- plogis(df$v2x_polyarchy + 0.01)
df$region <- factor(df$region)

df[!complete.cases(df),]

df$region <- factor(df$region)
e1 <- gdppc_grwt ~ lag(log(gdppc))

model_system <- list(
  build_model("deterministic",formula = gdppc ~ I(abs(lag(gdppc)*(1+gdppc_grwt)))),
  build_model("linear", formula = e1, boot = "resid"),
  build_model("exogen", formula = ~dem),
  build_model("exogen", formula = ~yjbest),
  build_model("exogen", formula = ~psecprop),
  build_model("exogen", formula = ~population),
  build_model("exogen", formula = ~region)
)

sys <- setup_system(
  models      = model_system,
  data        = df,
  train_start = 1970,
  test_start  = 2010,
  horizon     = 12,
  groupvar    = "gwcode",
  timevar     = "year",
  inner_sims  = 30,
  min_window  = 20
)

start <- Sys.time()
future::plan(future::multisession, workers = 8)
set.seed(42)
fit <- fit_system(sys, nsim = 128)
future::plan(future::sequential)
Sys.time() - start

my_coefs <- get_coefficients(fit)
my_coefs$window_size <- my_coefs$.window_end - my_coefs$.window_start

library(ggplot2)
ggplot(my_coefs, aes(x = window_size, y = p.value)) + geom_point()
ggplot(my_coefs, aes(x = .window_end, y = estimate)) + geom_point() + facet_wrap(~outcome + term) + geom_smooth()

# Coefficient-trajectory forecasting ----------------------------------------
# get_coefficients() shows where each coefficient sat across the random
# bootstrap windows. forecast_coefficients() turns that idea into a forward
# projection: it refits the linear model on a deterministic grid of expanding
# (or rolling) windows to get a clean coefficient time series, then projects
# each coefficient over the horizon as a random walk (the empirical analogue of
# a time-varying-parameter / state-space model). Use method = "drift" to let the
# central path continue the recent trend; "rw" (default) only fans out.
cf <- forecast_coefficients(fit, method = "rw")        # or window = "rolling"
plot_coefficient_forecast(cf)                          # observed path + forecast fan
plot_coefficient_forecast(fit, method = "drift")       # forecast directly from a fit

start <- Sys.time()
future::plan(future::multisession, workers = 8)
set.seed(42)
progressr::with_progress({
  res <- simulate_system(fit)
})
future::plan(future::sequential)
Sys.time() - start

sim_acc <- get_accuracy(res, "gdppc", df, test_start = 2010, by = "horizon")
plotsim(res, "gdppc", c(2, 20, 530), df) + ggplot2::facet_wrap(~gwcode, scales = "free_y")

# Sliding-window simulation -------------------------------------------------
# The random-window fit above (min_window) draws a random start AND end per
# draw, which confounds recency with sample size. A *sliding* fit instead
# refits each linear/glm/heterolm spec on a deterministic grid of window-end
# anchors (the same grid as forecast_coefficients()), so each fit is tied to a
# point in time. simulate_system() then decides, per forecast step, which window
# drives the coefficients via `window_policy`:
#   - "latest": always the most recent window (the forecast origin);
#   - "equal" : every window with equal probability;
#   - "decay" : recent windows dominate proximate horizons, flattening toward
#               uniform further out (tune with `decay`); or pass a custom
#               `weights` function(h, age) / horizon-by-window matrix.
future::plan(future::multisession, workers = 8)
set.seed(42)
fit_slide <- fit_system(sys, nsim = 128, window = "rolling", width = 20, step = 1)

# Inspect how the coefficients move across windows BEFORE picking a policy.
slide_coefs <- get_coefficients(fit_slide)
ggplot(slide_coefs[term != "(Intercept)"],
       aes(x = .window_end, y = estimate)) +
  geom_point(alpha = 0.2) + geom_smooth() +
  facet_wrap(~ outcome + term, scales = "free_y")

set.seed(42)
progressr::with_progress({
  res_latest <- simulate_system(fit_slide, window_policy = "latest")
  res_equal  <- simulate_system(fit_slide, window_policy = "equal")
  res_decay  <- simulate_system(fit_slide, window_policy = "decay", decay = 0.5)
})
future::plan(future::sequential)

acc_latest <- get_accuracy(res_latest, "gdppc", df, test_start = 2010, by = "horizon")
acc_equal  <- get_accuracy(res_equal,  "gdppc", df, test_start = 2010, by = "horizon")
acc_decay  <- get_accuracy(res_decay,  "gdppc", df, test_start = 2010, by = "horizon")

ggplot() +
  geom_line(data = acc_latest[horizon >= 1], aes(horizon, crps), colour = "black") +
  geom_line(data = acc_equal[horizon >= 1],  aes(horizon, crps), colour = "red") +
  geom_line(data = acc_decay[horizon >= 1],  aes(horizon, crps), colour = "blue") +
  labs(subtitle = "CRPS by horizon: latest (black) / equal (red) / decay (blue)")

lh_setup <- endogenr::setup_long_horizon(
  data      = df,
  formulas  = list("y_h|y" = lead_horizon(log(gdppc)) ~ factor(region) + log(gdppc),
                   "y_h|y+c+d" = lead_horizon(log(gdppc)) ~ factor(region) + log(gdppc) + yjbest + dem,
                   "y_h|y+c+d(region)" = lead_horizon(log(gdppc)) ~ factor(region)*dem + log(gdppc) + yjbest),
  horizons  = 1:20,
  groupvar  = "gwcode",
  timevar   = "year",
  test_start = 2004,
  boot      = "resid"
)

lh_res <- forecast_long_horizon(lh_setup, data = df, nsim = 10, inner_sims = 10)
lh_acc <- get_lh_accuracy(lh_res, df, scale = "native", inverse = exp)
lh_acc[,mean(crps),.(variant, horizon)]

library(ggplot2)
ggplot() +
  geom_line(data = sim_acc[horizon>=1], mapping = aes(x = horizon, y = crps), color = "black") +
  geom_line(data = lh_acc[variant == "y_h|y", .(crps = mean(crps)), horizon], mapping = aes(x = horizon, y = crps), color = "red") +
  geom_line(data = lh_acc[variant == "y_h|y+c+d", .(crps = mean(crps)), horizon], mapping = aes(x = horizon, y = crps), color = "blue") +
  geom_line(data = lh_acc[variant == "y_h|y+c+d(region)", .(crps = mean(crps)), horizon], mapping = aes(x = horizon, y = crps), color = "green")

compare_approaches(lh_acc,
                   sim_acc[horizon >= 1])[, .(crps = mean(crps),
                                              mae = mean(mae),
                                              winkler = mean(winkler)),
                                          .(approach, variant)]



# Experiment grid -----------------------------------------------------------
# Explore variations of the simulation pipeline (model system, train_start,
# and/or test_start) in one call. Each combination is a separate experiment;
# results stack into one long data.table that plugs into get_experiment_accuracy().

# Vary the gdppc_grwt equation: swap just the `linear` spec into the base system.
e1 <- gdppc_grwt ~ lag(log(gdppc))
e2 <- gdppc_grwt ~ lag(log(gdppc)) + lag(psecprop)
e3 <- gdppc_grwt ~ lag(log(gdppc)) + lag(psecprop) + lag(dem)
e4 <- gdppc_grwt ~ lag(log(gdppc)) + lag(psecprop) + lag(yjbest)
e5 <- gdppc_grwt ~ lag(log(gdppc)) + lag(psecprop) + lag(yjbest) + lag(dem)

e1r <- gdppc_grwt ~ region + lag(log(gdppc))
e2r <- gdppc_grwt ~ region + lag(log(gdppc)) + lag(psecprop)
e3r <- gdppc_grwt ~ region + lag(log(gdppc)) + lag(psecprop) + lag(dem)
e4r <- gdppc_grwt ~ region + lag(log(gdppc)) + lag(psecprop) + lag(yjbest)
e5r <- gdppc_grwt ~ region + lag(log(gdppc)) + lag(psecprop) + lag(yjbest) + lag(dem)


systems <- vary_model(model_system, list(
  e1 = build_model("linear", formula = e1, boot = "resid"),
  e2 = build_model("linear", formula = e2, boot = "resid"),
  d2 = build_model("univariate_fable", formula = dem ~ error("A") + trend("N") + season("N"), method = "ets")
))

# Parallelise ACROSS experiments (the default): set one plan; each experiment's
# inner fit/simulate runs sequentially on its worker. Cross 3 model variants
# with 2 forecast starts -> 6 experiments. Keep test_start + horizon - 1 within
# the data range (<= 2024).
start <- Sys.time()
future::plan(future::multisession, workers = 6)
set.seed(42)
progressr::with_progress({
  exp_res <- run_experiments(
    data        = df,
    models      = systems,
    train_start = 1970,
    test_start  = 2010,
    horizon     = 12,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 10,
    nsim        = 3,
    min_window  = 20
  )
})
future::plan(future::sequential)
Sys.time() - start

# One score block per experiment, each against its own test_start.
exp_acc <- get_experiment_accuracy(exp_res, "gdppc", df)
exp_acc[, .(crps = mean(crps), mae = mean(mae)),
        by = .(.experiment, model, test_start, horizon)]

library(ggplot2)
ggplot(exp_acc[, .(crps = mean(crps)), by = .(model, test_start, horizon)],
       aes(x = horizon, y = crps, colour = factor(model))) +
  geom_line() +
  facet_wrap(~test_start)
