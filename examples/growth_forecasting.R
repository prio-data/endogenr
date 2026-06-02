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
e1 <- gdppc_grwt ~ region + lag(dem) + lag(yjbest) + lag(log(gdppc)) + lag(psecprop)

model_system <- list(
  build_model("deterministic",formula = gdppc ~ I(abs(lag(gdppc)*(1+gdppc_grwt)))),
  build_model("linear", formula = e1, boot = "resid"),
  build_model("exogen", formula = ~yjbest),
  build_model("exogen", formula = ~dem),
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
ggplot(my_coefs, aes(x = .window_end, y = estimate)) + geom_point() + facet_wrap(~term)

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
e2 <- gdppc_grwt ~ region + lag(dem) + lag(log(gdppc))
e3 <- gdppc_grwt ~ region + lag(yjbest) + lag(log(gdppc)) + lag(psecprop)

systems <- vary_model(model_system, list(
  e1 = build_model("linear", formula = e1, boot = "resid"),
  e2 = build_model("linear", formula = e2, boot = "resid"),
  e3 = build_model("linear", formula = e3, boot = "resid")
))

# Parallelise ACROSS experiments (the default): set one plan; each experiment's
# inner fit/simulate runs sequentially on its worker. Cross 3 model variants
# with 2 forecast starts -> 6 experiments. Keep test_start + horizon - 1 within
# the data range (<= 2024).
start <- Sys.time()
future::plan(future::multisession, workers = 8)
set.seed(42)
progressr::with_progress({
  exp_res <- run_experiments(
    data        = df,
    models      = systems,
    train_start = 1970,
    test_start  = c(2006, 2010),
    horizon     = 12,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 30,
    nsim        = 64,
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
