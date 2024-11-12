source("R/functions.R")

#### test the select_col_per_row function ####
mm <- matrix(1:12, nrow = 3)
cols_selected <- c(1,3,4)
expected_output <- c(1, 8, 12)
select_col_per_row(mm, cols_selected) == expected_output

#### test the prediction interval function ####
df <- data.frame(x = 1:10, y = 2 + 3*(1:10) + rnorm(10))
m <- lm(y ~ x, df)
lmpred <- predict(m, se.fit = TRUE)
result <- replicate(10000, getpi(lmpred))
row_quantiles_apply <- function(mat, probs = c(0.025, 0.5, 0.975)) {
  apply(mat, 1, quantile, probs = probs)
}
test <- predict(m, interval = "prediction")[,c(2,1,3)]
result <- row_quantiles_apply(result) |> t()
all(abs(result - test) <= 0.1)

library(dplyr)
df <- endogenr::example_data
df <- tsibble::as_tsibble(df, key = "gwcode", index = "year")

create_panel_frame(gdppc_grwt ~ lag(gdppc_grwt), data = df)
formula <- gdppc_grwt ~ lag(gdppc_grwt)
# data <- df
# groupvar <- "gwcode"
# timevar <- "year"
# test_start <- 1990
# train_start <- 1970
# horizon <- 12
# inner_sims <- 10


e1 <- gdppc_grwt ~ lag(yjbest) + lag(gdppc_grwt) + lag(log(gdppc)) + lag(psecprop) + lag(zoo::rollmean(gdppc_grwt, k = 3, fill = NA, align = "right"))
c1 <- yjbest ~ lag(yjbest) + lag(log(gdppc)) + lag(log(population)) + lag(psecprop) + lag(dem) + lag(gdppc_grwt) + lag(zoo::rollmean(yjbest, k = 5, fill = NA, align = "right"))
d1 <- dem ~ lag(dem) + lag(gdppc_grwt) + lag(log(gdppc)) + lag(yjbest) + lag(psecprop) + lag(zoo::rollmean(dem, k = 3, fill = NA, align = "right"))

model_system <- list(
  build_model("deterministic",formula = gdppc ~ I(abs(lag(gdppc)*(1+gdppc_grwt)))),
  build_model("deterministic", formula = gdp ~ I(abs(gdppc*population))),
  build_model("parametric_distribution", formula = ~gdppc_grwt, distribution = "norm"),
  build_model("linear", formula = c1, boot = "resid"),
  build_model("linear", formula = d1, boot = "resid"),
  build_model("exogen", formula = ~psecprop),
  build_model("exogen", formula = ~population)
)

simulator_setup <- setup_simulator(models = model_system,
                data = df,
                train_start = 1970,
                test_start = 1990,
                horizon = 2,
                groupvar = "gwcode",
                timevar = "year",
                inner_sims = 2)

#simulator_setup$execution_order <- c("psecprop", "population", "gdppc_grwt", "yjbest", "dem", "gdppc", "gdp")

set.seed(42)
res <- simulate_endogenr(nsim = 2, simulator_setup = simulator_setup)

unique(res$.id) |> as.numeric()
unique(res$sim)


simdf$.id <- 1
simdf <- simdf |> dplyr::filter(year >= 1990)

sim_to_dist(res, "gdppc_grwt")

plotsim(res, "gdppc_grwt", c(2, 20, 530), df)

sapply(model_system, function(x) parse_formula(x)$outcome)


models <- lapply(model_system, build_model)

#process_independent_models(simdf, models)

m <- linearmodel(e1, boot = NULL, data = train)
m <- linearmodel(e1, boot = "resid", data = train)
pred <- m$fitted |> predict(interval = "prediction")
pred3 <- m$fitted |> predict(se.fit = T)
pred1 <- predict(m, simdf, 1990)
pred2 <- predict(m, simdf, 1990, what = "expectation")

microbenchmark::microbenchmark(predict(m, simdf, 1990, what = "pi"))

summary(pred1[["gdppc_grwt"]])
summary(pred2[["gdppc_grwt"]])


dt <- data.table::data.table(
  group = rep(c("A", "B"), each = 5),
  time = rep(1:5, times = 2),
  x = 1:10,
  y = 31:40
)
data.table::setkey(dt, group)
data.table::setindex(dt, time)
res <- create_panel_frame(y ~ lag(x), dt)$data
anyNA(res$lag_x) # this should be true (and that only happens when lag is dplyr::lag)


dm <- deterministicmodel(gdppc ~ I(abs(lag(gdppc)*(1+gdppc_grwt))))
predict(dm, 1989, simdf)

em <- exogenmodel(~psecprop, 2010, newdata = df, inner_sims = 10)
predict(em)

pm <- parametric_distribution_model(~gdppc_grwt, distribution = "norm", data = train)
predict(pm, data = simdf, test_start = 1990)







