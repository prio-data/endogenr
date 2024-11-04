
<!-- README.md is generated from README.Rmd. Please edit that file -->

# endogenr

<!-- badges: start -->
<!-- badges: end -->

The goal of endogenr is to make it easy to simulate dynamic systems from
regression models, mathematical equations, and exogenous inputs (either
based on a stochastic distribution, or given by some data). We assume a
panel-data structure, with two variables denoting the time and the space
dimensions.

The simulator works by identifying the dependency graph of the models
added to the system, and deriving the order of calculation based on this
graph.

It works with parallel simulation using the future::multisession
approach.

## Installation

You can install the development version of endogenr from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("prio-data/endogenr")
```

## Example

``` r
library(endogenr)

train_start = 1970; test_start = 2010; timevar = "year"; groupvar = "gwcode"
df <- example_data
train <- df[time >= train_start & time < test_start, env = list(time = timevar)]

e1 <- gdppc_grwt ~ dplyr::lag(yjbest) + dplyr::lag(gdppc_grwt) + dplyr::lag(log(gdppc)) + dplyr::lag(psecprop) + dplyr::lag(zoo::rollmean(gdppc_grwt, k = 3, fill = NA, align = "right"))
c1 <- yjbest ~ dplyr::lag(yjbest) + dplyr::lag(log(gdppc)) + dplyr::lag(log(population)) + dplyr::lag(psecprop) + dplyr::lag(dem) + dplyr::lag(gdppc_grwt) + dplyr::lag(zoo::rollmean(yjbest, k = 5, fill = NA, align = "right"))
d1 <- dem ~ dplyr::lag(dem) + dplyr::lag(gdppc_grwt) + dplyr::lag(log(gdppc)) + dplyr::lag(yjbest) + dplyr::lag(psecprop) + dplyr::lag(zoo::rollmean(dem, k = 3, fill = NA, align = "right"))

m1 <- StochasticStaticModel$new(~gdppc_grwt, "t_ls", fit_args = list(start = list(df = 1, mu = mean(train$gdppc_grwt), sigma = sd(train$gdppc_grwt))), censor = list("upper" = 0.8, "lower" = -0.8))
m2 <- DeterministicModel$new(gdppc ~ I(abs(dplyr::lag(gdppc)*(1+gdppc_grwt))))
m3 <- ExogenModel$new(~population, exogendata = df)
m4 <- DeterministicModel$new(gdp ~ I(abs(gdppc * population)))
m5 <- StochasticStaticModel$new(~round(yjbest), "nbinom", fit_args = list(discrete = TRUE, method = "mme"))
m6 <- LinearRegressionModel$new(d1, boot = T)
m7 <- ExogenModel$new(~psecprop, exogendata = df)

models <- list(m1, m2, m3, m4, m5, m6, m7)

simulator <- EndogenousSystem$new(
  df, 
  timevar = timevar, 
  groupvar = groupvar,
  train_start = train_start,
  test_start = test_start, 
  min_window = 10
)

models <- list(m1, m2, m3, m4, m5, m6, m7)
for(model in models){
  simulator$add_model(model)
}

simulator$simulate(horizon = 12, nsim = 2, parallel = F)
simulator$plot("gdppc_grwt", units = c(2, 101, 365, 710, 475))
```

<img src="man/figures/README-example-1.png" width="100%" />