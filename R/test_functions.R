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


simdf <- prepare_simulation_data(df, "gwcode", "year", 1970, 1990, 32)
train <- df |> dplyr::filter(year <= 1990)



e1 <- gdppc_grwt ~ lag(yjbest) + lag(gdppc_grwt) + lag(log(gdppc)) + lag(psecprop) + lag(zoo::rollmean(gdppc_grwt, k = 3, fill = NA, align = "right"))
c1 <- yjbest ~ lag(yjbest) + lag(log(gdppc)) + lag(log(population)) + lag(psecprop) + lag(dem) + lag(gdppc_grwt) + lag(zoo::rollmean(yjbest, k = 5, fill = NA, align = "right"))
d1 <- dem ~ lag(dem) + lag(gdppc_grwt) + lag(log(gdppc)) + lag(yjbest) + lag(psecprop) + lag(zoo::rollmean(dem, k = 3, fill = NA, align = "right"))

create_panel_frame(e1, train)

test_start <- 1990

model_system <- list(
  deterministicmodel(gdppc ~ I(abs(lag(gdppc)*(1+gdppc_grwt)))),
  deterministicmodel(gdp ~ I(abs(gdppc*population))),
  parametric_distribution_model(~gdppc_grwt, distribution = "norm", data = train),
  linearmodel(c1, boot = "resid", data = train),
  linearmodel(d1, boot = "resid", data = train),
  exogenmodel(~psecprop, impute_from = test_start, newdata = df),
  exogenmodel(~population, impute_from = test_start, newdata = df)
)

dependency_graph <-  igraph::make_empty_graph(directed = TRUE)
for(model in model_system){
  dependency_graph <- update_dependency_graph(model, dependency_graph)
}

execution_order <- get_execution_order(dependency_graph)


model_system[sapply(model_system, function(x) !x$independent)]

em <- exogenmodel(~psecprop, 1990, newdata = df)

simdf |> dplyr::rows_patch(em$fitted, by = c("gwcode", "year")) |> dplyr::filter(year == 1990)

start <- Sys.time()
simdf <- process_independent_models(simdf, model_system, test_start)
simdf <- process_dependent_models(simdf, model_system, test_start, 12, execution_order)
Sys.time() - start



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

em <- exogenmodel(~psecprop, 2010, newdata = df)
predict(em)

pm <- parametric_distribution_model(~gdppc_grwt, distribution = "norm", data = train)
predict(pm, data = simdf, test_start = 1990)







