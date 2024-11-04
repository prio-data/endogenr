#' Base Model Class
#'
#' @description
#' A base R6 Class that provides core functionality for modeling operations.
#' This class is not meant to be used directly but serves as a parent class
#' for specific model implementations.
#'
#' @details
#' The Model class implements common functionality such as model initialization,
#' panel data handling, and defines interface methods that child classes must implement.
#'
#' @export
Model <- R6::R6Class(
  "Model",
  public = list(
    #' @field name Name of the model outcome variable
    name = NULL,

    #' @field formula Model formula defining relationships between variables
    formula = NULL,

    #' @field boot Logical indicating whether to use bootstrap resampling
    boot = NULL,

    #' @field censor Optional censoring specification
    censor = NULL,

    #' @description
    #' Create a new model instance
    #' @param formula An object of class formula
    #' @param boot Logical indicating whether to use bootstrap resampling
    #' @param censor Optional censoring specification
    #' @return A new `Model` object
    initialize = function(formula, boot = FALSE, censor = NULL) {
      private$setup_model(formula, boot, censor)
    },

    #' @description
    #' Creates a model frame for panel data
    #' @param frm Formula to apply
    #' @param data Panel data as data.frame or data.table
    #' @return A data.table with the processed panel data
    panel_modelframe = function(frm, data) {
      private$validate_panel_inputs(frm, data)
      private$create_panel_frame(frm, data)
    },

    #' @description
    #' Virtual method for fitting the model
    #' @param data Training data
    #' @param time_parameters Optional time window parameters
    fit = function(data, time_parameters = NULL) {
      stop("'fit' method must be implemented in child classes")
    },

    #' @description
    #' Virtual method for making predictions
    #' @param data Data to make predictions for
    #' @param t Time point for prediction
    predict = function(data, t) {
      stop("'predict' method must be implemented in child classes")
    }
  ),

  private = list(
    setup_model = function(formula, boot, censor) {
      self$name <- rlang::f_lhs(formula) |> as.character()
      self$formula <- formula
      self$boot <- boot
      self$censor <- censor
    },

    validate_panel_inputs = function(frm, data) {
      if (!base::inherits(frm, "formula")) {
        stop("Input 'frm' must be a formula")
      }
      if (!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
      }
      return(data)
    },

    create_panel_frame = function(frm, data) {
      grp <- data.table::key(data)
      idx <- data.table::indices(data)

      frm <- stats::update(frm, paste(c(". ~ .", idx), collapse = "+"))

      # result <- data[, {
      #   mf <- stats::model.frame(frm, data = .SD, na.action = na.pass)
      #   as.list(mf)
      # }, by = grp]


      result <- data[, as.list(stats::model.frame(frm, data = .SD, na.action = NULL)), by = grp]

      data.table::setkeyv(result, grp)
      data.table::setindexv(result, idx)
      return(result)
    }
  )
)

#' Deterministic Model Class
#'
#' @description
#' An R6 Class representing a deterministic model where outcomes are calculated
#' directly from a formula without statistical fitting.
#'
#' @details
#' The DeterministicModel class is used when the relationship between variables
#' is known and can be expressed through a deterministic formula. The formula
#' must be wrapped in I() to ensure proper evaluation.
#'
#' @export
DeterministicModel <- R6::R6Class(
  "DeterministicModel",
  inherit = Model,
  public = list(
    #' @description
    #' Placeholder method for consistency with Model interface
    #' @param data Not used
    #' @param time_parameters Not used
    #' @return Character string indicating no fit needed
    fit = function(data, time_parameters = NULL) {
      "No fit needed for deterministic model."
    },

    #' @description
    #' Applies the deterministic formula to generate predictions
    #' @param data Data.table containing input variables
    #' @param t Time point for prediction
    #' @return Data.table with predictions for time t
    #' @examples
    #' \dontrun{
    #' model <- DeterministicModel$new(y ~ I(2 * x))
    #' pred <- model$predict(data, t = 1)
    #' }
    predict = function(data, t) {
      private$validate_prediction_inputs(data, t)
      result <- private$apply_deterministic_formula(data)
      return(private$filter_prediction_time(result, t))
    }
  ),

  private = list(
    validate_prediction_inputs = function(data, t) {
      if (!data.table::is.data.table(data)) {
        stop("Input 'data' must be a data.table")
      }
      if (is.null(t)) {
        stop("Time parameter 't' must be provided")
      }
    },

    apply_deterministic_formula = function(data) {
      res <- self$panel_modelframe(self$formula, data)
      y_star <- attr(stats::terms(self$formula)[1], "term.labels")

      if (!stringr::str_detect(y_star, "I\\(")) {
        stop("The formula to apply must be the first term and wrapped in I()")
      }

      grp <- data.table::key(data)
      idx <- data.table::indices(data)
      y <- as.character(self$formula[[2]])
      res[, outcome := new_outcome, env = list(outcome = y, new_outcome = y_star)]
      return(res[, .SD, .SDcols = c(grp, idx, y)])
    },

    filter_prediction_time = function(result, t) {
      idx <- data.table::indices(result)
      return(result[time == t, env = list(time = idx)])
    }
  )
)


#' Stochastic Static Model Class
#'
#' @description
#' An R6 Class representing a stochastic model where outcomes are generated
#' from a probability distribution fitted to the data.
#'
#' @details
#' The StochasticStaticModel fits a specified probability distribution to the data
#' and generates predictions by sampling from the fitted distribution. Supports
#' various distributions including normal, Cauchy, Gumbel, Student's t, gamma,
#' Poisson, and negative binomial.
#'
#' @export
StochasticStaticModel <- R6::R6Class(
  "StochasticStaticModel",
  inherit = Model,
  public = list(
    #' @field model The fitted distribution object
    model = NULL,

    #' @field distribution Character string specifying the probability distribution
    distribution = NULL,

    #' @field fit_args Additional arguments passed to fitdistrplus::fitdist()
    fit_args = NULL,

    #' @description
    #' Create a new stochastic static model instance
    #' @param formula Formula specifying the variable to model
    #' @param distribution Character string specifying the probability distribution
    #' @param fit_args List of additional arguments for fitdistrplus::fitdist()
    #' @param censor Optional censoring specification
    #' @return A new `StochasticStaticModel` object
    #' @examples
    #' \dontrun{
    #' model <- StochasticStaticModel$new(
    #'   formula = ~ returns,
    #'   distribution = "norm"
    #' )
    #' }
    initialize = function(formula, distribution, fit_args = NULL, censor = NULL) {
      private$validate_distribution(distribution)
      private$setup_stochastic_model(formula, distribution, fit_args, censor)
    },

    #' @description
    #' Fits the specified probability distribution to the data
    #' @param data Training data
    #' @param time_parameters Optional time window parameters (not used)
    fit = function(data, time_parameters = NULL) {
      self$model <- NULL
      model_data <- private$prepare_fit_data(data)
      self$model <- private$fit_distribution_with_fallback(model_data)
    },

    #' @description
    #' Generates predictions by sampling from the fitted distribution
    #' @param data Data to generate predictions for
    #' @param t Time point for prediction
    #' @return Data.table with predictions for time t
    predict = function(data, t) {
      private$validate_model_fit()
      return(private$generate_predictions(data, t))
    }
  ),

  private = list(
    supported_distributions = c("norm", "cauchy", "gumbel", "t_ls", "gamma", "pois", "nbinom"),

    validate_distribution = function(distribution) {
      if (!is.character(distribution)) {
        stop("Distribution must be specified as a character string")
      }
      if (!distribution %in% private$supported_distributions) {
        stop(sprintf("Unsupported distribution. Must be one of: %s",
                     paste(private$supported_distributions, collapse = ", ")))
      }
    },

    setup_stochastic_model = function(formula, distribution, fit_args, censor) {
      outcome_variable <- base::all.vars(formula)
      if(length(outcome_variable)!=1){
        stop("A stochastic model must have a formula with only one term on rhs and none on the lhs. Consider wrapping function with I()")
      }
      self$name <- outcome_variable
      self$formula <- formula
      self$distribution <- distribution
      self$fit_args <- fit_args
      self$censor <- censor
    },

    prepare_fit_data = function(data) {
      model_formula <- stats::update(self$formula, "~ . -1")
      input_data <- stats::model.matrix(model_formula, data) |>
        data.table::as.data.table()

      if (ncol(input_data) != 1) {
        stop("Formula must yield only one column. Consider wrapping function with I()")
      }

      return(unlist(input_data))
    },

    fit_distribution_with_fallback = function(data) {
      tryCatch({
        if (self$distribution == "t_ls") {
          t_env <- private$setup_t_distribution()
        }
        fit <- private$fit_distribution_internal(data)
        res <- private$create_distribution_object(fit)
        if (self$distribution == "t_ls") {
          on.exit(detach(t_env))
        }
        return(res)

      }, error = function(e) {
        warning("Error in distribution fitting: ", e$message, "\nUsing fallback distribution")
        return(private$create_fallback_distribution(data))
      })
    },

    fit_distribution_internal = function(data) {
      if (is.null(self$fit_args)) {
        fitdistrplus::fitdist(data, self$distribution)
      } else {
        args <- self$fit_args
        args$data <- data
        args$distr <- self$distribution
        do.call(fitdistrplus::fitdist, args)
      }
    },

    setup_t_distribution = function() {
      t_env <- new.env(parent = environment())
      t_env$dt_ls <- function(x, df=1, mu=0, sigma=1) 1/sigma * stats::dt((x - mu)/sigma, df)
      t_env$pt_ls <- function(q, df=1, mu=0, sigma=1) stats::pt((q - mu)/sigma, df)
      t_env$qt_ls <- function(p, df=1, mu=0, sigma=1) stats::qt(p, df)*sigma + mu
      t_env$rt_ls <- function(n, df=1, mu=0, sigma=1) stats::rt(n, df)*sigma + mu
      attach(t_env)
      return(t_env)
    },

    create_distribution_object = function(fit) {
      get_nbinom <- function(size, mu){
        prob <- size / (size + mu)
        distributional::dist_negative_binomial(size = size, prob = prob)
      }

      switch(self$distribution,
             "norm" = distributional::dist_normal(
               mean = fit$estimate["mean"],
               sd = fit$estimate["sd"]
             ),
             "cauchy" = distributional::dist_cauchy(
               location = fit$estimate["location"],
               scale = fit$estimate["scale"]
             ),
             "gumbel" = distributional::dist_gumbel(
               alpha = fit$estimate["alpha"],
               scale = fit$estimate["scale"]
             ),
             "gamma" = distributional::dist_gamma(
               shape = fit$estimate["shape"],
               rate = fit$estimate["rate"]
             ),
             "t_ls" = distributional::dist_student_t(
               df = fit$estimate["df"],
               mu = fit$estimate["mu"],
               sigma = fit$estimate["sigma"]
             ),
             "nbinom" = get_nbinom(fit$estimate["size"], fit$estimate["mu"]),
             "pois" = distributional::dist_poisson(
               lambda = fit$estimate["lambda"]
             ),
             stop("Unsupported distribution")
      )
    },

    validate_model_fit = function() {
      if (is.null(self$model)) {
        stop("Model must be fitted before prediction")
      }
    },

    generate_predictions = function(data, t) {
      idx <- data.table::indices(data)
      grp <- data.table::key(data)

      num_predictions <- data[, unique(group), env = list(group = grp)] |> length()

      pred_data <- data[time == t,
                        outcome := distributional::generate(self$model, num_predictions),
                        env = list(time = idx, outcome = self$name)]

      return(pred_data[time == t,
                       .(group, time, outcome),
                       env = list(time = idx, group = grp, outcome = self$name)])
    }
  )
)

#' Exogenous Model Class
#'
#' @description
#' An R6 Class representing a model that uses external (exogenous) data
#' to generate predictions.
#'
#' @details
#' The ExogenModel class handles cases where predictions are based on
#' external data sources. It requires the input data to be properly
#' grouped and sorted using data.table keys and indices.
#'
#' @export
ExogenModel <- R6::R6Class(
  "ExogenModel",
  inherit = Model,
  public = list(
    #' @field model Processed model data
    model = NULL,

    #' @field exogendata External data used for predictions
    exogendata = NULL,

    #' @field data_diagnostics List containing diagnostic information
    data_diagnostics = NULL,

    #' @description
    #' Create a new exogenous model instance
    #' @param formula Formula specifying how to process exogenous data
    #' @param exogendata Data.table containing exogenous variables
    #' @return A new `ExogenModel` object
    #' @examples
    #' \dontrun{
    #' data <- data.table::as.data.table(external_data)
    #' data.table::setkeyv(data, "group")
    #' data.table::setindexv(data, "time")
    #' model <- ExogenModel$new(
    #'   formula = ~ external_var,
    #'   exogendata = data
    #' )
    #' }
    initialize = function(formula, exogendata) {
      private$validate_inputs(formula, exogendata)
      self$name <- rlang::f_rhs(formula) |> as.character()
      self$formula <- formula
      self$exogendata <- exogendata
      self$data_diagnostics <- list()
    },

    #' @description
    #' Processes exogenous data for prediction
    #' @param data Not used (exogenous data provided at initialization)
    #' @param time_parameters Not used
    fit = function(data = NULL, time_parameters = NULL) {
      self$model <- NULL
      private$fit_exogenous_model()
      private$calculate_diagnostics()
    },

    #' @description
    #' Returns processed exogenous data
    #' @param data Not used
    #' @param t Not used
    #' @return Processed exogenous data
    predict = function(data, t = NULL) {
      private$validate_model()
      return(self$model)
    },

    #' @description
    #' Returns current state of the model
    #' @return List containing model summary and diagnostics
    get_model_state = function() {
      list(
        model_summary = private$summarize_model(),
        data_diagnostics = self$data_diagnostics
      )
    }
  ),

  private = list(
    validate_inputs = function(formula, exogendata) {
      if (!data.table::is.data.table(exogendata)) {
        stop("Exogen input data must be data.table")
      }
      if (is.null(data.table::key(exogendata))) {
        stop("Exogen input data must be a grouped data.table. See data.table::setkeyv()")
      }
      if (is.null(data.table::indices(exogendata))) {
        stop("Exogen input data must be a grouped and sorted (by time) data.table. See data.table::setindexv()")
      }
    },

    fit_exogenous_model = function() {
      grp <- data.table::key(self$exogendata)
      idx <- data.table::indices(self$exogendata)
      model_formula <- stats::update(self$formula, paste(c("~ . -1", grp, idx), collapse = " + "))

      tryCatch({
        self$model <- stats::model.matrix(model_formula, self$exogendata) |>
          data.table::as.data.table()
        data.table::setkeyv(self$model, grp)
        data.table::setindexv(self$model, idx)
      }, error = function(e) {
        stop(paste("Error fitting exogenous model:", e$message))
      })
    },

    calculate_diagnostics = function() {
      self$data_diagnostics <- list(
        n_observations = nrow(self$exogendata),
        n_groups = length(unique(self$exogendata[[data.table::key(self$exogendata)]])),
        time_range = base::range(self$exogendata[[data.table::indices(self$exogendata)]]),
        missing_values = colSums(is.na(self$exogendata)),
        timestamp = Sys.time()
      )
    },

    validate_model = function() {
      if (is.null(self$model)) {
        stop("Model has not been fitted. Call fit() before predict()")
      }
    },

    summarize_model = function() {
      list(
        formula = base::deparse(self$formula),
        model_dimensions = dim(self$model),
        group_key = data.table::key(self$model),
        time_index = data.table::indices(self$model)
      )
    }
  )
)


#' Linear Regression Model Class
#'
#' @description
#' An R6 Class implementing a linear regression model with optional
#' bootstrap capabilities for uncertainty estimation.
#'
#' @details
#' The LinearRegressionModel fits a linear regression to the data and can
#' use bootstrap resampling to generate predictions with uncertainty estimates.
#' It supports both residual and wild bootstrap methods.
#'
#' @export
LinearRegressionModel <- R6::R6Class(
  "LinearRegressionModel",
  inherit = Model,
  public = list(
    #' @field model The fitted linear model object
    model = NULL,

    #' @field new_formula Processed formula for fitting
    new_formula = NULL,

    #' @field coef_table Table of model coefficients and statistics
    coef_table = NULL,

    #' @field model_diagnostics Table of model fit diagnostics
    model_diagnostics = NULL,

    #' @field bootstrap_diagnostics List of bootstrap resampling diagnostics
    bootstrap_diagnostics = NULL,

    #' @description
    #' Create a new linear regression model instance
    #' @param formula Linear regression formula
    #' @param boot Whether to use bootstrap resampling
    #' @return A new `LinearRegressionModel` object
    #' @examples
    #' \dontrun{
    #' # Simple linear regression
    #' model <- LinearRegressionModel$new(
    #'   formula = y ~ x1 + x2,
    #'   boot = FALSE
    #' )
    #'
    #' # With bootstrap
    #' boot_model <- LinearRegressionModel$new(
    #'   formula = y ~ x1 + x2,
    #'   boot = TRUE
    #' )
    #' }
    initialize = function(formula, boot = FALSE) {
      super$initialize(formula = formula, boot = boot)
      private$initialize_diagnostics()
    },

    #' @description
    #' Fits the linear regression model to the data
    #' @param data Training data
    #' @param time_parameters Required if boot=TRUE, specifies bootstrap parameters
    fit = function(data, time_parameters = NULL) {
      self$model <- NULL
      private$validate_fit_inputs(data, time_parameters)
      private$prepare_model_formula(data)

      idx <- data.table::indices(data)

      if (!self$boot) {
        private$fit_standard_model(data)
      } else {
        private$fit_bootstrap_model(data, time_parameters)
      }
    },

    #' @description
    #' Generates predictions from the fitted model
    #' @param data Data to generate predictions for
    #' @param t Time point for prediction
    #' @return Data.table with predictions and uncertainty estimates
    predict = function(data, t) {
      private$validate_prediction_inputs(data, t)
      prediction_result <- private$generate_predictions(data, t)
      return(prediction_result)
    },

    #' @description
    #' Returns model diagnostics and performance metrics
    #' @return List containing model summary, diagnostics, and bootstrap info
    get_model_diagnostics = function() {
      list(
        model_summary = private$summarize_model(),
        diagnostics = self$model_diagnostics,
        bootstrap_info = self$bootstrap_diagnostics
      )
    }
  ),

  private = list(
    # Initialization
    initialize_diagnostics = function() {
      self$model_diagnostics <- data.table::data.table()
      self$bootstrap_diagnostics <- list()
    },

    # Validation methods
    validate_fit_inputs = function(data, time_parameters) {
      if (self$boot && is.null(time_parameters)) {
        stop("Time parameters required for bootstrap models")
      }
      if (!data.table::is.data.table(data)) {
        stop("Input data must be a data.table")
      }
    },

    validate_prediction_inputs = function(data, t) {
      if (is.null(self$model)) {
        stop("Model must be fitted before prediction")
      }
      if (is.null(t)) {
        stop("Time parameter 't' must be provided")
      }
    },

    # Model preparation
    prepare_model_formula = function(data) {
      grp <- data.table::key(data)
      idx <- data.table::indices(data)
      transformed_data <- self$panel_modelframe(self$formula, data) |>
        janitor::clean_names()

      varnames <- names(transformed_data)
      varnames <- varnames[!varnames %in% c(grp, idx)]
      self$new_formula <- paste(varnames[1], "~",
                                paste(varnames[2:length(varnames)], collapse = "+")) |>
        stats::formula()
    },

    # Model fitting methods
    fit_standard_model = function(data) {
      transformed_data <- private$transform_data(data)
      self$model <- private$fit_lm_with_diagnostics(self$new_formula, transformed_data)
    },

    fit_bootstrap_model = function(data, time_parameters) {
      transformed_data <- private$transform_data(data)
      self$model <- private$perform_bootstrap(
        self$new_formula, transformed_data, time_parameters
      )
    },

    fit_lm_with_diagnostics = function(formula, data) {
      idx <- data.table::indices(data)
      model <- stats::lm(formula, data = data)

      self$model_diagnostics <- data.table::rbindlist(list(self$model_diagnostics, broom::glance(model)))

      model_coefs <- data.table::data.table(
        broom::tidy(model),
        start = data[,min(time), env = list(time = idx)],
        end = data[,max(time), env = list(time = idx)]
      )
      self$coef_table <- data.table::rbindlist(list(self$coef_table, model_coefs))

      return(model)
    },

    # Bootstrap methods
    perform_bootstrap = function(formula, data, time_parameters) {
      grp <- data.table::key(data)
      idx <- data.table::indices(data)
      train <- data[time >= time_parameters$start &
                      time <= time_parameters$end,
                    env = list(time = idx)]

      data.table::setkeyv(train, grp)
      data.table::setindexv(train, idx)

      initial_model <- private$fit_lm_with_diagnostics(formula, train)
      bootstrap_data <- private$prepare_bootstrap_data(initial_model, train)
      bootstrap_result <- private$run_bootstrap(formula, bootstrap_data)

      return(bootstrap_result)
    },

    prepare_bootstrap_data = function(model, train) {
      model_dt <- data.table::data.table(
        row_names = names(model$residuals),
        residuals = model$residuals,
        fitted = model$fitted.values
      )

      train[, rownames := rownames(train)]
      train[model_dt, residuals := residuals, on = .(rownames = row_names)]
      train[model_dt, fitted := fitted, on = .(rownames = row_names)]

      return(train)
    },

    run_bootstrap = function(formula, data, type = "resid") {
      seed <- base::sample.int(.Machine$integer.max, 1)
      set.seed(seed)

      if (type == "resid") {
        resampled_residuals <- base::sample(data$residuals, size = nrow(data), replace = TRUE)
        y_star <- data$fitted + resampled_residuals
      } else if (type == "wild") {
        wild_residuals <- data$residuals * stats::rnorm(nrow(data))
        y_star <- data$fitted + wild_residuals
      } else {
        stop("Unknown bootstrap type")
      }

      # Store bootstrap diagnostics
      self$bootstrap_diagnostics[[length(self$bootstrap_diagnostics) + 1]] <- list(
        type = type,
        original_residuals = summary(data$residuals),
        resampled_residuals = if(type == "resid") summary(resampled_residuals) else summary(wild_residuals),
        timestamp = Sys.time()
      )

      response_var <- rlang::f_lhs(formula) |> as.character()
      new_data <- data.table::copy(data)
      data.table::set(new_data, j = response_var, value = y_star)

      return(private$fit_lm_with_diagnostics(formula, new_data))
    },

    # Prediction methods
    generate_predictions = function(data, t) {
      transformed_data <- private$transform_data(data)
      predictions <- stats::predict(self$model, newdata = transformed_data, se.fit = TRUE)

      return(private$process_predictions(predictions, transformed_data, t))
    },

    process_predictions = function(predictions, data, t) {
      se.PI <- sqrt(predictions$se.fit^2 + predictions$residual.scale^2)
      fit_uncertainty <- private$generate_prediction_intervals(predictions, se.PI)

      return(private$format_predictions(fit_uncertainty, data, t))
    },

    generate_prediction_intervals = function(predictions, se.PI) {
      tryCatch({
        predictions$fit + outer(se.PI, stats::rt(100, predictions$df))
      }, error = function(e) {
        warning("Error in generating prediction intervals, using normal approximation")
        predictions$fit + outer(se.PI, stats::rnorm(100))
      })
    },

    format_predictions = function(fit_uncertainty, data, t) {
      seed <- base::sample.int(.Machine$integer.max, 1)
      set.seed(seed)

      selector <- data.table::data.table(
        ridx = 1:dim(fit_uncertainty)[1],
        cidx = base::sample.int(100, dim(fit_uncertainty)[1], replace = TRUE)
      ) |> as.matrix()

      idx <- data.table::indices(data)
      data[, time := time + 1, env = list(time = idx)]

      outcome <- rlang::f_lhs(self$new_formula) |> as.character()
      data[, (outcome) := fit_uncertainty[selector]]

      grp <- data.table::key(data)
      keep <- c(grp, idx, outcome)
      result <- data[time == t, ..keep, env = list(time = idx)]
      data.table::setkeyv(result, grp)

      return(result)
    },

    # Helper methods
    transform_data = function(data) {
      grp <- data.table::key(data)
      idx <- data.table::indices(data)
      res <- self$panel_modelframe(self$formula, data) |> janitor::clean_names()
      data.table::setkeyv(res, grp)
      data.table::setindexv(res, idx)
      return(res)
    },

    summarize_model = function() {
      list(
        formula = base::deparse(self$new_formula),
        coefficients = broom::tidy(model),
        bootstrap_enabled = self$boot,
        latest_diagnostics = data.table::last(self$model_diagnostics)
      )
    }
  )
)
