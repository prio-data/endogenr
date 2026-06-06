# Experiment grid for endogenr system simulations -------------------------
#
# Run variations of the setup_system -> fit_system -> simulate_system pipeline
# (varying the model system, train_start, and/or test_start) and stack the
# results into one long data.table that plugs straight into the scoring
# helpers. Parallelism defaults to the experiment level (see run_experiments()).

# Internal: the outcome variable name a spec produces. Mirrors the rule used in
# setup_system(): independent families (exogen / parametric_distribution /
# univariate_fable) name their outcome on the (sole) formula, dependent
# families on the formula LHS.
.spec_outcome <- function(spec) {
  f <- spec$formula
  if (spec$type %in% get_independent_models()) {
    all.vars(f)[1L]
  } else {
    all.vars(rlang::f_lhs(f))[1L]
  }
}

# Internal: auto-name blank/missing entries `prefix1`, `prefix2`, ... and make
# the whole vector unique. Shared by vary_model() and .normalize_models(); same
# convention as longhorizon's .name_formulas().
.auto_name <- function(x, prefix) {
  nm <- names(x)
  if (is.null(nm)) nm <- rep("", length(x))
  blank <- is.na(nm) | !nzchar(nm)
  nm[blank] <- paste0(prefix, seq_along(x))[blank]
  names(x) <- make.unique(nm)
  x
}

# Internal: coerce `models` into a named list of full model systems. A single
# system (every element an endogenr_spec) is wrapped as `model1`; a named list
# of systems is auto-named and returned as-is.
.normalize_models <- function(models) {
  if (!is.list(models) || length(models) == 0L) {
    stop("`models` must be a model system or a named list of model systems.",
         call. = FALSE)
  }
  if (all(vapply(models, inherits, logical(1), "endogenr_spec"))) {
    return(list(model1 = models))
  }
  is_system <- vapply(models, function(m) {
    is.list(m) && length(m) > 0L &&
      all(vapply(m, inherits, logical(1), "endogenr_spec"))
  }, logical(1))
  if (!all(is_system)) {
    stop("`models` must be either a single model system (a list of ",
         "build_model() specs) or a named list of such systems.", call. = FALSE)
  }
  .auto_name(models, "model")
}

#' Build a window configuration for `run_experiments()`
#'
#' Creates a lightweight configuration object bundling the fit-side window
#' approach (passed to [fit_system()]) and the sim-side window policy (passed to
#' [simulate_system()]). Pass a single config to apply it uniformly across all
#' experiments, or a **named list** of configs to [run_experiments()] to cross
#' window approach as a fourth grid axis alongside `models`, `train_start`, and
#' `test_start`.
#'
#' @param window One of `"random"` (default, historical behaviour),
#'   `"rolling"` (fixed-width sliding window), or `"expanding"` (expanding
#'   window from `train_start`). Passed to [fit_system()].
#' @param width Integer or `NULL`. Rolling-window length; only used when
#'   `window = "rolling"`. Defaults to `min_window` (or 10) inside
#'   [fit_system()].
#' @param step Positive integer. Spacing between window-end anchors for sliding
#'   modes. Default `1L`. Ignored for `"random"`.
#' @param window_policy One of `"latest"` (default), `"equal"`, or `"decay"`.
#'   Controls which training window drives each forecast step in sliding modes;
#'   ignored for `"random"`. Passed to [simulate_system()].
#' @param decay Numeric in `(0, 1]`. Per-step recency retention for
#'   `window_policy = "decay"` at horizon 1. Default `0.5`.
#' @param weights Optional custom window weights overriding `window_policy`: a
#'   `function(h, age)` or a `horizon x n_window` matrix. Ignored for
#'   `"random"`.
#'
#' @return An `endogenr_window` object (a list with six fields).
#' @seealso [run_experiments()], [fit_system()], [simulate_system()]
#' @family experiments
#' @export
#'
#' @examples
#' # Default: random windows, latest policy (historical behaviour)
#' window_config()
#'
#' # Rolling window with decay policy
#' window_config("rolling", width = 15, step = 2,
#'               window_policy = "decay", decay = 0.7)
#'
#' # Named list as a grid axis in run_experiments()
#' \dontrun{
#' run_experiments(
#'   data = df, models = my_system,
#'   train_start = 1970, test_start = 2010, horizon = 10,
#'   groupvar = "gwcode", timevar = "year",
#'   inner_sims = 10, nsim = 20, min_window = 20,
#'   windows = list(
#'     rand = window_config(),
#'     roll = window_config("rolling", width = 20, window_policy = "latest")
#'   )
#' )
#' }
window_config <- function(window = c("random", "rolling", "expanding"),
                          width = NULL, step = 1L,
                          window_policy = c("latest", "equal", "decay"),
                          decay = 0.5, weights = NULL) {
  window        <- match.arg(window)
  window_policy <- match.arg(window_policy)
  step <- as.integer(step)
  if (is.na(step) || step < 1L) {
    stop("`step` must be a positive integer.", call. = FALSE)
  }
  if (!is.null(width)) {
    width <- as.integer(width)
    if (is.na(width) || width < 1L) {
      stop("`width` must be a positive integer.", call. = FALSE)
    }
  }
  if (!is.numeric(decay) || length(decay) != 1L || is.na(decay) ||
      decay <= 0 || decay > 1) {
    stop("`decay` must be a single number in (0, 1].", call. = FALSE)
  }
  if (!is.null(weights) && !is.function(weights) && !is.matrix(weights)) {
    stop("`weights` must be NULL, a function(h, age), or a numeric matrix.",
         call. = FALSE)
  }
  if (window == "random" && (window_policy != "latest" || !is.null(weights))) {
    warning("`window_policy` and `weights` are ignored when `window = \"random\"`.",
            call. = FALSE)
  }
  structure(
    list(window = window, width = width, step = step,
         window_policy = window_policy, decay = decay, weights = weights),
    class = "endogenr_window"
  )
}

# Internal: coerce `windows` into a named list of endogenr_window configs.
# NULL becomes a single default config; a single config is wrapped; a named
# list of configs is auto-named and returned as-is.
.normalize_windows <- function(windows) {
  if (is.null(windows)) {
    return(.auto_name(list(window_config()), "window"))
  }
  if (inherits(windows, "endogenr_window")) {
    return(.auto_name(list(windows), "window"))
  }
  if (!is.list(windows) || length(windows) == 0L) {
    stop("`windows` must be NULL, a window_config(), or a named list of ",
         "window_config() objects.", call. = FALSE)
  }
  is_wc <- vapply(windows, inherits, logical(1), "endogenr_window")
  if (!all(is_wc)) {
    stop("`windows` must be NULL, a window_config(), or a named list of ",
         "window_config() objects.", call. = FALSE)
  }
  .auto_name(windows, "window")
}

# Internal: build a compact, readable, unique label per experiment from the
# axes that actually vary. Falls back to the model name when nothing varies.
.experiment_labels <- function(grid, vary_model, vary_train, vary_test,
                               vary_window = FALSE) {
  parts <- vapply(seq_len(nrow(grid)), function(i) {
    comp <- character(0L)
    if (vary_model)  comp <- c(comp, grid$model[i])
    if (vary_train)  comp <- c(comp, paste0("train", grid$train_start[i]))
    if (vary_test)   comp <- c(comp, paste0("test", grid$test_start[i]))
    if (vary_window) comp <- c(comp, grid$window_cfg[i])
    if (length(comp) == 0L) comp <- grid$model[i]
    paste(comp, collapse = "_")
  }, character(1))
  make.unique(parts, sep = "_")
}


#' Build model-system variants by swapping specs into a base system
#'
#' Expands a single base model system into a named list of full systems, one per
#' variant, by replacing the base spec(s) whose outcome matches each override.
#' This is the ergonomic way to explore formula variations (e.g. several
#' candidate `gdppc_grwt` equations `e1`, `e2`, ...) without retyping the
#' unchanged specs each time. The result is exactly the named-list-of-systems
#' shape accepted by [run_experiments()].
#'
#' Matching is by outcome variable (see [build_model()] for how each model type
#' names its outcome): an override is spliced into the position of the base spec
#' producing the same outcome, preserving base order and leaving every other
#' spec untouched. An override whose outcome is **not** produced by `base` is an
#' error — `vary_model()` swaps specs, it never appends new ones (build the full
#' system by hand for that).
#'
#' @param base A model system: a list of specs from [build_model()].
#' @param variants A named list of overrides. Each entry is either a single
#'   [build_model()] spec or a list of specs; every spec replaces the base spec
#'   with the same outcome. Unnamed or partially named entries are auto-named
#'   `variant1`, `variant2`, ... and that name becomes the variant label.
#'
#' @return A named list of full model systems, suitable for the `models`
#'   argument of [run_experiments()].
#' @seealso [run_experiments()], [build_model()]
#' @family experiments
#' @export
#'
#' @examples
#' base <- list(
#'   build_model("deterministic",
#'               formula = gdppc ~ I(abs(lag(gdppc) * (1 + gdppc_grwt)))),
#'   build_model("linear", formula = gdppc_grwt ~ lag(log(gdppc)), boot = "resid"),
#'   build_model("exogen", formula = ~population)
#' )
#' systems <- vary_model(base, list(
#'   e1 = build_model("linear", formula = gdppc_grwt ~ lag(log(gdppc)), boot = "resid"),
#'   e2 = build_model("linear", formula = gdppc_grwt ~ lag(log(gdppc)) + lag(population),
#'                    boot = "resid")
#' ))
#' names(systems)
vary_model <- function(base, variants) {
  if (!is.list(base) || length(base) == 0L ||
      !all(vapply(base, inherits, logical(1), "endogenr_spec"))) {
    stop("`base` must be a model system: a non-empty list of build_model() specs.",
         call. = FALSE)
  }
  if (!is.list(variants) || length(variants) == 0L) {
    stop("`variants` must be a non-empty list of override specs.", call. = FALSE)
  }
  variants <- .auto_name(variants, "variant")
  base_outcomes <- vapply(base, .spec_outcome, character(1))

  lapply(stats::setNames(names(variants), names(variants)), function(v) {
    ov <- variants[[v]]
    overrides <- if (inherits(ov, "endogenr_spec")) list(ov) else ov
    if (!is.list(overrides) ||
        !all(vapply(overrides, inherits, logical(1), "endogenr_spec"))) {
      stop("Variant '", v, "' must be a build_model() spec or a list of specs.",
           call. = FALSE)
    }
    sys <- base
    for (spec in overrides) {
      oc  <- .spec_outcome(spec)
      pos <- which(base_outcomes == oc)
      if (length(pos) == 0L) {
        stop("Variant '", v, "': override outcome '", oc, "' is not produced ",
             "by any spec in `base`; vary_model() swaps specs, it does not ",
             "append new ones.", call. = FALSE)
      }
      sys[[pos[1L]]] <- spec
    }
    sys
  })
}


#' Run a grid of system-simulation experiments
#'
#' Runs the [setup_system()] -> [fit_system()] -> [simulate_system()] pipeline
#' once per experiment across a grid of variations and stacks the results into a
#' single long `data.table`. Vary the model system, the training start, and/or
#' the forecast start; every combination is a separate experiment with a
#' readable label. The output carries the same `panel_unit` / `panel_time`
#' attributes as [simulate_system()], so it plugs directly into
#' [get_experiment_accuracy()] and [plotsim()].
#'
#' @details
#' The experiment grid is the Cartesian product of up to four axes: the model
#' systems (`models`), `train_start`, `test_start`, and optionally `windows`
#' (see [window_config()]). An axis with a single value does not expand.
#' `horizon`, `inner_sims`, `min_window`, and `nsim` are shared across every
#' experiment.
#'
#' Each experiment is stamped with five leading columns: `.experiment` (a unique
#' label built only from the axes that vary), `model`, `train_start`,
#' `test_start`, and `window_cfg` (the name of the window configuration used).
#' The result also carries an `experiments` attribute (a metadata `data.table`)
#' recording the configuration of every experiment, which
#' [get_experiment_accuracy()] uses to recover per-experiment metadata.
#'
#' @section Parallelism:
#' Because `future` evaluates *nested* futures sequentially unless a multi-level
#' plan is set, you set **one** [future::plan()] and choose where it bites with
#' `parallel`:
#' \itemize{
#'   \item `"experiments"` (default): the experiment grid is mapped in parallel
#'     and each experiment's internal `fit`/`simulate` runs sequentially on its
#'     worker. Best when there are many experiments. Each worker builds its own
#'     simulation grid, so memory scales with the number of workers — keep that
#'     in mind for large panels.
#'   \item `"draws"`: experiments run sequentially and each [fit_system()] /
#'     [simulate_system()] parallelises over its `nsim` draws (the single-
#'     experiment behaviour). Best when there are few experiments but a large
#'     `nsim`.
#' }
#' Set the plan yourself, e.g. `future::plan(future::multisession, workers = 8)`.
#' Reproducibility uses `future.seed = TRUE`; seed once with [set.seed()] before
#' calling. The two schemes lay out their random streams differently, so they
#' are reproducible *within* a scheme but not bit-identical *across* schemes (the
#' same caveat as [fit_system()] versus the legacy single-call simulator).
#'
#' @param data A data.frame, data.table, or tsibble; see [setup_system()].
#' @param models Either a single model system (a list of [build_model()] specs)
#'   or a named list of such systems. [vary_model()] is the ergonomic way to
#'   build the named list by swapping individual specs.
#' @param train_start Integer scalar or vector of training starts to cross.
#' @param test_start Integer scalar or vector of forecast starts to cross.
#' @param horizon Integer. Forecast horizon, shared across experiments.
#' @param groupvar,timevar Character. Panel unit and time column names.
#' @param inner_sims Integer. Inner simulations per draw; see [setup_system()].
#' @param nsim Integer. Number of coefficient draws per experiment; see
#'   [fit_system()].
#' @param min_window Integer or `NULL`. Random-window refitting; see
#'   [setup_system()].
#' @param windows Either `NULL` (default; random window, latest policy),
#'   a single [window_config()] applied to every experiment, or a named list
#'   of [window_config()] objects that is crossed as a fourth grid axis. The
#'   config name appears as the `window_cfg` column.
#' @param globals Named list of user functions referenced by formulas, exported
#'   to parallel workers; see [setup_system()].
#' @param parallel One of `"experiments"` (default) or `"draws"`; see the
#'   Parallelism section.
#' @param keep_fits Logical. When `TRUE`, the per-experiment
#'   `endogenr_fitted_system` objects are attached to the result as a `fits`
#'   attribute (named by `.experiment`), so [get_coefficients()] /
#'   [plot_coefficients()] can be run per experiment. Default `FALSE` keeps the
#'   result lean.
#'
#' @return A `data.table` stacking every experiment's [simulate_system()]
#'   output, with leading columns `.experiment`, `model`, `train_start`,
#'   `test_start`, `window_cfg` and the usual simulation columns plus `.sim`.
#'   Carries `panel_unit`, `panel_time`, and `experiments` attributes (and
#'   `fits` when `keep_fits = TRUE`).
#' @seealso [vary_model()], [window_config()], [setup_system()], [fit_system()],
#'   [simulate_system()], [get_experiment_accuracy()]
#' @family experiments
#' @export
#'
#' @examples
#' \dontrun{
#' df <- endogenr::example_data
#' base <- list(
#'   build_model("deterministic",
#'               formula = gdppc ~ I(abs(lag(gdppc) * (1 + gdppc_grwt)))),
#'   build_model("linear", formula = gdppc_grwt ~ lag(log(gdppc)), boot = "resid"),
#'   build_model("exogen", formula = ~population)
#' )
#' systems <- vary_model(base, list(
#'   e1 = build_model("linear", formula = gdppc_grwt ~ lag(log(gdppc)), boot = "resid"),
#'   e2 = build_model("linear", formula = gdppc_grwt ~ lag(log(gdppc)) + lag(population),
#'                    boot = "resid")
#' ))
#'
#' future::plan(future::multisession, workers = 8)
#' set.seed(42)
#' res <- run_experiments(
#'   data = df, models = systems,
#'   train_start = 1970, test_start = c(2005, 2010), horizon = 12,
#'   groupvar = "gwcode", timevar = "year",
#'   inner_sims = 30, nsim = 64, min_window = 20
#' )
#' future::plan(future::sequential)
#'
#' acc <- get_experiment_accuracy(res, "gdppc", df)
#' }
run_experiments <- function(data, models, train_start, test_start, horizon,
                            groupvar, timevar, inner_sims, nsim = 1L,
                            min_window = NULL, windows = NULL, globals = NULL,
                            parallel = c("experiments", "draws"),
                            keep_fits = FALSE) {
  parallel <- match.arg(parallel)
  systems  <- .normalize_models(models)
  windows  <- .normalize_windows(windows)

  train_start <- as.integer(train_start)
  test_start  <- as.integer(test_start)
  if (anyNA(train_start) || length(train_start) == 0L ||
      anyNA(test_start) || length(test_start) == 0L) {
    stop("`train_start` and `test_start` must be non-empty integer-valued.",
         call. = FALSE)
  }
  nsim <- as.integer(nsim)

  # Cartesian grid over (model, train_start, test_start, window_cfg).
  grid <- expand.grid(
    .mi         = seq_along(systems),
    train_start = unique(train_start),
    test_start  = unique(test_start),
    .wi         = seq_along(windows),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  grid$model      <- names(systems)[grid$.mi]
  grid$window_cfg <- names(windows)[grid$.wi]
  grid$.experiment <- .experiment_labels(
    grid,
    vary_model  = length(systems) > 1L,
    vary_train  = length(unique(train_start)) > 1L,
    vary_test   = length(unique(test_start)) > 1L,
    vary_window = length(windows) > 1L
  )
  n_exp <- nrow(grid)

  # Per-experiment vectors, referenced by the worker below. Keeping them as
  # plain vectors (not the whole grid) makes the exported global set explicit.
  exp_idx        <- grid$.mi
  exp_train      <- grid$train_start
  exp_test       <- grid$test_start
  exp_label      <- grid$.experiment
  exp_wi         <- grid$.wi
  exp_window_cfg <- grid$window_cfg

  p <- if (requireNamespace("progressr", quietly = TRUE)) {
    progressr::progressor(steps = n_exp)
  } else {
    function(...) invisible(NULL)
  }

  # One experiment: build its own setup (hence its own simulation grid), fit,
  # and simulate. Errors are re-raised with the experiment label for context.
  run_one <- function(i) {
    out <- tryCatch({
      wc    <- windows[[exp_wi[i]]]
      setup <- setup_system(
        models = systems[[exp_idx[i]]], data = data,
        train_start = exp_train[i], test_start = exp_test[i], horizon = horizon,
        groupvar = groupvar, timevar = timevar,
        inner_sims = inner_sims, min_window = min_window, globals = globals
      )
      fitted <- fit_system(setup, nsim = nsim,
                           window = wc$window, width = wc$width, step = wc$step)
      sim <- if (wc$window == "random") {
        simulate_system(fitted)
      } else {
        simulate_system(fitted, window_policy = wc$window_policy,
                        decay = wc$decay, weights = wc$weights)
      }
      list(sim = sim, fit = if (isTRUE(keep_fits)) fitted else NULL)
    }, error = function(e) {
      stop("Experiment '", exp_label[i], "' failed: ", conditionMessage(e),
           call. = FALSE)
    })
    p()
    out
  }

  if (parallel == "experiments") {
    # Named-list globals are authoritative (no auto-detection), so enumerate
    # every name run_one references; spread user globals so formula NSE resolves
    # them on the workers. Inner fit/simulate futures run sequentially (nested).
    future_globals <- list(
      systems = systems, exp_idx = exp_idx, exp_train = exp_train,
      exp_test = exp_test, exp_label = exp_label, data = data,
      horizon = horizon, groupvar = groupvar, timevar = timevar,
      inner_sims = inner_sims, nsim = nsim, min_window = min_window,
      globals = globals, keep_fits = keep_fits, p = p,
      windows = windows, exp_wi = exp_wi
    )
    if (!is.null(globals)) {
      for (fn_name in names(globals)) future_globals[[fn_name]] <- globals[[fn_name]]
    }
    results <- future.apply::future_lapply(
      seq_len(n_exp), run_one,
      future.seed     = TRUE,
      future.globals  = future_globals,
      future.packages = c("endogenr", "data.table")
    )
  } else {
    results <- lapply(seq_len(n_exp), run_one)
  }

  # Stamp metadata onto each result and stack.
  sims <- vector("list", n_exp)
  for (i in seq_len(n_exp)) {
    s <- results[[i]]$sim
    s[, `:=`(.experiment = exp_label[i], model = grid$model[i],
             train_start = exp_train[i], test_start = exp_test[i],
             window_cfg = exp_window_cfg[i])]
    sims[[i]] <- s
  }
  out  <- data.table::rbindlist(sims, use.names = TRUE, fill = TRUE)
  lead <- c(".experiment", "model", "train_start", "test_start", "window_cfg")
  data.table::setcolorder(out, c(lead, setdiff(names(out), lead)))

  exp_meta <- data.table::data.table(
    .experiment   = exp_label, model = grid$model,
    train_start   = exp_train, test_start = exp_test,
    window_cfg    = exp_window_cfg,
    window        = vapply(exp_wi, function(i) windows[[i]]$window, character(1)),
    window_policy = vapply(exp_wi, function(i) windows[[i]]$window_policy, character(1)),
    horizon       = as.integer(horizon), inner_sims = as.integer(inner_sims),
    nsim          = nsim
  )
  data.table::setattr(out, "panel_unit", groupvar)
  data.table::setattr(out, "panel_time", timevar)
  data.table::setattr(out, "experiments", exp_meta)

  if (isTRUE(keep_fits)) {
    fits <- lapply(results, function(r) r$fit)
    names(fits) <- exp_label
    data.table::setattr(out, "fits", fits)
  }

  out[]
}


#' Score a grid of system-simulation experiments
#'
#' Computes CRPS, MAE, and Winkler scores for every experiment in a
#' [run_experiments()] result, returning one long `data.table` keyed by
#' `.experiment`. Each experiment is scored against `truth` with its **own**
#' `test_start`, so a grid that varies `test_start` is handled correctly (the
#' `horizon` column is `time - test_start + 1` per experiment). This is the
#' experiment-aware analogue of [get_accuracy()], which it calls under the hood.
#'
#' Rows before each experiment's `test_start` (the training period) are dropped
#' before scoring, so the reported horizons start at 1.
#'
#' @param experiment_results A `data.table` from [run_experiments()].
#' @param outcome Character. Name of the outcome variable to score.
#' @param truth A data.frame/data.table of observed values; must contain the
#'   group, time, and outcome columns.
#' @param by Character vector of grouping columns. Defaults to
#'   `c(<unit>, "horizon")`, matching [get_lh_accuracy()] / [compare_approaches()].
#' @param level Numeric. Coverage level for the Winkler score (default 50).
#' @param transform Function or `NULL`. Applied to both draws and truth before
#'   scoring; see [get_accuracy()].
#'
#' @return A `data.table` with leading columns `.experiment` and any present
#'   axis columns (`model`, `train_start`, `test_start`, `window_cfg`), the
#'   `by` columns, and `crps`, `mae`, `winkler` — one block of rows per
#'   experiment.
#' @seealso [run_experiments()], [get_accuracy()], [compare_approaches()]
#' @family experiments
#' @export
get_experiment_accuracy <- function(experiment_results, outcome, truth,
                                    by = NULL, level = 50, transform = NULL) {
  if (!data.table::is.data.table(experiment_results)) {
    experiment_results <- data.table::as.data.table(experiment_results)
  }
  if (!".experiment" %in% names(experiment_results)) {
    stop("`experiment_results` must carry an `.experiment` column ",
         "(produced by run_experiments()).", call. = FALSE)
  }

  ctx <- .infer_ctx(experiment_results, truth)
  if (is.null(ctx)) {
    stop("Could not infer panel context from `experiment_results` or `truth`. ",
         "Pass a result produced by run_experiments().", call. = FALSE)
  }
  unit_col <- ctx_unit(ctx)
  time_col <- ctx_time(ctx)
  if (is.null(by)) by <- c(unit_col, "horizon")

  # Per-experiment test_start: prefer the stamped metadata, fall back to the
  # column carried on the result.
  meta <- attr(experiment_results, "experiments", exact = TRUE)
  if (!is.null(meta)) {
    ts_by_exp <- stats::setNames(meta$test_start, meta$.experiment)
  } else if ("test_start" %in% names(experiment_results)) {
    u <- unique(experiment_results[, c(".experiment", "test_start"), with = FALSE])
    ts_by_exp <- stats::setNames(u$test_start, u$.experiment)
  } else {
    stop("Cannot determine per-experiment `test_start`: the result has no ",
         "`experiments` attribute and no `test_start` column.", call. = FALSE)
  }

  axis_cols <- intersect(c("model", "train_start", "test_start", "window_cfg"),
                         names(experiment_results))
  labels <- unique(as.character(experiment_results$.experiment))

  scored <- lapply(labels, function(lab) {
    sub <- experiment_results[experiment_results$.experiment == lab]
    ts  <- ts_by_exp[[lab]]
    if (is.null(ts) || is.na(ts)) {
      stop("No test_start recorded for experiment '", lab, "'.", call. = FALSE)
    }
    sub <- sub[sub[[time_col]] >= ts]  # forecast window only
    acc <- get_accuracy(sub, outcome = outcome, truth = truth, ctx = ctx,
                        level = level, test_start = ts, by = by,
                        transform = transform)
    acc[, .experiment := lab]
    for (col in axis_cols) {
      data.table::set(acc, j = col, value = sub[[col]][1L])
    }
    acc
  })

  out  <- data.table::rbindlist(scored, use.names = TRUE, fill = TRUE)
  lead <- c(".experiment", axis_cols)
  data.table::setcolorder(out, c(lead, setdiff(names(out), lead)))
  out[]
}
