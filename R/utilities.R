#' Draw a random training window
#'
#' Used by [fit_system()] to draw a random training window for `linear`,
#' `glm`, and `heterolm` models when `min_window` is set. The window is
#' anchored inside `[earliest_train_start, test_start - 1]`, with random
#' (possibly zero) padding on either side, and is guaranteed to be at
#' least `min_window` time steps long.
#'
#' @param earliest_train_start Integer. The earliest time step that may
#'   appear in the training window.
#' @param test_start Integer. The first time step in the forecast period;
#'   the training window always ends at `test_start - 1` or earlier.
#' @param min_window Integer or `NULL`. Minimum window length. When `NULL`
#'   the full range `[earliest_train_start, test_start - 1]` is returned
#'   unchanged.
#'
#' @return A list with `start`, `end`, and `window` (the latter `NULL`
#'   when `min_window` is `NULL`).
#' @family build
#' @export
#'
#' @examples
#' set.seed(1)
#' get_train_window(1970, 2010, min_window = 10)
get_train_window <- function(earliest_train_start, test_start, min_window = NULL){
  if(is.null(min_window)){
    return(list("start" = earliest_train_start, "end" = test_start -1, "window" = NULL))
  }
  full_range <- test_start - earliest_train_start
  if(min_window > full_range){
    stop("Train min_window must be smaller or equal to largest possible train-set.")
  }
  if(min_window < 1){
    stop("Train min_window must be 1 or larger")
  }

  start_increment <- sample.int(full_range-min_window, 1) - 1
  stop_decrement <- sample.int(full_range-min_window-start_increment, 1) - 1
  return(list("start" = earliest_train_start+start_increment, "end" = test_start - stop_decrement))
}

#' Inject a positional lag function into a formula's environment
#'
#' Replaces `lag()` in the formula's evaluation environment with a positional
#' shift function: `c(rep(NA, n), head(x, -n))`. This ensures `lag()` in model
#' formulas performs a within-group positional shift rather than `stats::lag()`.
#'
#' @param formula An R formula.
#' @return The formula with modified environment.
#' @family formula_helpers
#' @export
inject_positional_lag <- function(formula) {
  .positional_lag <- function(x, n = 1L) {
    c(rep(NA_real_, n), utils::head(x, -n))
  }
  formula_env <- new.env(parent = environment(formula))
  formula_env$lag <- .positional_lag
  environment(formula) <- formula_env
  formula
}

#' Per-group model.frame with optional key-column injection
#'
#' Evaluates `formula_eval` per group via `model.frame()`. When
#' `needs_keys = TRUE`, the scalar group-key values are recycled into the
#' per-group data frame before evaluation so that formula terms that
#' reference a key column (e.g. `factor(unit)`) can be evaluated.
#'
#' @param formula_eval A prepared formula (positional lag injected, index
#'   appended).
#' @param data A data.table.
#' @param all_keys Character vector of grouping columns (from [ctx_keys()]).
#' @param needs_keys Logical. When `TRUE`, recycle key values into the
#'   per-group frame.
#'
#' @return A data.table with `all_keys` + materialized columns.
#' @keywords internal
.model_frame_by_group <- function(formula_eval, data, all_keys, needs_keys) {
  if (needs_keys) {
    data[, {
      key_vals <- lapply(mget(all_keys), rep, .N)
      frame_data <- cbind(as.data.frame(.SD), key_vals)
      mf <- stats::model.frame(formula_eval, data = frame_data,
                               na.action = stats::na.pass)
      data.table::as.data.table(mf)
    }, by = all_keys]
  } else {
    data[, {
      mf <- stats::model.frame(formula_eval, data = .SD,
                               na.action = stats::na.pass)
      data.table::as.data.table(mf)
    }, by = all_keys]
  }
}

#' Rewrite formula term labels using a clean column-name mapping
#'
#' Walks the `factors` matrix of `terms(formula)` and replaces each raw
#' variable token (e.g. `"lag(log(gdppc))"`, `"factor(region)"`) with its
#' cleaned column name from `col_mapping`. Interaction structure (`:`) is
#' preserved: a term that was `g:x` in the raw formula becomes
#' `clean(g):clean(x)` in the output. Variables absent from `col_mapping`
#' pass through unchanged.
#'
#' @param formula An R formula.
#' @param col_mapping Named character vector mapping raw model.frame column
#'   names to cleaned names (as returned by [.build_mat_cache()]).
#'
#' @return Character vector of cleaned term labels (zero-length for
#'   intercept-only or empty formulas).
#' @keywords internal
.clean_term_labels <- function(formula, col_mapping) {
  tt  <- stats::terms(formula)
  fac <- attr(tt, "factors")

  # For intercept-only formulas, factors may be integer(0) without a dim
  # attribute (not a proper matrix); guard before calling ncol/apply.
  if (is.null(fac) || !is.matrix(fac) || ncol(fac) == 0L) return(character(0L))

  raw_vars   <- rownames(fac)
  clean_vars <- ifelse(raw_vars %in% names(col_mapping),
                       col_mapping[raw_vars],
                       raw_vars)
  names(clean_vars) <- raw_vars

  # Each column of `fac` is one term; contributing rows are those with entry != 0.
  # unname() strips the column-name attributes so callers get a plain character
  # vector (reformulate() does not care about names anyway).
  unname(apply(fac, 2L, function(col) {
    contributing <- raw_vars[col != 0L]
    paste(clean_vars[contributing], collapse = ":")
  }))
}

#' Materialize formula terms into a data.table
#'
#' Evaluates formula terms per group using `model.frame()` with positional lag
#' injection. Returns the materialized data.table with cleaned column names.
#'
#' For repeated calls with the same formula (e.g. inside a simulation loop),
#' pass `.mat_formula` and `.col_mapping` (cached from a prior call) to skip
#' the `update()` / `inject_positional_lag()` / `clean_names()` overhead.
#'
#' @param formula An R formula.
#' @param data A data.table (or data.frame/tsibble — will be coerced).
#' @param ctx A `panel_context` object.
#' @param .mat_formula Optional. A pre-prepared formula (already has index var
#'   appended and positional lag injected). When provided, `formula` is ignored
#'   for evaluation and used only as documentation.
#' @param .col_mapping Optional named character vector. Maps raw `model.frame()`
#'   column names to clean names. When provided, `janitor::clean_names()` is
#'   skipped and `setnames()` is used instead.
#'
#' @return A data.table with cleaned column names.
#' @family formula_helpers
#' @export
materialize_formula <- function(formula, data, ctx,
                                .mat_formula = NULL, .col_mapping = NULL) {
  all_keys <- ctx_keys(ctx)

  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(as.data.frame(data))
  }

  if (is.null(.mat_formula)) {
    idx <- ctx_time(ctx)
    .mat_formula <- stats::update(formula, paste(c(". ~ .", idx), collapse = "+"))
    .mat_formula <- inject_positional_lag(.mat_formula)
  }

  needs_keys <- any(all.vars(formula) %in% all_keys)

  eval_formula <- .mat_formula
  result <- .model_frame_by_group(eval_formula, data, all_keys, needs_keys)

  if (!is.null(.col_mapping)) {
    old <- intersect(names(.col_mapping), names(result))
    if (length(old) > 0L) {
      data.table::setnames(result, old, .col_mapping[old])
    }
  } else {
    result <- janitor::clean_names(result)
  }

  result
}

#' Derive a naive formula preserving term structure
#'
#' Reconstructs a formula suitable for `lm()` / `glm()` by mapping each raw
#' variable token in the original formula to its cleaned materialized column
#' name (via `col_mapping`), while preserving all term structure: interactions
#' (`:`), crossing (`*`), intercept suppression (`-1` / `0`), and higher-order
#' terms. This is the only correct way to derive the naive formula because
#' additive column-name reconstruction (the old approach) silently drops
#' interaction terms.
#'
#' @param formula The original R formula (before materialization).
#' @param col_mapping Named character vector. Maps raw `model.frame()` column
#'   names to cleaned names (as returned by [.build_mat_cache()]).
#' @param outcome Character. Cleaned outcome column name. When `NULL`
#'   (default), the LHS of `formula` is mapped through `col_mapping`.
#'
#' @return An R formula with cleaned variable names and full term structure.
#' @family formula_helpers
#' @export
derive_naive_formula <- function(formula, col_mapping, outcome = NULL) {
  tt <- stats::terms(formula)

  if (is.null(outcome)) {
    lhs_raw <- deparse(rlang::f_lhs(formula))
    outcome  <- if (lhs_raw %in% names(col_mapping)) col_mapping[[lhs_raw]] else lhs_raw
  }

  term_labels <- .clean_term_labels(formula, col_mapping)
  intercept   <- attr(tt, "intercept") == 1L

  if (length(term_labels) == 0L) {
    # Intercept-only formula (e.g. y ~ 1)
    return(stats::reformulate("1", response = outcome, intercept = TRUE))
  }

  stats::reformulate(term_labels, response = outcome, intercept = intercept)
}

#' Creates a panel data frame based on a formula and data.table
#'
#' Builds a materialization cache once via [.build_mat_cache()], then
#' materializes the formula and derives the naive formula from the original
#' term structure (preserving interactions, factor() calls, and intercept
#' suppression). Returns the materialized data together with the naive formula
#' and the cache components for reuse in `predict.*` methods.
#'
#' @param formula An R formula.
#' @param data A data.table (or data.frame/tsibble — will be coerced).
#' @param ctx A `panel_context` object.
#'
#' @return A list with `data` (data.table), `naive_formula`, `mat_formula`,
#'   and `col_mapping`.
#' @keywords internal
create_panel_frame <- function(formula, data, ctx) {
  cache        <- .build_mat_cache(formula, data, ctx)
  materialized <- materialize_formula(formula, data, ctx,
                                      .mat_formula = cache$mat_formula,
                                      .col_mapping = cache$col_mapping)
  naive_formula <- derive_naive_formula(formula, cache$col_mapping)
  list(data         = materialized,
       naive_formula = naive_formula,
       mat_formula   = cache$mat_formula,
       col_mapping   = cache$col_mapping)
}

#' Build a prepared formula and column-name mapping for fast materialization
#'
#' Pre-computes the formula (with index appended + positional lag injected) and
#' the raw → clean column name mapping. These can be passed to
#' [materialize_formula()] via `.mat_formula` and `.col_mapping` to skip
#' per-call overhead.
#'
#' @param formula An R formula.
#' @param data A data.table with at least a few rows (used to discover column
#'   names from `model.frame()`).
#' @param ctx A `panel_context` object.
#'
#' @return A list with `mat_formula` and `col_mapping`.
#' @keywords internal
.build_mat_cache <- function(formula, data, ctx) {
  idx      <- ctx_time(ctx)
  all_keys <- ctx_keys(ctx)

  mat_formula <- stats::update(formula, paste(c(". ~ .", idx), collapse = "+"))
  mat_formula <- inject_positional_lag(mat_formula)

  # Run model.frame on a small per-group sample to discover raw column names
  sample_dt  <- data[, utils::head(.SD, 2L), by = all_keys]
  needs_keys <- any(all.vars(formula) %in% all_keys)

  raw <- .model_frame_by_group(mat_formula, sample_dt, all_keys, needs_keys)

  raw_names   <- names(raw)
  cleaned     <- janitor::clean_names(raw)
  clean_names <- names(cleaned)

  # Only map names that actually changed
  changed     <- raw_names != clean_names
  col_mapping <- stats::setNames(clean_names[changed], raw_names[changed])

  list(mat_formula = mat_formula, col_mapping = col_mapping)
}

#' Pre-expand interaction/factor terms for heterolm via model.matrix
#' @keywords internal
#+
#+ `heterolm::hetero()` uses column-name lookups rather than R's formula
#+ algebra, so interaction terms (`a:b`) and factor dummy codes must be
#+ pre-computed as named columns. This helper runs `model.matrix` on the
#+ cleaned-label formula (which uses the already-materialized column names),
#+ adds any new columns to `mat_data` **in-place** (data.table reference
#+ semantics), and returns the flat column names for use in the heterolm
#+ formula string, plus the terms object and xlevels from the training
#+ model.frame so that predict time can reproduce the same basis (poly,
#+ factor contrasts, splines) coherently.
#+
#+ @param clean_labels Character vector of cleaned term labels (from
#+   `.clean_term_labels()`). May include interaction labels like `"a:b"`.
#+ @param mat_data A data.table, modified in-place by reference.
#+
#+ @return A list with:
#+   - `col_names`: flat cleaned column names (excluding intercept) to use in
#+     the heterolm formula (e.g. `c("g_b", "x", "g_b_x")`).
#+   - `terms_obj`: the `terms` object from the training `model.frame()`,
#+     carrying `predvars`/`xlevels` for coherent predict-time reconstruction.
#+   - `xlevels`: factor levels from the training data (for `xlev` at predict).
#+ @keywords internal
.hetero_expand_terms <- function(clean_labels, mat_data) {
  if (length(clean_labels) == 0L) {
    return(list(col_names = character(0L), terms_obj = NULL, xlevels = NULL))
  }

  rhs_form <- stats::reformulate(clean_labels)
  rhs_vars <- all.vars(rhs_form)

  # Identify complete rows using only the predictor variables
  avail_vars <- intersect(rhs_vars, names(mat_data))
  if (length(avail_vars) == 0L) {
    return(list(col_names = clean_labels, terms_obj = NULL, xlevels = NULL))
  }

  df_pred  <- as.data.frame(mat_data[, avail_vars, with = FALSE])
  complete <- which(stats::complete.cases(df_pred))
  if (length(complete) == 0L) {
    return(list(col_names = clean_labels, terms_obj = NULL, xlevels = NULL))
  }

  df_mm <- as.data.frame(mat_data[complete, ])

  # Build a model.frame to capture predvars (poly scaling, factor levels, etc.)
  # for coherent reconstruction at predict time.
  mf <- tryCatch(
    stats::model.frame(rhs_form, data = df_mm, na.action = stats::na.pass),
    error = function(e) NULL
  )
  if (is.null(mf)) {
    return(list(col_names = clean_labels, terms_obj = NULL, xlevels = NULL))
  }
  terms_obj <- attr(mf, "terms")
  xlevels   <- stats::.getXlevels(terms_obj, mf)

  mm <- tryCatch(
    stats::model.matrix(terms_obj, mf),
    error = function(e) NULL
  )
  if (is.null(mm)) {
    return(list(col_names = clean_labels, terms_obj = terms_obj, xlevels = xlevels))
  }

  raw_names   <- setdiff(colnames(mm), "(Intercept)")
  if (length(raw_names) == 0L) {
    return(list(col_names = character(0L), terms_obj = terms_obj, xlevels = xlevels))
  }
  clean_names <- janitor::make_clean_names(raw_names)

  # Add columns not yet in mat_data (NA-initialised, filled for complete rows)
  for (i in seq_along(raw_names)) {
    cn <- clean_names[[i]]
    if (!cn %in% names(mat_data)) {
      mat_data[, (cn) := NA_real_]
      data.table::set(mat_data, i = complete, j = cn,
                      value = as.numeric(mm[, raw_names[[i]]]))
    }
  }

  list(col_names = clean_names, terms_obj = terms_obj, xlevels = xlevels)
}

#+ Expand interaction/factor terms at predict time using the stored `terms_obj`
#+ and `xlevels` from fit time, ensuring poly/spline/factor bases are identical
#+ to those computed during training.
#+
#+ Adds any new flat columns to `mat_data` in-place (data.table reference
#+ semantics) and returns the clean column names (same contract as
#+ `.hetero_expand_terms()`).
#+ @keywords internal
.hetero_expand_from_terms <- function(terms_obj, xlevels, mat_data) {
  if (is.null(terms_obj)) return(invisible(NULL))
  rhs_vars <- all.vars(stats::reformulate(
    labels(terms_obj), intercept = FALSE
  ))
  avail_vars <- intersect(rhs_vars, names(mat_data))
  if (length(avail_vars) == 0L) return(invisible(NULL))

  df_pred  <- as.data.frame(mat_data[, avail_vars, with = FALSE])
  complete <- which(stats::complete.cases(df_pred))
  if (length(complete) == 0L) return(invisible(NULL))

  df_mm <- as.data.frame(mat_data[complete, ])
  mf <- tryCatch(
    stats::model.frame(terms_obj, data = df_mm,
                       xlev = xlevels, na.action = stats::na.pass),
    error = function(e) NULL
  )
  if (is.null(mf)) return(invisible(NULL))

  mm <- tryCatch(
    stats::model.matrix(terms_obj, mf, xlev = xlevels),
    error = function(e) NULL
  )
  if (is.null(mm)) return(invisible(NULL))

  raw_names   <- setdiff(colnames(mm), "(Intercept)")
  clean_names <- janitor::make_clean_names(raw_names)

  for (i in seq_along(raw_names)) {
    cn <- clean_names[[i]]
    if (!cn %in% names(mat_data)) {
      mat_data[, (cn) := NA_real_]
    }
    data.table::set(mat_data, i = complete, j = cn,
                    value = as.numeric(mm[, raw_names[[i]]]))
  }
  invisible(clean_names)
}