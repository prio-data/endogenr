#' Creates a panel data frame based on a formula and tsibble
#'
#' This correctly handles by-group functions in a formula, unlike stats::model.frame.
#'
#' @param formula
#' @param data
#'
#' @return
#' @export
#'
#' @examples
create_panel_frame <- function(formula, data) {
  get_naive_formula <- function(data) {
    # Get key and index variables from tsibble
    grp <- tsibble::key_vars(data)
    idx <- tsibble::index_var(data)

    # Get all variable names excluding key and index
    varnames <- names(data)
    varnames <- varnames[!varnames %in% c(grp, idx)]

    # Create new formula
    new_formula <- paste(varnames[1], "~", paste(varnames[2:length(varnames)], collapse = "+")) |>
      stats::formula()

    return(new_formula)
  }

  # Get key and index variables from tsibble
  grp <- tsibble::key_vars(data)
  idx <- tsibble::index_var(data)

  # Update formula to include index variable
  formula <- stats::update(formula, paste(c(". ~ .", idx), collapse = "+"))

  # Create model frame by group
  result <- data |>
    dplyr::group_by(!!!rlang::syms(grp)) |>
    dplyr::group_modify(~ {
      mf <- stats::model.frame(formula, data = ., na.action = NULL)
      tibble::as_tibble(mf)
    }) |>
    dplyr::ungroup()

  # Clean names
  result <- result |>
    janitor::clean_names() |>
    # Ensure tsibble structure is maintained
    tsibble::as_tsibble(
      key = grp,
      index = idx
    )

  # Get naive formula for the cleaned result
  naive_formula <- get_naive_formula(result)

  return(list(
    "data" = result,
    "naive_formula" = naive_formula
  ))
}
