df <- poldat::static_world |>
  dplyr::select(gwcode, year, gdppc, gdppc_grwt, population, best, psecprop, v2x_polyarchy) |>
  dplyr::group_by(gwcode) |>
  dplyr::arrange(year)

limited_filled_df <- df |> dplyr::filter(year >= 2018) |>
  tidyr::fill(gdppc, psecprop, population, best, v2x_polyarchy, .direction = "down")

example_data <- dplyr::bind_rows(
  df |> dplyr::filter(year < 2018),
  limited_filled_df
  ) |> dplyr::filter(year <= 2024)

usethis::use_data(example_data, overwrite = TRUE)
