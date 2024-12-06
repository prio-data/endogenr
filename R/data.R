#' example_data Panel-dataset for conflict research
#'
#' This is a slightly modified version of poldat::static_world, see <https://github.com/kvelleby/poldat>.
#'
#' @format ## `example_data`
#' A tsibble with 9010 rows and 11 columns:
#' \describe{
#'   \item{gwcode}{Slightly modified cShapes 2.0 gwcodes}
#'   \item{year}{Year}
#'   \item{gdppc}{GDP per capita sourced from PWT, WDI, Maddison, and WCDE}
#'   \item{gdppc_grwt}{Growth in GDPpc}
#'   \item{population}{Population sourced from WDI, WCDE, PWT, and Maddison}
#'   \item{yjbest}{Yao-Johnson transformed best estimate of battle-deaths from UCDP and PRIO}
#'   \item{dem}{Logit transformed liberal democracy score from V-Dem}
#'   \item{psecprop}{The proportion of the population with completed post-secondary education from WCDE}
#'   \item{best}{Best estimate of battle-deaths from UCDP and PRIO}
#'   \item{v2x_libdem}{Liberal democracy score from V-Dem}}
#'   \item{gdp}{Gross domestic product}
#' }
#' @source [PRIO Battle-deaths 3.1](https://www.prio.org/data/1)
#' @source [UCDP GED 23.1](https://ucdp.uu.se)
#' @source [cShapes 2.0](https://icr.ethz.ch/data/cshapes/)
#' @source [Varieties of Democracy 14](https://v-dem.net/)
#' @source [Penn World Tables 10.01](https://www.rug.nl/ggdc/productivity/pwt)
#' @source [Maddison Project 2020](https://www.rug.nl/ggdc/historicaldevelopment/maddison/releases/maddison-project-database-2020)
#' @source [World Development Indicators](https://data.worldbank.org/)
"example_data"
