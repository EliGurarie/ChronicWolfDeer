# =============================================================================
# Map plotting for expanded SIP output
#
#   mapResults_cwd() / mapResults_wtd()  <- user-facing
#     get_plot_params()                    <- internal: pick fill column + title
#     create_map_cwd() / create_map_wtd()  <- internal: build the ggplot
#
# Maps take expand_output() output, filter to one sim_year, join county_borders
# for geometry, and plot. Wolf presence outlines come from the wolf_presence
# column added by expand_output().
# =============================================================================


#' Fill column and title pieces for a sex/age combination (internal)
#'
#' @param sex One of `"all"`, `"male"`, `"female"`.
#' @param age One of `"all"`, `"adult"`, `"yearling"`.
#' @return A list with `fill_cwd`, `fill_wtd`, `count_filter`, and
#'   `title_suffix`.
#' @keywords internal
get_plot_params <- function(sex, age) {
  combo <- paste(sex, age, sep = "_")
  
  count_filter <- switch(combo,
                         "all_all"         = "wtd_count_sim",
                         "all_adult"       = "a_all_sim",
                         "all_yearling"    = "y_all_sim",
                         "male_all"        = "m_all_sim",
                         "female_all"      = "f_all_sim",
                         "male_adult"      = "ma_sim",
                         "female_adult"    = "fa_sim",
                         "male_yearling"   = "my_sim",
                         "female_yearling" = "fy_sim"
  )
  
  fill_cwd <- switch(combo,
                     "all_all"         = "i_all_sim_perc",
                     "all_adult"       = "ia_sim_perc",
                     "all_yearling"    = "iy_sim_perc",
                     "male_all"        = "im_sim_perc",
                     "female_all"      = "if_sim_perc",
                     "male_adult"      = "ima_sim_perc",
                     "female_adult"    = "ifa_sim_perc",
                     "male_yearling"   = "imy_sim_perc",
                     "female_yearling" = "ify_sim_perc"
  )
  
  fill_wtd <- switch(combo,
                     "all_all"         = "wtd_count_sim",
                     "all_adult"       = "a_all_sim",
                     "all_yearling"    = "y_all_sim",
                     "male_all"        = "m_all_sim",
                     "female_all"      = "f_all_sim",
                     "male_adult"      = "ma_sim",
                     "female_adult"    = "fa_sim",
                     "male_yearling"   = "my_sim",
                     "female_yearling" = "fy_sim"
  )
  
  title_suffix <- switch(combo,
                         "all_all"         = "total population",
                         "all_adult"       = "all adults",
                         "all_yearling"    = "all yearlings",
                         "male_all"        = "all males",
                         "female_all"      = "all females",
                         "male_adult"      = "all adult males",
                         "female_adult"    = "all adult females",
                         "male_yearling"   = "all yearling males",
                         "female_yearling" = "all yearling females"
  )
  
  list(fill_cwd = fill_cwd, fill_wtd = fill_wtd,
       count_filter = count_filter, title_suffix = title_suffix)
}


#' Build a CWD-prevalence county map (internal)
#' @keywords internal
create_map_cwd <- function(data_sf, fill_cwd, title_suffix, num_y) {
  ggplot2::ggplot(data_sf) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[[fill_cwd]])) +
    ggplot2::geom_sf(ggplot2::aes(col = .data$wolf_presence), fill = NA) +
    ggplot2::scale_fill_distiller(
      name = NULL, type = "div", palette = "RdYlGn", direction = -1,
      limits = c(0, 1), na.value = "transparent"
    ) +
    ggplot2::scale_color_manual(
      name = "Wolf presence",
      values = c("FALSE" = "lightgrey", "TRUE" = "darkseagreen4")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste0("Proportion of infected ones out of ",
                                 title_suffix, " after ", num_y, " years"))
}


#' Build an abundance county map (internal)
#' @keywords internal
create_map_wtd <- function(data_sf, fill_wtd, title_suffix, num_y) {
  ggplot2::ggplot(data_sf) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[[fill_wtd]])) +
    ggplot2::geom_sf(ggplot2::aes(col = .data$wolf_presence), fill = NA) +
    ggplot2::scale_fill_distiller(
      name = NULL, type = "div", palette = "YlGnBu", direction = -1,
      limits = c(0, 60000), na.value = "transparent"
    ) +
    ggplot2::scale_color_manual(
      name = "Wolf presence",
      values = c("FALSE" = "lightgrey", "TRUE" = "darkseagreen4")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste0("Estimated abundance of ",
                                 title_suffix, " after ", num_y, " years"))
}


# Internal: filter expanded output to one sim_year and attach county geometry.
attach_geometry <- function(output_df, sim_year) {
  one_year <- output_df[output_df$sim_year == sim_year, , drop = FALSE]
  geo <- merge(county_borders, one_year, by = "county")
  sf::st_as_sf(geo)
}


#' Map CWD prevalence by county
#'
#' Maps the simulated proportion of infected deer in a chosen sex/age group for
#' a single simulation year, with county outlines coloured by wolf presence.
#'
#' @param output_df Expanded output from [expand_output()].
#' @param sim_year Simulation year to map.
#' @param sex One of `"all"`, `"male"`, `"female"`. Default `"all"`.
#' @param age One of `"all"`, `"adult"`, `"yearling"`. Default `"all"`.
#'
#' @return A \pkg{ggplot2} object.
#' @seealso [mapResults_wtd()], [expand_output()]
#' @export
mapResults_cwd <- function(output_df, sim_year, sex = "all", age = "all") {
  pp <- get_plot_params(sex, age)
  data_sf <- attach_geometry(output_df, sim_year)
  create_map_cwd(data_sf, pp$fill_cwd, pp$title_suffix, num_y = sim_year)
}


#' Map deer abundance by county
#'
#' Maps simulated deer abundance for a chosen sex/age group for a single
#' simulation year, with county outlines coloured by wolf presence.
#'
#' @param output_df Expanded output from [expand_output()].
#' @param sim_year Simulation year to map.
#' @param sex One of `"all"`, `"male"`, `"female"`. Default `"all"`.
#' @param age One of `"all"`, `"adult"`, `"yearling"`. Default `"all"`.
#'
#' @return A \pkg{ggplot2} object.
#' @seealso [mapResults_cwd()], [expand_output()]
#' @export
mapResults_wtd <- function(output_df, sim_year, sex = "all", age = "all") {
  pp <- get_plot_params(sex, age)
  data_sf <- attach_geometry(output_df, sim_year)
  create_map_wtd(data_sf, pp$fill_wtd, pp$title_suffix, num_y = sim_year)
}