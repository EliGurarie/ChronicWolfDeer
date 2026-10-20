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
ggmapResults_cwd <- function(output_df, sim_year, sex = "all", age = "all") {
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
ggmapResults_wtd <- function(output_df, sim_year, sex = "all", age = "all") {
  pp <- get_plot_params(sex, age)
  data_sf <- attach_geometry(output_df, sim_year)
  create_ggmap_wtd(data_sf, pp$fill_wtd, pp$title_suffix, num_y = sim_year)
}


#' Build a CWD-prevalence county map (internal)
#' @keywords internal
create_ggmap_cwd <- function(data_sf, fill_cwd, title_suffix, num_y) {
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
create_ggmap_wtd <- function(data_sf, fill_wtd, title_suffix, num_y) {
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