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


map_base <- function(output_df, sim_year, column, breaks, palette, 
                   main, sub) {
  v <- output_df[output_df$sim_year == sim_year, ]
  v <- v[match(county_borders$county, v$county), ]
  cols <- grDevices::hcl.colors(length(breaks) - 1, palette, rev = TRUE)
  fill <- cols[findInterval(v[[column]], breaks, all.inside = TRUE)]
  fill[is.na(fill)] <- "grey85"
  g  <- sf::st_geometry(county_borders)
  bb <- basemap$bbox
  
  op <- graphics::par(mar = c(0, 0, 2, 0)); on.exit(graphics::par(op))
  plot(sf::st_geometry(basemap$states), col = "grey95", border = "grey60",
       xlim = bb[c("xmin", "xmax")], ylim = bb[c("ymin", "ymax")])
  plot(sf::st_geometry(basemap$lakes), col = "#DCEBF5", border = "#9DBBD3", add = TRUE)
  plot(g, col = fill, border = "white", lwd = 0.5, add = TRUE)
  plot(g[v$wolf_presence %in% TRUE], border = "darkseagreen4", lwd = 1.5, add = TRUE)
  plot(sf::st_geometry(basemap$river), col = "#2C7FB8", lwd = 3, add = TRUE)
  title(main, font.main = 1, line = 1)
  mtext(side = 3, sub)
  graphics::legend("bottomright", fill = rev(cols), border = NA, bty = "n", cex = 0.8,
                   legend = rev(paste0(breaks[-length(breaks)]*100, "%-", breaks[-1]*100, "%")))
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
#' @example examples/example_plot_results.R
#' @export
mapResults_cwd <- function(output_df, sim_year, sex = "all", age = "all",
                           title = TRUE, sub = TRUE) {
  pp <- get_plot_params(sex, age)
  map_base(output_df, sim_year, pp$fill_cwd,
           breaks = c(0, .005, .01, .02, .05, .1, .15, .2, .3, .5, 1), palette = "RdYlGn",
           main = ifelse(title, paste0("CWD Prevalence:  ", pp$title_suffix), ""), 
           sub  = paste0("year ", sim_year + min(output_df$year) + 1))
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
#' @seealso [mapResults_cwd()], [expand_output()]
#' @example examples/example_plot_results.R
#' @export
mapResults_wtd <- function(output_df, sim_year, sex = "all", age = "all") {
  pp <- get_plot_params(sex, age)
  map_base(output_df, sim_year, pp$fill_wtd,
           breaks = seq(0, 60000, by = 6000), palette = "YlGnBu",
           main = paste0("Abundance, ", pp$title_suffix),
           sub = paste0("year ", sim_year))
}


# Internal: filter expanded output to one sim_year and attach county geometry.
attach_geometry <- function(output_df, sim_year) {
  one_year <- output_df[output_df$sim_year == sim_year, , drop = FALSE]
  geo <- merge(county_borders, one_year, by = "county")
  sf::st_as_sf(geo)
}



# Internal: county-averaged timeline (mean +/- 2 SE across counties) by group.
create_timeline <- function(output_df, y_column, group, title, ylab) {
  s <- plyr::ddply(output_df, c("sim_year", group), function(x) {
    y <- x[[y_column]]; m <- mean(y, na.rm = TRUE)
    se <- stats::sd(y, na.rm = TRUE) / sqrt(sum(!is.na(y)))
    data.frame(mean = m, low = max(m - 2 * se, 0), high = m + 2 * se)
  })
  s[[group]] <- factor(s[[group]])
  ggplot2::ggplot(s, ggplot2::aes(.data$sim_year, .data$mean,
                                  colour = .data[[group]], fill = .data[[group]])) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$low, ymax = .data$high),
                         alpha = 0.3, colour = NA) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(title = title, x = "Simulation year", y = ylab,
                  colour = group, fill = group) +
    ggplot2::theme_minimal()
}


#' Plot county-averaged time series
#'
#' Mean (+/- 2 SE across counties) CWD prevalence or deer abundance for a
#' chosen sex/age group, by simulation year, with counties split by `group`.
#'
#' @param output_df Expanded output from [expand_output()].
#' @param sex One of `"all"`, `"male"`, `"female"`. Default `"all"`.
#' @param age One of `"all"`, `"adult"`, `"yearling"`. Default `"all"`.
#' @param group Column of `output_df` used to split counties. Default
#'   `"wolf_presence"`.
#' @return A \pkg{ggplot2} object.
#' @seealso [mapResults_cwd()], [mapResults_wtd()], [expand_output()]
#' @example examples/example_plot_results.R
#' @export
plotTimeSeries_cwd <- function(output_df, sex = "all", age = "all", group = "wolf_presence") {
  pp <- get_plot_params(sex, age)
  create_timeline(output_df, pp$fill_cwd, group,
                  paste("Proportion infected,", pp$title_suffix), "Proportion") +
    ggplot2::coord_cartesian(ylim = c(0, 1))
}

#' @rdname plotTimeSeries_cwd
#' @export
plotTimeSeries_count <- function(output_df, sex = "all", age = "all", group = "wolf_presence") {
  pp <- get_plot_params(sex, age)
  create_timeline(output_df, pp$count_filter, group,
                  paste("Abundance,", pp$title_suffix), "Deer count")
}


# Internal: classify counties by wolf density and assign per-county colours.
# "none" = zero/missing density; "some"/"lots" = median split of counties with
# wolves. Within a class, shade follows density (darkest = densest).
wolf_classes <- function(ramps) {
  w <- as.data.frame(wolf_density_by_county)[, c("county", "wolf_density")]
  pos <- w$wolf_density[w$wolf_density > 0 & !is.na(w$wolf_density)]
  w$class <- ifelse(is.na(w$wolf_density) | w$wolf_density == 0, "none",
                    ifelse(w$wolf_density <= stats::median(pos), "some", "lots"))
  w <- w[order(w$wolf_density, w$county), ]
  w$col <- stats::ave(w$class, w$class,
                      FUN = \(z) grDevices::colorRampPalette(ramps[[z[1]]])(length(z)))
  w
}
