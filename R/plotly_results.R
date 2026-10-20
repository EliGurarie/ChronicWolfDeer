#' Interactive observed vs simulated plots
#'
#' Displays a two-panel plotly figure comparing observed and simulated total
#' white-tailed deer abundance and CWD prevalence across counties and years.
#' Each county is a separate trace; points are sized by year; a 1:1 reference
#' line is shown in each panel.
#'
#' @param sim Output of [run_SIP()].
#' @export 

plotly_results <- function(sim){
  # helper: one scatter trace per county, observed x, simulated y
  make_plot <- function(df, x_col, y_col, ylab) {
    counties <- unique(df$county)
    p <- plot_ly()
    for (co in counties) {
      d <- df[df$county == co, ]
      d <- d[order(d$year), ]
      p <- add_trace(p,
                     data = d,
                     x = ~get(x_col), y = ~get(y_col),
                     type = "scatter", mode = "lines+markers",
                     name = co,
                     text = ~as.character(year),
                     hovertemplate = paste0(co, "<br>year: %{text}<br>obs: %{x}<br>sim: %{y}<extra></extra>"),
                     line = list(width = 1),
                     marker = list(size = ~(year - min(df$year) + 1) * 2),
                     showlegend = FALSE
      )
    }
    p |> layout(
      title = list(text = ylab),
      xaxis = list(title = "Observed"),
      yaxis = list(title = ylab),
      hovermode = "closest",
      shapes = list(list(
        type = "line",
        x0 = min(df[[x_col]], na.rm = TRUE), 
        y0 = min(df[[x_col]], na.rm = TRUE),
        x1 = max(df[[x_col]], na.rm = TRUE), 
        y1 = max(df[[x_col]], na.rm = TRUE),
        line = list(color = "red", width = 2, dash = "dash")
      ))
    )
  }
  
  plot_df <- merge(sim, full_data, by = c("year", "county")) |>
    mutate(N_sim = I_all_sim + S_all_sim,
           I_obs = prop_cwd * wtd_count,
           p_sim = I_all_sim / (I_all_sim + S_all_sim),
           N_sim_ma = IMA_sim + SMA_sim,
           N_sim_fa = IFA_sim + SFA_sim,
           N_sim_my = IMY_sim + SMY_sim,
           N_sim_fy = IFY_sim + SFY_sim)
  
  p_total    <- make_plot(plot_df, "wtd_count", "N_sim",    "Total WTD")
  p_prev     <- make_plot(plot_df, "prop_cwd",  "p_sim",    "")
  
  subplot(p_total, p_prev, nrows = 1, margin = 0.04, titleX = TRUE, titleY = TRUE) |>
    layout(yaxis  = list(scaleanchor = "x",  scaleratio = 1),
           yaxis2 = list(scaleanchor = "x2", scaleratio = 1),
           annotations = list(
             list(text = "Total WTD",      x = 0.2, y = 1.05, xref = "paper", yref = "paper", showarrow = FALSE, xanchor = "center"),
             list(text = "CWD Prevalence", x = 0.8, y = 1.05, xref = "paper", yref = "paper", showarrow = FALSE, xanchor = "center")
           ))
}


# Internal: default colour ramps (low to high wolf density) per wolf class.
wolf_ramps <- list(none = c("gold", "darkorange2"),
                   some = c("mediumpurple", "salmon"),
                   lots = c("steelblue", "midnightblue"))


# Internal: classify counties by wolf density and assign per-county colours.
# "none" = zero/missing density; "some"/"lots" = median split of counties with
# wolves. Within a class, shade follows density (darkest = densest).
wolf_classes <- function(ramps = wolf_ramps, cutoff = NULL) {
  w <- as.data.frame(wolf_density_by_county)[, c("county", "wolf_density")]
  x <- w$wolf_density
  none <- is.na(x) | x == 0
  if (is.null(cutoff)) cutoff <- stats::median(x[!none])
  w$wolf_class <- ifelse(none, "none", ifelse(x <= cutoff, "some", "lots"))
  scale01 <- \(v) if (diff(range(v)) > 0) (v - min(v)) / diff(range(v)) else rep(0.5, length(v))
  u <- stats::ave(ifelse(none, rank(w$county), x), w$wolf_class, FUN = scale01)
  w$wolf_col <- unname(mapply(\(cl, ui) grDevices::rgb(grDevices::colorRamp(ramps[[cl]])(ui),
                                                       maxColorValue = 255), w$wolf_class, u))
  w
}

# Internal: stored county colours, or recomputed if custom ramps are given.
wolf_colours <- function(ramps = NULL) {
  if (is.null(ramps)) sf::st_drop_geometry(county_borders)[, c("county", "wolf_class", "wolf_col")]
  else wolf_classes(ramps)
}

#' Interactive county-level time series
#'
#' Plots simulated CWD prevalence or deer abundance over time as one line per
#' county, for a chosen sex/age group. Counties are coloured by wolf density
#' class, each with its own colour ramp: no wolves, some wolves, lots of
#' wolves. Within a class, darker shades mean higher wolf density. "Some" vs
#' "lots" is a median split of the counties with wolves present in
#' [wolf_density_by_county]. Hovering shows the county, its wolf class, the
#' year and the value.
#'
#' By default the y-axis is on a square-root scale, which spreads out the many
#' counties with low values while (unlike a log scale) keeping zeros. Tick
#' labels are on the original scale.
#'
#' @param output_df Expanded output from [expand_output()].
#' @param sex One of `"all"`, `"male"`, `"female"`. Default `"all"`.
#' @param age One of `"all"`, `"adult"`, `"yearling"`. Default `"all"`.
#' @param variable `"prevalence"` (proportion infected) or `"abundance"`
#'   (deer count).
#' @param transform Y-axis scale: `"sqrt"` (default) or `"identity"`.
#' @param ramps Named list of two-colour ramps (low to high wolf density) for
#'   classes `none`, `some` and `lots`.
#' @param source Plotly source ID, used to retrieve hover events in Shiny via
#'   [plotly::event_data()]. Default `"timeseries"`.
#' @return A \pkg{plotly} object.
#' @seealso [map_wolf_classes()], [plotly_results()], [mapResults_cwd()],
#'   [expand_output()]
#' @example examples/example_plot_results.R
#' @export
plotly_timeseries <- function(output_df, sex = "all", age = "all",
                              variable = c("prevalence", "abundance"),
                              transform = c("sqrt", "identity"),
                              source = "timeseries") {
  variable  <- match.arg(variable)
  transform <- match.arg(transform)
  prev <- variable == "prevalence"
  
  pp <- get_plot_params(sex, age)
  w  <- sf::st_drop_geometry(county_borders)
  d  <- data.frame(county = output_df$county, year = output_df$year,
                   value = output_df[[if (prev) pp$fill_cwd else pp$fill_wtd]])
  d$class <- w$wolf_class[match(d$county, w$county)]
  d$label <- paste0(d$county, " (", d$class, " wolves)<br>", d$year, ": ",
                    if (prev) sprintf("%.3f", d$value) else format(round(d$value), big.mark = ","))
  d$y <- if (transform == "sqrt") sqrt(d$value) else d$value
  
  ylab  <- if (prev) "Prevalence" else "Deer count"
  yaxis <- list(title = ylab)
  if (transform == "sqrt") {
    ticks <- signif(pretty(sqrt(range(d$value, na.rm = TRUE)), n = 8)^2, 2)
    yaxis <- list(title = paste(ylab, "(square-root scale)"),
                  tickvals = sqrt(ticks),
                  ticktext = format(ticks, big.mark = ",", drop0trailing = TRUE, trim = TRUE))
  }
  
  plotly::plot_ly(d, x = ~year, y = ~y, color = ~county,
                  colors = stats::setNames(w$wolf_col, w$county),
                  type = "scatter", mode = "lines",
                  customdata = ~county, source = source,
                  text = ~label, hovertemplate = "%{text}<extra></extra>") |>
    plotly::layout(showlegend = FALSE,
                   title = paste(if (prev) "CWD prevalence," else "Abundance,", pp$title_suffix),
                   xaxis = list(title = "Year"), yaxis = yaxis)
}


#' Map of counties coloured by wolf density class
#'
#' Base-graphics map of Wisconsin counties using the same per-county colours
#' as [plotly_timeseries()], serving as its legend. Optionally highlights one
#' or more counties (e.g. the one hovered over in the time series) by fading
#' all others and outlining the highlighted ones.
#'
#' @param highlight Character vector of county names to highlight, or `NULL`.
#' @param ramps Named list of two-colour ramps for classes `none`, `some`,
#'   `lots`; see [plotly_timeseries()].
#' @return Invisibly, the county class table (`county`, `wolf_density`,
#'   `class`, `col`).
#' @seealso [plotly_timeseries()]
#' @examples
#' map_wolf_classes()
#' map_wolf_classes(highlight = "Iowa")
#' 
#' @export
map_wolf_classes <- function(highlight = NULL) {
  cb  <- county_borders
  g   <- sf::st_geometry(cb)
  hit <- cb$county %in% highlight
  col <- cb$wolf_col
  if (any(hit)) col[!hit] <- grDevices::adjustcolor(col[!hit], alpha.f = 0.8)
  
#  op <- par(mar = c(0, 0, 2, 0)); on.exit(par(op))
  plot(g, col = col, border = "white", lwd = 0.5)
  if (any(hit)) {
    plot(g[hit], border = "black", lwd = 3, add = TRUE)
    text(sf::st_coordinates(sf::st_centroid(g[hit])),
                   cb$county[hit], adj = 0.5, font = 2)
  }
  title("Wolf densities", font = 1)
  # legend swatch: colour of the median-density county in each class
  o   <- order(cb$wolf_density, cb$county)
  mid <- tapply(cb$wolf_col[o], cb$wolf_class[o], \(z) z[ceiling(length(z) / 2)])
  legend("bottomleft", bty = "n", border = NA,
                   fill = mid[c("none", "some", "lots")],
                   legend = c("No wolves", "Some wolves (< 1/km^2)", "Many wolves (> 1/km^2)"))
  
  invisible(sf::st_drop_geometry(cb)[, c("county", "wolf_density", "wolf_class", "wolf_col")])
}

utils::globalVariables(c("county", "year", "y", "label"))