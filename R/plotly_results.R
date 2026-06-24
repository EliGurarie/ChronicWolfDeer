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
    layout(annotations = list(
      list(text = "Total WTD",      x = 0.2, y = 1.05, xref = "paper", yref = "paper", showarrow = FALSE, xanchor = "center"),
      list(text = "CWD Prevalence", x = 0.8, y = 1.05, xref = "paper", yref = "paper", showarrow = FALSE, xanchor = "center")
    ))
}
