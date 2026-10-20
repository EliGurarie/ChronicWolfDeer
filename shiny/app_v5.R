a <- lapply(c("shiny", "bslib", "sf", "ChronicWolfDeer", "plotly"), library, character.only = TRUE)


p1 <- list(beta_ffy = 0.30, beta_mf = 1.35, beta_mm = 2.49,
           theta_m = 0.18, theta_f = 0.50, tau = 2.7,
           alpha = 2, xi = 2, eta = 3.5, omega = 10,
           kappa_NF_y = 10, kappa_EF_y = 6.5)

p_best <- list(beta_ffy = 0.30, beta_mf = 1.35, beta_mm = 2.49,
               theta_m = 0.18, theta_f = 0.50, tau = 2.7,
               alpha = 2, xi = 2, eta = 3.5, omega = 10,
               kappa_NF_y = 30, kappa_EF_y = 20)


# Globals
{
  
  
  wolf_ramps <- list(none = c("gold", "orange"),
                     some = c("salmon", "mediumpurple"),
                     lots = c("steelblue", "midnightblue"))
  
  county_borders <- merge(county_borders[, c("county", "geometry")],
                          ChronicWolfDeer:::wolf_classes(wolf_ramps))
  
  county_borders <- subset(county_borders, county %in% unique(full_data$county))
  M_distance <- county_borders |> st_centroid() |> st_distance() |> 
    units::set_units("km") |> units::drop_units()
  
  row.names(M_distance) <- colnames(M_distance) <- county_borders$county
  
  # ---- Model fit: best-fit baseline, log-likelihood helper, sensitivity ranges ----
  kap_best <- unlist(p_best[c("kappa_NF_y", "kappa_EF_y")])
  
  # log-likelihood of a parameter list; calc_logLik_SIP rebuilds beta_mf/mm from the C multipliers
  ll_of <- function(p, dls) {
    p <- c(p, C_betamf = p$beta_mf / p$beta_ffy, C_betamm = p$beta_mm / p$beta_ffy)
    calc_logLik_SIP_all(p, p, character(0), dls)
  }
  
  fit_best <- run_SIP(p_best, n_years = 9, prep_data(full_data, 2014, kap_best, M_distance))
  ll_best  <- ll_of(p_best, prep_data_all(full_data, kap_best, M_distance))
  
  sens_range <- list(tau = c(1, 5), beta_ffy = c(0.05, 1), beta_mf = c(0.2, 5), beta_mm = c(0.5, 10),
                     theta_m = c(0.05, 1), theta_f = c(0.05, 1), alpha = c(0, 20), xi = c(0.5, 10),
                     eta = c(0.5, 10), omega = c(0, 60), kappa_NF_y = c(1, 40), kappa_EF_y = c(1, 40))
  
  groups <- c("Total" = "all_all", "Adult females" = "female_adult", "Adult males" = "male_adult")
  
  # Sub-tabs: one map per group, plus a time series with its wolf-class legend map
  map_tabs <- function(prefix) do.call(navset_pill, c(
    unname(Map(\(lab, k) nav_panel(lab, plotOutput(paste0(prefix, "_", k), height = "725px")),
               names(groups), groups)),
    list(nav_panel("Time series",
                   radioButtons(paste0("ts_group_", prefix), NULL, groups, inline = TRUE),
                   layout_columns(col_widths = c(8, 4),
                                  plotly::plotlyOutput(paste0("ts_", prefix), height = "725px"),
                                  plotOutput(paste0("ts_map_", prefix), height = "450px"))),
         id = paste0("sub_", prefix))))
  
  # Year slider is hidden on time-series sub-tabs and on the model-fit tab
  ts_showing <- "(input.main == 'Infection' && input.sub_cwd == 'Time series') ||
                 (input.main == 'Abundance' && input.sub_wtd == 'Time series') ||
                 input.main == 'Model fit' || input.main == 'About' || input.main == 'Contact'"
  
  sens_par_choices <- list(
    Demography = c("Recruitment multiplier (\u03C4)" = "tau"),
    Infection  = c("Baseline transmission, F\u2013F & yearlings (\u03B2_FFY)" = "beta_ffy",
                   "Male\u2013female transmission (\u03B2_MF)"               = "beta_mf",
                   "Male\u2013male transmission (\u03B2_MM)"                 = "beta_mm",
                   "CWD survival multiplier, adult males (\u03B8_M)"         = "theta_m",
                   "CWD survival multiplier, adult females (\u03B8_F)"       = "theta_f"),
    Predation  = c("Selection for infected (\u03B1)"    = "alpha",
                   "Selection for males (\u03BE)"       = "xi",
                   "Selection for yearlings (\u03B7)"   = "eta",
                   "Deer killed / wolf / year (\u03C9)" = "omega"),
    Dispersal  = c("Dispersal scale, North, km (\u03BA_N)" = "kappa_NF_y",
                   "Dispersal scale, South, km (\u03BA_S)" = "kappa_EF_y"))
  
  sens_label <- setNames(names(unlist(unname(sens_par_choices))), unlist(sens_par_choices))
  
}




ui_SIP <- fluidPage(
  theme = bs_theme(version = 5, bootswatch = "sandstone"),
  
  titlePanel("Can Ma'iingan Keep Chronic Wasting Disease Away?"),
  hr(),
  sidebarLayout(
    
    sidebarPanel(width = 3,
                 numericInput("Tmax", "Simulation length (years)", value = 50, min = 2, max = 300, step = 10),
                 actionButton("do", "Run simulation", class = "btn-primary w-100"),
                 
                 accordion(open = "Predation",
                           accordion_panel("Demography",
                                                   input_switch("fix_demo", "Fix at estimated values", value = TRUE),
                                                   sliderInput("tau", HTML("Recruitment multiplier (&tau;)"), 0.5, 5, p_best$tau, step = 0.05),
                                                   sliderInput("fty_male_modif",   HTML("Fawn-to-yearling survival multiplier, males (&phi;<sub>FY,M</sub>)"),   0.1, 3, 1, step = 0.05),
                                                   sliderInput("fty_female_modif", HTML("Fawn-to-yearling survival multiplier, females (&phi;<sub>FY,F</sub>)"), 0.1, 3, 1, step = 0.05)
                           ),
                           accordion_panel("Infection",
                                           input_switch("fix_inf", "Fix at estimated values", value = TRUE),
                                           sliderInput("beta_ffy", HTML("Baseline transmission, F&ndash;F &amp; yearlings (&beta;<sub>FFY</sub>)"), 0.01, 3, p_best$beta_ffy, step = 0.01),
                                           sliderInput("K1", HTML("Male&ndash;female multiplier (K<sub>1</sub> = &beta;<sub>MF</sub> / &beta;<sub>FFY</sub>)"), 1, 30, round(p_best$beta_mf / p_best$beta_ffy, 1), step = 0.5),
                                           sliderInput("K2", HTML("Male&ndash;male multiplier (K<sub>2</sub> = &beta;<sub>MM</sub> / &beta;<sub>FFY</sub>)"),   1, 50, round(p_best$beta_mm / p_best$beta_ffy, 1), step = 0.5),
                                           sliderInput("theta_m", HTML("CWD survival multiplier, adult males (&theta;<sub>M</sub>)"),   0.01, 1, p_best$theta_m, step = 0.01),
                                           sliderInput("theta_f", HTML("CWD survival multiplier, adult females (&theta;<sub>F</sub>)"), 0.01, 1, p_best$theta_f, step = 0.01)
                           ),
                           accordion_panel("Dispersal",
                                            sliderInput("kappa_NF_y", HTML("Dispersal scale, North, km (&kappa;<sub>N</sub>)"), 0.5, 40, p_best$kappa_NF_y, step = 0.5),
                                            sliderInput("kappa_EF_y", HTML("Dispersal scale, South, km (&kappa;<sub>S</sub>)"), 0.5, 40, p_best$kappa_EF_y, step = 0.5)
                            ),
                            accordion_panel("Predation",
                                            input_switch("wolf_switch", "Remove wolves", value = FALSE),
                                            conditionalPanel("!input.wolf_switch",
                                                             sliderInput("alpha", HTML("Selection for infected (&alpha;)"),    0, 20, p_best$alpha, step = 0.5),
                                                             sliderInput("xi",    HTML("Selection for males (&xi;)"),          0.5, 10, p_best$xi, step = 0.5),
                                                             sliderInput("eta",   HTML("Selection for yearlings (&eta;)"),     0.5, 10, p_best$eta, step = 0.5),
                                                             sliderInput("omega", HTML("Deer killed / wolf / year (&omega;)"), 0, 60, p_best$omega, step = 1)
      )
    )
  )
), 

mainPanel(width = 9,
          navset_tab(id = "main",
                     nav_panel("Infection", map_tabs("cwd")),
                     nav_panel("Abundance", map_tabs("wtd")),
                     nav_panel("Model fit", navset_pill(id = "sub_fit",
                                                        nav_panel("Current settings",
                                                                  h5(textOutput("fit_ll")),
                                                                  plotly::plotlyOutput("fit_current", height = "650px", width = "1200px")),
                                                        nav_panel("Best fit", plotly::plotlyOutput("fit_best", height = "650px", width = "1200px")),
                                                        nav_panel("Sensitivity",
                                                                  layout_columns(col_widths = c(3, 2, 2, 3, 2),
                                                                                 selectInput("sens_par", "Parameter", sens_par_choices),
                                                                                 numericInput("sens_lo", "Low",  sens_range[[1]][1]),
                                                                                 numericInput("sens_hi", "High", sens_range[[1]][2]),
                                                                                 sliderInput("sens_n", "Grid points", 5, 30, 15, step = 1),
                                                                                 actionButton("do_sens", "Run sweep", class = "btn-primary mt-4")),
                                                                  plotly::plotlyOutput("sens_plot", height = "550px")))),
                     nav_panel("About",
                               div(style = "max-width: 900px; margin-top: 30px; font-size: 17px;",
                                   h2("About this app"),
                                   p("This app explores how wolf predation, and in particular wolves' selection for
                            deer infected with chronic wasting disease (CWD), might affect the spread of CWD
                            and the abundance of white-tailed deer in Wisconsin."),
                                   p("It runs a county-level, annual metapopulation model of deer demography and CWD,
                            with eight sex (male/female), age (yearling/adult) and infection
                            (susceptible/infected) classes. Each year the model steps through
                            (1) recruitment and survival, including harvest, other mortality and wolf predation,
                            (2) CWD transmission, and (3) dispersal of yearling males between counties.
                            Deer and CWD parameters were fitted to Wisconsin DNR data from 2014\u20132022;
                            wolf densities come from DNR wolf maps."),
                                   h4("How to use it"),
                                   tags$ol(
                                     tags$li("Set parameters in the sidebar. Demography and infection are fixed at their
                                    fitted values unless unlocked; predation and dispersal can be varied freely."),
                                     tags$li("Click ", strong("Run simulation"), "."),
                                     tags$li(strong("Infection"), " and ", strong("Abundance"), ": maps by year (use the
                                    slider or play button) and county time series."),
                                     tags$li(strong("Model fit"), ": observed vs.µ simulated data for 2014\u20132022, and
                                    how the fit changes with each parameter.")),
                                   h4("Caveats"),
                                   p("The model is exploratory. Wolf selection for infected deer (\u03B1) has not been
                            measured, so results describe what ", em("could"), " happen under different
                            assumptions, not predictions."),
                                   h4("Contact"),
                                   HTML("<p><strong>Ekaterina Khadonova</strong><br>
                               Email: <a href='mailto:ekhadonova@gmail.com'>ekhadonova@gmail.com</a><br>
                               LinkedIn: <a href='https://www.linkedin.com/in/ekaterina-khadonova/' target='_blank'>linkedin.com/in/ekaterina-khadonova/</a></p>"),
                                   HTML("<p><strong>Elie Gurarie</strong>, SUNY College of Environmental Science and Forestry<br>
                               Email: <a href='mailto:egurarie@esf.edu'>egurarie@esf.edu</a></p>"),
                                   div(style = "margin-top: 40px; display: flex; justify-content: center; align-items: center; gap: 50px;",
                                       tags$img(src = "esf-logo.png",    height = "100px"),
                                       tags$img(src = "uw-logo.png",     height = "100px"),
                                       tags$img(src = "glifwc-logo.png", height = "120px",
                                                style = "background: #5a5a5a; border-radius: 50%; padding: 6px;"))))          ),
          tags$style(HTML("
        .slider-animate-container { text-align: left; margin-top: .25rem; }
        .slider-animate-button { font-size: 2.5rem; color: var(--bs-primary); opacity: 1; }")),
          conditionalPanel(paste0("!(", ts_showing, ")"),
                           sliderInput("year", "Year", min = 1, max = 50, value = 50, step = 1, width = "100%",µ
                                       animate = animationOptions(interval = 250, loop = FALSE))
          )
) 
) 
) 

server_SIP <- function(input, output, session) {
  
  kap <- reactive(c(kappa_NF_y = input$kappa_NF_y, kappa_EF_y = input$kappa_EF_y))
  
  dl <- reactive(prep_data(full_data, year1 = 2022, M_distance = M_distance, kappas = kap()))
  observeEvent(input$sens_par, {
    updateNumericInput(session, "sens_lo", value = sens_range[[input$sens_par]][1])
    updateNumericInput(session, "sens_hi", value = sens_range[[input$sens_par]][2])
  })
  
  
  observeEvent(input$fix_inf, if (input$fix_inf) {
    updateSliderInput(session, "beta_ffy", value = p_best$beta_ffy)
    updateSliderInput(session, "K1", value = round(p_best$beta_mf / p_best$beta_ffy, 1))
    updateSliderInput(session, "K2", value = round(p_best$beta_mm / p_best$beta_ffy, 1))
    updateSliderInput(session, "theta_m", value = p_best$theta_m)
    updateSliderInput(session, "theta_f", value = p_best$theta_f)
  })
  observeEvent(input$fix_demo, if (input$fix_demo) {
    updateSliderInput(session, "tau", value = p_best$tau)
    updateSliderInput(session, "fty_male_modif", value = 1)
    updateSliderInput(session, "fty_female_modif", value = 1)
  })
  
  # current slider settings -> parameter list
  pars <- reactive({
    demo <- if (input$fix_demo) list() else list(tau = input$tau)
    inf  <- if (input$fix_inf) list() else
      list(theta_m = input$theta_m, theta_f = input$theta_f,
           beta_ffy = input$beta_ffy,
           beta_mf  = input$beta_ffy * input$K1,
           beta_mm  = input$beta_ffy * input$K2)
    modifyList(p_best, c(list(alpha = input$alpha, xi = input$xi,
                              eta = input$eta, omega = input$omega), demo, inf))
  })
  
  # fawn-to-yearling and wolf-removal modifications to a data_list
  mod_data <- function(d) {
    fty <- if (input$fix_demo) c(1, 1) else c(input$fty_male_modif, input$fty_female_modif)
    d$data_year1 <- within(d$data_year1, {
      fty_survival_male   <- pmin(fty_survival_male   * fty[1], 1)
      fty_survival_female <- pmin(fty_survival_female * fty[2], 1)
      if (input$wolf_switch) wolf_abundance <- 0
    })
    d
  }
  
  # forecast from the last data year
  sim <- eventReactive(input$do, run_SIP(pars(), n_years = input$Tmax, data_list = mod_data(dl())))
  
  # fit period 2014-2022: free-running simulation + one-year-ahead log-likelihood
  fit <- eventReactive(input$do, list(
    sim = run_SIP(pars(), n_years = 9, mod_data(prep_data(full_data, 2014, kap(), M_distance))),
    ll  = ll_of(pars(), lapply(prep_data_all(full_data, kap(), M_distance), mod_data))))
  
  out <- reactive(expand_output(sim()))
  
  observeEvent(out(), {
    n <- max(out()$sim_year)
    updateSliderInput(session, "year", max = n, value = n)
  })
  
  yr <- reactive(min(input$year, max(out()$sim_year)))
  
  # Maps: one output per variable x group, frames cached per run
  lapply(groups, \(k) {
    sa <- strsplit(k, "_")[[1]]
    output[[paste0("cwd_", k)]] <- renderPlot(mapResults_cwd(out(), yr(), sa[1], sa[2]), res = 96) |>
      bindCache(input$do, yr(), "cwd", k, cache = "session")
    output[[paste0("wtd_", k)]] <- renderPlot(mapResults_wtd(out(), yr(), sa[1], sa[2]), res = 96) |>
      bindCache(input$do, yr(), "wtd", k, cache = "session")
  })
  
  # Time series + hover-linked legend map, for prevalence (cwd) and abundance (wtd)
  lapply(c(cwd = "prevalence", wtd = "abundance"), \(v) {
    prefix <- if (v == "prevalence") "cwd" else "wtd"
    src    <- paste0("ts_", prefix)
    
    output[[src]] <- plotly::renderPlotly({
      sa <- strsplit(input[[paste0("ts_group_", prefix)]], "_")[[1]]
      plotly_timeseries(out(), sa[1], sa[2], variable = v,
                        transform = if (v == "prevalence") "sqrt" else "identity",
                        source = src) |>
        plotly::event_register("plotly_hover")
    })
    
    output[[paste0("ts_map_", prefix)]] <- renderPlot({
      req(out())
      map_wolf_classes(plotly::event_data("plotly_hover", source = src)$customdata[[1]])
    }, res = 96)
    outputOptions(output, paste0("ts_map_", prefix), priority = -1)
  })
  
  # Model fit: observed vs simulated, current settings and best fit
  output$fit_current <- plotly::renderPlotly(plotly_results(fit()$sim))
  output$fit_best    <- plotly::renderPlotly(plotly_results(fit_best))
  output$fit_ll <- renderText(sprintf("log-likelihood: %.1f   (best fit %.1f;  \u0394 = %+.1f)",
                                      fit()$ll, ll_best, fit()$ll - ll_best))
  
  # Sensitivity: one-at-a-time sweep of a parameter, others at current settings
  sens <- eventReactive(input$do_sens, {
    par <- input$sens_par; p0 <- pars(); kp <- kap()
    req(input$sens_lo < input$sens_hi)
    x <- seq(input$sens_lo, input$sens_hi, length = input$sens_n)
    dls <- lapply(prep_data_all(full_data, kp, M_distance), mod_data)
    ll  <- withProgress(message = paste("Sweeping", par), value = 0,
                        vapply(seq_along(x), \(i) {
                          incProgress(1 / length(x))
                          v <- x[i]
                          if (par %in% names(kp)) {
                            kp[par] <- v
                            ll_of(p0, lapply(prep_data_all(full_data, kp, M_distance), mod_data))
                          } else ll_of(modifyList(p0, setNames(list(v), par)), dls)
                        }, numeric(1)))
    list(par = par, d = data.frame(x, dll = ll - ll_best),
         now = if (par %in% names(kp)) kp[[par]] else p0[[par]], best = p_best[[par]])
  })
  
  output$sens_plot <- plotly::renderPlotly({
    s  <- sens()
    vl <- \(x, col) list(type = "line", x0 = x, x1 = x, yref = "paper", y0 = 0, y1 = 1,
                         line = list(color = col, dash = "dash"))
    plotly::plot_ly(s$d, x = ~x, y = ~dll, type = "scatter", mode = "lines+markers") |>
      plotly::layout(xaxis = list(title = s$par),
                     yaxis = list(title = "\u0394 log-likelihood vs best fit"),
                     shapes = list(vl(s$best, "red"), vl(s$now, "grey"),
                                   list(type = "line", xref = "paper", x0 = 0, x1 = 1,
                                        y0 = -1.92, y1 = -1.92,
                                        line = list(color = "black", dash = "dot"))))
  })
} 

shinyApp(ui_SIP, server_SIP)