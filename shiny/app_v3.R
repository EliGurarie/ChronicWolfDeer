# write("options(shiny.autoload.r=FALSE)", ".Rprofile", append = TRUE)

a <- lapply(c("shiny", "bslib", "shinythemes", "sf", "ChronicWolfDeer"), 
            library, character = TRUE)

p_best <- list(beta_ffy = 0.30, beta_mf = 1.35, beta_mm = 2.49,
               theta_m = 0.18, theta_f = 0.50, tau = 2.7,
               alpha = 2, xi = 2, eta = 3.5, omega = 10,
               kappa_NF_y = 10, kappa_EF_y = 6.5)


wolf_ramps <- list(none = c("gold", "darkorange2"),
                    some = c("mediumpurple", "salmon"),
                    lots = c("steelblue", "midnightblue"))

county_borders <- subset(county_borders, county %in% unique(full_data$county))
M_distance <- county_borders |> st_centroid() |> st_distance() |> units::drop_units()
row.names(M_distance) <- colnames(M_distance) <- county_borders$county

groups <- c("Total" = "all_all", "Adult females" = "female_adult", "Adult males" = "male_adult")

map_tabs <- function(prefix, ...) 
  do.call(navset_pill, 
          c(unname(Map(\(lab, k)
                       nav_panel(lab, plotOutput(paste0(prefix, "_", k), height = "725px")), names(groups), groups)), list(...)))


ui_SIP <- fluidPage(
  # theme = shinytheme("sandstone"),  # Try different themes like "cosmo", "journal", etc.
  shinyjs::useShinyjs(),
  theme = bs_theme(version = 5, bootswatch = "sandstone"),
  
  titlePanel("Modeling chronic wasting disease spread in white-tailed deer"),
  hr(),
  sidebarLayout(
    
    sidebarPanel(
      numericInput("Tmax", "Simulation length (years)", value = 50, min = 2, max = 300, step = 10),
      
      actionButton("do", "Run simulation", class = "btn-primary w-100"), 
      
      accordion(open = "Predation",
                accordion_panel("Demography",
                                input_switch("fix_demo", "Fix at estimated values", value = TRUE),
                                conditionalPanel("!input.fix_demo",
                                                 sliderInput("tau", "Recruitment multiplier (tau)", 0.5, 5, p_best$tau, step = 0.05),
                                                 # TODO: revisit -- fitted with unmodified fty?
                                                 sliderInput("fty_male_modif", "Fawn-to-yearling survival multiplier (males)", 0.1, 3, 1, step = 0.05),
                                                 sliderInput("fty_female_modif", "Fawn-to-yearling survival multiplier (females)", 0.1, 3, 1, step = 0.05)
                                )
                ),
                accordion_panel("Infection",
                                input_switch("fix_inf", "Fix at estimated values", value = TRUE),
                                conditionalPanel("!input.fix_inf",
                                                 sliderInput("beta_modif", "Transmission multiplier", 0.1, 10, 1, step = 0.05),
                                                 sliderInput("theta_m", "Infected adult male survival multiplier", 0.01, 1, p_best$theta_m, step = 0.01),
                                                 sliderInput("theta_f", "Infected adult female survival multiplier", 0.01, 1, p_best$theta_f, step = 0.01)
                                )
                ),
                accordion_panel("Dispersal",
                                sliderInput("kappa_NF_y", "Male yearling dispersal scale, North (km)", 0.5, 40, p_best$kappa_NF_y, step = 0.5),
                                sliderInput("kappa_EF_y", "Male yearling dispersal scale, South (km)", 0.5, 40, p_best$kappa_EF_y, step = 0.5)
                ),
                accordion_panel("Predation",
                                input_switch("wolf_switch", "Remove wolves", value = FALSE),
                                conditionalPanel("!input.wolf_switch",
                                                 sliderInput("alpha", "Selection for infected (alpha)", 0, 20, p_best$alpha, step = 0.5),
                                                 sliderInput("xi", "Selection for males (xi)", 0.5, 10, p_best$xi, step = 0.5),
                                                 sliderInput("eta", "Selection for yearlings (eta)", 0.5, 10, p_best$eta, step = 0.5),
                                                 sliderInput("omega", "Deer killed / wolf / year (omega)", 0, 60, p_best$omega, step = 1)
                                )
                )
      )
    ),
    
    
    mainPanel(
      navset_tab(
        nav_panel("Infection", map_tabs("cwd",
                                        nav_panel("Time series",
                                                  radioButtons("ts_group", NULL, groups, inline = TRUE),
                                                  layout_columns(col_widths = c(8, 4),
                                                                 plotly::plotlyOutput("cwd_ts", height = "725px"),
                                                                 plotOutput("ts_map", height = "450px"))))),
        nav_panel("Abundance", map_tabs("wtd"))
      ),
      tags$style(HTML("
        .slider-animate-container { text-align: left; margin-top: .25rem; }
        .slider-animate-button { font-size: 2.5rem; color: var(--bs-primary); opacity: 1; }")),
      sliderInput("year", "Year", min = 1, max = 50, value = 50, step = 1, width = "100%",
                  animate = animationOptions(interval = 250, loop = FALSE))
    )
    
  ) # sidebarPayout()
) # fluidPage()

# Define server logic
server_SIP <- function(input, output, session) {
  
  dl <- reactive(prep_data(full_data, year1 = 2022, M_distance = M_distance,
                           kappas = c(kappa_NF_y = input$kappa_NF_y, kappa_EF_y = input$kappa_EF_y)))
  
  sim <- eventReactive(input$do, {
    d <- dl()
    fty <- if (input$fix_demo) c(1, 1) else c(input$fty_male_modif, input$fty_female_modif)
    d$data_year1 <- within(d$data_year1, {
      fty_survival_male   <- pmin(fty_survival_male   * fty[1], 1)
      fty_survival_female <- pmin(fty_survival_female * fty[2], 1)
      if (input$wolf_switch) wolf_abundance <- 0
    })
    demo <- if (input$fix_demo) list() else list(tau = input$tau)
    inf  <- if (input$fix_inf)  list() else
      list(theta_m = input$theta_m, theta_f = input$theta_f,
           beta_ffy = p_best$beta_ffy * input$beta_modif,
           beta_mf  = p_best$beta_mf  * input$beta_modif,
           beta_mm  = p_best$beta_mm  * input$beta_modif)
    p <- modifyList(p_best, c(list(alpha = input$alpha, xi = input$xi,
                                   eta = input$eta, omega = input$omega), demo, inf))
    run_SIP(p, n_years = input$Tmax, data_list = d)
  })
  
  out <- reactive(expand_output(sim()))
  
  observeEvent(out(), {
    n <- max(out()$sim_year)
    updateSliderInput(session, "year", max = n, value = n)
  })
  
  yr <- reactive(min(input$year, max(out()$sim_year)))
  
  output$cwd_ts <- plotly::renderPlotly({
    sa <- strsplit(input$ts_group, "_")[[1]]
    plotly_timeseries(out(), sa[1], sa[2], source = "ts", ramps = wolf_ramps) |>
      plotly::event_register("plotly_hover")
  })
  
  output$ts_map <- renderPlot({
    req(out())
    map_wolf_classes(plotly::event_data("plotly_hover", source = "ts")$customdata[[1]],
                      ramps = wolf_ramps)
  }, res = 96)
  outputOptions(output, "ts_map", priority = -1)
  
  
  lapply(groups, \(k) {
    sa <- strsplit(k, "_")[[1]]
    output[[paste0("cwd_", k)]] <- renderPlot(mapResults_cwd(out(), yr(), sa[1], sa[2]), res = 96) |>
      bindCache(input$do, yr(), "cwd", k, cache = "session")
    output[[paste0("wtd_", k)]] <- renderPlot(mapResults_wtd(out(), yr(), sa[1], sa[2]), res = 96) |>
      bindCache(input$do, yr(), "wtd", k, cache = "session")
  })
} # server_SIP



# Run the application
# and finally...
shinyApp(ui_SIP, server_SIP)
