library(ChronicWolfDeer)

## ---- Layer 1: one annual step (internal engine) -------------------------
## step_SIP_matrix() advances a single year. Not normally called directly;
## shown here to illustrate the bottom of the stack.

## ---- Layer 2: several years on toy data ---------------------------------
## run_SIP_matrix() iterates the step. Also internal; runModel_sim() wraps it.

## ---- Layer 3: runModel_sim() -- the function you call -------------------

# Toy data
  toy <- make_toy_sip_data(n_counties = 5, n_years = 4, seed = 1)
  
  sim <- runModel_sim(toy$parameters, Tmx = 4,
                      initial_data = toy$init_data,
                      M_distance   = toy$M_distance,
                      init_year    = 2014)
  head(sim)

# Load real data

  data(full_data)
  data("county_borders")

  
# removes Menominee
  
  counties <- full_data$county |> unique()
  county_borders <- subset(county_borders, county %in% counties)
  
# obtains county distance matrix
  require(sf)
  
  M_distance <- county_borders |> st_centroid() |> st_distance() |> units::drop_units() 
  row.names(M_distance) <- colnames(M_distance) <- county_borders$county
  
# computes forward using 
  p0 <- list(alpha = 2, xi = 1.5, eta = 1.5, omega = 5, theta_m = 0.8, 
       theta_f = 0.85, beta_ffy = 0.02, beta_mf = 0.1, beta_mm = 0.6, 
       kappa_NF_y = 20, kappa_EF_y = 15, sigma_abund = 0.1)
  
  sim_realdata <- runModel_sim(p0, Tmx = 4,
                               initial_data = subset(full_data, year == 2014),
                               M_distance   = M_distance,
                               init_year    = 2014)
  
  
  
  ## ---- Layer 4: expand and map -------------------------------------------
  ## expand_output() adds aggregate counts, proportions, and wolf_presence.
  ## mapResults_*() filter to one sim_year, attach county geometry, and plot.
  
  sim_expanded <- expand_sim(sim_realdata)
  
  mapResults_cwd(sim_expanded, sim_year = 4)                 # total prevalence
  mapResults_cwd(sim_expanded, sim_year = 4, sex = "male", age = "adult")
  mapResults_wtd(sim_expanded, sim_year = 4)                 # total abundance
  
  