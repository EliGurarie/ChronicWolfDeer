library(ChronicWolfDeer)

data(full_data)
data(county_borders)

# county-centroid distance matrix (drops Menominee to match full_data)
require(sf)
county_borders <- subset(county_borders, county %in% unique(full_data$county))
M_distance <- county_borders |> st_centroid() |> st_distance() |> units::drop_units()
row.names(M_distance) <- colnames(M_distance) <- county_borders$county

# baseline parameters
beta0  <- 4.2e-3
betas  <- beta0 * c(beta_mama = 30,
                    beta_myfa = 5, beta_famy = 5, beta_fama = 5, beta_mafa = 5,
                    beta_fyfa = 1, beta_fafy = 1, beta_fafa = 1)
kappas <- c(kappa_NF_y = 10, kappa_EF_y = 6.5)

p0 <- list(
  beta_mm  = betas[["beta_mama"]],
  beta_mf  = betas[["beta_mafa"]],
  beta_ffy = betas[["beta_fafa"]],
  theta_m  = 0.13,   # CWD mortality penalty - m
  theta_f  = 0.35,   # CWD mortality penalty - f
  alpha    = 2,      # selectivity infected
  xi       = 2,      # selectivity males
  eta      = 3.5,    # selectivity yearlings
  omega    = 25,     # deer killed/ wolf / year
  kappa_NF_y = kappas[["kappa_NF_y"]],   # dispersal forest
  kappa_EF_y = kappas[["kappa_EF_y"]]    # dispersal farm
)

# kappas are fixed -> dispersal kernel built once
  data_list <- prep_data(full_data, year1 = 2014, kappas, M_distance)

# log-likelihood at the baseline (overlay nothing)
  calc_logLik_SIP(p0, p0, pars_to_fit = character(0), data_list)

# evaluate a beta_ffy profile
  beta_profile <- seq(0,.4, length = 200)
  ll <- sapply(beta_profile, function(b)
    calc_logLik_SIP(list(beta_ffy = b), p0, "beta_ffy", data_list))
  plot(beta_profile, ll)


  
# Collective sample
  
  data_lists <- prep_data_all(full_data, kappas, M_distance)
  
  # thousands of times inside the Sobol/grid campaign
  calc_logLik_SIP_all(list(beta_ffy = 6e-6, alpha = 3), p0,
                      pars_to_fit = c("beta_ffy", "alpha"), data_lists)
  
