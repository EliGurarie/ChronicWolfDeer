library(sf)
library(ChronicWolfDeer)


data(full_data)
data(county_borders)

# match county sets (drops Menominee)
counties       <- unique(full_data$county)
county_borders <- subset(county_borders, county %in% counties)

# county-centroid distance matrix
require(sf)
M_distance <- county_borders |> st_centroid() |> st_distance() |> units::drop_units()
row.names(M_distance) <- colnames(M_distance) <- county_borders$county

p0 <- list(alpha = 2, xi = 1.5, eta = 1.5, 
           omega = 5, theta_m = 0.8,
           theta_f = 0.85, beta_ffy = 0.02, 
           beta_mf = 0.1, beta_mm = 0.6,
           tau = 2, 
           kappa_NF_y = 20, kappa_EF_y = 15, sigma_abund = 0.1)

data_list    <- prep_data(full_data, year1 = 2014,
                          kappas = p0[c("kappa_NF_y", "kappa_EF_y")], M_distance)

sim_realdata <- run_SIP(p0, n_years = 100, data_list)

# plot observed vs. modeled

plotly_results(sim_realdata)
