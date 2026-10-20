library(sf); library(ChronicWolfDeer)
county_borders <- subset(county_borders, county %in% unique(full_data$county))
M_distance <- county_borders |> st_centroid() |> st_distance() |> units::drop_units()
row.names(M_distance) <- colnames(M_distance) <- county_borders$county

data_list <- prep_data(full_data, year1 = 2022,
                       kappas = c(kappa_NF_y = 10, kappa_EF_y = 6.5), M_distance)
p <- list(beta_ffy = 0.30, beta_mf = 1.35, beta_mm = 2.49,
          theta_m = 0.18, theta_f = 0.50, tau = 2.7,
          alpha = 2, xi = 2, eta = 3.5, omega = 10)
out <- expand_output(run_SIP(p, n_years = 20, data_list))

mapResults_cwd(out, sim_year = 20)

# quicky animation (use x11() or other device - not in Rstudio IDE)
x11()
for(i in 1:20) mapResults_cwd(out, sim_year = i)

# abundance
mapResults_wtd(out, sim_year = 20, sex = "male", age = "yearling")


# plotly time series

plotly_timeseries(out)
plotly_timeseries(out, sex = "female", age = "adult")
