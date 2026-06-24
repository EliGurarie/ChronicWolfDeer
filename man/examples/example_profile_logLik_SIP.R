library(ChronicWolfDeer)
library(future.apply)
plan(multisession, workers = max(1, floor(availableCores() * 0.8)))

data(full_data)
data(county_borders)
require(sf)
county_borders <- subset(county_borders, county %in% unique(full_data$county))
M_distance <- county_borders |> st_centroid() |> st_distance() |> units::drop_units()
row.names(M_distance) <- colnames(M_distance) <- county_borders$county

beta0  <- .01
p0 <- list(beta_mm = beta0 * 30, beta_mf = beta0 * 5, beta_ffy = beta0,
           theta_m = 0.13, theta_f = 0.35,
           alpha = 2, xi = 2, eta = 3.5, omega = 25,
           kappa_NF_y = 10, kappa_EF_y = 6.5)
kappas <- c(kappa_NF_y = 10, kappa_EF_y = 6.5)

data_lists <- prep_data_all(full_data, kappas, M_distance)

# beta_ffy: expect a clear peak (data see transmission where CWD exists)
prof_beta <- profile_logLik_SIP("beta_ffy", beta0 * 2^seq(-3, 6, length = 25),
                                p0, data_lists)

with(prof_beta, plot(value, loglik, type = "b", log = "x"))

# xi, eta: expect near-flat (little wolf/CWD overlap in the data)
prof_xi  <- profile_logLik_SIP("xi",  seq(0.5, 6, length = 25), p0, data_lists)
prof_eta <- profile_logLik_SIP("eta", seq(0.1, 15, length = 25), p0, data_lists)

par(mfrow = c(3,3))
with(prof_xi,  plot(value, loglik, type = "b"))
with(prof_eta, plot(value, loglik, type = "b"))
