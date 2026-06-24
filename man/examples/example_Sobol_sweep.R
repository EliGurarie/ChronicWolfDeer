unloadNamespace("ChronicWolfDeer")
library(ChronicWolfDeer)
library(randtoolbox)
library(future.apply)

data(full_data); data(county_borders)
require(sf)
county_borders <- subset(county_borders, county %in% unique(full_data$county))
M_distance <- county_borders |> st_centroid() |> st_distance() |> units::drop_units()
row.names(M_distance) <- colnames(M_distance) <- county_borders$county

beta0 <- .01

p0 <- c(beta_ffy = beta0,
        theta_m = 0.13, theta_f = 0.35,
        alpha = 2, xi = 2, eta = 3.5, omega = 25,
        kappa_NF_y = 10, kappa_EF_y = 6.5)

kappas <- c(kappa_NF_y = 10, kappa_EF_y = 6.5)
data_lists <- prep_data_all(full_data, kappas, M_distance)

# search box: lower, upper, log (TRUE -> sample geometrically)

bounds <- data.frame(
  row.names = c("beta_ffy", "alpha", "xi",  "eta", "theta_m", "theta_f", "omega"),
  lower     = c(beta0/8,    0.5,     0.5,   0.5,   0.05,      0.05,      1),
  upper     = c(beta0*8,    5,       6,     6,     0.95,      0.95,      60),
  log       = c(TRUE,       FALSE,   FALSE, FALSE, FALSE,     FALSE,    FALSE)
)

plan(multisession, workers = 15)
system.time(
  ll_profile <- performSobolSweep(bounds, 
                                  pars_to_fit = c("beta_ffy", "theta_m", "theta_f", "omega"),
                                  p0, data_lists, N = 2000)
)

require(ggplot2)
require(plyr)
ggplot(ll_profile |> subset(loglik > -700) |> 
         mutate(omega_factor = cut(omega, seq(0,60,5))), 
       aes(theta_m, loglik, 
           col = log(beta_ffy))) + 
  geom_point() + 
  facet_wrap(.~omega_factor)