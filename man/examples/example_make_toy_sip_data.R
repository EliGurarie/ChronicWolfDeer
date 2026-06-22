library(ChronicWolfDeerFunction)

toy <- make_toy_sip_data(n_counties = 4, n_years = 3, seed = 42)
names(toy)
str(toy$init_data, list.len = 8)
toy$parameters
