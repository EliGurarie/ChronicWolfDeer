#' Build a small synthetic data set for examples and tests
#'
#' Generates a self-consistent toy data set with the column layout the model
#' expects, so that the documented examples and tests run without the real
#' Wisconsin data. Values are plausible but arbitrary.
#'
#' @param n_counties Number of fake counties. Default 5.
#' @param n_years Number of calendar years in `full_data`, starting at
#'   `init_year`. Default 5.
#' @param init_year First calendar year. Default 2014.
#' @param seed Optional integer seed for reproducibility. Default 1.
#'
#' @return A list with:
#'   * `init_data`  - one row per county (the `init_year` slice), suitable for
#'     [run_SIP_matrix()] and as `initial_data` for [runModel_sim()].
#'   * `full_data`  - `n_counties * n_years` rows.
#'   * `M_distance` - symmetric `n_counties`-by-`n_counties` distance matrix
#'     (km), ordered to match `init_data$county`.
#'   * `parameters` - a raw-scale parameter list with all entries
#'     [step_SIP_matrix()] requires.
#'
#' @export
#' @example man/examples/example_make_toy_sip_data.R
#' 
make_toy_sip_data <- function(n_counties = 5, n_years = 5,
                              init_year = 2014, seed = 1) {
  if (!is.null(seed)) set.seed(seed)

  counties <- paste0("County_", LETTERS[seq_len(n_counties)])

  # roughly half the counties are "northern" (label contains "North")
  region <- ifelse(seq_len(n_counties) <= ceiling(n_counties / 2),
                   "Northern Forest", "Eastern Farmland")

  # county centroids (km) and pairwise distance matrix
  xy <- cbind(stats::runif(n_counties, 0, 300),
              stats::runif(n_counties, 0, 300))
  M_distance <- as.matrix(stats::dist(xy))
  dimnames(M_distance) <- list(counties, counties)

  one_year <- function(year) {
    data.frame(
      county = counties,
      year   = year,
      # seed counts
      S.M.A0 = round(stats::runif(n_counties, 6000, 10000)),
      S.F.A0 = round(stats::runif(n_counties, 7000, 11000)),
      S.M.Y0 = round(stats::runif(n_counties, 1000, 2000)),
      S.F.Y0 = round(stats::runif(n_counties, 1500, 2800)),
      I.M.A0 = round(stats::runif(n_counties, 50, 400)),
      I.F.A0 = round(stats::runif(n_counties, 10, 120)),
      I.M.Y0 = round(stats::runif(n_counties, 10, 180)),
      I.F.Y0 = round(stats::runif(n_counties, 2, 60)),
      # observed summaries used by the fit wrapper
      wtd_count = round(stats::runif(n_counties, 18000, 28000)),
      prop_cwd  = round(stats::runif(n_counties, 0.005, 0.05), 4),
      # demography / mortality rates
      mean_mu_male_harv    = round(stats::runif(n_counties, 0.20, 0.30), 3),
      mean_mu_fem_harv     = round(stats::runif(n_counties, 0.25, 0.32), 3),
      mean_wtd_fdr         = round(stats::runif(n_counties, 0.85, 1.05), 3),
      nh_mu_adult_human    = 0.0222,
      nh_mu_yearling_human = 0.0425,
      nh_mu_adult_other    = 0.0335,
      nh_mu_yearling_other = 0.0635,
      fty_survival_male    = round(stats::runif(n_counties, 0.25, 0.35), 3),
      fty_survival_female  = round(stats::runif(n_counties, 0.38, 0.45), 3),
      # context
      dnr_region        = region,
      K                 = round(stats::runif(n_counties, 35000, 50000)),
      wolf_presence     = stats::runif(n_counties) > 0.3,
      wolf_abundance    = round(stats::runif(n_counties, 0, 40), 1),
      wolf_density_rank = sample(c("low", "medium", "high"), n_counties, replace = TRUE),
      stringsAsFactors  = FALSE
    )
  }

  full_data <- do.call(rbind, lapply(seq(init_year, init_year + n_years - 1), one_year))
  init_data <- full_data[full_data$year == init_year, , drop = FALSE]
  rownames(init_data) <- NULL

  parameters <- list(
    alpha = 2, xi = 1.5, eta = 1.5,        # predation preferences
    omega = 5,                             # deer killed per wolf per year
    theta_m = 0.80, theta_f = 0.85,        # infected-adult survival multipliers
    beta_ffy = 0.02, beta_mf = 0.10, beta_mm = 0.60,
    kappa_NF_y = 20, kappa_EF_y = 15,      # dispersal kernel scales
    sigma_abund = 0.1
  )

  list(init_data = init_data, full_data = full_data,
       M_distance = M_distance, parameters = parameters)
}
