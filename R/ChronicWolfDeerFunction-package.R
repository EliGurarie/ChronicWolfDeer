#' ChronicWolfDeerFunction: Wolf-CWD SIP model for Wisconsin deer
#'
#' Tools to simulate white-tailed deer demography and chronic wasting
#' disease (CWD) dynamics across Wisconsin counties using a discrete-time,
#' sex-age structured, susceptible-infected-perished (SIP) meta-population
#' model.
#'
#' The typical workflow is:
#'
#' 1. Assemble a county-year data frame (see [make_toy_sip_data()] for the
#'    required column layout).
#' 2. Iterate the model with [run_SIP_matrix()] (single seed year) or the
#'    higher-level wrappers [runModel_sim()].
#' 3. Reshape the matrix output with [output_transform()].
#'
#' The one-step engine [step_SIP_matrix()] is exposed for testing and for
#' users who want to drive the projection manually.
#'
#' @keywords internal
"_PACKAGE"

## Imports used across the package -------------------------------------------
#' @importFrom magrittr %>%
#' @importFrom rlang .data sym
#' @importFrom stats plogis qlogis runif dist
NULL

# Note: `!!` (quasiquotation) is handled by the tidy-eval data mask inside
# dplyr verbs and does not require an explicit import.
NULL

## Quiet R CMD check NOTEs about non-standard-evaluation column names that
## are created or referenced inside dplyr/tidyr pipelines.
utils::globalVariables(c(
  # input / observed columns
  "year", "county", "wtd_count", "prop_cwd", "n_analyzed", "n_cwd_positive",
  "K", "dnr_region", "wolf_presence", "wolf_density_rank",
  "fty_survival_male", "fty_survival_female",
  "S.M.A0", "S.F.A0", "S.M.Y0", "S.F.Y0",
  "I.M.A0", "I.F.A0", "I.M.Y0", "I.F.Y0",
  # pivot helpers / simulated columns
  "pop_type", "sim_year", "rowid",
  "S.M.A", "S.F.A", "S.M.Y", "S.F.Y", "I.M.A", "I.F.A", "I.M.Y", "I.F.Y",
  "ma_sim", "fa_sim", "my_sim", "fy_sim",
  "m_all_sim", "f_all_sim", "a_all_sim", "y_all_sim",
  "i_all_sim", "s_all_sim", "wtd_count_sim",
  # package constant
  "sip_classes"
))
