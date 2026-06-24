#' Profile the SIP log-likelihood over one parameter
#'
#' Sweeps a single parameter across `values`, holding all others at `p_base`,
#' and returns the log-likelihood summed over all year-pairs at each value. Use
#' to see whether the data constrain a parameter (curved profile with a peak)
#' or not (flat). Evaluations are independent and run in parallel via
#' [future.apply::future_sapply()].
#'
#' @param par Name of the parameter to profile.
#' @param values Numeric vector of values to evaluate.
#' @param p_base Baseline parameter list.
#' @param data_lists Output of [prep_data_all()].
#' @param parallel Evaluate in parallel via [future.apply::future_sapply()]?
#'   Requires a `future::plan()`. Default `FALSE`.
#' @return A data frame with columns `value` and `loglik`.
#' @export
profile_logLik_SIP <- function(par, values, p_base, data_lists, w = 0.5,
                               parallel = FALSE) {
  f <- function(v) calc_logLik_SIP_all(setNames(list(v), par), p_base, par, data_lists, w)
  loglik <- if (parallel) future.apply::future_sapply(values, f, future.seed = NULL)
  else sapply(values, f)
  data.frame(value = values, loglik = loglik)
}

