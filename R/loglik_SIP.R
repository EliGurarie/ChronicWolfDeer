#' SIP model likelihood
#'
#' Pseudo-likelihood comparing a one-year-ahead simulation to observed data on
#' two axes: population growth (`log lambda`, normal) and CWD prevalence
#' (binomial on the tested sample). Three functions:
#'
#' * `build_comp()` assembles the per-county comparison matrix from a 2-year
#'   simulation and the matching observed rows.
#' * `loglik_SIP()` scores that matrix.
#' * `calc_logLik_SIP()` is the end-to-end wrapper: overlay fitted parameters,
#'   run the model 2 years, build the comparison, return the log-likelihood.
#'
#' @param sim,obs Simulation output ([run_SIP()]) and observed panel, each
#'   covering exactly the two years compared.
#' @param comp Matrix from `build_comp()`.
#' @param zeta SD of observed `log lambda` (fixed; empirical ~0.16).
#' @param w Weight on the growth term vs prevalence (0-1). The only knob.
#' @param eps Clipping bound for simulated prevalence.
#' @param p_fit Named natural-scale vector of fitted parameter values.
#' @param p_base Baseline parameter list.
#' @param pars_to_fit Names in `p_fit` that overlay `p_base`.
#' @param data_list Output of [prep_data()].
#' @return `build_comp`: a numeric matrix. `loglik_SIP` / `calc_logLik_SIP`:
#'   a scalar log-likelihood (`calc_logLik_SIP` returns `-1e10` on failure).
#' @name loglik_SIP
NULL


#' @rdname loglik_SIP
#' @export
build_comp <- function(sim, obs) {
  sim1 <- subset(sim, year == min(year))
  sim2 <- subset(sim, year == max(year))
  obs1 <- subset(obs, year == min(year))
  obs2 <- subset(obs, year == max(year))

  cbind(
    loglambda_obs  = log(obs2$wtd_count / obs1$wtd_count),
    loglambda_sim  = log((sim2$S_all_sim + sim2$I_all_sim) / (sim1$S_all_sim + sim1$I_all_sim)),
    n_analyzed     = obs2$n_analyzed,
    n_cwd_positive = obs2$n_cwd_positive,
    prop_cwd_sim   = sim2$I_all_sim / (sim2$S_all_sim + sim2$I_all_sim),
    N_prev         = obs1$wtd_count
  )
}


#' @rdname loglik_SIP
#' @export
loglik_SIP <- function(comp, zeta = 0.16, w = 0.5, eps = 1e-6) {
  p_clip <- pmin(pmax(comp[, "prop_cwd_sim"], eps), 1 - eps)
  na <- comp[, "n_analyzed"]; np <- comp[, "n_cwd_positive"]

  ll_prev <- ifelse(!is.na(na) & na > 0 & !is.na(np) & np <= na,
                    dbinom(np, size = na, prob = p_clip, log = TRUE),
                    ifelse(!is.na(na) & na > 0, -1e6, 0))

  ll_lam <- dnorm(comp[, "loglambda_obs"], mean = comp[, "loglambda_sim"],
                  sd = zeta, log = TRUE)

  w * sum(ll_lam) + (1 - w) * sum(ll_prev)
}


#' @rdname loglik_SIP
#' @export
calc_logLik_SIP <- function(p_fit, p_base, pars_to_fit, data_list, w = 0.5) {
  p <- p_base
  p[pars_to_fit] <- p_fit[pars_to_fit]
  p$beta_mf <- p$beta_ffy * 5
  p$beta_mm <- p$beta_ffy * 30

  sim <- tryCatch(run_SIP(p, n_years = 2, data_list),
                  error = function(e) NULL)
  if (is.null(sim)) return(-1e10)

  yrs <- c(data_list$year1, data_list$year1 + 1)
  build_comp(sim, subset(data_list$full_data, year %in% yrs)) |> loglik_SIP(w = w)
}


' Prepare data lists for all consecutive year-pairs
#'
#' Builds one [prep_data()] list per seed year from `year_start` to
#' `year_end - 1`. The dispersal kernel is identical across pairs (kappas and
#' counties are fixed), but each list carries a different `data_year1`. Call
#' this once before the fitting campaign.
#'
#' @param full_data Full county-year panel.
#' @param kappas Named vector with `kappa_NF_y` and `kappa_EF_y`.
#' @param M_distance Square county-by-county distance matrix.
#' @param year_start First seed year. Default 2014.
#' @param year_end Last observed year. Default 2022.
#' @return A named list of `data_list` objects, one per seed year.
#' @export
prep_data_all <- function(full_data, kappas, M_distance,
                          year_start = 2014, year_end = 2022) {
  years <- year_start:(year_end - 1)
  setNames(lapply(years, prep_data,
                  full_data = full_data, kappas = kappas, M_distance = M_distance),
           years)
}


#' Log-likelihood summed across all consecutive year-pairs
#'
#' Calls [calc_logLik_SIP()] for each seed year in `data_lists` and returns
#' the unweighted sum. Use [prep_data_all()] to build `data_lists` once before
#' fitting.
#'
#' @param p_fit Named natural-scale vector of fitted parameter values.
#' @param p_base Baseline parameter list.
#' @param pars_to_fit Names of parameters in `p_fit` that overlay `p_base`.
#' @param data_lists Output of [prep_data_all()].
#' @param w Weight on the growth term vs prevalence (0-1). Default 0.5.
#' @return Scalar total log-likelihood (`-1e10` if any simulation fails).
#' @export
calc_logLik_SIP_all <- function(p_fit, p_base, pars_to_fit, data_lists, w = 0.5) {
  sum(sapply(data_lists, calc_logLik_SIP,
             p_fit = p_fit, p_base = p_base, pars_to_fit = pars_to_fit, w = w))
}


#' Set up a Sobol' design matrix
#'
#' @param bounds Data frame with `lower`, `upper`, `log` columns and parameter
#'   names as row names.
#' @param N Number of Sobol' points.
#' @return Numeric matrix of `N` rows scaled to the parameter ranges.
#' @export
setupSobol <- function(bounds, N) {
  u  <- randtoolbox::sobol(N, dim = nrow(bounds))
  lo <- ifelse(bounds$log, log(bounds$lower), bounds$lower)
  hi <- ifelse(bounds$log, log(bounds$upper), bounds$upper)
  X  <- sweep(sweep(u, 2, hi - lo, "*"), 2, lo, "+")
  X[, bounds$log] <- exp(X[, bounds$log])
  colnames(X) <- rownames(bounds)
  X
}


#' Run a parallelised Sobol' likelihood sweep
#'
#' Builds a Sobol' design over `bounds`, evaluates [calc_logLik_SIP_all()] at
#' each point for the parameters named in `pars_to_fit` (which must be a subset
#' of `rownames(bounds)`), and returns the design plus log-likelihoods sorted
#' best-first.
#'
#' @param bounds Data frame with `lower`, `upper`, `log` and parameter names as
#'   row names. May cover more parameters than `pars_to_fit`.
#' @param pars_to_fit Character vector of parameter names to fit (subset of
#'   `rownames(bounds)`).
#' @param p0 Baseline parameter list.
#' @param data_lists Output of [prep_data_all()].
#' @param N Number of Sobol' points.
#' @param w Weight on growth vs prevalence term. Default 0.5.
#' @param parallel Use [future.apply::future_sapply()]? Default `TRUE`.
#' @return Data frame of `N` rows: fitted parameter columns + `loglik`, sorted
#'   descending by `loglik`.
#' @example examples/example_Sobol_sweep.R
#' @export
performSobolSweep <- function(bounds, pars_to_fit, p0, data_lists, N,
                              w = 0.5, parallel = TRUE) {
  X <- setupSobol(subset(bounds, rownames(bounds) %in% pars_to_fit), N)
  
  f <- function(i)
    calc_logLik_SIP_all(setNames(as.list(X[i, ]), pars_to_fit),
                        p0, pars_to_fit, data_lists, w)
  
  ll <- if (parallel) future.apply::future_sapply(seq_len(N), f, future.seed = NULL)
  else sapply(seq_len(N), f)
  
  res <- data.frame(X, loglik = ll)
  res[order(-res$loglik), ]
}