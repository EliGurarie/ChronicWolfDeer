#' Null-coalescing operator
#'
#' Returns `a` unless it is `NULL`, in which case it returns `b`. Used to
#' supply defaults for optional parameter-list entries.
#'
#' @param a,b Values; `a` is returned unless `NULL`.
#'
#' @return `a` if non-`NULL`, otherwise `b`.
#' @keywords internal
#' @examples
#' ChronicWolfDeerFunction:::`%||%`(NULL, 5)
#' ChronicWolfDeerFunction:::`%||%`(2, 5)
`%||%` <- function(a, b) if (!is.null(a)) a else b


#' Recover a parameter on the raw (natural) scale
#'
#' Looks up a single base parameter in a named list/vector, accepting the value
#' on the raw scale or on a transformed scale, and returns it back-transformed
#' to the raw scale. Transformed aliases are recognised by name prefix:
#'
#' * `log_<base>`    -> back-transformed with [base::exp()]
#' * `qlogis_<base>` -> back-transformed with [stats::plogis()]
#'
#' The lookup order is raw, then `log_`, then `qlogis_`.
#'
#' @param params A named list or named numeric vector of parameters.
#' @param base Character scalar; the base parameter name (e.g. `"beta_ffy"`).
#'
#' @return The parameter value on the raw scale, or `NULL` if no raw or
#'   recognised transformed alias is present.
#' @keywords internal
#' @examples
#' p <- list(beta_ffy = 0.02, log_alpha = log(2), qlogis_eta = qlogis(0.6))
#' ChronicWolfDeerFunction:::get_raw_param(p, "beta_ffy") # 0.02
#' ChronicWolfDeerFunction:::get_raw_param(p, "alpha")    # 2
#' ChronicWolfDeerFunction:::get_raw_param(p, "eta")      # 0.6
#' ChronicWolfDeerFunction:::get_raw_param(p, "missing")  # NULL
get_raw_param <- function(params, base) {
  if (is.null(names(params))) stop("Parameters must be a named list or vector.")

  # 1) raw, as-is
  if (base %in% names(params)) {
    return(params[[base]])
  }
  # 2) log-transformed -> exp
  logn <- paste0("log_", base)
  if (logn %in% names(params)) {
    return(exp(params[[logn]]))
  }
  # 3) logit-transformed -> plogis
  logitn <- paste0("qlogis_", base)
  if (logitn %in% names(params)) {
    return(plogis(params[[logitn]]))
  }
  # not found in any form
  NULL
}


#' Coerce a parameter object to a raw, model-ready list
#'
#' Converts a mixed-scale parameter object (raw and/or `log_`/`qlogis_`
#' transformed entries) into a flat, raw-scale named list with canonical names,
#' suitable to pass as the `p` argument of [step_SIP_matrix()].
#'
#' Base names handled explicitly (back-transformed if supplied on a transformed
#' scale) are `beta_ffy`, `alpha`, `xi`, `eta`, and `sigma_abund`. Any other
#' fields (e.g. `kappa_NF_y`, `kappa_EF_y`, `omega`, `theta_m`, `theta_f`,
#' `beta_mf`, `beta_mm`) are carried forward untouched.
#'
#' If absent, `beta_mf` and `beta_mm` are derived from `beta_ffy` using the
#' rule-of-thumb multipliers 5 and 30, respectively.
#'
#' @section Note:
#' Only the names in `base_names` are back-transformed. Transformed aliases of
#' the pass-through parameters (for example `qlogis_theta_m`) are **not**
#' recognised and would be carried forward under their transformed name. Extend
#' `base_names` if you intend to fit those on a transformed scale. This helper
#' is provided as a standalone utility; the simulation wrappers
#' ([run_SIP()]) currently expect raw-scale input
#' and do not call it automatically.
#'
#' @param params A named list or named numeric vector of parameters, possibly
#'   containing `log_`/`qlogis_` transformed entries.
#'
#' @return A named list of raw-scale parameters with `NULL` entries dropped.
#' @export
#' @examples
#' fit_scale <- list(
#'   log_beta_ffy = log(0.02),
#'   log_alpha    = log(2),
#'   qlogis_eta   = qlogis(0.6),
#'   xi           = 1.5,
#'   sigma_abund  = 0.1,
#'   theta_m      = 0.8, theta_f = 0.85,
#'   kappa_NF_y   = 20,  kappa_EF_y = 15,
#'   omega        = 5
#' )
#' coerce_to_raw_params(fit_scale)
coerce_to_raw_params <- function(params) {
  base_names <- c("beta_ffy", "alpha", "xi", "eta", "sigma_abund")

  out <- lapply(base_names, function(nm) get_raw_param(params, nm))
  names(out) <- base_names
  out <- as.list(out)

  # carry forward extras (kappas, omegas, thetas, derived betas) untouched,
  # dropping any transformed aliases of the base names we already resolved.
  drop_aliases <- unlist(lapply(
    base_names, function(x) c(x, paste0("log_", x), paste0("qlogis_", x))
  ))
  extras <- setdiff(names(params), drop_aliases)
  for (nm in extras) out[[nm]] <- params[[nm]]

  # derive cross-sex / male betas if not supplied
  if (is.null(out$beta_mf) && !is.null(out$beta_ffy)) out$beta_mf <- out$beta_ffy * 5
  if (is.null(out$beta_mm) && !is.null(out$beta_ffy)) out$beta_mm <- out$beta_ffy * 30

  if (!is.null(out$sigma_abund) && out$sigma_abund <= 0) {
    stop("sigma_abund must be > 0 on raw scale. Got: ", out$sigma_abund)
  }

  Filter(Negate(is.null), out)
}




#' Obtain parameters from fitted set
#'
#' Takes a baseline parameter list and one row of fitted values (e.g. the best
#' row of a Sobol' sweep), overlays the fitted parameters, derives `beta_mf`
#' and `beta_mm` from their multipliers, and returns a parameter list ready for
#' [run_SIP()].
#'
#' @param p0 Baseline parameter list.
#' @param fit One-row data frame (or named vector) of fitted values. Diagnostic
#'   columns (`loglik`, `delta_ll`) and the `C_*` multipliers are dropped after
#'   use.
#' @return A parameter list for [run_SIP()].
#' @export
format_parameters <- function(p0, fit) {
  fit <- as.list(fit)
  p <- modifyList(as.list(p0), fit)
  p$beta_mf <- p$beta_ffy * p$C_betamf
  p$beta_mm <- p$beta_ffy * p$C_betamm
  p[c("C_betamf", "C_betamm", "loglik", "delta_ll")] <- NULL
  p
}

