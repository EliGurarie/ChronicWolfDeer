#' SIP model parameters
#'
#' Parameters of the sex-age-structured CWD meta-population model, passed to
#' [run_SIP()] as a named list (`p`). Predation parameters act through a
#' softmax allocation of each county's kill budget across the eight classes
#' (susceptible/infected x male/female x adult/yearling).
#'
#' @name parameters
#'
#' @section Wolf predation:
#' * `alpha`: wolf selectivity for infected deer relative to susceptible.
#'   `alpha > 1` means infected deer are taken disproportionately.
#' * `xi`: wolf preference for males relative to females. `xi > 1` means
#'   males are taken disproportionately.
#' * `eta`: wolf preference for yearlings relative to adults. `eta > 1` means
#'   yearlings are taken disproportionately.
#' * `omega`: kill budget, deer killed per wolf per year. Total kill in a
#'   county is `omega * n.wolf`, allocated across classes by `alpha`, `xi`,
#'   and `eta`.
#'
#' @section CWD survival:
#' * `theta_m`, `theta_f`: survival multipliers for infected adult males and
#'   females. Infected adult survival is `theta * phi_healthy`, where
#'   `phi_healthy` is the harvest/non-harvest/predation survival of a
#'   susceptible. Values below 1 mean CWD mortality is additive to all other
#'   sources. Yearlings have no theta multiplier.
#'
#' @section Transmission:
#' Transmission is frequency-dependent. Each beta is an absolute rate filling
#' cells of a 4x4 matrix (target x source; classes `m.a`, `f.a`, `m.y`,
#' `f.y`) that is multiplied by the infected proportion of the county
#' population. Per-class infection probability is `1 - exp(-rate)`.
#' * `beta_mm`: male-male transmission (`m.a`<->`m.a`, `m.a`<->`m.y`
#'   excluding `m.y`->`m.y`).
#' * `beta_mf`: male-female transmission (`f.a`<->`m.a`, `f.a`<->`m.y`).
#' * `beta_ffy`: baseline for all remaining contacts: `f.a`<->`f.a`,
#'   `f.a`<->`f.y`, `f.y`<->`f.y`, `m.y`<->`m.y`, `m.y`<->`f.y`.
#'
#' Female yearlings and adult males (`f.y`<->`m.a`) do not transmit to each
#' other
#'
#' @section Demography:
#' * `tau`: multiplier on per-capita recruitment. Recruitment is
#'   `(mean_wtd_fdr / 2) * tau`, so `tau = 1` recovers the raw fawn-to-doe
#'   ratio. Calibrated to match observed population growth.
#'
#' @section Dispersal:
#' * `kappa_NF_y`, `kappa_EF_y`: scale of the male-yearling dispersal kernel
#'   (km) in the Northern Forest and Eastern Forest regions, respectively.
#'
#' @seealso [run_SIP()]
NULL