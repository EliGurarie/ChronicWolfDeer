# =============================================================================
# run_model: the SIP model run engine
#
#   run_SIP()           <- user-facing: simulate forward, return tidy data frame
#     run_SIP_matrix()    <- internal: iterate the annual step over n_years
#       step_SIP_matrix()   <- internal: advance one year (the actual mechanics)
#
# Fixed inputs (data_year1, dispersal_matrix, counties, classes) come prebuilt
# from prep_data(); the engine never rebuilds the dispersal kernel.
# =============================================================================


#' Iterate the SIP model over multiple years (internal engine)
#'
#' @param parameters Raw-scale parameter list.
#' @param n_years Number of years to simulate (>= 2); year 1 is the seed.
#' @param data_list Output of [prep_data()].
#' @param ... Forwarded to [step_SIP_matrix()].
#' @return A list of eight counties-by-`n_years` matrices plus `data`.
#' @keywords internal
run_SIP_matrix <- function(parameters, n_years, data_list, ...) {
  d        <- data_list$data_year1
  counties <- data_list$counties
  n.sites  <- length(counties)

  mk <- function() matrix(NA_real_, nrow = n.sites, ncol = n_years)
  S.M.A <- mk(); S.F.A <- mk(); S.M.Y <- mk(); S.F.Y <- mk()
  I.M.A <- mk(); I.F.A <- mk(); I.M.Y <- mk(); I.F.Y <- mk()

  S.M.A[, 1] <- d$S.M.A0; S.F.A[, 1] <- d$S.F.A0
  S.M.Y[, 1] <- d$S.M.Y0; S.F.Y[, 1] <- d$S.F.Y0
  I.M.A[, 1] <- d$I.M.A0; I.F.A[, 1] <- d$I.F.A0
  I.M.Y[, 1] <- d$I.M.Y0; I.F.Y[, 1] <- d$I.F.Y0

  for (j in 2:n_years) {
    nextstep <- step_SIP_matrix(
      s.m.a = S.M.A[, j - 1], i.m.a = I.M.A[, j - 1],
      s.f.a = S.F.A[, j - 1], i.f.a = I.F.A[, j - 1],
      s.m.y = S.M.Y[, j - 1], i.m.y = I.M.Y[, j - 1],
      s.f.y = S.F.Y[, j - 1], i.f.y = I.F.Y[, j - 1],
      mean_mu_male_harv = d$mean_mu_male_harv,
      mean_mu_fem_harv  = d$mean_mu_fem_harv,
      mean_wtd_fdr      = d$mean_wtd_fdr,
      nh_mu_adult_other    = d$nh_mu_adult_other,
      nh_mu_yearling_other = d$nh_mu_yearling_other,
      nh_mu_adult_human    = d$nh_mu_adult_human,
      nh_mu_yearling_human = d$nh_mu_yearling_human,
      phi_fty_male   = d$fty_survival_male,
      phi_fty_female = d$fty_survival_female,
      p                = parameters,
      dispersal_matrix = data_list$dispersal_matrix,
      W                = d$wolf_presence,
      n.wolf           = d$wolf_abundance,
      K                = d$K,
      counties         = counties,
      classes          = data_list$classes,
      ...
    )

    I.M.A[, j] <- nextstep$I.M.A; I.F.A[, j] <- nextstep$I.F.A
    I.M.Y[, j] <- nextstep$I.M.Y; I.F.Y[, j] <- nextstep$I.F.Y
    S.M.A[, j] <- nextstep$S.M.A; S.F.A[, j] <- nextstep$S.F.A
    S.M.Y[, j] <- nextstep$S.M.Y; S.F.Y[, j] <- nextstep$S.F.Y
  }

  list(S.M.A = S.M.A, S.F.A = S.F.A, I.M.A = I.M.A, I.F.A = I.F.A,
       S.M.Y = S.M.Y, S.F.Y = S.F.Y, I.M.Y = I.M.Y, I.F.Y = I.F.Y,
       data = d)
}


#' Advance the SIP model by one annual step (internal engine)
#'
#' Projects every county forward one year: (1) density-dependent recruitment and
#' survival via a per-county Lefkovitch projection, (2) frequency-dependent CWD
#' transmission, (3) dispersal of male yearlings between counties.
#'
#' @param s.m.a,s.f.a,s.m.y,s.f.y,i.m.a,i.f.a,i.m.y,i.f.y Per-county counts.
#' @param mean_mu_male_harv,mean_mu_fem_harv Harvest mortality rates.
#' @param mean_wtd_fdr Fawn-to-doe ratio (recruitment = `mean_wtd_fdr / 2`).
#' @param nh_mu_adult_other,nh_mu_yearling_other,nh_mu_adult_human,nh_mu_yearling_human
#'   Non-harvest mortality rates.
#' @param phi_fty_male,phi_fty_female Fawn-to-yearling survival probabilities.
#' @param K Carrying capacities. `W` Wolf presence. `n.wolf` Wolf abundance.
#' @param p Raw-scale parameter list (`alpha`, `xi`, `eta`, `omega`, `theta_m`,
#'   `theta_f`, `beta_mm`, `beta_mf`, `beta_ffy`).
#' @param dispersal_matrix Prebuilt row-stochastic dispersal matrix.
#' @param counties County names (order of all vector arguments).
#' @param classes The eight class labels; defaults to [sip_classes].
#' @return A named list of eight per-county vectors (`S.M.A` ... `I.F.Y`).
#' @keywords internal
step_SIP_matrix <- function(s.m.a, s.f.a, i.m.a, i.f.a,
                            s.m.y, s.f.y, i.m.y, i.f.y,
                            mean_mu_male_harv, mean_mu_fem_harv, mean_wtd_fdr,
                            nh_mu_adult_other, nh_mu_yearling_other,
                            nh_mu_adult_human, nh_mu_yearling_human,
                            phi_fty_male, phi_fty_female,
                            K, W, n.wolf, p,
                            dispersal_matrix, counties,
                            classes = sip_classes) {

  n.classes <- length(classes)
  n.sites   <- length(counties)

  is_inf  <- grepl("i", classes)
  is_male <- grepl("m", classes)
  is_yrlg <- grepl("y", classes)

  # Step 1. Population Dynamics ----
  PopMatrix0 <- cbind(s.m.a, s.f.a, s.m.y, s.f.y, i.m.a, i.f.a, i.m.y, i.f.y)
  dimnames(PopMatrix0) <- list(counties, classes)

  LeslieArray <- array(0, dim = c(n.classes, n.classes, n.sites),
                       dimnames = list(classes, classes, counties))

  rho <- (mean_wtd_fdr / 2)

  ## 1.3 Predation softmax weights (vectorized)
  pref <- exp(is_inf * log(p$alpha) + is_male * log(p$xi) + is_yrlg * log(p$eta))

  weights <- PopMatrix0 * matrix(pref, n.sites, n.classes, byrow = TRUE)
  prop_pred_class <- weights / (rowSums(weights) + 1e-12)

  max.tot  <- ifelse(as.logical(W), p$omega * n.wolf, 0)
  capped   <- pmin(max.tot, rowSums(PopMatrix0) * 0.75)

  prop.predated.byclass.all <- pmin((capped * prop_pred_class) / pmax(PopMatrix0, 1e-12), 0.75)
  prop.predated.byclass.all[!is.finite(prop.predated.byclass.all)] <- 0
  dimnames(prop.predated.byclass.all) <- list(counties, classes)

  ## 1.6 Per-class survival via multiplicative hazards (vectorized)
  male_yearling_prop   <- 0.48
  female_yearling_prop <- 0.29

  phi_matrix <- matrix(NA_real_, n.sites, n.classes,
                       dimnames = list(counties, classes))

  for (cl_idx in seq_len(n.classes)) {
    adult <- !is_yrlg[cl_idx]
    male  <-  is_male[cl_idx]
    inf   <-  is_inf[cl_idx]

    if (adult) {
      harvest_rate <- if (male) (1 - male_yearling_prop) * mean_mu_male_harv
      else      (1 - female_yearling_prop) * mean_mu_fem_harv
      nonharvest_rate <- nh_mu_adult_human + nh_mu_adult_other
    } else {
      harvest_rate <- if (male) male_yearling_prop * mean_mu_male_harv
      else      female_yearling_prop * mean_mu_fem_harv
      nonharvest_rate <- nh_mu_yearling_human + nh_mu_yearling_other
    }
    wolf_rate <- prop.predated.byclass.all[, cl_idx]

    harvest_rate    <- pmin(pmax(harvest_rate, 0), 1)
    nonharvest_rate <- pmin(pmax(nonharvest_rate, 0), 1)
    wolf_rate       <- pmin(pmax(wolf_rate, 0), 1)

    S_base <- (1 - harvest_rate) * (1 - nonharvest_rate) * (1 - wolf_rate)
    if (inf && adult) S_base <- (if (male) p$theta_m else p$theta_f) * S_base

    phi_matrix[, cl_idx] <- pmin(pmax(S_base, 0), 1)
  }

  phi_matrix[is.na(phi_matrix)] <- 0
  phi_matrix[is.infinite(phi_matrix)] <- 1

  ## 1.7 Fill the Lefkovitch entries
  LeslieArray["s.m.y", "s.f.a", ] <- pmin(pmax(rho * phi_fty_male, 0.00001), 0.99)
  LeslieArray["s.m.a", "s.m.y", ] <- phi_matrix[, "s.m.y"]
  LeslieArray["s.m.a", "s.m.a", ] <- phi_matrix[, "s.m.a"]
  LeslieArray["s.f.y", "s.f.a", ] <- pmin(pmax(rho * phi_fty_female, 0.00001), 0.99)
  LeslieArray["s.f.a", "s.f.y", ] <- phi_matrix[, "s.f.y"]
  LeslieArray["s.f.a", "s.f.a", ] <- phi_matrix[, "s.f.a"]

  LeslieArray["s.m.y", "i.f.a", ] <- pmin(pmax(rho * phi_fty_male, 0.00001), 0.99)
  LeslieArray["s.f.y", "i.f.a", ] <- pmin(pmax(rho * phi_fty_female, 0.00001), 0.99)
  LeslieArray["i.m.y", "i.f.a", ] <- 0
  LeslieArray["i.m.a", "i.m.y", ] <- phi_matrix[, "i.m.y"]
  LeslieArray["i.m.a", "i.m.a", ] <- phi_matrix[, "i.m.a"]
  LeslieArray["i.f.y", "i.f.a", ] <- 0
  LeslieArray["i.f.a", "i.f.y", ] <- phi_matrix[, "i.f.y"]
  LeslieArray["i.f.a", "i.f.a", ] <- phi_matrix[, "i.f.a"]

  ## 1.8 Beverton-Holt density correction, then project
  PopMatrix.postLeslie <- matrix(NA_real_, nrow = n.sites, ncol = n.classes,
                                 dimnames = list(counties, classes))

  for (c_idx in seq_len(n.sites)) {
    county <- counties[c_idx]
    pop0 <- as.numeric(PopMatrix0[county, classes])
    delta <- 1 / (1 + sum(pop0) / K[c_idx])

    Lc <- LeslieArray[, , county]
    Lc["s.m.y", "s.f.a"] <- Lc["s.m.y", "s.f.a"] * delta
    Lc["s.f.y", "s.f.a"] <- Lc["s.f.y", "s.f.a"] * delta
    Lc["i.m.y", "i.f.a"] <- Lc["i.m.y", "i.f.a"] * delta
    Lc["i.f.y", "i.f.a"] <- Lc["i.f.y", "i.f.a"] * delta

    PopMatrix.postLeslie[county, ] <- as.vector(Lc %*% pop0)
  }

  PopMatrix.postLeslie <- pmax(PopMatrix.postLeslie, 0)
  rownames(PopMatrix.postLeslie) <- counties
  colnames(PopMatrix.postLeslie) <- classes

  # Step 2. Infection ----
  sa_classes <- c("m.a", "f.a", "m.y", "f.y")
  TransmissionMatrix <- matrix(0, length(sa_classes), length(sa_classes),
                               dimnames = list(sa_classes, sa_classes))
  TransmissionMatrix["m.a", "m.a"] <- p$beta_mm
  TransmissionMatrix["m.a", "m.y"] <- p$beta_mm
  TransmissionMatrix["m.y", "m.a"] <- p$beta_mm
  TransmissionMatrix["f.a", "m.a"] <- p$beta_mf
  TransmissionMatrix["m.a", "f.a"] <- p$beta_mf
  TransmissionMatrix["f.a", "m.y"] <- p$beta_mf
  TransmissionMatrix["m.y", "f.a"] <- p$beta_mf
  TransmissionMatrix["f.a", "f.a"] <- p$beta_ffy
  TransmissionMatrix["f.y", "f.a"] <- p$beta_ffy
  TransmissionMatrix["f.a", "f.y"] <- p$beta_ffy
  TransmissionMatrix["m.y", "m.y"] <- p$beta_ffy
  TransmissionMatrix["f.y", "m.y"] <- p$beta_ffy
  TransmissionMatrix["m.y", "f.y"] <- p$beta_ffy
  TransmissionMatrix["f.y", "f.y"] <- p$beta_ffy

  susc <- c("s.m.a", "s.f.a", "s.m.y", "s.f.y")
  infs <- c("i.m.a", "i.f.a", "i.m.y", "i.f.y")

  InfectionShuffle <- array(0, dim = c(n.classes, n.classes, n.sites),
                            dimnames = list(classes, classes, counties))

  for (c in counties) {
    pop_vec <- PopMatrix.postLeslie[c, ]
    old_infected <- sum(pop_vec[infs])
    old_pop <- old_infected + sum(pop_vec[susc])

    if (old_pop > 0) {
      I_prop <- pop_vec[infs] / old_pop
      P.infected <- as.vector(TransmissionMatrix %*% I_prop)
    } else {
      P.infected <- rep(0, length(susc))
    }
    if (any(P.infected < 0 | P.infected > 1)) {
      warning("Invalid P.infected values; clamping applied")
      P.infected <- pmin(pmax(P.infected, 0), 1)
    }

    S_vec <- pop_vec[susc]
    new_infected <- pmin(P.infected * S_vec, S_vec)

    allowed_new_infections <- max(0, 0.99 * old_pop - old_infected)
    if (sum(new_infected) > allowed_new_infections && sum(new_infected) > 0) {
      new_infected <- new_infected * (allowed_new_infections / sum(new_infected))
    }

    InfectionShuffle[, , c] <- 0
    InfectionShuffle[infs, infs, c] <- diag(1, length(infs))
    for (k in seq_along(susc)) {
      p_inf <- if (S_vec[k] > 0) pmin(pmax(new_infected[k] / S_vec[k], 0), 1) else 0
      InfectionShuffle[susc[k], infs[k], c] <- p_inf
      InfectionShuffle[susc[k], susc[k], c] <- 1 - p_inf
    }
  }

  PopMatrix.postInfection <- sapply(counties, function(c) {
    PopMatrix.postLeslie[c, ] %*% InfectionShuffle[, , c]
  }) |> t()
  rownames(PopMatrix.postInfection) <- counties
  colnames(PopMatrix.postInfection) <- classes
  PopMatrix.postInfection <- pmax(PopMatrix.postInfection, 0)

  for (c in counties) {
    pop_c <- PopMatrix.postInfection[c, ]
    N_c <- sum(pop_c); I_c <- sum(pop_c[infs])
    if (N_c > 0 && I_c / N_c > 0.50) {
      scale <- (I_c - 0.50 * N_c) / I_c
      for (k in seq_along(infs)) {
        move <- scale * pop_c[infs[k]]
        PopMatrix.postInfection[c, infs[k]] <- pop_c[infs[k]] - move
        PopMatrix.postInfection[c, susc[k]] <- pop_c[susc[k]] + move
      }
    }
  }

  # Step 3. Dispersal (male yearlings only) ----
  PopMatrix.postDispersal <- PopMatrix.postInfection
  PopMatrix.postDispersal[, "s.m.y"] <-
    as.numeric(PopMatrix.postInfection[, "s.m.y"] %*% dispersal_matrix)
  PopMatrix.postDispersal[, "i.m.y"] <-
    as.numeric(PopMatrix.postInfection[, "i.m.y"] %*% dispersal_matrix)
  PopMatrix.postDispersal <- pmax(PopMatrix.postDispersal, 0)

  list(
    S.M.A = PopMatrix.postDispersal[, "s.m.a"],
    S.F.A = PopMatrix.postDispersal[, "s.f.a"],
    I.M.A = PopMatrix.postDispersal[, "i.m.a"],
    I.F.A = PopMatrix.postDispersal[, "i.f.a"],
    S.M.Y = PopMatrix.postDispersal[, "s.m.y"],
    S.F.Y = PopMatrix.postDispersal[, "s.f.y"],
    I.M.Y = PopMatrix.postDispersal[, "i.m.y"],
    I.F.Y = PopMatrix.postDispersal[, "i.f.y"]
  )
}


#' Run the Wolf-CWD SIP model
#'
#' Simulates deer demography and CWD dynamics for `n_years`, starting from the
#' seed year in `data_list` (from [prep_data()]). Returns one tidy row per
#' county per simulated year.
#'
#' @param parameters Raw-scale parameter list. Required: `alpha`, `xi`, `eta`,
#'   `omega`, `theta_m`, `theta_f`, `beta_mm`, `beta_mf`, `beta_ffy`.
#' @param n_years Number of years to simulate (>= 2). Year 1 is the seed.
#' @param data_list Output of [prep_data()].
#' @param ... Advanced: forwarded to the engine.
#' @return A data frame: eight class counts (`SMA_sim` ... `IFY_sim`),
#'   `I_all_sim`, `S_all_sim`, `county`, `sim_year`, `year`. `NULL` on failure.
#' @export
#' @example examples/example_run_SIP.R
run_SIP <- function(parameters, n_years, data_list, ...) {
  if (!is.numeric(n_years) || n_years < 2)
    stop("run_SIP: n_years must be an integer >= 2.")

  myfit <- tryCatch(
    run_SIP_matrix(parameters = parameters, n_years = n_years,
                   data_list = data_list, ...),
    error = function(e) {
      warning("run_SIP: run_SIP_matrix failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(myfit)) return(NULL)

  result_list <- vector("list", length = n_years)
  for (t in 1:n_years) {
    SMA <- myfit$S.M.A[, t]; SFA <- myfit$S.F.A[, t]
    SMY <- myfit$S.M.Y[, t]; SFY <- myfit$S.F.Y[, t]
    IMA <- myfit$I.M.A[, t]; IFA <- myfit$I.F.A[, t]
    IMY <- myfit$I.M.Y[, t]; IFY <- myfit$I.F.Y[, t]

    result_list[[t]] <- data.frame(
      SMA_sim = SMA, SFA_sim = SFA, SMY_sim = SMY, SFY_sim = SFY,
      IMA_sim = IMA, IFA_sim = IFA, IMY_sim = IMY, IFY_sim = IFY,
      I_all_sim = IMA + IFA + IMY + IFY,
      S_all_sim = SMA + SFA + SMY + SFY,
      county   = data_list$data_year1$county,
      sim_year = t,
      year     = data_list$year1 + t - 1,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, result_list)
}
