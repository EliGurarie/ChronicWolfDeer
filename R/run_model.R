# =============================================================================
# run_model: the SIP model run engine
#
#   runModel_sim()    <- user-facing: simulate forward, return tidy data frame
#     run_SIP_matrix()    <- internal: iterate the annual step over Tmax years
#       step_SIP_matrix()   <- internal: advance one year (the actual mechanics)
#
# Only runModel_sim() is exported. The two engine functions are internal.
# Plumbing arguments (M_distance, counties, classes) are forwarded down the
# stack via `...`; the per-county data columns are unpacked from init_data
# inside run_SIP_matrix().
# =============================================================================


#' Iterate the SIP model over multiple years (internal engine)
#'
#' Seeds the eight class matrices from `init_data` (year 1) and advances them
#' with [step_SIP_matrix()] for `Tmax` years. County order is taken from
#' `init_data$county`. Per-county data columns are unpacked here and handed to
#' the annual step.
#'
#' The dispersal kernel is constant across years, so it is built once (via
#' [calculate_prob_dispersal()]) and reused every step rather than rebuilt.
#'
#' @param init_data One row per county, with the seed columns (`S.M.A0` ...
#'   `I.F.Y0`), rate columns, and context columns (`county`, `dnr_region`,
#'   `wolf_presence`, `wolf_abundance`, `K`).
#' @param parameters Raw-scale parameter list (passed as `p`; must include
#'   `kappa_NF_y` and `kappa_EF_y` for the dispersal kernel).
#' @param Tmax Number of years to simulate (>= 2); year 1 is the seed.
#' @param M_distance Square county-by-county distance matrix, ordered to match
#'   `init_data$county`.
#' @param ... Forwarded to [step_SIP_matrix()] (e.g. `classes`).
#'
#' @return A list of eight counties-by-`Tmax` matrices plus `data` (the input).
#' @keywords internal
run_SIP_matrix <- function(init_data, parameters, Tmax, M_distance, ...) {
  if (!is.data.frame(init_data) || nrow(init_data) == 0)
    stop("run_SIP_matrix: init_data must be a non-empty data frame.")
  if (!is.numeric(Tmax) || Tmax < 2)
    stop("run_SIP_matrix: Tmax must be an integer >= 2.")
  
  counties <- as.character(init_data$county)
  n.sites  <- length(counties)
  
  # Dispersal kernel: constant across years, build once.
  dispersal_matrix <- calculate_prob_dispersal(
    counties, M_distance, parameters$kappa_NF_y, parameters$kappa_EF_y,
    init_data$dnr_region
  )
  dispersal_matrix[dispersal_matrix < 0] <- 0
  dispersal_matrix <- dispersal_matrix / rowSums(dispersal_matrix)
  
  mk <- function() matrix(NA_real_, nrow = n.sites, ncol = Tmax)
  S.M.A <- mk(); S.F.A <- mk(); S.M.Y <- mk(); S.F.Y <- mk()
  I.M.A <- mk(); I.F.A <- mk(); I.M.Y <- mk(); I.F.Y <- mk()
  
  # seed year 1
  S.M.A[, 1] <- init_data$S.M.A0; S.F.A[, 1] <- init_data$S.F.A0
  S.M.Y[, 1] <- init_data$S.M.Y0; S.F.Y[, 1] <- init_data$S.F.Y0
  I.M.A[, 1] <- init_data$I.M.A0; I.F.A[, 1] <- init_data$I.F.A0
  I.M.Y[, 1] <- init_data$I.M.Y0; I.F.Y[, 1] <- init_data$I.F.Y0
  
  for (j in 2:Tmax) {
    nextstep <- step_SIP_matrix(
      s.m.a = S.M.A[, j - 1], i.m.a = I.M.A[, j - 1],
      s.f.a = S.F.A[, j - 1], i.f.a = I.F.A[, j - 1],
      s.m.y = S.M.Y[, j - 1], i.m.y = I.M.Y[, j - 1],
      s.f.y = S.F.Y[, j - 1], i.f.y = I.F.Y[, j - 1],
      mean_mu_male_harv = init_data$mean_mu_male_harv,
      mean_mu_fem_harv  = init_data$mean_mu_fem_harv,
      mean_wtd_fdr      = init_data$mean_wtd_fdr,
      nh_mu_adult_other    = init_data$nh_mu_adult_other,
      nh_mu_yearling_other = init_data$nh_mu_yearling_other,
      nh_mu_adult_human    = init_data$nh_mu_adult_human,
      nh_mu_yearling_human = init_data$nh_mu_yearling_human,
      phi_fty_male   = init_data$fty_survival_male,
      phi_fty_female = init_data$fty_survival_female,
      p                = parameters,
      dispersal_matrix = dispersal_matrix,
      W                = init_data$wolf_presence,
      n.wolf           = init_data$wolf_abundance,
      K                = init_data$K,
      counties         = counties,
      ...
    )
    
    I.M.A[, j] <- nextstep$I.M.A; I.F.A[, j] <- nextstep$I.F.A
    I.M.Y[, j] <- nextstep$I.M.Y; I.F.Y[, j] <- nextstep$I.F.Y
    S.M.A[, j] <- nextstep$S.M.A; S.F.A[, j] <- nextstep$S.F.A
    S.M.Y[, j] <- nextstep$S.M.Y; S.F.Y[, j] <- nextstep$S.F.Y
  }
  
  list(S.M.A = S.M.A, S.F.A = S.F.A, I.M.A = I.M.A, I.F.A = I.F.A,
       S.M.Y = S.M.Y, S.F.Y = S.F.Y, I.M.Y = I.M.Y, I.F.Y = I.F.Y,
       data = init_data)
}



#' Advance the SIP model by one annual step (internal engine)
#'
#' Projects every county forward one year: (1) density-dependent recruitment and
#' survival via a per-county Lefkovitch projection, (2) frequency-dependent CWD
#' transmission, (3) dispersal of male yearlings between counties. All vector
#' arguments are indexed by county and share the order given in `counties`.
#'
#' The dispersal kernel is constant across years, so it is built once by the
#' caller ([run_SIP_matrix()]) and passed in via `dispersal_matrix` rather than
#' rebuilt every step.
#'
#' @param s.m.a,s.f.a,s.m.y,s.f.y,i.m.a,i.f.a,i.m.y,i.f.y Per-county counts in
#'   each of the eight susceptible/infected x sex x age classes.
#' @param mean_mu_male_harv,mean_mu_fem_harv Harvest mortality rates.
#' @param mean_wtd_fdr Fawn-to-doe ratio (recruitment = `mean_wtd_fdr / 2`).
#' @param nh_mu_adult_other,nh_mu_yearling_other,nh_mu_adult_human,nh_mu_yearling_human
#'   Non-harvest mortality rates (other/natural and anthropogenic).
#' @param phi_fty_male,phi_fty_female Fawn-to-yearling survival probabilities.
#' @param K Carrying capacities (Beverton-Holt density correction).
#' @param W Wolf-presence flags. `n.wolf` Wolf abundances.
#' @param p Raw-scale parameter list (`alpha`, `xi`, `eta`, `omega`, `theta_m`,
#'   `theta_f`, `beta_mm`, `beta_mf`, `beta_ffy`).
#' @param dispersal_matrix Row-stochastic county-by-county dispersal matrix from
#'   [calculate_prob_dispersal()], ordered to match `counties`.
#' @param counties County names giving the order of all vector arguments.
#' @param classes The eight class labels; defaults to [sip_classes].
#'
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
  
  # Fixed class properties (constant for the canonical class ordering).
  is_inf  <- grepl("i", classes)
  is_male <- grepl("m", classes)
  is_yrlg <- grepl("y", classes)
  
  # Step 1. Population Dynamics ----
  
  ## 1.1 Starting population matrix
  PopMatrix0 <- cbind(s.m.a, s.f.a, s.m.y, s.f.y, i.m.a, i.f.a, i.m.y, i.f.y)
  dimnames(PopMatrix0) <- list(counties, classes)
  
  ## 1.2 Leslie/Lefkovitch array (one slice per county)
  LeslieArray <- array(0, dim = c(n.classes, n.classes, n.sites),
                       dimnames = list(classes, classes, counties))
  
  rho <- (mean_wtd_fdr / 2)
  
  ## 1.3 Predation softmax weights (vectorized over classes and counties) ----
  # Per-class utility is fixed; weight = exp(utility) * abundance per county.
  class_utils <- is_inf * log(p$alpha) + is_male * log(p$xi) + is_yrlg * log(p$eta)
  pref <- exp(class_utils)                                  # length n.classes
  
  weights <- PopMatrix0 * matrix(pref, n.sites, n.classes, byrow = TRUE)
  row_w   <- rowSums(weights) + 1e-12
  prop_pred_class <- weights / row_w                         # n.sites x n.classes
  
  # total kill budget per county (0 where no wolves), capped at 75% of pop
  pop_tot  <- rowSums(PopMatrix0)
  max.tot  <- ifelse(as.logical(W), p$omega * n.wolf, 0)
  capped   <- pmin(max.tot, pop_tot * 0.75)
  
  pop_safe <- pmax(PopMatrix0, 1e-12)
  prop.predated.byclass.all <- pmin((capped * prop_pred_class) / pop_safe, 0.75)
  prop.predated.byclass.all[!is.finite(prop.predated.byclass.all)] <- 0
  dimnames(prop.predated.byclass.all) <- list(counties, classes)
  
  ## 1.6 Per-class survival via multiplicative hazards (vectorized) ----
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
  
  ## 1.7 Fill the Lefkovitch entries per county
  # Recruitment (adult female -> yearling), survival (yearling -> adult, adult -> adult)
  LeslieArray["s.m.y", "s.f.a", ] <- pmin(pmax(rho * phi_fty_male, 0.00001), 0.99)
  LeslieArray["s.m.a", "s.m.y", ] <- phi_matrix[, "s.m.y"]
  LeslieArray["s.m.a", "s.m.a", ] <- phi_matrix[, "s.m.a"]
  LeslieArray["s.f.y", "s.f.a", ] <- pmin(pmax(rho * phi_fty_female, 0.00001), 0.99)
  LeslieArray["s.f.a", "s.f.y", ] <- phi_matrix[, "s.f.y"]
  LeslieArray["s.f.a", "s.f.a", ] <- phi_matrix[, "s.f.a"]
  
  # Infected mothers bear susceptible fawns (no vertical transmission)
  LeslieArray["s.m.y", "i.f.a", ] <- pmin(pmax(rho * phi_fty_male, 0.00001), 0.99)
  LeslieArray["s.f.y", "i.f.a", ] <- pmin(pmax(rho * phi_fty_female, 0.00001), 0.99)
  LeslieArray["i.m.y", "i.f.a", ] <- 0
  LeslieArray["i.m.a", "i.m.y", ] <- phi_matrix[, "i.m.y"]
  LeslieArray["i.m.a", "i.m.a", ] <- phi_matrix[, "i.m.a"]
  LeslieArray["i.f.y", "i.f.a", ] <- 0
  LeslieArray["i.f.a", "i.f.y", ] <- phi_matrix[, "i.f.y"]
  LeslieArray["i.f.a", "i.f.a", ] <- phi_matrix[, "i.f.a"]
  
  ## 1.8 Beverton-Holt density correction on recruitment, then project
  PopMatrix.postLeslie <- matrix(NA_real_, nrow = n.sites, ncol = n.classes,
                                 dimnames = list(counties, classes))
  
  for (c_idx in seq_len(n.sites)) {
    county <- counties[c_idx]
    pop0 <- as.numeric(PopMatrix0[county, classes])
    N0 <- sum(pop0)
    
    delta <- 1 / (1 + N0 / K[c_idx])  # Beverton-Holt: smooth, never negative
    
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
  
  ## 2.1 Inter-class transmission matrix (rows = target, cols = source)
  sa_classes <- c("m.a", "f.a", "m.y", "f.y")
  TransmissionMatrix <- matrix(0, length(sa_classes), length(sa_classes),
                               dimnames = list(sa_classes, sa_classes))
  
  # male-male
  TransmissionMatrix["m.a", "m.a"] <- p$beta_mm
  TransmissionMatrix["m.a", "m.y"] <- p$beta_mm
  TransmissionMatrix["m.y", "m.a"] <- p$beta_mm
  # male-female (adult females with any males)
  TransmissionMatrix["f.a", "m.a"] <- p$beta_mf
  TransmissionMatrix["m.a", "f.a"] <- p$beta_mf
  TransmissionMatrix["f.a", "m.y"] <- p$beta_mf
  TransmissionMatrix["m.y", "f.a"] <- p$beta_mf
  # female / yearling background
  TransmissionMatrix["f.a", "f.a"] <- p$beta_ffy
  TransmissionMatrix["f.y", "f.a"] <- p$beta_ffy
  TransmissionMatrix["f.a", "f.y"] <- p$beta_ffy
  TransmissionMatrix["m.y", "m.y"] <- p$beta_ffy
  TransmissionMatrix["f.y", "m.y"] <- p$beta_ffy
  TransmissionMatrix["m.y", "f.y"] <- p$beta_ffy
  TransmissionMatrix["f.y", "f.y"] <- p$beta_ffy
  
  ## 2.2 Susceptible / infected class names
  susc <- c("s.m.a", "s.f.a", "s.m.y", "s.f.y")
  infs <- c("i.m.a", "i.f.a", "i.m.y", "i.f.y")
  
  ## 2.3 Per-county susceptible -> infected transition array
  InfectionShuffle <- array(0, dim = c(n.classes, n.classes, n.sites),
                            dimnames = list(classes, classes, counties))
  
  ## 2.4 Frequency-dependent force of infection per county
  for (c in counties) {
    pop_vec <- PopMatrix.postLeslie[c, ]
    
    old_infected    <- sum(pop_vec[infs])
    old_susceptible <- sum(pop_vec[susc])
    old_pop <- old_infected + old_susceptible
    
    if (old_pop > 0) {
      I_vec <- c("m.a" = pop_vec["i.m.a"], "f.a" = pop_vec["i.f.a"],
                 "m.y" = pop_vec["i.m.y"], "f.y" = pop_vec["i.f.y"])
      I_prop <- I_vec / old_pop
      P.infected <- as.vector(TransmissionMatrix %*% I_prop)  # ordered like sa_classes
    } else {
      P.infected <- rep(0, length(susc))
    }
    
    if (any(P.infected < 0 | P.infected > 1)) {
      warning("Invalid P.infected values; clamping applied")
      P.infected <- pmin(pmax(P.infected, 0), 1)
    }
    
    S_vec <- pop_vec[susc]
    new_infected <- pmin(P.infected * S_vec, S_vec)
    
    # keep prevalence <= 99% within the step
    allowed_new_infections <- max(0, 0.99 * old_pop - old_infected)
    if (sum(new_infected) > allowed_new_infections && sum(new_infected) > 0) {
      new_infected <- new_infected * (allowed_new_infections / sum(new_infected))
    }
    
    InfectionShuffle[, , c] <- 0
    InfectionShuffle[infs, infs, c] <- diag(1, length(infs))
    
    s_to_i <- setNames(infs, susc)
    for (k in seq_along(susc)) {
      s_cls <- susc[k]
      i_cls <- s_to_i[[s_cls]]
      S_count <- S_vec[k]
      p_inf <- if (S_count > 0) new_infected[k] / S_count else 0
      p_inf <- pmin(pmax(p_inf, 0), 1)
      InfectionShuffle[s_cls, i_cls, c] <- p_inf
      InfectionShuffle[s_cls, s_cls, c] <- 1 - p_inf
    }
  }
  
  ## 2.5 Apply the transition array (compartment-preserving)
  PopMatrix.postInfection <- sapply(counties, function(c) {
    PopMatrix.postLeslie[c, ] %*% InfectionShuffle[, , c]
  }) |> t()
  rownames(PopMatrix.postInfection) <- counties
  colnames(PopMatrix.postInfection) <- classes
  PopMatrix.postInfection <- pmax(PopMatrix.postInfection, 0)
  
  # post-step prevalence cap at 50%
  for (c in counties) {
    pop_c <- PopMatrix.postInfection[c, ]
    N_c <- sum(pop_c)
    I_c <- sum(pop_c[infs])
    if (N_c > 0 && I_c / N_c > 0.50) {
      excess <- I_c - 0.50 * N_c
      scale  <- excess / I_c
      for (k in seq_along(infs)) {
        move <- scale * pop_c[infs[k]]
        PopMatrix.postInfection[c, infs[k]] <- pop_c[infs[k]] - move
        PopMatrix.postInfection[c, susc[k]] <- pop_c[susc[k]] + move
      }
    }
  }
  
  # Step 3. Dispersal (male yearlings only) ----
  # dispersal_matrix is precomputed and constant across years (see run_SIP_matrix).
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
#' Simulates white-tailed deer demography and chronic wasting disease dynamics
#' across Wisconsin counties for `Tmx` years, starting from `initial_data`.
#' Returns one tidy row per county per simulated year.
#'
#' This is the user-facing entry point for model runs. It drives an internal
#' engine (`run_SIP_matrix()` -> `step_SIP_matrix()`) that applies, each year,
#' density-dependent recruitment and survival, frequency-dependent transmission,
#' and male-yearling dispersal between counties.
#'
#' @param parameters Raw-scale parameter list. Required entries: `alpha`, `xi`,
#'   `eta` (wolf preferences), `omega` (deer per wolf per year), `theta_m`,
#'   `theta_f` (infected-adult survival multipliers), `beta_mm`, `beta_mf`,
#'   `beta_ffy` (transmission rates), `kappa_NF_y`, `kappa_EF_y` (dispersal
#'   kernel scales for North / South).
#' @param Tmx Number of years to simulate (>= 2). Year 1 holds the seed values.
#' @param initial_data One row per county with the seed counts (`S.M.A0` ...
#'   `I.F.Y0`), mortality/recruitment rates, and context columns (`county`,
#'   `dnr_region`, `wolf_presence`, `wolf_abundance`, `K`). The `full_data`
#'   dataset (subset to one year) has this layout.
#' @param M_distance Square county-by-county distance matrix, ordered to match
#'   `initial_data$county`. Build it from `county_borders` centroids (see
#'   examples).
#' @param init_year Calendar year of the seed row; output `year` is
#'   `init_year + sim_year`. Default 2014.
#' @param ... Advanced: forwarded to the engine (e.g. `classes`).
#'
#' @return A data frame with one row per county-year: the eight simulated class
#'   counts (`SMA_sim` ... `IFY_sim`), `I_all_sim`, `S_all_sim`, `county`,
#'   `sim_year`, and `year`. `NULL` if the simulation failed.
#'
#' @export
#' @example man/examples/example_runModel_sim.R
runModel_sim <- function(parameters, Tmx, initial_data, M_distance,
                         init_year = 2014, ...) {
  if (!is.data.frame(initial_data) || nrow(initial_data) == 0)
    stop("runModel_sim: initial_data must be provided and non-empty.")
  if (!is.numeric(Tmx) || Tmx < 1)
    stop("runModel_sim: Tmx must be a positive integer.")
  
  myfit <- tryCatch(
    run_SIP_matrix(init_data = initial_data, parameters = parameters,
                   Tmax = Tmx, M_distance = M_distance, ...),
    error = function(e) {
      warning("runModel_sim: run_SIP_matrix failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(myfit)) return(NULL)
  
  result_list <- vector("list", length = Tmx)
  for (t in seq_len(Tmx)) {
    SMA <- myfit$S.M.A[, t]; SFA <- myfit$S.F.A[, t]
    SMY <- myfit$S.M.Y[, t]; SFY <- myfit$S.F.Y[, t]
    IMA <- myfit$I.M.A[, t]; IFA <- myfit$I.F.A[, t]
    IMY <- myfit$I.M.Y[, t]; IFY <- myfit$I.F.Y[, t]
    
    result_list[[t]] <- data.frame(
      SMA_sim = SMA, SFA_sim = SFA, SMY_sim = SMY, SFY_sim = SFY,
      IMA_sim = IMA, IFA_sim = IFA, IMY_sim = IMY, IFY_sim = IFY,
      I_all_sim = IMA + IFA + IMY + IFY,
      S_all_sim = SMA + SFA + SMY + SFY,
      county   = initial_data$county,
      sim_year = t,
      year     = init_year + t,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, result_list)
}
