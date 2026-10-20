#' Pre-assemble fixed model inputs for fast repeated fitting
#'
#' Builds everything that does not change across parameter sets once, so the
#' likelihood loop can call [run_SIP()] / [calc_logLik_SIP()] thousands of times
#' without rebuilding the dispersal kernel. Kappas are fixed here (not fitted).
#'
#' @param full_data Full county-year panel (carried through for the likelihood).
#' @param year1 Seed year.
#' @param kappas Named vector/list with `kappa_NF_y` and `kappa_EF_y`.
#' @param M_distance Square county-by-county distance matrix.
#' @return A list with `data_year1`, `dispersal_matrix`, `year1`, `counties`,
#'   `classes`, `full_data`.
#' @export
prep_data <- function(full_data, year1, kappas, M_distance) {
  
  if (stats::median(M_distance[M_distance > 0]) > 1000)
    warning("M_distance looks like metres; kappa is in km.")
  
  data_year1 <- subset(full_data, year == year1)

  dispersal_matrix <- calculate_prob_dispersal(
    as.character(data_year1$county), M_distance,
    kappas[["kappa_NF_y"]], kappas[["kappa_EF_y"]], data_year1$dnr_region
  )
  dispersal_matrix[dispersal_matrix < 0] <- 0
  dispersal_matrix <- dispersal_matrix / rowSums(dispersal_matrix)

  list(data_year1       = data_year1,
       dispersal_matrix = dispersal_matrix,
       year1            = year1,
       counties         = as.character(data_year1$county),
       classes          = sip_classes,
       full_data        = full_data)
}
