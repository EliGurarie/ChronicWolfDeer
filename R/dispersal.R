#' Build the inter-county dispersal Markov matrix
#'
#' Constructs the row-stochastic county-to-county transition matrix used to
#' redistribute dispersing male yearlings. Dispersal follows an exponential
#' kernel with a region-specific mean distance \eqn{\kappa}: the probability of
#' moving from county \eqn{i} to county \eqn{j} is
#' \eqn{\exp(-d_{ij} / \kappa_i)} for \eqn{d_{ij} > 0}. The diagonal (staying
#' put) is set to `1 - rowSums(off-diagonal)`, floored at 0, and each row is
#' then normalised to sum to 1.
#'
#' @param counties Character vector of county names, in the same order as the
#'   rows/columns of `distances`.
#' @param distances Numeric county-by-county distance matrix (e.g. centroid
#'   distances in km). Must be square with `length(counties)` rows.
#' @param kappa_north,kappa_south Mean dispersal distances (kernel scale) for
#'   northern and southern regions, respectively.
#' @param region_vector Character vector (length `length(counties)`) of region
#'   labels. A county is treated as northern when its label contains the
#'   substring `"North"` (case-sensitive); otherwise it is southern.
#'
#' @return A row-stochastic numeric matrix with `length(counties)` rows and
#'   columns, row names set to `counties`.
#'
#' @section Note:
#' Region assignment keys on the literal substring `"North"`. Confirm your
#' `dnr_region` labels match this convention (e.g. a label like
#' `"West Central Region"` is treated as southern).
#'
#' @keywords internal
#' @examples
#' counties <- c("Aville", "Bborough", "Ctown")
#' d <- matrix(c(0, 10, 20, 10, 0, 15, 20, 15, 0), 3, 3)
#' regions <- c("Northern Forest", "Eastern Farmland", "Eastern Farmland")
#' ChronicWolfDeerFunction:::calculate_prob_dispersal(
#'   counties, d, kappa_north = 20, kappa_south = 15, region_vector = regions
#' )
#' @export
calculate_prob_dispersal <- function(counties, distances,
                                     kappa_north, kappa_south,
                                     region_vector) {
  prob.dispersal <- sapply(seq_along(counties), function(j) {
    # smaller kappa = more local dispersal; larger = flatter tail
    kappa_value <- ifelse(grepl("North", region_vector[j]), kappa_north, kappa_south)
    ifelse(distances[, j] > 0, exp(-distances[, j] / kappa_value), 0)
  }) |> t()

  rownames(prob.dispersal) <- counties

  diag(prob.dispersal) <- pmax(1 - rowSums(prob.dispersal), 0)
  prob.dispersal <- apply(prob.dispersal, 1, function(x) x / sum(x)) |> t()
  prob.dispersal
}
