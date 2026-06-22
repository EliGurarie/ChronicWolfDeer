#' SIP sex-age-infection class labels
#'
#' The fixed ordering of the eight model classes used throughout the package:
#' susceptible/infected (`s`/`i`) x male/female (`m`/`f`) x adult/yearling
#' (`a`/`y`). This ordering defines the column order of every population matrix
#' and the row/column order of the per-county Lefkovitch and infection arrays,
#' so it must not be reordered.
#'
#' @details
#' Kept as an internal package object. It can later be promoted to a documented
#' dataset via `usethis::use_data(sip_classes, internal = TRUE)` without
#' changing any calling code, because the engine takes it through the
#' `classes` argument (defaulted to this object).
#'
#' @format A length-8 character vector.
#' @keywords internal
"sip_classes"

# Defined here so the package works before any sysdata.rda is built.
# (Promote to data/ or sysdata.rda later if desired.)
sip_classes <- c("s.m.a", "s.f.a", "s.m.y", "s.f.y",
                 "i.m.a", "i.f.a", "i.m.y", "i.f.y")
