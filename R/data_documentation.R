#' Wisconsin white-tailed deer and CWD surveillance data
#'
#' County-year panel data for Wisconsin white-tailed deer populations,
#' chronic wasting disease (CWD) surveillance, harvest, mortality, wolf
#' presence, and derived model inputs. One row per county per year,
#' spanning 2014--2022 across 71 Wisconsin counties.
#'
#' @format A data frame with 639 rows and 33 variables:
#'
#' **Identifiers**
#' - `county` — County name.
#' - `year` — Calendar year.
#'
#' **Observed deer abundance and CWD**
#' - `wtd_count` — Estimated total white-tailed deer population (Sex-Age-Kill
#'   model; Storm, WDNR).
#' - `prop_cwd` — Apparent CWD prevalence: positive tests / samples analyzed.
#'   Counties with no samples are set to zero (WDNR public CWD testing data).
#' - `n_analyzed` — Number of deer tested for CWD.
#' - `n_cwd_positive` — Number of CWD-positive deer detected.
#'
#' **Initial compartment counts** (derived from `wtd_count` and `prop_cwd`)
#' - `S.M.A0`, `S.F.A0`, `S.M.Y0`, `S.F.Y0` — Susceptible adult male/female
#'   and yearling male/female counts.
#' - `I.M.A0`, `I.F.A0`, `I.M.Y0`, `I.F.Y0` — Infected adult male/female and
#'   yearling male/female counts.
#'
#' **Spatial and regional context**
#' - `dnr_region` — Wisconsin DNR administrative region (drives mortality and
#'   dispersal region assignments).
#' - `K` — County carrying capacity (deer), derived from region-specific
#'   maximum densities (Zouwen & Warnke 1995) multiplied by county area; set
#'   to observed population if population exceeds the literature estimate.
#' - `area` — County area (km²).
#' - `shape_area` — County polygon area from shapefile (m²).
#' - `x.km`, `y.km` — County centroid easting and northing (km).
#' - `Z` — County centroid as a complex number (`x.km + y.km * 1i`); used for
#'   inter-county distance calculations.
#'
#' **Wolf covariates**
#' - `wolf_presence` — Logical; `TRUE` if wolves are present in the county
#'   (WDNR wolf density maps).
#' - `wolf_abundance` — Estimated wolf abundance derived from WDNR wolf
#'   density maps and county area (see appendix).
#' - `wolf_density_rank` — Categorical wolf density: `"low"`, `"medium"`,
#'   or `"high"`.
#'
#' **Harvest mortality** (averaged 2014--2022; WDNR harvest data, Storm)
#' - `mean_mu_male_harv` — Mean annual harvest mortality rate for males
#'   (antlered harvest). Yearling bucks comprise ~50% of harvested males,
#'   allocated equally between `MA` and `MY`.
#' - `mean_mu_fem_harv` — Mean annual harvest mortality rate for females
#'   (antlerless harvest). Yearlings comprise ~30% of harvested females:
#'   `mu_FA = 2/3 * TH_antlerless`, `mu_FY = 1/3 * TH_antlerless`.
#'
#' **Non-harvest mortality** (Norton et al. 2021, tab. 3; age- and
#' region-specific; does not include CWD or predation)
#' - `nh_mu_adult_human` — Human-associated adult mortality (roadkill,
#'   poaching).
#' - `nh_mu_yearling_human` — Human-associated yearling mortality.
#' - `nh_mu_adult_other` — Other natural adult mortality (starvation, disease
#'   excluding CWD).
#' - `nh_mu_yearling_other` — Other natural yearling mortality.
#'
#' **Recruitment inputs**
#' - `fty_survival_male`, `fty_survival_female` — Fawn-to-yearling survival
#'   probabilities by sex and region, approximated from daily hazard data
#'   (Warbington et al. 2017, fig. 4) over 110 days using
#'   \eqn{\phi_t = \prod_{d=1}^{t}(1 - h_d)}.
#' - `mean_wtd_fdr` — Mean fawn-to-doe ratio (WDNR, Storm); county-specific,
#'   averaged across years.
#'
#' @source
#' - **Deer abundance, harvest, fawn-to-doe ratio**: Wisconsin DNR Sex-Age-Kill
#'   model (Storm, WDNR; <https://apps.dnr.wi.gov/deermetrics/>).
#' - **CWD prevalence**: WDNR public CWD testing database
#'   (Wisconsin DNR; <https://dnr.wisconsin.gov/topic/WildlifeHabitat/cwd.html>).
#' - **Non-harvest mortality**: Norton et al. (2021).
#' - **Fawn-to-yearling survival**: Warbington et al. (2017).
#' - **Carrying capacity**: Zouwen & Warnke (1995).
#' - **Wolf abundance**: Wisconsin DNR wolf density maps (WDNR 2023).
"full_data"


#' SIP sex-age-infection class labels
#'
#' Fixed ordering of the eight compartment classes used throughout the model.
#' Defines the column order of every population matrix and the row/column order
#' of per-county Lefkovitch and infection arrays. Must not be reordered.
#'
#' @format A length-8 character vector:
#'
#' - `s.m.a` — Susceptible adult male
#' - `s.f.a` — Susceptible adult female
#' - `s.m.y` — Susceptible yearling male
#' - `s.f.y` — Susceptible yearling female
#' - `i.m.a` — Infected adult male
#' - `i.f.a` — Infected adult female
#' - `i.m.y` — Infected yearling male
#' - `i.f.y` — Infected yearling female
#'
#' @examples
#' data(sip_classes)
#' sip_classes
"sip_classes"


#' Wisconsin county borders
#'
#' Polygon geometries for 72 Wisconsin counties, used for mapping model outputs.
#'
#' @format An `sf` object with 72 features and 2 variables:
#' - `county` — County name.
#' - `geometry` — `MULTIPOLYGON` geometry (NAD83(HARN) / Wisconsin Transverse
#'   Mercator, EPSG:3071).
#'
#' @source Wisconsin DNR county boundary shapefile.
"county_borders"


#' Wolf density estimates by Wisconsin county
#'
#' Per-county wolf density (wolves per unit area) with upper and lower bounds.
#'
#' @format A tibble with 72 rows and 4 variables:
#' - `county` — County name.
#' - `wolf_density` — Estimated wolf density (`NaN` where unavailable).
#' - `wolf_density_up` — Upper bound.
#' - `wolf_density_low` — Lower bound.
#'
#' @source Wisconsin DNR wolf monitoring data.
"wolf_density_by_county"