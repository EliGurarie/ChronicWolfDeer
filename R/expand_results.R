#' Expand simulation output into per-class counts, proportions, and wolf presence
#'
#' Takes the output of [runModel_sim()] (one row per county per simulated year,
#' with the eight `*_sim` class counts) and adds the aggregate counts and
#' within-group proportions used for plotting and exploration, plus a
#' `wolf_presence` flag joined from [wolf_density_by_county].
#'
#' This is the bridge between a raw model run and any downstream plotting,
#' summarising, or app display: `runModel_sim() |> expand_output() |> ...`.
#'
#' @param sim A data frame from [runModel_sim()], containing `SMA_sim` ...
#'   `IFY_sim`, `county`, `sim_year`, and `year`.
#'
#' @return `sim` with renamed dotted class columns (`S.M.A` ... `I.F.Y`), all
#'   aggregate counts (e.g. `ma_sim`, `i_all_sim`, `wtd_count_sim`), within-group
#'   proportions (e.g. `ima_sim_perc`, `i_all_sim_perc`), and a logical
#'   `wolf_presence` column.
#'
#' @seealso [mapResults_cwd()], [mapResults_wtd()]
#' @export
expand_sim <- function(sim) {
  d <- sim
  # rename the eight simulated counts to the dotted class names
  names(d)[match(c("SMA_sim","SFA_sim","SMY_sim","SFY_sim",
                   "IMA_sim","IFA_sim","IMY_sim","IFY_sim"), names(d))] <-
    c("S.M.A","S.F.A","S.M.Y","S.F.Y","I.M.A","I.F.A","I.M.Y","I.F.Y")
  
  d <- dplyr::mutate(
    d,
    wtd_count_sim = S.M.A + S.F.A + S.M.Y + S.F.Y + I.M.A + I.F.A + I.M.Y + I.F.Y,
    ma_sim = S.M.A + I.M.A,
    fa_sim = S.F.A + I.F.A,
    my_sim = S.M.Y + I.M.Y,
    fy_sim = S.F.Y + I.F.Y,
    
    im_sim = I.M.A + I.M.Y,
    if_sim = I.F.A + I.F.Y,
    iy_sim = I.M.Y + I.F.Y,
    ia_sim = I.M.A + I.F.A,
    
    sm_sim = S.M.A + S.M.Y,
    sf_sim = S.F.A + S.F.Y,
    sa_sim = S.M.A + S.F.A,
    sy_sim = S.M.Y + S.F.Y,
    
    m_all_sim = S.M.A + I.M.A + S.M.Y + I.M.Y,
    f_all_sim = S.F.A + I.F.A + S.F.Y + I.F.Y,
    a_all_sim = S.M.A + S.F.A + I.M.A + I.F.A,
    y_all_sim = S.M.Y + S.F.Y + I.M.Y + I.F.Y,
    i_all_sim = I.M.A + I.F.A + I.M.Y + I.F.Y,
    s_all_sim = S.M.A + S.F.A + S.M.Y + S.F.Y,
    
    sma_sim_perc = ifelse(ma_sim > 0, S.M.A / ma_sim, 0),
    sfa_sim_perc = ifelse(fa_sim > 0, S.F.A / fa_sim, 0),
    ima_sim_perc = ifelse(ma_sim > 0, I.M.A / ma_sim, 0),
    ifa_sim_perc = ifelse(fa_sim > 0, I.F.A / fa_sim, 0),
    smy_sim_perc = ifelse(my_sim > 0, S.M.Y / my_sim, 0),
    sfy_sim_perc = ifelse(fy_sim > 0, S.F.Y / fy_sim, 0),
    imy_sim_perc = ifelse(my_sim > 0, I.M.Y / my_sim, 0),
    ify_sim_perc = ifelse(fy_sim > 0, I.F.Y / fy_sim, 0),
    
    ma_sim_perc = ifelse(a_all_sim > 0, ma_sim / a_all_sim, 0),
    fa_sim_perc = ifelse(a_all_sim > 0, fa_sim / a_all_sim, 0),
    my_sim_perc = ifelse(y_all_sim > 0, my_sim / y_all_sim, 0),
    fy_sim_perc = ifelse(y_all_sim > 0, fy_sim / y_all_sim, 0),
    
    im_sim_perc = ifelse(m_all_sim > 0, im_sim / m_all_sim, 0),
    if_sim_perc = ifelse(f_all_sim > 0, if_sim / f_all_sim, 0),
    ia_sim_perc = ifelse(a_all_sim > 0, ia_sim / a_all_sim, 0),
    iy_sim_perc = ifelse(y_all_sim > 0, iy_sim / y_all_sim, 0),
    
    sm_sim_perc = ifelse(m_all_sim > 0, sm_sim / m_all_sim, 0),
    sf_sim_perc = ifelse(f_all_sim > 0, sf_sim / f_all_sim, 0),
    sa_sim_perc = ifelse(a_all_sim > 0, sa_sim / a_all_sim, 0),
    sy_sim_perc = ifelse(y_all_sim > 0, sy_sim / y_all_sim, 0),
    
    m_all_sim_perc = ifelse(wtd_count_sim > 0, m_all_sim / wtd_count_sim, 0),
    f_all_sim_perc = ifelse(wtd_count_sim > 0, f_all_sim / wtd_count_sim, 0),
    a_all_sim_perc = ifelse(wtd_count_sim > 0, a_all_sim / wtd_count_sim, 0),
    y_all_sim_perc = ifelse(wtd_count_sim > 0, y_all_sim / wtd_count_sim, 0),
    i_all_sim_perc = ifelse(wtd_count_sim > 0, i_all_sim / wtd_count_sim, 0),
    s_all_sim_perc = ifelse(wtd_count_sim > 0, s_all_sim / wtd_count_sim, 0)
  )
  
  # wolf presence: density present and non-zero
  wp <- data.frame(
    county = wolf_density_by_county$county,
    wolf_presence = !is.na(wolf_density_by_county$wolf_density) &
      wolf_density_by_county$wolf_density > 0
  )
  dplyr::left_join(d, wp, by = "county")
}