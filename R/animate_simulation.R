#' Animated county maps of simulation output
#'
#' Writes a GIF of [mapResults_cwd()] or [mapResults_wtd()] maps, one frame per
#' simulation year. Given a named list of outputs, frames show the scenarios
#' side by side on a common colour scale, labelled beneath each panel.
#'
#' @param outs Output of [expand_output()], or a (named) list of them.
#' @param file Path of the GIF to write.
#' @param variable `"cwd"` (prevalence) or `"wtd"` (abundance).
#' @param sex,age Passed to the map function; see [mapResults_cwd()].
#' @param years Simulation years to animate. Default: all.
#' @param delay Seconds per frame.
#' @param width,height Frame size in pixels.
#' @return The file path, invisibly.
#' @export
animate_simulation <- function(outs, file, variable = c("cwd", "wtd"), sex = "all", age = "all",
                     years = NULL, delay = 0.25, width = 650 * length(outs), height = 750,
                     loop = FALSE, hold = 0, ...) {
  variable <- match.arg(variable)
  if (is.data.frame(outs)) outs <- list(outs)
  if (is.null(years)) years <- sort(unique(outs[[1]]$sim_year))
  years <- c(years, rep(utils::tail(years, 1), round(hold / delay)))
  mapf <- if (variable == "cwd") mapResults_cwd else mapResults_wtd
  lab  <- !is.null(names(outs))
  gifski::save_gif({
    for (y in years) {
      graphics::par(mfrow = c(1, length(outs)), oma = c(if (lab) 3 else 0, 0, 0, 0))
      for (o in outs) mapf(o, y, sex, age, ...)
      if (lab) graphics::mtext(names(outs), side = 1, outer = TRUE, line = 1, cex = 1.6, font = 2,
                               at = (seq_along(outs) - 0.5) / length(outs))
    }
  }, gif_file = file, width = width, height = height, delay = delay, loop = loop, res = 96)
  invisible(file)
}