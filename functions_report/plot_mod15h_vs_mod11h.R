################################################################################
#    plot_mod15h_vs_mod11h: creates a plot for a pdf report                    #
#                                                                              #
#    Copyright (C) 2025 Björn C. Rall (https://orcid.org/0000-0002-3191-8389)  #
#                                                                              #
#    This program is free software: you can redistribute it and/or modify      #
#    it under the terms of the GNU General Public License as published by      #
#    the Free Software Foundation, either version 3 of the License, or         #
#    (at your option) any later version.                                       #
#                                                                              #
#    This program is distributed in the hope that it will be useful,           #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of            #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
#    GNU General Public License for more details.                              #
#                                                                              #
#    You should have received a copy of the GNU General Public License         #
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.    #
################################################################################
#' 
#' @description
#'     find the description including parameters here:
#'         https://github.com/b-c-r/CRITTERcode/blob/main/README.md
#'     
#'     find further details including the full statistics here:
#'         https://github.com/b-c-r/CRITTERstatistics/blob/main/statisticsReport.pdf
#'     
#'     if you are interested in the full scientific paper follow:
#'         https://doi.org/10.1101/2025.02.22.639633
#'     
#'     if you use this code, please cite:
#'         Rall et al. (2025): Habitat complexity reduces feeding strength of
#'         freshwater predators (CRITTER) - Code. Zenodo.
#'         https://doi.org/10.5281/zenodo.14894598
#'
#' @return Creates a nice table for a pdf report.
#'
#' @examples
#'
#' # find an example in our report:
#' # https://github.com/b-c-r/CRITTERstatistics/blob/main/statisticsReport.Rmd
#' 

plot_mod15h_vs_mod11h <- function(
    model_fit_15h,                                                              # the mle2 fit object
    model_fit_11h,                                                              # the mle2 fit object
    include_habitat_pics = T,                                                   # include the habitat pictograms, default = True
    pic_x1 =  c( 1.500,  1.600,  2.500,  2.600),                                # lower (left) x values for the 4 habitat pictures, the vector has values for 4 pictograms
    pic_x2 =  c( 1.900,  2.000,  2.900,  3.000),                                # upper (right) x values for the 4 habitat pictures, the vector has values for 4 pictograms
    pic_y1a = c(15.000, 30.000, 27.500, 14.000),                                # lower (left) y values (plot a) for the 4 habitat pictures, the vector has values for 4 pictograms
    pic_y2a = c(17.500, 32.500, 30.000, 16.500),                                # upper (right) y values (plot a) for the 4 habitat pictures, the vector has values for 4 pictograms  
    pic_y1b = c( 0.050,  0.025,  0.030,  0.060),                                # lower (left) y values (plot b) for the 4 habitat pictures, the vector has values for 4 pictograms
    pic_y2b = c( 0.055,  0.030,  0.035,  0.065),                                # upper (right) y values (plot b) for the 4 habitat pictures, the vector has values for 4 pictograms  
    ci_reps = 10000,                                                            # number of samples for the confidence interval lines
    ci_levels = c(0.025, 0.975),                                                # lower and upper confidence limits
    x_res = 1000,                                                               # number of x values for regression line
    ylim = c(13, 40),
    pch = 16,
    cex = 0.5,
    ci_col = "lightgrey",
    journal_style = FALSE                                                       # should the style according to journal guidelines used? Default: FALSE
){
  
  ##############################################################################
  # Calculate the f_max | t_h vs. ring number line
  ##############################################################################
  
  n_rings_sim <- seq(
    0,
    max(model_fit_11h@data$n_rings),
    length.out = x_res
  ) # creates x values for the regression line  
  
  f_max_rings <- 1 / 10^(model_fit_11h@coef[1] + model_fit_11h@coef[2] * n_rings_sim)
  t_h_rings <- 10^(model_fit_11h@coef[1] + model_fit_11h@coef[2] * n_rings_sim)
  
  ##############################################################################
  # Simulation of confidence bands
  ##############################################################################
  
  ci_samples_11 <- MASS::mvrnorm(
    ci_reps,
    mu = model_fit_11h@coef,
    Sigma = bbmle::vcov(model_fit_11h)
  ) # samples from the variance - co-variance matrix parameters
  
  ci_samples_15 <- MASS::mvrnorm(
    ci_reps,
    mu = model_fit_15h@coef,
    Sigma = bbmle::vcov(model_fit_15h)
  ) # samples from the variance - co-variance matrix parameters
  
  ##############################################################################
  
  regline_f_max <- foreach::foreach(
    i = 1:nrow(ci_samples_11),
    .combine = "cbind"
  ) %do% {
    1 / 10^(ci_samples_11[i,1] + ci_samples_11[i,2] * n_rings_sim)
  }
  
  regline_t_h <- foreach::foreach(
    i = 1:nrow(ci_samples_11),
    .combine = "cbind"
  ) %do% {
    10^(ci_samples_11[i,1] + ci_samples_11[i,2] * n_rings_sim)
  }
  
  ##############################################################################

  ci_lines_f_max <- foreach::foreach(
    i = 1:nrow(regline_f_max),
    .combine = "rbind"
  ) %do% {
    quantile(regline_f_max[i,], probs = ci_levels, names = FALSE)
  } # selects the results for the confidence intervals for model 11 f_max
  
  ci_lines_t_h <- foreach::foreach(
    i = 1:nrow(regline_t_h),
    .combine = "rbind"
  ) %do% {
    quantile(regline_t_h[i,], probs = ci_levels, names = FALSE)
  } # selects the results for the confidence intervals for model 11 t_h
  
  ##############################################################################
  # Plot Settings
  ##############################################################################
  
  par(
    mfrow=c(1,2),
    oma = c(4.1,.1,.1,.1),
    mar = c(.25,4.25,.25,.25),
    las = 1
  )

  ##############################################################################
  # Plot (a)
  ##############################################################################
  
  plot(
    c(0,1.975,2.025,2.975,3.025),
    1/model_fit_15h@coef[1:5],
    ylim = ylim,
    pch = pch,
    cex = cex,
    xlab = "",
    ylab = expression(F[max]),
    type = "n",
    bty = ifelse(journal_style, "l", "o"),
    tck = ifelse(journal_style, -0.01, NA)
  ) # setup the plot
  
  polygon(
    c(n_rings_sim, rev(n_rings_sim)),
    c(ci_lines_f_max[,1], rev(ci_lines_f_max[,2])),
    col = ci_col,
    border = ci_col
  ) # adds a polygon to the plot with the CI limits as borders (model 11h)
  
  # adds best fit lines (model 11h):
  lines(n_rings_sim, f_max_rings, lwd = 2, col = "white")
  lines(n_rings_sim, f_max_rings, lwd = 1, col = "black")
  
  # add CI for model 15h
  for(i in 1:5){
    arrows(
      c(0,1.975,2.025,2.975,3.025)[i],
      1/10^quantile(ci_samples_15[,i], ci_levels[1]),
      c(0,1.975,2.025,2.975,3.025)[i],
      1/10^quantile(ci_samples_15[,i], ci_levels[2]),
      length = 0.05,
      angle = 90,
      code = 3
    )
  }
  
  # adds point estimates model 15:
  points(
    c(0,1.975,2.025,2.975,3.025),
    1/10^model_fit_15h@coef[1:5],
    pch = pch,
    cex = cex
  )
  mtext(
    expression(paste("(a) ", F[max], sep = "")),
    adj = .05,
    line = -1.5
  ) # adds plot letter and information
  
  # add pictogram
  
  gh_path <- "https://raw.githubusercontent.com/b-c-r/CRITTERdata/refs/heads/main/"
  f_path <- "pictures/"
  
  for(i in 1:4){
    if(include_habitat_pics){
      png_url <- paste(gh_path, f_path, "habitat_complexity_level_0", i, ".png", sep ="")
      graphics::rasterImage(
        png::readPNG(RCurl::getURLContent(png_url)),                            # see https://stackoverflow.com/questions/12888120/loading-png-files-directly-from-url
        pic_x1[i],
        pic_y1a[i],
        pic_x2[i],
        pic_y2a[i]
      )
    }
  }
  
  ##############################################################################
  # Plot (b)
  ##############################################################################
  
  plot(
    c(0,1.975,2.025,2.975,3.025),
    model_fit_15h@coef[1:5],
    ylim = 1/rev(ylim),
    pch = pch,
    cex = cex,
    xlab = "",
    ylab = expression(T[h]),
    type = "n",
    bty = ifelse(journal_style, "l", "o"),
    tck = ifelse(journal_style, -0.01, NA)
  ) # setup the plot
  
  polygon(
    c(n_rings_sim, rev(n_rings_sim)),
    c(ci_lines_t_h[,1], rev(ci_lines_t_h[,2])),
    col = ci_col,
    border = ci_col
  ) # adds a polygon to the plot with the CI limits as borders (model 11h)
  
  # adds best fit lines (model 11h):
  lines(n_rings_sim, t_h_rings, lwd = 2, col = "white")
  lines(n_rings_sim, t_h_rings, lwd = 1, col = "black")
  
  # add CI for model 15h
  for(i in 1:5){
    arrows(
      c(0,1.975,2.025,2.975,3.025)[i],
      10^quantile(ci_samples_15[,i], ci_levels[1]),
      c(0,1.975,2.025,2.975,3.025)[i],
      10^quantile(ci_samples_15[,i], ci_levels[2]),
      length = 0.05,
      angle = 90,
      code = 3
    )
  }
  
  # adds point estimates model 15:
  points(
    c(0,1.975,2.025,2.975,3.025),
    10^model_fit_15h@coef[1:5],
    pch = pch,
    cex = cex
  )
  mtext(
    expression(paste("(b) ", T[h], sep = "")),
    adj = .05,
    line = -1.5
  ) # adds plot letter and information
  
  # add pictogram
  
  for(i in 1:4){
    if(include_habitat_pics){
      png_url <- paste(gh_path, f_path, "habitat_complexity_level_0", i, ".png", sep ="")
      graphics::rasterImage(
        png::readPNG(RCurl::getURLContent(png_url)),                            # see https://stackoverflow.com/questions/12888120/loading-png-files-directly-from-url
        pic_x1[i],
        pic_y1b[i],
        pic_x2[i],
        pic_y2b[i]
      )
    }
  }

  ##############################################################################
  # axis text
  ##############################################################################
  
  title(
    xlab = "number of rings",
    outer = TRUE
  )
  
}

