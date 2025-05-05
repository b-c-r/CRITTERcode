################################################################################
#    plot_mod15h: creates a plot for a pdf report                              #
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

plot_mod15h <- function(
    model_fit,                                                                  # the mle2 fit object
    include_habitat_pics = T,                                                   # include the habitat pictograms, default = True
    pic_x1 = rep(140.0, 4),                                                     # lower (left) x values for the 4 habitat pictures, the vector has values for 4 pictograms
    pic_x2 = rep(180.0, 4),                                                     # upper (right) x values for the 4 habitat pictures, the vector has values for 4 pictograms
    pic_y1 = rep( 30.0, 4),                                                     # lower (left) y values for the 4 habitat pictures, the vector has values for 4 pictograms
    pic_y2 = rep( 40.0, 4),                                                     # upper (right) y values for the 4 habitat pictures, the vector has values for 4 pictograms  
    ci_reps = 10000,                                                            # number of samples for the confidence interval lines
    ci_levels = c(0.025, 0.975),                                                # lower and upper confidence limits
    x_res = 1000,                                                               # number of x values for regression line
    ylim = c(0, 40),
    pch = 16,
    cex = 0.5,
    ci_col = "lightgrey",
    no_threads = max(
      c(1,2,5,10,20,40,50)[c(1,2,5,10,20,40,50) <= parallel::detectCores()]
      ),                                                                        # number of threads that should be used for simulation
    export_functions_to_workers = "rrpe_sim",
    journal_style = FALSE                                                       # should the style according to journal guidelines used? Default: FALSE
){
  
  ##############################################################################
  # Simulation of best fits
  ##############################################################################
  
  n_in_sim <- seq(
    0,
    max(model_fit@data$n_initial),
    length.out = x_res
  ) # creates x values for the regression line  
  
  best_fit_c0 <- rrpe_sim(
    fr_style = "Holling",
    n_initial = n_in_sim,
    p = 1,
    t_h = 10^model_fit@coef[1],
    a = 10^(model_fit@coef[6] + model_fit@coef[7] * 0)
  ) # simulates regression line for complexity level 0
  
  best_fit_c1 <- rrpe_sim(
    fr_style = "Holling",
    n_initial = n_in_sim,
    p = 1,
    t_h = 10^model_fit@coef[2],
    a = 10^(model_fit@coef[6] + model_fit@coef[7] * 2)
  ) # simulates regression line for complexity level 1
  
  best_fit_c2 <- rrpe_sim(
    fr_style = "Holling",
    n_initial = n_in_sim,
    p = 1,
    t_h = 10^model_fit@coef[3],
    a = 10^(model_fit@coef[6] + model_fit@coef[7] * 2)
  ) # simulates regression line for complexity level 2
  
  best_fit_c3 <- rrpe_sim(
    fr_style = "Holling",
    n_initial = n_in_sim,
    p = 1,
    t_h = 10^model_fit@coef[4],
    a = 10^(model_fit@coef[6] + model_fit@coef[7] * 3)
  ) # simulates regression line for complexity level 3
  
  best_fit_c4 <- rrpe_sim(
    fr_style = "Holling",
    n_initial = n_in_sim,
    p = 1,
    t_h = 10^model_fit@coef[5],
    a = 10^(model_fit@coef[6] + model_fit@coef[7] * 3)
  ) # simulates regression line for complexity level 4
  
  ##############################################################################
  # Simulation of confidence bands
  ##############################################################################
  
  ci_samples <- MASS::mvrnorm(
    ci_reps,
    mu = model_fit@coef,
    Sigma = bbmle::vcov(model_fit)
  ) # samples from the variance - co-variance matrix parameters

  ## setup the cluster for parallel computing on the local machine
  cl <- parallel::makeCluster(no_threads)
  doParallel::registerDoParallel(cl)
  
  regline_c0 <- foreach::foreach(
    i = 1:nrow(ci_samples),
    .export = export_functions_to_workers,
    .combine = "cbind"
  ) %dopar% {
    
    rrpe_sim(
      fr_style = "Holling",
      n_initial = n_in_sim,
      p = 1,
      t_h = 10^ci_samples[i,1],
      a = 10^(ci_samples[i,6] + ci_samples[i,7] * 0)
    )[,2]
  } # simulates all regression lines from the samples for complexity level 0
  
  regline_c1 <- foreach::foreach(
    i = 1:nrow(ci_samples),
    .export = export_functions_to_workers,
    .combine = "cbind"
  ) %dopar% {

    rrpe_sim(
      fr_style = "Holling",
      n_initial = n_in_sim,
      p = 1,
      t_h = 10^ci_samples[i,2],
      a = 10^(ci_samples[i,6] + ci_samples[i,7] * 2)
    )[,2]
  } # simulates all regression lines from the samples for complexity level 1
  
  regline_c2 <- foreach::foreach(
    i = 1:nrow(ci_samples),
    .export = export_functions_to_workers,
    .combine = "cbind"
  ) %dopar% {

    rrpe_sim(
      fr_style = "Holling",
      n_initial = n_in_sim,
      p = 1,
      t_h = 10^ci_samples[i,3],
      a = 10^(ci_samples[i,6] + ci_samples[i,7] * 2)
    )[,2]
  } # simulates all regression lines from the samples for complexity level 2
  
  regline_c3 <- foreach::foreach(
    i = 1:nrow(ci_samples),
    .export = export_functions_to_workers,
    .combine = "cbind"
  ) %dopar% {

    rrpe_sim(
      fr_style = "Holling",
      n_initial = n_in_sim,
      p = 1,
      t_h = 10^ci_samples[i,4],
      a = 10^(ci_samples[i,6] + ci_samples[i,7] * 3)
    )[,2]
  } # simulates all regression lines from the samples for complexity level 3
  
  regline_c4 <- foreach::foreach(
    i = 1:nrow(ci_samples),
    .export = export_functions_to_workers,
    .combine = "cbind"
  ) %dopar% {

    rrpe_sim(
      fr_style = "Holling",
      n_initial = n_in_sim,
      p = 1,
      t_h = 10^ci_samples[i,5],
      a = 10^(ci_samples[i,6] + ci_samples[i,7] * 3)
    )[,2]
  } # simulates all regression lines from the samples for complexity level 4

  parallel::stopCluster(cl)
  
  ##############################################################################

  ci_lines_c0 <- foreach::foreach(
    i = 1:nrow(regline_c0),
    .combine = "rbind"
  ) %do% {
    quantile(regline_c0[i,], probs = ci_levels, names = FALSE)
  } # selects the results for the confidence intervals for complexity level 0
  
  ci_lines_c1 <- foreach::foreach(
    i = 1:nrow(regline_c1),
    .combine = "rbind"
  ) %do% {
    quantile(regline_c1[i,], probs = ci_levels, names = FALSE)
  } # selects the results for the confidence intervals for complexity level 1
  
  ci_lines_c2 <- foreach::foreach(
    i = 1:nrow(regline_c2),
    .combine = "rbind"
  ) %do% {
    quantile(regline_c2[i,], probs = ci_levels, names = FALSE)
  } # selects the results for the confidence intervals for complexity level 2  
  
  ci_lines_c3 <- foreach::foreach(
    i = 1:nrow(regline_c3),
    .combine = "rbind"
  ) %do% {
    quantile(regline_c3[i,], probs = ci_levels, names = FALSE)
  } # selects the results for the confidence intervals for complexity level 3  
  
  ci_lines_c4 <- foreach::foreach(
    i = 1:nrow(regline_c4),
    .combine = "rbind"
  ) %do% {
    quantile(regline_c4[i,], probs = ci_levels, names = FALSE)
  } # selects the results for the confidence intervals for complexity level 4 
  
  ##############################################################################
  # Plot Settings
  ##############################################################################
  
  par(
    mfrow=c(2,3),
    oma = c(4.1,3.75,.1,.1),
    mar = c(.25,.25,.25,.25),
    las = 1
  )

  ##############################################################################
  # Plot (a)
  ##############################################################################
  
  plot(
    model_fit@data[[2]][model_fit@data[[4]] == 0],
    model_fit@data[[1]][model_fit@data[[4]] == 0],
    xlim = c(0, max(model_fit@data$n_initial)),
    ylim = ylim,
    pch = pch,
    cex = cex,
    xlab = "",
    ylab = "",
    type = "n",
    bty = ifelse(journal_style, "l", "o"),
    tck = ifelse(journal_style, -0.01, NA)
  ) # setup the plot
  
  polygon(
    c(n_in_sim, rev(n_in_sim)),
    c(ci_lines_c0[,1], rev(ci_lines_c0[,2])),
    col = ci_col,
    border = ci_col
  ) # adds a polygon to the plot with the CI limits as borders
  
  # adds best fit lines:
  lines(best_fit_c0$n_initial, best_fit_c0$n_eaten, lwd = 2, col = "white")
  lines(best_fit_c0$n_initial, best_fit_c0$n_eaten, lwd = 1, col = "black")
  
  # adds data:
  points(
    model_fit@data[[2]][model_fit@data[[4]] == 0],
    model_fit@data[[1]][model_fit@data[[4]] == 0],
    pch = pch,
    cex = cex
  )
  mtext(
    "(a) complexity level 0",
    adj = .05,
    line = -1.5
  ) # adds plot letter and information

  ##############################################################################
  # Plot (b)
  ##############################################################################
  
  plot(
    model_fit@data[[2]][model_fit@data[[4]] == 1],
    model_fit@data[[1]][model_fit@data[[4]] == 1],
    xlim = c(0, max(model_fit@data$n_initial)),
    ylim = ylim,
    pch = pch,
    cex = cex,
    xaxt = "n",
    yaxt = "n",
    ylab = "",
    xlab = "",
    type = "n",
    bty = ifelse(journal_style, "l", "o"),
    tck = ifelse(journal_style, -0.01, NA)
  )  # setup the plot
  
  polygon(
    c(n_in_sim, rev(n_in_sim)),
    c(ci_lines_c1[,1], rev(ci_lines_c1[,2])),
    col = ci_col,
    border = ci_col
  ) # adds a polygon to the plot with the CI limits as borders
  
  # adds best fit lines:
  lines(best_fit_c1$n_initial, best_fit_c1$n_eaten, lwd = 2, col = "white")
  lines(best_fit_c1$n_initial, best_fit_c1$n_eaten, lwd = 1, col = "black")
  
  # adds data:
  points(
    model_fit@data[[2]][model_fit@data[[4]] == 1],
    model_fit@data[[1]][model_fit@data[[4]] == 1],
    pch = pch,
    cex = cex
  )
  mtext(
    "(b) complexity level 1",
    adj = .05,
    line = -1.5
  ) # adds plot letter and information
  
  # add pictogram
  
  gh_path <- "https://raw.githubusercontent.com/b-c-r/CRITTERdata/refs/heads/main/"
  f_path <- "pictures/"
  
  if(include_habitat_pics){
    png_url <- paste(gh_path, f_path, "habitat_complexity_level_01.png", sep ="")
    graphics::rasterImage(
      png::readPNG(RCurl::getURLContent(png_url)),                              # see https://stackoverflow.com/questions/12888120/loading-png-files-directly-from-url
      pic_x1[1],
      pic_y1[1],
      pic_x2[1],
      pic_y2[1]
    )
  }
  
  ##############################################################################
  # Plot (c)
  ##############################################################################
  
  plot(
    model_fit@data[[2]][model_fit@data[[4]] == 3],
    model_fit@data[[1]][model_fit@data[[4]] == 3],
    xlim = c(0, max(model_fit@data$n_initial)),
    ylim = ylim,
    pch = pch,
    cex = cex,
    xaxt = "n",
    yaxt = "n",
    ylab = "",
    xlab = "",
    type = "n",
    bty = ifelse(journal_style, "l", "o"),
    tck = ifelse(journal_style, -0.01, NA)
  )  # setup the plot
  
  polygon(
    c(n_in_sim, rev(n_in_sim)),
    c(ci_lines_c3[,1], rev(ci_lines_c3[,2])),
    col = ci_col,
    border = ci_col
  ) # adds a polygon to the plot with the CI limits as borders
  
  # adds best fit lines:
  lines(best_fit_c3$n_initial, best_fit_c3$n_eaten, lwd = 2, col = "white")
  lines(best_fit_c3$n_initial, best_fit_c3$n_eaten, lwd = 1, col = "black")
  
  # adds data:
  points(
    model_fit@data[[2]][model_fit@data[[4]] == 3],
    model_fit@data[[1]][model_fit@data[[4]] == 3],
    pch = pch,
    cex = cex
  )
  mtext(
    "(c) complexity level 3",
    adj = .05,
    line = -1.5
  ) # adds plot letter and information
  
  # add pictogram
  
  gh_path <- "https://raw.githubusercontent.com/b-c-r/CRITTERdata/refs/heads/main/"
  f_path <- "pictures/"
  
  if(include_habitat_pics){
    png_url <- paste(gh_path, f_path, "habitat_complexity_level_03.png", sep ="")
    graphics::rasterImage(
      png::readPNG(RCurl::getURLContent(png_url)),                              # see https://stackoverflow.com/questions/12888120/loading-png-files-directly-from-url
      pic_x1[3],
      pic_y1[3],
      pic_x2[3],
      pic_y2[3]
    )
  }
  
  ##############################################################################
  # Empty Plot
  ##############################################################################
  
  plot(
    0,
    0,
    ylim = ylim,
    pch = pch,
    cex = cex,
    yaxt = "n",
    xaxt = "n",
    xlab = "",
    ylab = "",
    type = "n",
    bty = "n"
  )
  
  
  ##############################################################################
  # Plot (d)
  ##############################################################################
  
  plot(
    model_fit@data[[2]][model_fit@data[[4]] == 2],
    model_fit@data[[1]][model_fit@data[[4]] == 2],
    xlim = c(0, max(model_fit@data$n_initial)),
    ylim = ylim,
    pch = pch,
    cex = cex,
    xlab = "",
    ylab = "",
    type = "n",
    bty = ifelse(journal_style, "l", "o"),
    tck = ifelse(journal_style, -0.01, NA)
  ) # setup the plot
  
  polygon(
    c(n_in_sim, rev(n_in_sim)),
    c(ci_lines_c2[,1], rev(ci_lines_c2[,2])),
    col = ci_col,
    border = ci_col
  ) # adds a polygon to the plot with the CI limits as borders
  
  # adds best fit lines:
  lines(best_fit_c2$n_initial, best_fit_c2$n_eaten, lwd = 2, col = "white")
  lines(best_fit_c2$n_initial, best_fit_c2$n_eaten, lwd = 1, col = "black")
  
  # adds data:
  points(
    model_fit@data[[2]][model_fit@data[[4]] == 2],
    model_fit@data[[1]][model_fit@data[[4]] == 2],
    pch = pch,
    cex = cex
  )
  mtext(
    "(d) complexity level 2",
    adj = .05,
    line = -1.5
  ) # adds plot letter and information
  
  # add pictogram
  
  gh_path <- "https://raw.githubusercontent.com/b-c-r/CRITTERdata/refs/heads/main/"
  f_path <- "pictures/"
  
  if(include_habitat_pics){
    png_url <- paste(gh_path, f_path, "habitat_complexity_level_02.png", sep ="")
    graphics::rasterImage(
      png::readPNG(RCurl::getURLContent(png_url)),                              # see https://stackoverflow.com/questions/12888120/loading-png-files-directly-from-url
      pic_x1[2],
      pic_y1[2],
      pic_x2[2],
      pic_y2[2]
    )
  }
  
  ##############################################################################
  # Plot (e)
  ##############################################################################
  
  plot(
    model_fit@data[[2]][model_fit@data[[4]] == 4],
    model_fit@data[[1]][model_fit@data[[4]] == 4],
    xlim = c(0, max(model_fit@data$n_initial)),
    ylim = ylim,
    yaxt = "n",
    pch = pch,
    cex = cex,
    xlab = "",
    ylab = "",
    type = "n",
    bty = ifelse(journal_style, "l", "o"),
    tck = ifelse(journal_style, -0.01, NA)
  ) # setup the plot
  
  polygon(
    c(n_in_sim, rev(n_in_sim)),
    c(ci_lines_c4[,1], rev(ci_lines_c4[,2])),
    col = ci_col,
    border = ci_col
  ) # adds a polygon to the plot with the CI limits as borders
  
  # adds best fit lines:
  lines(best_fit_c4$n_initial, best_fit_c4$n_eaten, lwd = 2, col = "white")
  lines(best_fit_c4$n_initial, best_fit_c4$n_eaten, lwd = 1, col = "black")
  
  # adds data:
  points(
    model_fit@data[[2]][model_fit@data[[4]] == 4],
    model_fit@data[[1]][model_fit@data[[4]] == 4],
    pch = pch,
    cex = cex
  )
  mtext(
    "(e) complexity level 4",
    adj = .05,
    line = -1.5
  ) # adds plot letter and information
  
  # add pictogram
  
  gh_path <- "https://raw.githubusercontent.com/b-c-r/CRITTERdata/refs/heads/main/"
  f_path <- "pictures/"
  
  if(include_habitat_pics){
    png_url <- paste(gh_path, f_path, "habitat_complexity_level_04.png", sep ="")
    graphics::rasterImage(
      png::readPNG(RCurl::getURLContent(png_url)),                              # see https://stackoverflow.com/questions/12888120/loading-png-files-directly-from-url
      pic_x1[4],
      pic_y1[4],
      pic_x2[4],
      pic_y2[4]
    )
  }
  
  ##############################################################################
  # axis text
  ##############################################################################
  
  title(
    xlab = "initial prey items in arena",
    ylab = "prey eaten per day",
    outer = TRUE
  )
  
}

