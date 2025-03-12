################################################################################
#    plot_mod05r: creates a plot for a pdf report                              #
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
#'     if you prefer to download a pdf, including the full statistics, follow:
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
#' @return Creates a nice plot for a pdf report.
#'
#' @examples
#' 
#' # see code in the statistical report (https://github.com/b-c-r/CRITTERstatistics)
#' 

plot_mod05r <- function(
    model_fit,                                                                  # the mle2 fit object
    ci_reps = 1000,                                                             # number of samples for the confidence interval lines
    ci_levels = c(0.025, 0.975),                                                # lower and upper confidence limits
    x_res = 1000,                                                               # number of x values for regression line
    ylim = c(0, 25),
    pch = 16,
    cex = 0.5,
    ci_col = "lightgrey",
    no_threads = 10                                                             # number of threads that should be used for simulation
){
  
  ##############################################################################
  # Simulation of best fits
  ##############################################################################
  
  n_in_sim <- seq(
    0,
    max(model_fit@data$n_initial),
    length.out = x_res
  ) # creates x values for the regression line  
  
  best_fit_hab0 <- rrpe_sim(
    fr_style = "Real",
    n_initial = n_in_sim,
    p = 1,
    f_max = 10^model_fit@coef[1],
    n_half = 10^model_fit@coef[3]
  ) # simulates regression line for habitat absent
  
  best_fit_hab1 <- rrpe_sim(
    fr_style = "Real",
    n_initial = n_in_sim,
    p = 1,
    f_max = 10^model_fit@coef[2],
    n_half = 10^model_fit@coef[3]
  ) # simulates regression line for habitat present
  
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
  
  regline_hab0 <- foreach::foreach(
    i = 1:nrow(ci_samples),
    .combine = "cbind"
  ) %dopar% {
    
    gh_path <- "https://raw.githubusercontent.com/b-c-r/CRITTERcode/refs/heads/main/"
    f_path <- "functions_habitat_statistics/"
    source(paste(gh_path, f_path, "rrpe_sim.R", sep = "")) 
    
    rrpe_sim(
      fr_style = "Real",
      n_initial = n_in_sim,
      p = 1,
      f_max = 10^ci_samples[i,1],
      n_half = 10^ci_samples[i,3]
    )[,2]
  } # simulates all regression lines from the samples (case no habitat)
  
  regline_hab1 <- foreach::foreach(
    i = 1:nrow(ci_samples),
    .combine = "cbind"
  ) %dopar% {
    
    gh_path <- "https://raw.githubusercontent.com/b-c-r/CRITTERcode/refs/heads/main/"
    f_path <- "functions_habitat_statistics/"
    source(paste(gh_path, f_path, "rrpe_sim.R", sep = "")) 
    
    rrpe_sim(
      fr_style = "Real",
      n_initial = n_in_sim,
      p = 1,
      f_max = 10^ci_samples[i,2],
      n_half = 10^ci_samples[i,3]
    )[,2]
  } # simulates all regression lines from the samples (case with habitat)
  
  parallel::stopCluster(cl)
  
  ##############################################################################
  ci_lines_hab0 <- foreach::foreach(
    i = 1:nrow(regline_hab0),
    .combine = "rbind"
  ) %do% {
    quantile(regline_hab0[i,], probs = ci_levels, names = FALSE)
  } # selects the results for the confidence intervals (case no habitat)
  
  ci_lines_hab1 <- foreach::foreach(
    i = 1:nrow(regline_hab1),
    .combine = "rbind"
  ) %do% {
    quantile(regline_hab1[i,], probs = ci_levels, names = FALSE)
  } # selects the results for the confidence intervals (case with habitat)
  
  
  ##############################################################################
  # Plot Settings
  ##############################################################################
  
  par(
    mfrow=c(1,2),
    oma = c(4.1,3.75,.1,.1),
    mar = c(.25,.25,.25,.25),
    las = 1
  )

  ##############################################################################
  # Plot (a)
  ##############################################################################
  
  plot(
    model_fit@data[[2]][model_fit@data[[3]] == 0],
    model_fit@data[[1]][model_fit@data[[3]] == 0],
    ylim = ylim,
    pch = pch,
    cex = cex,
    xlab = "",
    ylab = "",
    type = "n"
  ) # setup the plot
  
  polygon(
    c(n_in_sim, rev(n_in_sim)),
    c(ci_lines_hab0[,1], rev(ci_lines_hab0[,2])),
    col = ci_col,
    border = ci_col
  ) # adds a polygon to the plot with the CI limits as borders
  
  # adds best fit lines:
  lines(best_fit_hab0$n_initial, best_fit_hab0$n_eaten, lwd = 2, col = "white")
  lines(best_fit_hab0$n_initial, best_fit_hab0$n_eaten, lwd = 1, col = "black")
  
  # adds data:
  points(
    model_fit@data[[2]][model_fit@data[[3]] == 0],
    model_fit@data[[1]][model_fit@data[[3]] == 0],
    pch = pch,
    cex = cex
  )
  mtext(
    "(a) habitat absent",
    adj = .05,
    line = -1.5
  ) # adds plot letter and information
  
  ##############################################################################
  # Plot (b)
  ##############################################################################
  
  plot(
    model_fit@data[[2]][model_fit@data[[3]] > 0],
    model_fit@data[[1]][model_fit@data[[3]] > 0],
    ylim = ylim,
    pch = pch,
    cex = cex,
    yaxt = "n",
    ylab = "",
    xlab = "",
    type = "n"
  )  # setup the plot
  
  polygon(
    c(n_in_sim, rev(n_in_sim)),
    c(ci_lines_hab1[,1], rev(ci_lines_hab1[,2])),
    col = ci_col,
    border = ci_col
  ) # adds a polygon to the plot with the CI limits as borders
  
  # adds best fit lines:
  lines(best_fit_hab1$n_initial, best_fit_hab1$n_eaten, lwd = 2, col = "white")
  lines(best_fit_hab1$n_initial, best_fit_hab1$n_eaten, lwd = 1, col = "black")
  
  # adds data:
  points(
    model_fit@data[[2]][model_fit@data[[3]] > 0],
    model_fit@data[[1]][model_fit@data[[3]] > 0],
    pch = pch,
    cex = cex
  )
  mtext(
    "(b) habitat present",
    adj = .05,
    line = -1.5
  ) # adds plot letter and information
  
  ##############################################################################
  # axis text
  ##############################################################################
  
  title(
    xlab = "initial prey items in arena",
    ylab = "prey eaten per day",
    outer = TRUE
  )
  
}

