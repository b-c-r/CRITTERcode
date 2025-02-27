################################################################################
#   rrpe_fit_mod12r: fitting the RRPE to data                                  #
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
#' @description `rrpe_fit_mod12r` allows for fitting the Rogers' Random Predator
#'     Equation function (Royama 1971, Rogers 1972) using the Michaelis-Menten
#'     version of the Type II functional response (Real 1977). You only need to
#'     provide the number of resources eaten and the initial resource density.
#'     The function does the rest. But see the parameters section and the
#'     example for more possibilities. E.g., if your trials ran a day, but you
#'     want to have the maximum feeding rate on an hourly basis, you can enter
#'     t_end = 24.
#'      
#'     Required packages and their dependencies to be installed:
#'     - `emdbook` (Bolker 2023)
#'     - `dplyr` (Wickham et al. 2023)
#'     - `foreach` (Microsoft and Weston 2022)
#'     - `lhs` (Carnell 2024)
#'     Required packages to be attached:
#'     - `dplyr` (Wickham et al. 2023)
#'     - `foreach` (Microsoft and Weston 2022)
#'
#' @references Bolker (2023) emdbook: support functions and data for "Ecological
#'     models and data". Version 1.3.13.
#'     https://CRAN.R-project.org/package=emdbook
#' @references Carnell (2024) lhs: latin hypercube samples. Version 1.2.0.
#'     https://CRAN.R-project.org/package=lhs
#' @references Microsoft and Weston (2022) foreach: provides foreach looping
#'     construct. Version 1.5.2. https://CRAN.R-project.org/package=foreach
#' @references Real (1977) The kinetics of functional response. Am Nat 111, 289-
#'     300. https://doi.org/10.1086/283161
#' @references Rogers (1972) Random search and insect population models. J Anim
#'     Ecol 41, 369-383. https://doi.org/10.2307/3474
#' @references Royama (1971) A comparative study of models for predation and
#'     parasitism. Res Popul Ecol 13, 1-91. https://doi.org/10.1007/BF02511547
#' @references Wickham et al. (2023) dplyr: a grammar of data manipulation.
#'     1.1.4. https://CRAN.R-project.org/package=dplyr
#'
#' @include rrpe_sim.R
#' @include rrpe_nll_mod12r.R
#' @include rrpe_scan_mod12r.R
#'
#' @param n_eaten integer (or float); the prey items that were eaten throughout
#'     the experimental trial. A vector.
#' @param n_initial integer or float; a vector of initial prey densities.
#' @param n_rings number of ring structures (0, 2, or 3), a single integer value.
#' @param complexity level of complexity (0-4), a single integer value.
#' @param p integer or float; a single value of a fixed predator density.
#'     The default value is 1.
#' @param t_end integer or float; the time were the feeding ends. A single
#'     value; default = 1 (e.g. 1 day).
#' @param no_lhs_samples a single integer value; the number of random
#'     latin hypercube samplings.
#' @param range_multiplier The multipliers with which the current best
#'     parameters should be multiplied for the validation random latin hypercube
#'     sampling.
#' @param rel_f_max_range These two values are multiplied by the largest feeding
#'      value of n_eaten to set the initial boundaries to seek for a
#'      reasonable starting value of f_max, default = c(0.6, 0.95)
#' @param rel_n_half_range These two values are multiplied by the largest
#'      starting density value, n_initial, to set the initial boundaries to seek
#'      for a reasonable starting value of n_half, default = c(0.2, 0.8)
#' @param slope_range the range of initial slopes tested for the log-linear
#'      relationship with number of rings. It applies for any FR variable if
#'      fitted to number of rings. The intercept is calculated using complexity
#'      level = 0. Default is c(-0.05, 0.05).
#' @param witer_max How many fits should be performed without convergence?
#' @param mle2_tol The tolerance of a single mle2 fit.
#' @param val_tol The tolerance of the validation.
#' @param set_seed set seed for better reproducibility? default = TRUE.
#' @param seed_value seed value, default = 123.
#'
#' @examples
#' 
#' rm(list=ls())
#' 
#' library("foreach")
#' library("dplyr")
#' 
#' gh_path <- "https://raw.githubusercontent.com/b-c-r/CRITTERcode/refs/heads/main/functions_habitat_statistics/"
#' 
#' source(paste(gh_path, "rrpe_sim.R", sep = ""))
#' source(paste(gh_path, "rrpe_nll_mod12r.R", sep = ""))
#' source(paste(gh_path, "rrpe_scan_mod12r.R", sep = ""))
#' source(paste(gh_path, "rrpe_fit_mod12r.R", sep = ""))
#'
#' ################################################################################
#' ## Data
#' fr_data <- read.csv(
#'   "https://raw.githubusercontent.com/b-c-r/CRITTERdata/refs/heads/main/critter_data.csv"
#'   )
#'   
#' fr_data_ie <- subset(fr_data, predator == "Ischnura elegans")
#' fr_data_ng <- subset(fr_data, predator == "Notonecta glauca")
#' 
#' ################################################################################
#' ## Fit Ischnura elegans
#' mod12r_fit_ie <- rrpe_fit_mod12r(
#'   n_eaten = fr_data_ie$n_eaten,
#'   n_initial = fr_data_ie$n_initial,
#'   n_rings = fr_data_ie$ring_count,
#'   complexity = fr_data_ie$complexity_level
#' )
#'
#' bbmle::summary(mod12r_fit_ie)
#' BIC(mod12r_fit_ie)
#' 
#' ################################################################################
#' ## Fit Notonecta glauca
#' 
#' mod12r_fit_ng <- rrpe_fit_mod12r(
#'   n_eaten = fr_data_ng$n_eaten,
#'   n_initial = fr_data_ng$n_initial,
#'   n_rings = fr_data_ng$ring_count,
#'   complexity = fr_data_ng$complexity_level
#' )
#' 
#' bbmle::summary(mod12r_fit_ng)
#' BIC(mod12r_fit_ng)
#' 
#' 

rrpe_fit_mod12r <- function(
    n_eaten,
    n_initial,
    n_rings,
    complexity,
    p = 1,
    t_end = 1,
    no_lhs_samples = 1000,
    range_multiplier = c(1.0001, 1.001, 1.1, 1.2, 1.3),
    rel_f_max_range = c(0.6, 0.95),
    rel_n_half_range = c(0.2, 0.8),
    slope_range = c(-0.1, 0.1),
    witer_max = 25,
    mle2_tol = 1e-12,
    val_tol = 6,
    set_seed = TRUE,
    seed_value = 123
){
  
  if(set_seed) set.seed(seed_value) # set the seed to assure reproducible
  
  f_max_range_intercept  <- rel_f_max_range *  max(n_eaten[complexity == 0])/t_end
  n_half_range_0 <- rel_n_half_range * max(n_initial[complexity == 0])
  n_half_range_1 <- rel_n_half_range * max(n_initial[complexity == 1])
  n_half_range_2 <- rel_n_half_range * max(n_initial[complexity == 2])
  n_half_range_3 <- rel_n_half_range * max(n_initial[complexity == 3])
  n_half_range_4 <- rel_n_half_range * max(n_initial[complexity == 4])
  
  # initial lhs sampling
  initial_guess <- rrpe_scan_mod12r(
    n_eaten = n_eaten,
    n_initial = n_initial,
    n_rings = n_rings,
    complexity  = complexity,
    p = p,
    f_max_range_intercept_log10 = log10(f_max_range_intercept),
    f_max_range_slope_log10     = slope_range,
    n_half_range_0_log10        = log10(n_half_range_0),
    n_half_range_1_log10        = log10(n_half_range_1),
    n_half_range_2_log10        = log10(n_half_range_2),
    n_half_range_3_log10        = log10(n_half_range_3),
    n_half_range_4_log10        = log10(n_half_range_4),
    t_end = t_end,
    no_lhs_samples = no_lhs_samples
  )
  
  fit <- list()
  nll_fits <- c()

  fit[[1]] <- bbmle::mle2(
    minuslogl = rrpe_nll_mod12r,
    start = list(
      f_max_intercept_log10  =  initial_guess$f_max_intercept_log10,
      f_max_slope_log10  =  initial_guess$f_max_slope_log10,
      n_half_0_log10 = initial_guess$n_half_0_log10,
      n_half_1_log10 = initial_guess$n_half_1_log10,
      n_half_2_log10 = initial_guess$n_half_2_log10,
      n_half_3_log10 = initial_guess$n_half_3_log10,
      n_half_4_log10 = initial_guess$n_half_4_log10
    ),
    fixed = list(
      t_end = t_end,
      p = p
    ),
    data = list(
      n_eaten = n_eaten,
      n_initial = n_initial,
      n_rings = n_rings,
      complexity  = complexity
    ),
    control = list(reltol = mle2_tol)
  )
  nll_fits[1] <- round(bbmle::logLik(fit[[1]])[][1], val_tol)

  message("#########################################")
  message(paste("parameter values after fit 1:", sep = ""))
  message(paste("f_max_intercept_log10: ", round(bbmle::coef(fit[[1]])[2],2)))
  message(paste("f_max_slope_log10: ",     round(bbmle::coef(fit[[1]])[3],2)))
  message(paste("n_half_0_log10: ",        round(bbmle::coef(fit[[1]])[4],2)))
  message(paste("n_half_1_log10: ",        round(bbmle::coef(fit[[1]])[5],2)))
  message(paste("n_half_2_log10: ",        round(bbmle::coef(fit[[1]])[6],2)))
  message(paste("n_half_3_log10: ",        round(bbmle::coef(fit[[1]])[7],2)))
  message(paste("n_half_4_log10: ",        round(bbmle::coef(fit[[1]])[8],2)))
  message(paste("nll: ", nll_fits[1]))
  message("")
  message("start validating the result")
  message("")
  message("#########################################")
  message("")

  witer <- 2
  while(witer <= witer_max){

    start_parms <- foreach::foreach(
      i = 1:length(range_multiplier),
      .combine = "rbind") %do% {

        f_max_range_intercept <- c(10^bbmle::coef(fit[[witer-1]])[2] / range_multiplier[i], 10^bbmle::coef(fit[[witer-1]])[2] * range_multiplier[i])
        f_max_range_slope      <- c(
          bbmle::coef(fit[[witer-1]])[3] - range_multiplier[i] * abs(bbmle::coef(fit[[witer-1]])[3]),
          bbmle::coef(fit[[witer-1]])[3] + range_multiplier[i] * abs(bbmle::coef(fit[[witer-1]])[3])
        )
        n_half_range_0 <- c(10^bbmle::coef(fit[[witer-1]])[4] / range_multiplier[i], 10^bbmle::coef(fit[[witer-1]])[4] * range_multiplier[i])
        n_half_range_1 <- c(10^bbmle::coef(fit[[witer-1]])[5] / range_multiplier[i], 10^bbmle::coef(fit[[witer-1]])[5] * range_multiplier[i])
        n_half_range_2 <- c(10^bbmle::coef(fit[[witer-1]])[6] / range_multiplier[i], 10^bbmle::coef(fit[[witer-1]])[6] * range_multiplier[i])
        n_half_range_3 <- c(10^bbmle::coef(fit[[witer-1]])[7] / range_multiplier[i], 10^bbmle::coef(fit[[witer-1]])[7] * range_multiplier[i])
        n_half_range_4 <- c(10^bbmle::coef(fit[[witer-1]])[8] / range_multiplier[i], 10^bbmle::coef(fit[[witer-1]])[8] * range_multiplier[i])
        
        rrpe_scan_mod12r(
          n_eaten = n_eaten,
          n_initial = n_initial,
          n_rings = n_rings,
          complexity  = complexity,
          p = p,
          f_max_range_intercept_log10  = log10(f_max_range_intercept),
          f_max_range_slope_log10      = f_max_range_slope,
          n_half_range_0_log10         = log10(n_half_range_0),
          n_half_range_1_log10         = log10(n_half_range_1),
          n_half_range_2_log10         = log10(n_half_range_2),
          n_half_range_3_log10         = log10(n_half_range_3),
          n_half_range_4_log10         = log10(n_half_range_4),
          t_end = t_end,
          no_lhs_samples = round(no_lhs_samples/length(range_multiplier))
        )
      }
    
    start_parms <- start_parms %>%
      dplyr::filter(nll == min(nll))
    
    fit[[witer]] <- bbmle::mle2(
      minuslogl = rrpe_nll_mod12r,
      start = list(
        f_max_intercept_log10  =  start_parms$f_max_intercept_log10,
        f_max_slope_log10  =  start_parms$f_max_slope_log10,
        n_half_0_log10 = start_parms$n_half_0_log10,
        n_half_1_log10 = start_parms$n_half_1_log10,
        n_half_2_log10 = start_parms$n_half_2_log10,
        n_half_3_log10 = start_parms$n_half_3_log10,
        n_half_4_log10 = start_parms$n_half_4_log10
      ),
      fixed = list(
        t_end = t_end,
        p = p
      ),
      data = list(
        n_eaten = n_eaten,
        n_initial = n_initial,
        n_rings = n_rings,
        complexity  = complexity
      ),
      control = list(reltol = mle2_tol)
    )

    nll_fits[witer] <- round(bbmle::logLik(fit[[witer]])[][1], val_tol)

    message("#########################################")
    message(paste("parameter values after fit ", witer, ":", sep = ""))
    message(paste("f_max_intercept_log10: ",  round(bbmle::coef(fit[[witer]])[2],2)))
    message(paste("f_max_slope_log10: ",      round(bbmle::coef(fit[[witer]])[3],2)))
    message(paste("n_half_0_log10: ",         round(bbmle::coef(fit[[witer]])[4],2)))
    message(paste("n_half_1_log10: ",         round(bbmle::coef(fit[[witer]])[5],2)))
    message(paste("n_half_2_log10: ",         round(bbmle::coef(fit[[witer]])[6],2)))
    message(paste("n_half_3_log10: ",         round(bbmle::coef(fit[[witer]])[7],2)))
    message(paste("n_half_4_log10: ",         round(bbmle::coef(fit[[witer]])[8],2)))
    message(paste("nll: ", nll_fits[witer]))
    message("#########################################")
    message("")

    if(witer >= 3){
      if(length(unique(nll_fits[(witer-2):witer])) == 1) break
    }
    witer <- witer+1

  }

  if(witer < witer_max){
    return(fit[[witer]])
  } else{
    ret_i <- witer[nll_fits == min(nll_fits)]
    message("no stable fit found, returning best fit")
    return(fit[[ret_i]])
  }
}
