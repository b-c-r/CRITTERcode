################################################################################
#   gen_fr_fit: fitting the generalized functional response model to data      #
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
#' @description `gen_fr_fit` automatically fits the generalized functional
#'     response model (Real 1977, Rosenbaum and Rall 2018) to data. In the
#'     simplest case, you only need to provide the number of resources eaten
#'     (count data) and the initial resource density (count data): the function
#'     does the rest, including initial parameter value guessing. See the
#'     parameters section and the code example for more options. If your
#'     experiment ran a day, but you want to have the maximum feeding rate on an
#'     hourly basis, you can enter t_end = 24. See also the description of
#'     `gen_fr_nll` and `gen_fr_parms_scan` for further information.
#'
#'     Required packages and their dependencies to be installed:
#'       - `bbmle` (Bolker et al 2023)
#'       - `dplyr` (Wickham et al. 2023)
#'       - `foreach` (Microsoft and Weston 2022)
#'       - `lhs` (Carnell 2024)
#'       - `odin` (FitzJohn and Jombart 2024)
#'     Required packages to be attached:
#'       - `dplyr` (Wickham et al. 2023)
#'       - `foreach` (Microsoft and Weston 2022)
#'
#' @references Bolker (2023) bbmle: tools for general Maximum Likelihood
#'     Estimation. Version 1.0.25.1. https://doi.org/10.32614/CRAN.package.bbmle
#' @references Carnell (2024) lhs: latin hypercube samples. Version 1.2.0.
#'     https://CRAN.R-project.org/package=lhs
#' @references FitzJohn and Jombart (2024) odin: ODE generation and integration.
#'     Ver. 1.2.6. https://doi.org/10.32614/CRAN.package.odin
#'     see also: https://github.com/mrc-ide/odin
#' @references Microsoft and Weston (2022) foreach: provides foreach looping
#'     construct. Version 1.5.2. https://CRAN.R-project.org/package=foreach
#' @references Real (1977) The kinetics of functional response. Am Nat 111, 289-
#'     300. https://doi.org/10.1086/283161
#' @references Rosenbaum and Rall (2018) Fitting functional responses: Direct
#'     parameter estimation by simulating differential equations. Methods Ecol
#'     Evol 9, 2076-2090. https://doi.org/10.1111/2041-210X.13039
#' @references Wickham et al. (2023) dplyr: a grammar of data manipulation.
#'     1.1.4. https://CRAN.R-project.org/package=dplyr
#' 
#' @include gen_fr_compile.R
#' @include gen_fr_sim.R
#' @include gen_fr_nll.R
#' @include gen_fr_parms_scan.R
#' 
#' @param n_eaten integer (or float); the prey items that were eaten throughout
#'     the experimental trial. A vector.
#' @param n_initial integer (or float); the initial prey density. A vector of
#'     the same length as n_eaten.
#' @param p The predator density. A single value
#' @param t_start integer or float; the time were the feeding starts. A single
#'     value; default = 0.
#' @param t_end integer or float; the time were the feeding ends. A single
#'     value; default = 1 (e.g. 1 day).
#' @param t_length integer or float; the number of time steps that should be
#'     generated. The more time steps, the more precise the simulation. A single
#'     value; default = 1000.
#' @param penalty a penalty that is added to the nll if the value of q is below
#'     q_low or above q_up. The default= 1000. Equation:
#'         if(q < q_low) nll + penalty*(q-q_low)^2
#'         if(q > q_up) nll + penalty*(q-q_up)^2
#' @param q_low lower soft boundary of q, default = 0 (Type II FR).
#' @param q_up upper soft boundary of q, default = 1 (Type III FR).
#' @param no_lhs_sample a single integer value; the number of random latin
#'     hypercube samplings. Default = 1000.
#' @param range_multiplier The multipliers with which the current best
#'     parameters should be multiplied for the validation random latin hypercube
#'     sampling. Default = c(1.0001, 1.001, 1.1, 1.5, 2).
#' @param witer_max How many fits should be performed without convergence?
#'     Default = 25.
#' @param val_tol The tolerance of the validation. Default = 6 (decimal places).
#' @param mle2_tol The tolerance of a single mle2 fit. Default = 1e-12.
#' @param maxit the maximum number of iterations of a single mle2-fit. Default
#'     here is 5000 (the mle2 default is 500, but the generalized functional
#'     response model requires more iterations)
#'
#' @examples
#' 
#' 
#' library("foreach")
#' library("dplyr")
#' 
#' source(here::here("functions_gen_fr", "gen_fr_compile.R"))
#' source(here::here("functions_gen_fr", "gen_fr_sim.R"))
#' source(here::here("functions_gen_fr", "gen_fr_nll.R"))
#' source(here::here("functions_gen_fr", "gen_fr_parms_scan.R"))
#' source(here::here("functions_gen_fr", "gen_fr_fit.R"))
#' 
#' fr_data <- read.csv("data/fr_data.csv")
#' treats <- sort(unique(fr_data$treatment))
#' fr_data_treat <- subset(fr_data, treatment == treats[1])
#' 
#' gen_fr_fit(
#'   n_eaten = fr_data_treat$n_eaten,
#'   n_initial = fr_data_treat$n_initial
#' )
#' 

gen_fr_fit <- function(
    n_eaten,
    n_initial,
    p = 1,
    t_start = 0,
    t_end = 1,
    t_length = 1000,
    penalty = 1000,
    q_low = 0,
    q_up = 1,
    no_lhs_samples = 1000,
    range_multiplier = c(1.0001, 1.001, 1.1, 1.5, 2),
    witer_max = 25,
    val_tol = 6,
    mle2_tol = 1e-12,
    maxit = 5000
){
  
  gen_fr_compile()
  
  # initial lhs sampling
  initial_guess <- gen_fr_parms_scan(
    n_eaten = n_eaten,
    n_initial = n_initial,
    p = p,
    f_max_range_log10 = log10(c(1, max(n_eaten))/t_end),
    n_half_range_log10 = log10(c(1, max(n_initial))),
    q_range = c(q_low, q_up),
    t_start = t_start,
    t_end = t_end,
    t_length = t_length,
    penalty = penalty,
    q_low = q_low,
    q_up = q_up,
    no_lhs_samples = no_lhs_samples
    )

  fit <- list()
  nll_fits <- c()

  fit[[1]] <- bbmle::mle2(
    minuslogl = gen_fr_nll,
    start = list(
      f_max_log10 = initial_guess$f_max_log10,
      n_half_log10 = initial_guess$n_half_log10,
      q = initial_guess$q
    ),
    fixed = list(
      t_start = t_start,
      t_end = t_end,
      t_length = t_length,
      p = p,
      penalty = penalty,
      q_low = q_low,
      q_up = q_up
    ),
    data = list(
      n_eaten = n_eaten,
      n_initial = n_initial
    ),
    control = list(reltol = mle2_tol, maxit = maxit)
  )
  nll_fits[1] <- round(bbmle::logLik(fit[[1]])[][1], val_tol)


  message("#########################################")
  message(paste("parameter values after fit 1:", sep = ""))
  message(paste("f_max_log10: ", round(bbmle::coef(fit[[1]])[2],2)))
  message(paste("n_half_log10: ", round(bbmle::coef(fit[[1]])[3],2)))
  message(paste("q: ", round(bbmle::coef(fit[[1]])[4],2)))
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
        f_max_range <- c(10^bbmle::coef(fit[[witer-1]])[2]/range_multiplier[i], 10^bbmle::coef(fit[[witer-1]])[2]*range_multiplier[i])
        n_half_range <- c(10^bbmle::coef(fit[[witer-1]])[3]/range_multiplier[i], 10^bbmle::coef(fit[[witer-1]])[3]*range_multiplier[i])
        h_range <- c( (1+bbmle::coef(fit[[witer-1]])[4])  /  range_multiplier[i], (1+bbmle::coef(fit[[witer-1]])[4])  *  range_multiplier[i] )
        q_range <- h_range-1
        q_range[q_range < q_low] <- q_low
        q_range[q_range > q_up] <- q_up

        gen_fr_parms_scan(
          n_eaten = n_eaten,
          n_initial = n_initial,
          p = p,
          f_max_range_log10 = log10(f_max_range),
          n_half_range_log10 = log10(n_half_range),
          q_range = q_range,
          t_start = t_start,
          t_end = t_end,
          t_length = t_length,
          penalty = penalty,
          q_low = q_low,
          q_up = q_up,
          no_lhs_samples = ceiling(no_lhs_samples/length(range_multiplier))
        )
      }

    start_parms <- start_parms %>%
      dplyr::filter(nll == min(nll))

    fit[[witer]] <- bbmle::mle2(
      minuslogl = gen_fr_nll,
      start = list(
        f_max_log10 = start_parms$f_max_log10,
        n_half_log10 = start_parms$n_half_log10,
        q = start_parms$q
      ),
      fixed = list(
        t_start = t_start,
        t_end = t_end,
        t_length = t_length,
        p = p,
        penalty = penalty,
        q_low = q_low,
        q_up = q_up
      ),
      data = list(
        n_eaten = n_eaten,
        n_initial = n_initial
      ),
      control = list(reltol = mle2_tol, maxit = maxit)
    )

    nll_fits[witer] <- round(bbmle::logLik(fit[[witer]])[][1], val_tol)


    message("#########################################")
    message(paste("parameter values after fit ", witer, ":", sep = ""))
    message(paste("f_max_log10: ", round(bbmle::coef(fit[[witer]])[2],2)))
    message(paste("n_half_log10: ", round(bbmle::coef(fit[[witer]])[3],2)))
    message(paste("q: ", round(bbmle::coef(fit[[witer]])[4],2)))
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
