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
#' @include gen_fr_compile.R
#' @include gen_fr_sim.R
#' @include gen_fr_nll.R
#' @include gen_fr_parms_scan.R
#' 
#' @return returns a mle2 fit result.
#' 
#' @examples
#' 
#' # find an executable example here:
#' # https://github.com/b-c-r/CRITTERcode/examples_habitat_statistics/examples_habitat_statistics/gen_fr_fit_examples.R
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
    rel_f_max_range = c(0.6, 0.95),
    rel_n_half_range = c(0.2, 0.8),
    witer_max = 25,
    val_tol = 6,
    mle2_tol = 1e-12,
    maxit = 5000
){
  
  gen_fr_compile()
  
  f_max_range  <- rel_f_max_range *  max(n_eaten)/t_end
  n_half_range <- rel_n_half_range * max(n_initial)
  
  # initial lhs sampling
  initial_guess <- gen_fr_parms_scan(
    n_eaten = n_eaten,
    n_initial = n_initial,
    p = p,
    f_max_range_log10 = log10(f_max_range),
    n_half_range_log10 = log10(n_half_range),
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
