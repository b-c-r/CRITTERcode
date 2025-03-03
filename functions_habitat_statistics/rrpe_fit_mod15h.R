################################################################################
#   rrpe_fit_mod15h: fitting the RRPE to data                                  #
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
#'      please also consider citing for the underlying method:
#'         Bolker (2008) Ecological models and data in R, Princeton University Press,
#'         Princeton, New Jersey.
#'         https://math.mcmaster.ca/~bolker/emdbook/index.html
#' 
#' @include rrpe_sim.R
#' 
#' @return Returns a single negative log-likelihood value.
#' 
#' @examples
#' 
#' # find an executable example here:
#' # https://github.com/b-c-r/CRITTERcode/examples_habitat_statistics/examples_habitat_statistics/mod15h_examples.R
#' 

rrpe_fit_mod15h <- function(
    n_eaten,
    n_initial,
    n_rings,
    complexity,
    p = 1,
    t_end = 1,
    no_lhs_samples = 1000,
    range_multiplier = c(1.0001, 1.001, 1.1, 1.5, 2),
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
  
  f_max_range_0  <- rel_f_max_range *  max(n_eaten[complexity == 0])/t_end
  f_max_range_1  <- rel_f_max_range *  max(n_eaten[complexity == 1])/t_end
  f_max_range_2  <- rel_f_max_range *  max(n_eaten[complexity == 2])/t_end
  f_max_range_3  <- rel_f_max_range *  max(n_eaten[complexity == 3])/t_end
  f_max_range_4  <- rel_f_max_range *  max(n_eaten[complexity == 4])/t_end
  n_half_range_0 <- rel_n_half_range * max(n_initial[complexity == 0])
  
  # initial lhs sampling
  initial_guess <- rrpe_scan_mod15h(
    n_eaten = n_eaten,
    n_initial = n_initial,
    n_rings = n_rings,
    complexity  = complexity,
    p = p,
    t_h_range_0_log10 = log10(1/rev(f_max_range_0)),
    t_h_range_1_log10 = log10(1/rev(f_max_range_1)),
    t_h_range_2_log10 = log10(1/rev(f_max_range_2)),
    t_h_range_3_log10 = log10(1/rev(f_max_range_3)),
    t_h_range_4_log10 = log10(1/rev(f_max_range_4)),
    a_range_intercept_log10 = log10(c(f_max_range_0[1]/n_half_range_0[2], f_max_range_0[2]/n_half_range_0[1])),
    a_range_slope_log10 = slope_range,
    t_end = t_end,
    no_lhs_samples = no_lhs_samples
  )
  
  fit <- list()
  nll_fits <- c()

  fit[[1]] <- bbmle::mle2(
    minuslogl = rrpe_nll_mod15h,
    start = list(
      t_h_0_log10  =  initial_guess$t_h_0_log10,
      t_h_1_log10  =  initial_guess$t_h_1_log10,
      t_h_2_log10  =  initial_guess$t_h_2_log10,
      t_h_3_log10  =  initial_guess$t_h_3_log10,
      t_h_4_log10  =  initial_guess$t_h_4_log10,
      a_intercept_log10 = initial_guess$a_intercept_log10,
      a_slope_log10 = initial_guess$a_slope_log10
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
  message(paste("t_h_0_log10: ", round(bbmle::coef(fit[[1]])[2],2)))
  message(paste("t_h_1_log10: ", round(bbmle::coef(fit[[1]])[3],2)))
  message(paste("t_h_2_log10: ", round(bbmle::coef(fit[[1]])[4],2)))
  message(paste("t_h_3_log10: ", round(bbmle::coef(fit[[1]])[5],2)))
  message(paste("t_h_4_log10: ", round(bbmle::coef(fit[[1]])[6],2)))
  message(paste("a_intercept_log10: ", round(bbmle::coef(fit[[1]])[7],2)))
  message(paste("a_slope_log10: "    , round(bbmle::coef(fit[[1]])[8],2)))
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
        t_h_range_0  <- c(10^bbmle::coef(fit[[witer-1]])[2]  / range_multiplier[i], 10^bbmle::coef(fit[[witer-1]])[2]  * range_multiplier[i])
        t_h_range_1  <- c(10^bbmle::coef(fit[[witer-1]])[3]  / range_multiplier[i], 10^bbmle::coef(fit[[witer-1]])[3]  * range_multiplier[i])
        t_h_range_2  <- c(10^bbmle::coef(fit[[witer-1]])[4]  / range_multiplier[i], 10^bbmle::coef(fit[[witer-1]])[4]  * range_multiplier[i])
        t_h_range_3  <- c(10^bbmle::coef(fit[[witer-1]])[5]  / range_multiplier[i], 10^bbmle::coef(fit[[witer-1]])[5]  * range_multiplier[i])
        t_h_range_4  <- c(10^bbmle::coef(fit[[witer-1]])[6]  / range_multiplier[i], 10^bbmle::coef(fit[[witer-1]])[6]  * range_multiplier[i])
        a_range_intercept <- c(10^bbmle::coef(fit[[witer-1]])[7] / range_multiplier[i], 10^bbmle::coef(fit[[witer-1]])[7] * range_multiplier[i])
        a_range_slope     <- c(
          bbmle::coef(fit[[witer-1]])[8] - range_multiplier[i] * abs(bbmle::coef(fit[[witer-1]])[8]),
          bbmle::coef(fit[[witer-1]])[8] + range_multiplier[i] * abs(bbmle::coef(fit[[witer-1]])[8])
        )
        
        rrpe_scan_mod15h(
          n_eaten = n_eaten,
          n_initial = n_initial,
          n_rings = n_rings,
          complexity  = complexity,
          p = p,
          t_h_range_0_log10  =  log10(t_h_range_0),
          t_h_range_1_log10  =  log10(t_h_range_1),
          t_h_range_2_log10  =  log10(t_h_range_2),
          t_h_range_3_log10  =  log10(t_h_range_3),
          t_h_range_4_log10  =  log10(t_h_range_4),
          a_range_intercept_log10 = log10(a_range_intercept),
          a_range_slope_log10 = a_range_slope,
          t_end = t_end,
          no_lhs_samples = round(no_lhs_samples/length(range_multiplier))
        )
      }

    start_parms <- start_parms %>%
      dplyr::filter(nll == min(nll))

    fit[[witer]] <- bbmle::mle2(
      minuslogl = rrpe_nll_mod15h,
      start = list(
        t_h_0_log10  =  start_parms$t_h_0_log10,
        t_h_1_log10  =  start_parms$t_h_1_log10,
        t_h_2_log10  =  start_parms$t_h_2_log10,
        t_h_3_log10  =  start_parms$t_h_3_log10,
        t_h_4_log10  =  start_parms$t_h_4_log10,
        a_intercept_log10 = start_parms$a_intercept_log10,
        a_slope_log10 = start_parms$a_slope_log10
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
    message(paste("t_h_0_log10: ",  round(bbmle::coef(fit[[witer]])[2],2)))
    message(paste("t_h_1_log10: ",  round(bbmle::coef(fit[[witer]])[3],2)))
    message(paste("t_h_2_log10: ",  round(bbmle::coef(fit[[witer]])[4],2)))
    message(paste("t_h_3_log10: ",  round(bbmle::coef(fit[[witer]])[5],2)))
    message(paste("t_h_4_log10: ",  round(bbmle::coef(fit[[witer]])[6],2)))
    message(paste("a_intercept_log10: ", round(bbmle::coef(fit[[witer]])[7],2)))
    message(paste("a_slope_log10: "    , round(bbmle::coef(fit[[witer]])[8],2)))
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
