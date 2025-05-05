################################################################################
#   gen_fr_parms_scan: scans for reasonable starting parameters before fitting #
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
#'     
#' @return Returns a data frame with a single row of parameter values.
#' 
#' @examples
#' 
#' # find an executable example here:
#' # https://github.com/b-c-r/CRITTERcode/examples_habitat_statistics/examples_habitat_statistics/gen_fr_parms_scan_examples.R
#' 

gen_fr_parms_scan <- function(
    n_eaten,
    n_initial,
    p = 1,
    f_max_range_log10,
    n_half_range_log10,
    q_range = c(0,1),
    t_start = 0,
    t_end = 1,
    t_length = 1000,
    penalty = 1000,
    q_low = 0,
    q_up = 1,
    no_lhs_samples = 1000
){
  
  lhsvals <- lhs::randomLHS(no_lhs_samples, 3)
  
  f_max_range_log10  <- log10((lhsvals[,1] * (10^f_max_range_log10[2]  - 10^f_max_range_log10[1] )) + 10^f_max_range_log10[1])
  n_half_range_log10 <- log10((lhsvals[,2] * (10^n_half_range_log10[2] - 10^n_half_range_log10[1])) + 10^n_half_range_log10[1])
  q_range            <- (lhsvals[,3] * (q_range[2]-q_range[1])) + q_range[1]
  
  ## calculate nlls
  nlls <- foreach::foreach(
    i = 1:no_lhs_samples,
    .combine = "c") %do% {
      
      nll <- gen_fr_nll(
        n_eaten = n_eaten,
        n_initial = n_initial,
        p = p,
        f_max_log10 = f_max_range_log10[i],
        n_half_log10 = n_half_range_log10[i],
        q = q_range[i],
        t_start = t_start,
        t_end = t_end ,
        t_length = t_length,
        penalty = penalty,
        q_low = q_low,
        q_up = q_up
      )
      
      return(nll)
    }
  
  sel_parms <- data.frame(
    f_max_log10 = f_max_range_log10[nlls == min(nlls)],
    n_half_log10 = n_half_range_log10[nlls == min(nlls)],
    q = q_range[nlls == min(nlls)],
    nll = nlls[nlls == min(nlls)]
  )
  
  return(sel_parms)
}
