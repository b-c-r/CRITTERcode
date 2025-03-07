################################################################################
#   mod06r_rrpe_scan: scans for reasonable starting parameters                 #
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
#' @include mod06r_rrpe_nll.R
#' 
#' @return Returns a single negative log-likelihood value.
#' 
#' @examples
#' 
#' # find an executable example here:
#' # https://github.com/b-c-r/CRITTERcode/examples_habitat_statistics/examples_habitat_statistics/mod06r_examples.R
#' 

mod06r_rrpe_scan <- function(
    n_eaten,
    n_initial,
    complexity,
    p = 1,
    f_max_hab0_log10_range,
    f_max_hab1_log10_range,
    n_half_hab0_log10_range,
    n_half_hab1_log10_range,
    t_end = 1,
    no_lhs_samples = 1000
){
  
  # generation of 
  lhsvals <- lhs::randomLHS(no_lhs_samples, 4)
  
  f_max_hab0_log10_range  <- log10((lhsvals[, 1] * ( 10^f_max_hab0_log10_range[2] -  10^f_max_hab0_log10_range[1])) +  10^f_max_hab0_log10_range[1])
  f_max_hab1_log10_range  <- log10((lhsvals[, 2] * ( 10^f_max_hab1_log10_range[2] -  10^f_max_hab1_log10_range[1])) +  10^f_max_hab1_log10_range[1])
  n_half_hab0_log10_range <- log10((lhsvals[, 3] * (10^n_half_hab0_log10_range[2] - 10^n_half_hab0_log10_range[1])) + 10^n_half_hab0_log10_range[1])
  n_half_hab1_log10_range <- log10((lhsvals[, 4] * (10^n_half_hab1_log10_range[2] - 10^n_half_hab1_log10_range[1])) + 10^n_half_hab1_log10_range[1])
  
  ## calculate nlls
  nlls <- foreach::foreach(
    i = 1:no_lhs_samples,
    .combine = "c") %do% {
      
      nll <- mod06r_rrpe_nll(
          n_eaten = n_eaten,
          n_initial = n_initial,
          complexity = complexity,
          p = p,
          f_max_hab0_log10 = f_max_hab0_log10_range[i],
          f_max_hab1_log10 = f_max_hab1_log10_range[i],
          n_half_hab0_log10 = n_half_hab0_log10_range[i],
          n_half_hab1_log10 = n_half_hab1_log10_range[i],
          t_end = t_end
        )
      
    }
  
  sel_parms <- data.frame(
    f_max_hab0_log10  = f_max_hab0_log10_range[nlls == min(nlls)],
    f_max_hab1_log10  = f_max_hab1_log10_range[nlls == min(nlls)],
    n_half_hab0_log10 = n_half_hab0_log10_range[nlls == min(nlls)],
    n_half_hab1_log10 = n_half_hab1_log10_range[nlls == min(nlls)],
    nll = nlls[nlls == min(nlls)]
  )
  
  return(sel_parms)
}
