################################################################################
#   rrpe_parms_scan_mod12h: scans for reasonable starting parameters           #
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
#' @include rrpe_sim.R
#' @include rrpe_nll_mod12h.R
#' 
#' @return returns a data frame with a single row of parameter values (lowest nll).
#' 
#' @examples
#' 
#' # find an executable example here:
#' # https://github.com/b-c-r/CRITTERcode/examples_habitat_statistics/examples_habitat_statistics/mod12h_examples.R
#' 

rrpe_scan_mod12h <- function(
    n_eaten,
    n_initial,
    n_rings,
    complexity,
    p = 1,
    t_h_range_intercept_log10,
    t_h_range_slope_log10,
    a_range_0_log10,
    a_range_1_log10,
    a_range_2_log10,
    a_range_3_log10,
    a_range_4_log10,
    t_end = 1,
    no_lhs_samples = 1000
){
  
  # generation of 
  lhsvals <- lhs::randomLHS(no_lhs_samples, 7)
  
  t_h_range_intercept_log10 <- log10((lhsvals[, 1] * (10^t_h_range_intercept_log10[2] - 10^t_h_range_intercept_log10[1])) + 10^t_h_range_intercept_log10[1])
  t_h_range_slope_log10     <-       (lhsvals[, 2] * (       t_h_range_slope_log10[2] -        t_h_range_slope_log10[1])) +        t_h_range_slope_log10[1]
  a_range_0_log10           <- log10((lhsvals[, 3] * (          10^a_range_0_log10[2] -           10^a_range_0_log10[1])) +           10^a_range_0_log10[1])
  a_range_1_log10           <- log10((lhsvals[, 4] * (          10^a_range_1_log10[2] -           10^a_range_1_log10[1])) +           10^a_range_1_log10[1])
  a_range_2_log10           <- log10((lhsvals[, 5] * (          10^a_range_2_log10[2] -           10^a_range_2_log10[1])) +           10^a_range_2_log10[1])
  a_range_3_log10           <- log10((lhsvals[, 6] * (          10^a_range_3_log10[2] -           10^a_range_3_log10[1])) +           10^a_range_3_log10[1])
  a_range_4_log10           <- log10((lhsvals[, 7] * (          10^a_range_4_log10[2] -           10^a_range_4_log10[1])) +           10^a_range_4_log10[1])
  
  ## calculate nlls
  nlls <- foreach::foreach(
    i = 1:no_lhs_samples,
    .combine = "c") %do% {
      
      nll <- rrpe_nll_mod12h(
          n_eaten = n_eaten,
          n_initial = n_initial,
          n_rings = n_rings,
          complexity = complexity,
          p = p,
          t_h_intercept_log10 = t_h_range_intercept_log10[i],
          t_h_slope_log10 = t_h_range_slope_log10[i],
          a_0_log10 = a_range_0_log10[i],
          a_1_log10 = a_range_1_log10[i],
          a_2_log10 = a_range_2_log10[i],
          a_3_log10 = a_range_3_log10[i],
          a_4_log10 = a_range_4_log10[i],
          t_end = t_end
        )
      
    }
  
  sel_parms <- data.frame(
    t_h_intercept_log10  =  t_h_range_intercept_log10[nlls == min(nlls)],
    t_h_slope_log10  =  t_h_range_slope_log10[nlls == min(nlls)],
    a_0_log10 = a_range_0_log10[nlls == min(nlls)],
    a_1_log10 = a_range_1_log10[nlls == min(nlls)],
    a_2_log10 = a_range_2_log10[nlls == min(nlls)],
    a_3_log10 = a_range_3_log10[nlls == min(nlls)],
    a_4_log10 = a_range_4_log10[nlls == min(nlls)],
    nll = nlls[nlls == min(nlls)]
  )
  
  return(sel_parms)
}
