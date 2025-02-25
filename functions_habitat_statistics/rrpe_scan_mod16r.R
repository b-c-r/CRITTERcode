################################################################################
#   rrpe_parms_scan_mod16r: scans for reasonable starting parameters           #
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
#' @description `rrpe_parms_scan_mod16r` creates Latin hypercube samples for the
#'     functional response parameters in a reasonable range and calculates the
#'     according negative log-likelihood values. It returns the parameter values
#'     with the lowest negative log likelihood of these samples. Non-linear
#'     maximum likelihood fitting procedures require starting parameters,
#'     generally based on an educated guess (e.g., Bolker 2008). Moreover,
#'     these fits may end up in local best fits, and users should re-fit the
#'     data using different starting parameters (Bolker 2008). To overcome
#'     manually eyeballing as well as re-shuffling the starting parameters,
#'     Jager and Ashauer (2018) suggested creating samples in a reasonable
#'     parameter range using and choosing the starting parameters (from the
#'     lowest nll value) from these samples. To reduce the number of required
#'     samples by keeping the variance of parameter values as wide as possible,
#'     we use Latin hypercube sampling. Moreover, the function requires the
#'     model parameters on log-scale, as this transformation (1) accelerates the
#'     fitting procedure and (2) prevents biologically irrelevant negative
#'     estimations that would crash the fitting algorithm.
#'     `rrpe_parms_scan_mod16r` requires the lhs package (Carnell 2024).
#' 
#'     Required packages and their dependencies to be installed:
#'       - `emdbook` (Bolker 2023)
#'       - `foreach` (Microsoft and Weston 2022)
#'       - `lhs` (Carnell 2024)
#'     Required packages to be attached:
#'       - `foreach` (Microsoft and Weston 2022)
#' 
#' @references Bolker (2008) Ecological models and data in R, Princeton
#'     University Press, Princeton, New Jersey.
#'     https://math.mcmaster.ca/~bolker/emdbook/index.html
#' @references Bolker (2023) emdbook: support functions and data for "Ecological
#'     models and data". Version 1.3.13.
#'     https://CRAN.R-project.org/package=emdbook
#' @references Carnell (2024) lhs: latin hypercube samples. Version 1.2.0.
#'     https://CRAN.R-project.org/package=lhs
#' @references Microsoft and Weston (2022) foreach: provides foreach looping
#'     construct. Version 1.5.2. https://doi.org/10.32614/CRAN.package.foreach
#' @references Jager and Ashauer (2018) Modelling survival under chemical stress
#'     Leanpub. https://leanpub.com/guts_book
#'
#' @include rrpe_sim.R
#' @include rrpe_nll_mod16r.R
#' 
#' @param n_eaten integer (or float); the prey items that were eaten throughout the experimental trial. A vector.
#' @param n_initial integer or float; a vector of initial prey densities.
#' @param complexity level of complexity (0-4), a single integer value.
#' @param p integer or float; a single value of a fixed predator density. The def
#' @param f_max_range_range_0_log10 A range for f_max for the respective complexity level, two values.
#' @param f_max_range_range_1_log10 A range for f_max for the respective complexity level, two values.
#' @param f_max_range_range_2_log10 A range for f_max for the respective complexity level, two values.
#' @param f_max_range_range_3_log10 A range for f_max for the respective complexity level, two values.
#' @param f_max_range_range_4_log10 A range for f_max for the respective complexity level, two values.
#' @param n_half_range_range_0_log10 A range for n_half for the respective complexity level, two values.
#' @param n_half_range_range_1_log10 A range for n_half for the respective complexity level, two values.
#' @param n_half_range_range_2_log10 A range for n_half for the respective complexity level, two values.
#' @param n_half_range_range_3_log10 A range for n_half for the respective complexity level, two values.
#' @param n_half_range_range_4_log10 A range for n_half for the respective complexity level, two values.
#' @param t_end integer or float; the time were the feeding ends. A single value; default = 1 (e.g. 1 day).
#' @param no_lhs_samples a single integer value; the number of random latin hypercube samplings.
#' 
#' @return Returns a data frame with a single row of parameter values.
#' 
#' @examples
#' 
#' rm(list=ls())
#' 
#' library("foreach")
#' 
#' gh_path <- "https://raw.githubusercontent.com/b-c-r/CRITTERcode/refs/heads/main/functions_habitat_statistics/"
#' 
#' source(paste(gh_path, "rrpe_sim.R", sep = ""))
#' source(paste(gh_path, "rrpe_nll_mod16h.R", sep = ""))
#' source(paste(gh_path, "rrpe_scan_mod16r.R", sep = ""))
#' 
#' fr_data <- read.csv("https://raw.githubusercontent.com/b-c-r/CRITTERdata/refs/heads/main/critter_data.csv")
#' fr_data_ie <- subset(fr_data, predator == "Ischnura elegans")
#' 
#' rrpe_scan_mod16r(
#'   n_eaten = fr_data_ie$n_eaten,
#'   n_initial = fr_data_ie$n_initial,
#'   complexity  = fr_data_ie$complexity_level,
#'   p = 1,
#'   f_max_range_0_log10 = log10(c(1, max(fr_data_ie$n_eaten[fr_data_ie$complexity_level == 0]))),
#'   f_max_range_1_log10 = log10(c(1, max(fr_data_ie$n_eaten[fr_data_ie$complexity_level == 1]))),
#'   f_max_range_2_log10 = log10(c(1, max(fr_data_ie$n_eaten[fr_data_ie$complexity_level == 2]))),
#'   f_max_range_3_log10 = log10(c(1, max(fr_data_ie$n_eaten[fr_data_ie$complexity_level == 3]))),
#'   f_max_range_4_log10 = log10(c(1, max(fr_data_ie$n_eaten[fr_data_ie$complexity_level == 4]))),
#'   n_half_range_0_log10 = log10(c(1, max(fr_data_ie$n_initial[fr_data_ie$complexity_level == 0]))),
#'   n_half_range_1_log10 = log10(c(1, max(fr_data_ie$n_initial[fr_data_ie$complexity_level == 1]))),
#'   n_half_range_2_log10 = log10(c(1, max(fr_data_ie$n_initial[fr_data_ie$complexity_level == 2]))),
#'   n_half_range_3_log10 = log10(c(1, max(fr_data_ie$n_initial[fr_data_ie$complexity_level == 3]))),
#'   n_half_range_4_log10 = log10(c(1, max(fr_data_ie$n_initial[fr_data_ie$complexity_level == 4]))),
#'   t_end = 1,
#'   no_lhs_samples = 100
#' )
#' 
#' #############################################################################
#' 
#' fr_data_ng <- subset(fr_data, predator == "Notonecta glauca")
#' 
#' rrpe_scan_mod16r(
#'   n_eaten = fr_data_ng$n_eaten,
#'   n_initial = fr_data_ng$n_initial,
#'   complexity  = fr_data_ng$complexity_level,
#'   p = 1,
#'   f_max_range_0_log10 = log10(c(1, max(fr_data_ng$n_eaten[fr_data_ng$complexity_level == 0]))),
#'   f_max_range_1_log10 = log10(c(1, max(fr_data_ng$n_eaten[fr_data_ng$complexity_level == 1]))),
#'   f_max_range_2_log10 = log10(c(1, max(fr_data_ng$n_eaten[fr_data_ng$complexity_level == 2]))),
#'   f_max_range_3_log10 = log10(c(1, max(fr_data_ng$n_eaten[fr_data_ng$complexity_level == 3]))),
#'   f_max_range_4_log10 = log10(c(1, max(fr_data_ng$n_eaten[fr_data_ng$complexity_level == 4]))),
#'   n_half_range_0_log10 = log10(c(1, max(fr_data_ng$n_initial[fr_data_ng$complexity_level == 0]))),
#'   n_half_range_1_log10 = log10(c(1, max(fr_data_ng$n_initial[fr_data_ng$complexity_level == 1]))),
#'   n_half_range_2_log10 = log10(c(1, max(fr_data_ng$n_initial[fr_data_ng$complexity_level == 2]))),
#'   n_half_range_3_log10 = log10(c(1, max(fr_data_ng$n_initial[fr_data_ng$complexity_level == 3]))),
#'   n_half_range_4_log10 = log10(c(1, max(fr_data_ng$n_initial[fr_data_ng$complexity_level == 4]))),
#'   t_end = 1,
#'   no_lhs_samples = 100
#' )
#' 

rrpe_scan_mod16r <- function(
    n_eaten,
    n_initial,
    complexity,
    p = 1,
    f_max_range_0_log10,
    f_max_range_1_log10,
    f_max_range_2_log10,
    f_max_range_3_log10,
    f_max_range_4_log10,
    n_half_range_0_log10,
    n_half_range_1_log10,
    n_half_range_2_log10,
    n_half_range_3_log10,
    n_half_range_4_log10,
    t_end = 1,
    no_lhs_samples = 1000
){
  
  # generation of 
  lhsvals <- lhs::randomLHS(no_lhs_samples, 10)
  
  f_max_range_0_log10   <- (lhsvals[, 1] * ( f_max_range_0_log10[2] -  f_max_range_0_log10[1])) +  f_max_range_0_log10[1]
  f_max_range_1_log10   <- (lhsvals[, 2] * ( f_max_range_1_log10[2] -  f_max_range_1_log10[1])) +  f_max_range_1_log10[1]
  f_max_range_2_log10   <- (lhsvals[, 3] * ( f_max_range_2_log10[2] -  f_max_range_2_log10[1])) +  f_max_range_2_log10[1]
  f_max_range_3_log10   <- (lhsvals[, 4] * ( f_max_range_3_log10[2] -  f_max_range_3_log10[1])) +  f_max_range_3_log10[1]
  f_max_range_4_log10   <- (lhsvals[, 5] * ( f_max_range_4_log10[2] -  f_max_range_4_log10[1])) +  f_max_range_4_log10[1]
  n_half_range_0_log10  <- (lhsvals[, 6] * (n_half_range_0_log10[2] - n_half_range_0_log10[1])) + n_half_range_0_log10[1]
  n_half_range_1_log10  <- (lhsvals[, 7] * (n_half_range_1_log10[2] - n_half_range_1_log10[1])) + n_half_range_1_log10[1]
  n_half_range_2_log10  <- (lhsvals[, 8] * (n_half_range_2_log10[2] - n_half_range_2_log10[1])) + n_half_range_2_log10[1]
  n_half_range_3_log10  <- (lhsvals[, 9] * (n_half_range_3_log10[2] - n_half_range_3_log10[1])) + n_half_range_3_log10[1]
  n_half_range_4_log10  <- (lhsvals[,10] * (n_half_range_4_log10[2] - n_half_range_4_log10[1])) + n_half_range_4_log10[1]
  
  ## calculate nlls
  nlls <- foreach::foreach(
    i = 1:no_lhs_samples,
    .combine = "c") %do% {
      
      nll <- rrpe_nll_mod16r(
          n_eaten = n_eaten,
          n_initial = n_initial,
          complexity = complexity,
          p = p,
          f_max_0_log10 = f_max_range_0_log10[i],
          f_max_1_log10 = f_max_range_1_log10[i],
          f_max_2_log10 = f_max_range_2_log10[i],
          f_max_3_log10 = f_max_range_3_log10[i],
          f_max_4_log10 = f_max_range_4_log10[i],
          n_half_0_log10 = n_half_range_0_log10[i],
          n_half_1_log10 = n_half_range_1_log10[i],
          n_half_2_log10 = n_half_range_2_log10[i],
          n_half_3_log10 = n_half_range_3_log10[i],
          n_half_4_log10 = n_half_range_4_log10[i],
          t_end = t_end
        )
      
    }
  
  sel_parms <- data.frame(
    f_max_0_log10  =  f_max_range_0_log10[nlls == min(nlls)],
    f_max_1_log10  =  f_max_range_1_log10[nlls == min(nlls)],
    f_max_2_log10  =  f_max_range_2_log10[nlls == min(nlls)],
    f_max_3_log10  =  f_max_range_3_log10[nlls == min(nlls)],
    f_max_4_log10  =  f_max_range_4_log10[nlls == min(nlls)],
    n_half_0_log10 = n_half_range_0_log10[nlls == min(nlls)],
    n_half_1_log10 = n_half_range_1_log10[nlls == min(nlls)],
    n_half_2_log10 = n_half_range_2_log10[nlls == min(nlls)],
    n_half_3_log10 = n_half_range_3_log10[nlls == min(nlls)],
    n_half_4_log10 = n_half_range_4_log10[nlls == min(nlls)],
    nll = nlls[nlls == min(nlls)]
  )
  
  return(sel_parms)
}
