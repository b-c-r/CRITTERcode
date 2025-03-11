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
#' @description `gen_fr_parms_scan` creates Latin hypercube samples for the
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
#'     it is recommended to use Latin hypercube sampling. `gen_fr_parms_scan`
#'     requires the lhs package (Carnell 2024). See also the description of
#'     `gen_fr_nll` for further information.
#' 
#'     Required packages and their dependencies to be installed:
#'       - `foreach` (Microsoft and Weston 2022)
#'       - `lhs` (Carnell 2024)
#'       - `odin` (FitzJohn and Jombart 2024)
#'     Required packages to be attached:
#'       - `foreach` (Microsoft and Weston 2022)
#' 
#' @references Bolker (2008) Ecological models and data in R, Princeton
#'     University Press, Princeton, New Jersey.
#'     https://math.mcmaster.ca/~bolker/emdbook/index.html
#' @references Carnell (2024) lhs: latin hypercube samples. Version 1.2.0.
#'     https://CRAN.R-project.org/package=lhs
#' @references FitzJohn and Jombart (2024) odin: ODE generation and integration.
#'     Ver. 1.2.6. https://doi.org/10.32614/CRAN.package.odin
#'     see also: https://github.com/mrc-ide/odin
#' @references Microsoft and Weston (2022) foreach: provides foreach looping
#'     construct. Version 1.5.2. https://doi.org/10.32614/CRAN.package.foreach
#' @references Jager and Ashauer (2018) Modelling survival under chemical stress
#'     Leanpub. https://leanpub.com/guts_book
#' 
#' @include gen_fr_compile.R
#' @include gen_fr_sim.R
#' @include gen_fr_nll.R
#' 
#' @param n_eaten integer (or float); the prey items that were eaten throughout
#'     the experimental trial. A vector.
#' @param n_initial integer (or float); the initial prey density. A vector of
#'     the same length as n_eaten.
#' @param p The predator density. A single value
#' @param f_max_range_log10 float; a range (2 values) of the log10 of the
#'     maximum feeding rate.
#' @param n_half_range_log10 float; a range (2 values) of the log10 of the half
#'     saturation density.
#' @param q_range float; shape parameter, a range (2 values). A strict type II
#'     functional has q = 0, a strict type III functional response has q = 1.
#'     The values should match the values q_low and q_up below.
#'     Default is c(0,1).
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
#'     
#' @return Returns a data frame with a single row of parameter values.
#' 
#' @examples
#' 
#' 
#' library("foreach")
#' 
#' source(here::here("functions_gen_fr", "gen_fr_compile.R"))
#' source(here::here("functions_gen_fr", "gen_fr_sim.R"))
#' source(here::here("functions_gen_fr", "gen_fr_nll.R"))
#' source(here::here("functions_gen_fr", "gen_fr_parms_scan.R"))
#' 
#' gen_fr_compile()
#' 
#' fr_data <- read.csv("data/fr_data.csv")
#' treats <- sort(unique(fr_data$treatment))
#' fr_data_treat <- subset(fr_data, treatment == treats[1])
#' 
#' gen_fr_parms_scan(
#'   n_eaten = fr_data_treat$n_eaten,
#'   n_initial = fr_data_treat$n_initial,
#'   p = 1,
#'   f_max_range_log10 = log10(c(1, max(fr_data_treat$n_eaten))),
#'   n_half_range_log10 = log10(c(1, max(fr_data_treat$n_initial))),
#'   q_range = c(0, 1),
#'   t_start = 0,
#'   t_end = 1,
#'   t_length = 100,
#'   penalty = 1000,
#'   q_low = 0,
#'   q_up = 1,
#'   no_lhs_samples = 1000
#' )
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
