################################################################################
#   gen_fr_fit_all: fits the gen. FR model to all 10 treatments in parallel    #
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
#' @description `gen_fr_fit_all` fits the generalized functional response model
#'     (Real 1977, Rosenbaum and Rall 2018) by running `gen_fr_fit` for all
#'     treatments in parallel.
#'     
#'     Required packages and their dependencies to be installed:
#'       - `bbmle` (Bolker et al 2023)
#'       - `doParallel` (Microsoft and Weston 2022)
#'       - `dplyr` (Wickham et al. 2023)
#'       - `foreach` (Microsoft and Weston 2022)
#'       - `lhs` (Carnell 2024)
#'       - `odin` (FitzJohn and Jombart 2024)
#'       - `purrr` (Wickham and Henry 2025)
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
#' @references Microsoft and Weston (2022) doParallel: foreach parallel adaptor
#'     for the 'parallel' package. Version 1.0.17.
#'     https://CRAN.R-project.org/package=doParallel
#' @references Microsoft and Weston (2022) foreach: provides foreach looping
#'     construct. Version 1.5.2. https://CRAN.R-project.org/package=foreach
#' @references Real (1977) The kinetics of functional response. Am Nat 111, 289-
#'     300. https://doi.org/10.1086/283161
#' @references Rosenbaum and Rall (2018) Fitting functional responses: Direct
#'     parameter estimation by simulating differential equations. Methods Ecol
#'     Evol 9, 2076-2090. https://doi.org/10.1111/2041-210X.13039
#' @references Wickham et al. (2023) dplyr: a grammar of data manipulation.
#'     1.1.4. https://CRAN.R-project.org/package=dplyr
#' @references Wickham et al. (2023) purrr: functional programming tools.
#'     Version 1.0.4. https://doi.org/10.32614/CRAN.package.purrr 
#'
#' @include gen_fr_compile.R
#' @include gen_fr_sim.R
#' @include gen_fr_nll.R
#' @include gen_fr_parms_scan.R
#' @include gen_fr_fit.R
#'
#' @param data the input data - should be pipeable in the tidyverse (not tested)
#' @param name_initial the column name of the initial resource density. The data
#'     must be integer; the prey items at the start of the  experimental trial.
#' @param name_eaten the column name of the resource items eaten. The data
#'     must be integer; the prey items that were eaten throughout the
#'     experimental trial.
#' @param name_treatments the column name of the treatments column.
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
#' @param no_threads Number of threads that should be used to compute the tests
#'     in parallel. Default is to ten, as 10 tests should be computed. Most
#'     modern computers, even laptops should be able to run 10 thread in
#'     parallel. if you're unsure run parallel::detectCores().
#' 

gen_fr_fit_all <- function(
    data,
    name_initial,
    name_eaten,
    name_treatments,
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
    maxit = 5000,
    no_threads = max(c(1, 2, 5, 10)[c(1, 2, 5, 10) <= max_cores])
    ){
  
  # select the required data:
  data_internal <- data %>%
    dplyr::select(all_of(c(name_initial, name_eaten, name_treatments))) %>%
    purrr::set_names(c("n_initial", "n_eaten", "treatment"))
  
  # extract unique treatments for the for-loop and sort data for ordered output:
  treats <- sort(unique(data_internal$treatment))
  
  cl <- parallel::makeCluster(no_threads)
  doParallel::registerDoParallel(cl)
  
  result <- foreach::foreach(
    i = 1:length(treats),
    .packages = c("foreach", "dplyr"),
    .export = ls(globalenv())
    ) %dopar% {

    # subset the data:
    data_i <- subset(data_internal, treatment == treats[i])
    data_orig_i <- subset(data, treatment == treats[i])
    
    # run the fitting function:
    fr_res <- gen_fr_fit(
      n_eaten = data_i$n_eaten,
      n_initial = data_i$n_initial,
      p = p,
      t_start = t_start,
      t_end = t_end,
      t_length = t_length,
      penalty = penalty,
      q_low = q_low,
      q_up = q_up,
      no_lhs_samples = no_lhs_samples,
      range_multiplier = range_multiplier,
      witer_max = witer_max,
      val_tol = val_tol,
      mle2_tol = mle2_tol,
      maxit = maxit
    )
    
    # create the output:
    out_table <- list(
      treatment = treats[i],
      q_test_results = fr_res#,
      # data_orig = data_orig_i
    )
    
    return(out_table)
  }
  parallel::stopCluster(cl)
  return(result)
}
