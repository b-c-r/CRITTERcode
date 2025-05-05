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
#' @include gen_fr_fit.R
#' 
#' @return returns a list of fitting results and data.
#'
#' # find an example in our report:
#' # https://github.com/b-c-r/CRITTERstatistics/blob/main/statisticsReport.Rmd
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
    no_threads = max(c(1, 2, 5, 10)[c(1, 2, 5, 10) <= parallel::detectCores()])
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
    i = 1:length(treats)#,
    #.packages = export_packages_to_workers,
    #.export = export_functions_to_workers
    ) %dopar% {
      
      library("dplyr")
      library("foreach")
      
      # import the functions from GitHub:
      gh_path <- "https://raw.githubusercontent.com/b-c-r/CRITTERcode/refs/heads/main/" 
      f_path <- "functions_type_statistics/"
      source(paste(gh_path, f_path, "gen_fr_compile.R", sep = ""))
      source(paste(gh_path, f_path, "gen_fr_sim.R", sep = ""))
      source(paste(gh_path, f_path, "gen_fr_nll.R", sep = ""))
      source(paste(gh_path, f_path, "gen_fr_parms_scan.R", sep = ""))
      source(paste(gh_path, f_path, "gen_fr_fit.R", sep = ""))

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
        q_test_results = fr_res,
        data_orig = data_orig_i
      )
      
      return(out_table)
  }
  parallel::stopCluster(cl)
  return(result)
}
