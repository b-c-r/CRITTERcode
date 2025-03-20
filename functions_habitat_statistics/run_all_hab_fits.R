################################################################################
#   run_all_hab_fits: fits all habitat models to all data                      #
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
#' @return Returns a single negative log-likelihood value.
#' 
#' @examples
#' 
#' # find an executable example here:
#' # https://github.com/b-c-r/CRITTERcode/examples_habitat_statistics/examples_habitat_statistics/mod15r_examples.R
#' 

run_all_hab_fits <- function(
    x,
    model_numbers = rep(1:16, 4),
    style_letters = c(rep("h",16), rep("r",16), rep("h",16), rep("r",16)),
    predator_spec = c(rep("Ischnura elegans", 32), rep("Notonecta glauca", 32)),
    no_threads = max(2^(0:6)[2^(0:6) <= parallel::detectCores()]),
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
  
  library("foreach")
  library("dplyr")
  
  cl <- parallel::makeCluster(no_threads)
  doParallel::registerDoParallel(cl)
  
  result <- foreach::foreach(
    i = 1:length(model_numbers),
    .packages = c("foreach", "dplyr"),
    .export = ls(globalenv())
  ) %dopar% {
    
    if(model_numbers[i] < 10){
      func_name <- paste("mod0", model_numbers[i], style_letters[i], "_rrpe_fit", sep ="")
    } else {
      func_name <- paste("mod", model_numbers[i], style_letters[i], "_rrpe_fit", sep ="")
    }
    
    fr_data_i <- subset(x, predator == predator_spec[i])
    
    out <- do.call(
      what = func_name,
      args = list(
        n_eaten = fr_data_i$n_eaten,
        n_initial = fr_data_i$n_initial,
        n_rings = fr_data_i$ring_count,
        complexity = fr_data_i$complexity_level,
        p = p,
        t_end = t_end,
        no_lhs_samples = no_lhs_samples,
        range_multiplier = range_multiplier,
        rel_f_max_range = rel_f_max_range,
        rel_n_half_range = rel_n_half_range,
        slope_range = slope_range,
        witer_max = witer_max,
        mle2_tol = mle2_tol,
        val_tol = val_tol,
        set_seed = set_seed,
        seed_value = seed_value
      )
    )
    
    # create the output:
    out_table <- list(
      model_name = func_name,
      predator_spec = predator_spec[i],
      mle2_object = out
    )
    
    return(out_table)
  }
  
  parallel::stopCluster(cl)
  return(result)
  
}
