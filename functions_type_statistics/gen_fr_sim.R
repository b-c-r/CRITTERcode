################################################################################
#    gen_fr_sim: simulates a generalized functional response time series       #
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
#'
#' @return Returns a data frame with n_initial and the according n_eaten.
#' 
#' @examples
#' 
#' # find an executable example here:
#' # https://github.com/b-c-r/CRITTERcode/examples_habitat_statistics/examples_habitat_statistics/gen_fr_sim_examples.R
#' 

gen_fr_sim <- function(
    n_initial,
    p = 1,
    f_max,
    n_half,
    q,
    t_start = 0,
    t_end = 1,
    t_length = 1000
){
  
  # foreach loop across all initial prey densities:
  n_eaten <- foreach::foreach(
    i = 1:length(n_initial),
    .combine = c
  ) %do% {
    
    # assign values to model parameters
    gen_fr_model_assigned <- gen_fr_model$new(
      n_initial = n_initial[i],
      f_max = f_max,
      n_half = n_half,
      q = q,
      p = p
    )
    
    # set time steps to be computed
    tt <- seq(
      t_start,
      t_end,
      length.out = t_length
    )
    
    # calculate the remaining prey items/density
    remaining_prey <- gen_fr_model_assigned$run(tt)[[t_length,2]]
    
    # calculate the number of prey eaten (start - final prey density):
    n_initial[i] - remaining_prey
  }
  
  # return the data frame
  data.frame(
    n_initial = n_initial,
    n_eaten = n_eaten
  )
}
