################################################################################
#   gen_fr_nll: estimates the negative log-likelihood of the gen. FR model     #
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
#'
#' @return Returns a single negative log-likelihood value.
#'
#' @examples
#' 
#' # find an executable example here:
#' # https://github.com/b-c-r/CRITTERcode/examples_habitat_statistics/examples_habitat_statistics/gen_fr_nll_examples.R
#' 

gen_fr_nll <- function(
    n_eaten,
    n_initial,
    p = 1,
    f_max_log10,
    n_half_log10,
    q,
    t_start = 0,
    t_end = 1,
    t_length = 1000,
    penalty = 1000,
    q_low = 0,
    q_up = 1
){
  
  nll <- Inf
  
  try({
    eaten_simulated <- gen_fr_sim(
      n_initial = n_initial,
      p = p,
      f_max = 10^f_max_log10,
      n_half = 10^n_half_log10,
      q = q,
      t_start = t_start,
      t_end = t_end,
      t_length = t_length
    )
    
    lls <- dbinom(
      x = n_eaten,
      size = n_initial,
      prob = eaten_simulated$n_eaten/n_initial,
      log = T
    )
    
    nll <- -1*sum(lls)
    
  }, silent = TRUE)
  
  if(nll == Inf){return(nll)}
  
  # penalty if q < 0 or q > 1 (changeable in function header!)
  if(q < q_low){
    nll <- nll + penalty*(q-q_low)^2
  } else{
    if(q >= q_up){
      nll <- nll + penalty*(q-q_up)^2
    } else{
      nll <- nll
    }
  }
  
  return(nll)
}
