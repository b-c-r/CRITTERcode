################################################################################
#   mod03r_rrpe_nll: estimates the negative log-likelihood                     #
#                                                                              #
#    Copyright (C) 2025                                                        #
#       Björn C. Rall (https://orcid.org/0000-0002-3191-8389)                  #
#       Mireia Aranbarri (https://orcid.org/0009-0001-3506-0914)               #
#       Lorea Flores (https://orcid.org/0000-0002-0082-4072)                   #
#       Ioar de Guzmán (https://orcid.org/0000-0001-8894-8477)                 #
#       Aitor Larrañaga (https://orcid.org/0000-0002-0185-9154)                #
#       Julia Reiss (https://orcid.org/0000-0002-3740-0046)                    #
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
#' 
#' @return Returns a single negative log-likelihood value.
#' 
#' @examples
#' 
#' # find an executable example here:
#' # https://github.com/b-c-r/CRITTERcode/examples_habitat_statistics/examples_habitat_statistics/mod03r_examples.R
#' 

mod03r_rrpe_nll <- function(
    n_eaten,
    n_initial,
    n_rings,
    complexity,
    p = 1,
    f_max_log10,
    n_half_intercept_log10,
    n_half_slope,
    t_end = 1
){
  
  nll <- Inf
  
  try({
    eaten_simulated <- foreach::foreach(
      i = 1:length(n_initial),
      .combine = "rbind"
    ) %do% {
      
      rrpe_sim(
        fr_style = "Real",
        n_initial = n_initial[i],
        p = p,
        f_max = 10^(f_max_log10),
        n_half = 10^(n_half_intercept_log10 + n_half_slope * n_rings[i]),
        t_end = t_end
      )
    }
    
    lls <- dbinom(
      x = n_eaten,
      size = n_initial,
      prob = eaten_simulated$n_eaten/n_initial,
      log = T
    )
    
    nll <- -1*sum(lls)
    
  }, silent = TRUE)
  
  return(nll)
}
