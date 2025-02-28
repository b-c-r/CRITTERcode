################################################################################
#    rrpe_sim: simulates a type II RRPE functional response                    #
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
#' @return returns a data frame with simulated n_initial and n_eaten.
#' 
#' @examples
#'  
#' # find executable examples here:
#' # https://github.com/b-c-r/CRITTERcode/examples_habitat_statistics/examples_habitat_statistics/rrpe_sim_examples.R
#' 

rrpe_sim <- function(
    fr_style = "Real",
    n_initial,
    p = 1,
    a = 0,
    t_h = 0,
    f_max = 0,
    n_half = 0,
    t_end = 1
    ){
  
  if(fr_style == "Real"){

    remaining_prey <- n_half * emdbook::lambertW( n_initial/n_half * exp((n_initial - f_max * p * t_end)/n_half))
    
    else{
      if(fr_style == "Holling"){
        
        remaining_prey <- emdbook::lambertW(a * t_h * n_initial * exp(a * (t_h * n_initial - p * t_end))) / (a * t_h)

      } else{
        
        stop("please provide a valid functional response style: Real or Holling")
        
      }
    }
  }
  
  n_eaten <- n_initial - remaining_prey
  
  data.frame(
    n_initial = n_initial,
    n_eaten = n_eaten
  )

}

 