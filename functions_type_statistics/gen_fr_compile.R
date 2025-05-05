################################################################################
#    gen_fr_compile: compiles a functional response ODE using odin             #
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
#' @return the function compiles the functional response model and returns a
#'     function called `gen_fr_model` that is needed to simulate the time series
#'     of decaying resource items.
#'    
#' @examples
#' 
#' # find an executable example here:
#' # https://github.com/b-c-r/CRITTERcode/examples_habitat_statistics/examples_habitat_statistics/gen_fr_compile_examples.R
#' 

gen_fr_compile <- function(){
  gen_fr_model <<- odin::odin({
    deriv(n) <- -f_max * n^(1+q) / (n_half^(1+q) + n^(1+q)) * p                 # the ODE of the generalized functional response
    initial(n) <- n_initial                                                     # assign initial prey density
    
    # all parameters must be filled with a value, "user()" means that this
    # information can be added later and is not hard coded
    n_initial <- user()
    f_max <- user()
    n_half <- user()
    q <- user()
    p <- user()
  })
}
