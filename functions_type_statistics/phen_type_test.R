################################################################################
#    phen_type_test: phenomenological test for functional response type        #
#                                                                              #
#    Copyright (C) 2025                                                        #
#       Björn C. Rall (https://orcid.org/0000-0002-3191-8389)                  #
#       Mireia Aranbarri (https://orcid.org/0009-0001-3506-0914)               #
#       Ioar de Guzmán (https://orcid.org/0000-0001-8894-8477)                 #
#       Aitor Larrañaga (https://orcid.org/0000-0002-0185-9154)                #
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
#' @return returns a list of frair test results
#'
#' @examples
#' 
#' # find an executable example here:
#' # https://github.com/b-c-r/CRITTERcode/examples_habitat_statistics/examples_habitat_statistics/phen_type_test_examples.R
#' 

phen_type_test <- function(
    data,
    name_initial,
    name_eaten,
    name_treatments
    ){
  
  # select the required data:
  data_internal <- data %>%
    dplyr::select(all_of(c(name_initial, name_eaten, name_treatments))) %>%
    purrr::set_names(c("n_initial", "n_eaten", "treatment"))
  
  # extract unique treatments for the for-loop and sort data for ordered output:
  treats <- sort(unique(data_internal$treatment))
  
  # repeat the frair_test() across treatments and save it with data attached:
  foreach::foreach(i = 1:length(treats)) %do% {
    data_i <- subset(data_internal, treatment == treats[i])
    data_orig_i <- subset(data, treatment == treats[i])
    table <- list(
      treatment = treats[i],
      frair_test_results = frair::frair_test(n_eaten ~ n_initial, data = data_i),
      data_orig = data_orig_i
    )
  }
}
