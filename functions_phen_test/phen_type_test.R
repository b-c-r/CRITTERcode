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
#' @description `phen_type_test` is a wrapper around `frair_test` Pritchard et
#'     al. (2017a, b) to test for the type of the functional response using the
#'     proportion of prey eaten as response to the initial prey density. If the
#'     proportion of prey eaten increases at low prey densities and decreases
#'     after reaching a maximum, there is evidence for a type III functional
#'     response (see Juliano 2001). Proportion data is fitted using a binomial
#'     GLM, but see Crawley (2012), chapter 16 for an introduction to this
#'     topic.
#'     
#'     Required packages and their dependencies to be installed:
#'     - `dplyr` (Wickham et al. 2023)
#'     - `purrr` (Wickham and Henry 2025)
#'     - `foreach` (Microsoft and Weston 2022)
#'     - `frair` (Pritchard et al. 2017b)
#'     Required packages to be attached:
#'     - `dplyr` (Wickham et al. 2023)
#'     - `foreach` (Microsoft and Weston 2022)
#' 
#' @references Crawley (2012) The R Book. Wiley, Chichester, West Sussex, UK
#' @references Juliano (2001) Nonlinear curve fitting: predation and functional
#'     response curves. In: Scheiner, Gurevitch (eds.). Design and analysis of
#'     ecological experiments. 178–196
#' @references Microsoft and Weston (2022) foreach: provides foreach looping
#'     construct. Version 1.5.2. https://doi.org/10.32614/CRAN.package.foreach
#' @references Pritchard et al. (2017a) frair: an R package for fitting and
#'     comparing consumer functional responses. Methods Ecol Evol 8, 1528-1534.
#'     https://doi.org/10.1111/2041-210X.12784
#' @references Pritchard et al. (2017b) frair: tools for functional response
#'     analysis. Version 0.5.100. https://cran.r-project.org/web/packages/frair/
#' @references Wickham et al. (2023) dplyr: a grammar of data manipulation.
#'     Version 1.1.4. https://doi.org/10.32614/CRAN.package.dplyr  
#' @references Wickham et al. (2023) purrr: functional programming tools.
#'     Version 1.0.4. https://doi.org/10.32614/CRAN.package.purrr 
#' 
#' @param data the input data
#' @param name_initial the column name
#' @param name_eaten
#' @param name_treatments
#'
#' @examples
#' 
#' library("foreach")
#' library("dplyr")
#' source(here::here("functions_phen_test", "phen_type_test.R"))
#' 
#' fr_data <- data.frame(
#'   n_start =   rep(c(1, 2, 4, 8, 16), 4),
#'   n_eaten =   c(0,0,2,7,8,1,2,4,6,7,1,1,3,4,4,1,1,4,4,4),
#'   treatment = c(rep("A", 5), rep("B", 5), rep("A", 5), rep("B", 5)),
#'   predator = c(rep("Avicularia avicularia", 10), rep("Cicindela gallica", 10))
#' )
#' 
#' phen_type_test(
#'   data = fr_data,
#'   name_initial = "n_start",
#'   name_eaten = "n_eaten",
#'   name_treatments = "treatment"
#' )
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
