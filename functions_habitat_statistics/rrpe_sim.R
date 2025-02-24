################################################################################
#    rrpe_sim: simulates a type II RRPE functional response (Real-style)       #
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
#' @description `rrpe_sim` runs the Real-style Type II functional (Real 1977,
#'     1979). The feeding rate, F, is determined by the model parameters maximum
#'     feeding rate, F_max, and half saturation density, N_half:
#'         
#'     F = F_max * N / (N_half + N),
#'     
#'     where N is the resource density. However, this function requires a
#'     constant resource density, which is not given in most functional response
#'     experiments. To take the temporal decline of the resource into account, 
#'     we apply Rogers' Random Equation (Royama 1971, Rogers 1972):
#'     
#'     N_eaten = N_initial * (1 - exp(F_max / N_half * (N_eaten / F_max - P * T_end))),
#'     
#'     where N_eaten are the eaten resource items, N_initial is the initial
#'     resource density, P is the predator density, and T_end is the
#'     experimental duration (time). This function contains N_eaten on both
#'     sides of the equation and is only solvable iteratively (Juliano 2001,
#'     Vonesh and Bolker 2005). Bolker (2008) solved this issue by using the
#'     Lambert W function (Corless et al. 1996):
#'     
#'     N_eaten = N_initial - lambertW(1 / N_half * N_initial * exp(-f_max / N_half * (P * T_end - N_initial / F_max))) / (1 / N_half)
#'     
#'     We apply this function to compute the type II functional response in its
#'     Michaelis-Menten version (Real style).
#'     
#'     Required packages and their dependencies to be installed:
#'     - `emdbook` (Bolker et al. 2023)
#'     Required packages to be attached:
#'     - none
#' 
#' @references Bolker (2008) Ecological models and data in R, Princeton
#'     University Press, Princeton, New Jersey.
#'     https://math.mcmaster.ca/~bolker/emdbook/index.html
#' @references Bolker (2023) emdbook: support functions and data for "Ecological
#'     models and data". Version 1.3.13.
#'     https://CRAN.R-project.org/package=emdbook
#' @references Corless et al. (1996) On the LambertW function. Adv Comput Math
#'     5, 329-359. https://doi.org/10.1007/BF02124750
#' @references Juliano (2001) Nonlinear curve fitting: Predation and functional
#'     response curves. In: Scheiner and Gurevitch: Design and analysis of
#'     ecological experiments. Chapman and Hall, London, 2nd Edition.
#' @references Real (1977) The kinetics of functional response. Am Nat 111, 289-
#'     300. https://doi.org/10.1086/283161
#' @references Real (1979) Ecological determinants of functional response.
#'     Ecology 60, 481-485. https://doi.org/10.2307/1936067
#' @references Rogers (1972) Random search and insect population models. J Anim
#'     Ecol 41, 369-383. https://doi.org/10.2307/3474
#' @references Royama (1971) A comparative study of models for predation and
#'     parasitism. Res Popul Ecol 13, 1-91. https://doi.org/10.1007/BF02511547
#' @references Vonesh and Bolker (2005) Compensatory larval responses shift
#'     trade-offs associated with predator-induced hatching plasticity. Ecology
#'     86, 1580-1591. https://doi.org/10.1890/04-0535
#'
#' @param n_initial integer or float; a vector of initial prey densities.
#' @param p integer or float; a single value of a fixed predator density.
#'     The default value is 1.
#' @param f_max float; maximum feeding rate, a single value.
#' @param n_half float; half saturation density, a single value.
#' @param t_end integer or float; the time were the feeding ends. A single
#'     value; default = 1 (e.g. 1 day).
#' 
#' @return Returns a data frame with n_initial and n_eaten.
#' 
#' @examples
#'  
#' out <- rrpe_sim(
#'   n_initial = 1:100,
#'   p = 1,
#'   f_max = 10,
#'   n_half = 10
#' )
#' 
#' plot(out$n_initial, out$n_eaten, type = "l")
#'

rrpe_sim <- function(
    n_initial,
    p = 1,
    f_max,
    n_half,
    t_end = 1
    ){
  
  remaining_prey <- emdbook::lambertW(
    1/n_half * n_initial * exp(
      -f_max/n_half * (p * t_end - n_initial/f_max))
    ) / (1/n_half)
  
  n_eaten <- n_initial - remaining_prey
  
  data.frame(
    n_initial = n_initial,
    n_eaten = n_eaten
    )
}

 