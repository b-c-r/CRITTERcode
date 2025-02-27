################################################################################
#   rrpe_nll_mod13h: estimates the negative log-likelihood of mod13h           #
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
#'
#' @description `rrpe_nll_mod13h` calculates the negative log likelihood of
#'     using experimental functional response data (eaten resource items as a
#'     function of resource density) using the Michaelis-Menten Type
#'     II functional response (Real 1977) from the `rrpe_sim` function, but
#'     transform the values to the Holling type II (Holling 1959): T_h and a.
#'     We calculated the likelihood by assuming a binomial distribution of the
#'     depended data, i.e., whether a resource item can be eaten or not at the
#'     end of the experimental trial. See (Bolker 2008), chapter 8 for details.
#'     Moreover, the function requires the model parameters on log-scale, as
#'     this transformation (1) accelerates the fitting procedure and (2)
#'     prevents biologically irrelevant negative estimations that would crash
#'     the fitting algorithm.
#'     
#'     Required packages and their dependencies to be installed:
#'         - `emdbook` (Bolker 2023)
#'     
#'     Required packages to be attached:
#'         - none
#'     
#' @references Bolker (2008) Ecological models and data in R, Princeton
#'     University Press, Princeton, New Jersey.
#'     https://math.mcmaster.ca/~bolker/emdbook/index.html
#' @references Bolker (2023) emdbook: support functions and data for "Ecological
#'     models and data". Version 1.3.13.
#'     https://CRAN.R-project.org/package=emdbook
#' @references Holling (1959) Some characteristics of simple types of predation
#'     and parasitism. Can Entomol 91, 385-398.
#'     https://doi.org/10.4039/Ent91385-7
#' @references Real (1977) The kinetics of functional response. Am Nat 111, 289-
#'     300. https://doi.org/10.1086/283161
#'     
#' @include rrpe_sim.R
#'
#' @param n_eaten integer (or float); the prey items that were eaten throughout the experimental trial. A vector.
#' @param n_initial integer or float; a vector of initial prey densities.
#' @param complexity level of complexity (0-4), a single integer value.
#' @param p integer or float; a single value of a fixed predator density. The default value is 1.
#' @param t_h_0_log10 t_h for the respective complexity level, single value.
#' @param t_h_1_log10 t_h for the respective complexity level, single value.
#' @param t_h_2_log10 t_h for the respective complexity level, single value.
#' @param t_h_3_log10 t_h for the respective complexity level, single value.
#' @param t_h_4_log10 t_h for the respective complexity level, single value.
#' @param a_log10 a, single value.
#' @param t_end integer or float; the time were the feeding ends. A single value; default = 1 (e.g. 1 day).
#'
#'     
#' @return Returns a single negative log-likelihood value.
#' 
#' @examples
#' 
#' rm(list=ls())
#' 
#' library("foreach")
#' 
#' gh_path <- "https://raw.githubusercontent.com/b-c-r/CRITTERcode/refs/heads/main/functions_habitat_statistics/"
#' 
#' source(paste(gh_path, "rrpe_sim.R", sep = ""))
#' source(paste(gh_path, "rrpe_nll_mod13h.R", sep = ""))
#' 
#' fr_data <- read.csv("https://raw.githubusercontent.com/b-c-r/CRITTERdata/refs/heads/main/critter_data.csv")
#' fr_data_ie <- subset(fr_data, predator == "Ischnura elegans")
#' 
#' fit_ie_mod13h <- bbmle::mle2(
#'   minuslogl = rrpe_nll_mod13h,
#'   start = list(
#'     t_h_0_log10  =    log10(max(fr_data_ie$n_eaten[fr_data_ie$complexity_level == 0])),
#'     t_h_1_log10  =    log10(max(fr_data_ie$n_eaten[fr_data_ie$complexity_level == 1])),
#'     t_h_2_log10  =    log10(max(fr_data_ie$n_eaten[fr_data_ie$complexity_level == 2])),
#'     t_h_3_log10  =    log10(max(fr_data_ie$n_eaten[fr_data_ie$complexity_level == 3])),
#'     t_h_4_log10  =    log10(max(fr_data_ie$n_eaten[fr_data_ie$complexity_level == 4])),
#'     a_log10 = log10(mean(fr_data_ie$n_initial[fr_data_ie$complexity_level == 0]))
#'   ),
#'   fixed = list(
#'     t_end = 1,
#'     p = 1
#'   ),
#'   data = list(
#'     n_eaten = fr_data_ie$n_eaten,
#'     n_initial = fr_data_ie$n_initial,
#'     complexity = fr_data_ie$complexity_level
#'   ),
#'   control = list(reltol = 1e-12),
#'   
#' )
#' 
#' bbmle::summary(fit_ie_mod13h)
#' 
#' #############################################################################
#' 
#' fr_data_ng <- subset(fr_data, predator == "Notonecta glauca")
#' 
#' fit_ng_mod13h <- bbmle::mle2(
#'   minuslogl = rrpe_nll_mod13h,
#'   start = list(
#'     t_h_0_log10  =    log10(max(fr_data_ng$n_eaten[fr_data_ng$complexity_level == 0])),
#'     t_h_1_log10  =    log10(max(fr_data_ng$n_eaten[fr_data_ng$complexity_level == 1])),
#'     t_h_2_log10  =    log10(max(fr_data_ng$n_eaten[fr_data_ng$complexity_level == 2])),
#'     t_h_3_log10  =    log10(max(fr_data_ng$n_eaten[fr_data_ng$complexity_level == 3])),
#'     t_h_4_log10  =    log10(max(fr_data_ng$n_eaten[fr_data_ng$complexity_level == 4])),
#'     a_log10 = log10(mean(fr_data_ng$n_initial[fr_data_ng$complexity_level == 0]))
#'   ),
#'   fixed = list(
#'     t_end = 1,
#'     p = 1
#'   ),
#'   data = list(
#'     n_eaten = fr_data_ng$n_eaten,
#'     n_initial = fr_data_ng$n_initial,
#'     complexity = fr_data_ng$complexity_level
#'   ),
#'   control = list(reltol = 1e-12),
#'   
#' )
#' 
#' bbmle::summary(fit_ng_mod13h)
#' 

rrpe_nll_mod13h <- function(
    n_eaten,
    n_initial,
    complexity,
    p = 1,
    t_h_0_log10,
    t_h_1_log10,
    t_h_2_log10,
    t_h_3_log10,
    t_h_4_log10,
    a_log10,
    t_end = 1
){
  
  eaten_simulated <- foreach::foreach(
    i = 1:length(n_initial),
    .combine = "rbind"
  ) %do% {
    
    if(complexity[i] == 0){
      t_h_log10 <- t_h_0_log10
    }
    if(complexity[i] == 1){
      t_h_log10 <- t_h_1_log10
    }
    if(complexity[i] == 2){
      t_h_log10 <- t_h_2_log10
    }
    if(complexity[i] == 3){
      t_h_log10 <- t_h_3_log10
    }
    if(complexity[i] == 4){
      t_h_log10 <- t_h_4_log10
    }
    
    rrpe_sim(
      n_initial = n_initial[i],
      p = p,
      f_max = 1 / 10^t_h_log10,
      n_half = 1 / (10^a_log10 * 10^t_h_log10),
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
  
  return(nll)
}
