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
#' @description `gen_fr_sim` simulates time series across different initial prey
#'     densities. It returns only the number of initial prey items and the
#'     number of prey eaten at the end of the time series, mimicking common
#'     functional response laboratory experiments. The underlying functional
#'     response model is the generalized functional response model (Real 1977,
#'     1979). Note that not integers, but floating numbers are returned by
#'     `gen_fr_sim`. `gen_fr_sim` depends on `gen_fr_compile`. Please find more
#'     details there.
#'     
#'     Required packages and their dependencies to be installed:
#'       - `odin` (FitzJohn and Jombart 2024)
#'       - `foreach` (Microsoft and Weston 2022)
#'     Required packages to be attached:
#'       - `foreach` (Microsoft and Weston 2022)
#'       
#' 
#' @references FitzJohn and Jombart (2024) odin: ODE generation and integration.
#'     Ver. 1.2.6. https://doi.org/10.32614/CRAN.package.odin
#'     see also: https://github.com/mrc-ide/odin
#' @references Microsoft and Weston (2022) foreach: provides foreach looping
#'     construct. Version 1.5.2. https://doi.org/10.32614/CRAN.package.foreach
#' @references Real (1977) The kinetics of functional response. Am Nat 111, 289-
#'     300. https://doi.org/10.1086/283161
#' @references Real (1979) Ecological determinants of functional response.
#'     Ecology 60, 481-485. https://doi.org/10.2307/1936067     
#' 
#' @include gen_fr_compile.R
#'
#' @param n_initial integer or float; a vector of initial prey densities.
#' @param p integer or float; a single value of a fixed predator density. The
#'     default value is 1.
#' @param f_max float; maximum feeding rate, a single value.
#' @param n_half float; half saturation density, a single value.
#' @param q float; shape parameter, a single value. A strict type II functional
#'     has q = 0, a strict type III functional response has q = 1.
#' @param t_start integer or float; the time were the feeding starts. A single
#'     value; default = 0.
#' @param t_end integer or float; the time were the feeding ends. A single
#'     value; default = 1 (e.g. 1 day).
#' @param t_length integer or float; the number of time steps that should be
#'     generated. The more time steps, the more precise the simulation. A single
#'     value; default = 1000.
#' 
#' @return Returns a data frame with n_initial and the according n_eaten.
#' 
#' @examples
#' 
#' library("foreach")
#' source(here::here("functions_gen_fr", "gen_fr_compile.R"))
#' source(here::here("functions_gen_fr", "gen_fr_sim.R"))
#' 
#' # compile the functional response ODE in C
#' gen_fr_compile()
#' 
#' # simulate the functional response:
#' out <- gen_fr_sim(
#'  n_initial = 1:100,                                                          # vector of initial prey densities
#'  p = 1,                                                                      # fixed predator density
#'  f_max = 10,                                                                 # maximum feeding rate
#'  n_half = 25,                                                                # half saturation density
#'  q = 1                                                                       # shape parameter (1 = s-shaped)
#' )
#' 
#' # plot the functional response
#' plot(
#'   out$n_initial,
#'   out$n_eaten,
#'   type = "l"
#' )
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
