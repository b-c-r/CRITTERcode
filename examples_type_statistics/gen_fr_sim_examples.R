################################################################################
#    Example code for gen_fr_sim                                               #
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
#
# find a description including here:
#     https://github.com/b-c-r/CRITTERcode/blob/main/README.md
#     
# find further details including the full statistics here:
#     https://github.com/b-c-r/CRITTERstatistics/blob/main/statisticsReport.pdf
#     
# if you are interested in the full scientific paper follow:
#     https://doi.org/10.1101/2025.02.22.639633
#     
# if you use this code, please cite:
#     Rall et al. (2025): Habitat complexity reduces feeding strength of
#         freshwater predators (CRITTER) - Code. Zenodo.
#         https://doi.org/10.5281/zenodo.14894598
#
################################################################################
## Setup

# please install following packages, if not already done:
# install.packages("foreach")
# install.packages("odin")

# empty environment:
rm(list=ls())

# attach packages
library("foreach")

# import the functions from GitHub:
gh_path <- "https://raw.githubusercontent.com/b-c-r/CRITTERcode/refs/heads/main/"
f_path <- "functions_type_statistics/"

source(paste(gh_path, f_path, "gen_fr_compile.R", sep = ""))
source(paste(gh_path, f_path, "gen_fr_sim.R", sep = ""))

################################################################################

gen_fr_compile()                                                                # compile the ODE

# simulate the functional response:
out <- gen_fr_sim(
  n_initial = 1:100,                                                            # vector of initial prey densities
  p = 1,                                                                        # fixed predator density
  f_max = 10,                                                                   # maximum feeding rate
  n_half = 25,                                                                  # half saturation density
  q = 1                                                                         # shape parameter (1 = s-shaped)
)
 
# plot the functional response
plot(
  out$n_initial,
  out$n_eaten,
  type = "l"
)

