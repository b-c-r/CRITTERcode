################################################################################
#    Example code for gen_fr_compile                                           #
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
# install.packages("odin")

# empty environment:
rm(list=ls())

# import the functions from GitHub:
gh_path <- "https://raw.githubusercontent.com/b-c-r/CRITTERcode/refs/heads/main/"
f_path <- "functions_type_statistics/"

source(paste(gh_path, f_path, "gen_fr_compile.R", sep = ""))

# compiling the model in C using odin
gen_fr_compile()

# assign parameter values
gen_fr_model_assigned <- gen_fr_model$new(
 n_initial = 5,                                                                 # initial prey number (or density)
 f_max = 18,                                                                    # maximum feeding rate
 n_half = 3,                                                                    # half saturation density
 q = 0,                                                                         # shape parameter (q = 0 is a type II functional response)
 p = 1                                                                          # predator number (or density, fixed)
)

# defining time steps to compute, the more the better but slower
tt <- seq(
 0,                                                                             # starts at time zero                                          
 1,                                                                             # ends at time = 1, e.g. one day
 length.out = 1000                                                              # computes 1000 steps
)

# simulate the time series of decaying prey
out <- gen_fr_model_assigned$run(tt)

# plot results
plot(
 out[,1],
 out[,2],
 type = "l"
)

