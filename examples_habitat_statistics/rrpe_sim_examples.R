################################################################################
#   Example code to fit mod12h to data                                         #
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
# if you prefer to download a pdf, including the full statistics, follow:
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
# please also consider citing for the underlying method:
#     Bolker (2008) Ecological models and data in R, Princeton University Press,
#         Princeton, New Jersey.
#         https://math.mcmaster.ca/~bolker/emdbook/index.html
#
################################################################################
## Setup

# please install following packages, if not already done:

# install.packages("emdbook")


# empty environment:
rm(list=ls())


# import the functions from GitHub:
gh_path <- "https://raw.githubusercontent.com/b-c-r/CRITTERcode/refs/heads/main/"
f_path <- "functions_habitat_statistics/"

source(paste(gh_path, f_path, "rrpe_sim.R", sep = ""))                          # simulates a feeding functional response based on parameters and initial prey density

################################################################################
## example for Holling-style

out_h <- rrpe_sim(
  style = "Holling",
  n_initial = 1:100,                                                            # x-axis values (initial prey density)
  p = 1,                                                                        # predator density - single value
  a = .5,                                                                       # attack rate - single value
  t_h = .1                                                                      # handling time - single value
)

plot(out_h$n_initial, out_h$n_eaten, type = "l")                                # plots the output

################################################################################
## example for Real-style

out <- rrpe_sim(
  style = "Real",
  n_initial = 1:100,                                                             # x-axis values (initial prey density)
  p = 1,                                                                         # predator density - single value
  f_max = 10,                                                                    # maximum feeding rate - single value
  n_half = 20                                                                    # half saturation density - single value
)

plot(out$n_initial, out$n_eaten, type = "l")                                    # plots the output
