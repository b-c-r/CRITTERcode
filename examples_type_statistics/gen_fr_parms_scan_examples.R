################################################################################
#    Example code for gen_fr_parms_scan                                        #
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
# install.packages("lhs")
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
source(paste(gh_path, f_path, "gen_fr_nll.R", sep = ""))
source(paste(gh_path, f_path, "gen_fr_parms_scan.R", sep = ""))

################################################################################
## Data

# import data from GitHub repository:
fr_data <- read.csv(
  "https://raw.githubusercontent.com/b-c-r/CRITTERdata/refs/heads/main/critter_data.csv"
)

# select a predator and complexity level
fr_data_ie <- subset(fr_data, predator == "Ischnura elegans" & complexity_level == 4)

gen_fr_compile()                                                                # compile the ODE

# scan for nlls:
gen_fr_parms_scan(
 n_eaten = fr_data_ie$n_eaten,
 n_initial = fr_data_ie$n_initial,
 p = 1,
 f_max_range_log10 = log10(c(1, max(fr_data_ie$n_eaten))),
 n_half_range_log10 = log10(c(1, max(fr_data_ie$n_initial))),
 q_range = c(0, 1),
 t_start = 0,
 t_end = 1,
 t_length = 100,
 penalty = 1000,
 q_low = 0,
 q_up = 1,
 no_lhs_samples = 1000
)
