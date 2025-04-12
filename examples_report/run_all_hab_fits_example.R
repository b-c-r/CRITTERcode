################################################################################
#   Example code for run_all_hab_fits                                          #
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
################################################################################
## Setup

# please install following packages, if not already done:

# install.packages("bbmle")
# install.packages("dplyr")
# install.packages("emdbook")
# install.packages("foreach")
# install.packages("lhs")


# empty environment:
rm(list=ls())

# attach libraries
library("dplyr")
library("foreach")

# import the functions from GitHub:
hab_stats_files <- supportR::github_ls(
  "https://github.com/b-c-r/CRITTERcode",
  "functions_habitat_statistics"
)$name

gh_path <- "https://raw.githubusercontent.com/b-c-r/CRITTERcode/refs/heads/main/"
f_path <- "functions_habitat_statistics/"
for(i in 1:length(hab_stats_files)){
  source(paste(gh_path, f_path, hab_stats_files[i], sep = ""))
}

report_files <- supportR::github_ls(
  "https://github.com/b-c-r/CRITTERcode",
  "functions_report"
)$name

f_path <- "functions_report/"
for(i in 1:length(report_files)){
  source(paste(gh_path, f_path, report_files[i], sep = ""))
}

################################################################################
## Data

# import data from GitHub repository:
fr_data <- read.csv(
  "https://raw.githubusercontent.com/b-c-r/CRITTERdata/refs/heads/main/critter_data.csv"
)

# select a predator
fr_data_ie <- subset(fr_data, predator == "Ischnura elegans")

################################################################################
## example for: mod15h_rrpe_nll.R

hab_results <- run_all_hab_fits(
  x = fr_data,
  model_numbers = 15,
  style_letters = "h",
  predator_spec = "Ischnura elegans",
  no_threads = 1
)

################################################################################
## example for: plot_mod15h.R

plot_mod15h(
  model_fit = hab_results[[1]],                                                 # the mle2 fit object
  include_habitat_pics = T,                                                     # include the habitat pictograms, default = True
  pic_x1 = rep(100.0, 4),                                                       # lower (left) x values for the 4 habitat pictures, the vector has values for 4 pictograms
  pic_x2 = rep(120.0, 4),                                                       # upper (right) x values for the 4 habitat pictures, the vector has values for 4 pictograms
  pic_y1 = rep( 22.0, 4),                                                       # lower (left) y values for the 4 habitat pictures, the vector has values for 4 pictograms
  pic_y2 = rep( 25.0, 4),                                                       # upper (right) y values for the 4 habitat pictures, the vector has values for 4 pictograms  
  ylim = c(0, 25),
  ci_reps = 100,                                                                # number of samples for the confidence interval lines, default for publication is 10000
  ci_levels = c(0.025, 0.975),                                                  # lower and upper confidence limits
  x_res = 100,                                                                  # number of x values for regression line (200+ for a smooth shape, default = 1000)
)
