################################################################################
#   Example code to fit mod05r to data                                         #
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
gh_path <- "https://raw.githubusercontent.com/b-c-r/CRITTERcode/refs/heads/main/"
f_path <- "functions_habitat_statistics/"

source(paste(gh_path, f_path, "rrpe_sim.R", sep = ""))                          # simulates a feeding functional response based on parameters and initial prey density
source(paste(gh_path, f_path, "mod05r_rrpe_nll.R", sep = ""))                   # calculates negative log likelihood
source(paste(gh_path, f_path, "mod05r_rrpe_scan.R", sep = ""))                  # calculates a set negative log likelihoods (nll) of random parameters and returns parameters from lowest nll
source(paste(gh_path, f_path, "mod05r_rrpe_fit.R", sep = ""))                   # fits functional response model to data
source(paste(gh_path, "functions_report/", "plot_mod5r.R", sep = ""))            # creates a nice plot

################################################################################
## Data

# import data from GitHub repository:
fr_data <- read.csv(
  "https://raw.githubusercontent.com/b-c-r/CRITTERdata/refs/heads/main/critter_data.csv"
)

# select a predator
fr_data_ie <- subset(fr_data, predator == "Ischnura elegans")

################################################################################
## example for: mod05r_rrpe_nll.R

fit_ie_mod05r <- bbmle::mle2(                                                   # mle2 - the maximum likelihood estimator from the bbmle package
  minuslogl = mod05r_rrpe_nll,                                                  # mle2 requires the nll function for optimization
  start = list(
    f_max_hab0_log10 = 1,                                                       # start value for maximum feeding rate (log10 scale) when habitat is absent
    f_max_hab1_log10 = 1,                                                       # start value for maximum feeding rate (log10 scale) when habitat is present
    n_half_log10 = 1                                                             # start value for half saturation density (log10 scale)
  ),
  fixed = list(
    t_end = 1,                                                                  # end time of the experiment, here 1 day
    p = 1                                                                       # number of predators in the experiment: here 1 predator per vessel
  ),
  data = list(
    n_eaten = fr_data_ie$n_eaten,                                               # data: number of prey eaten, as integer
    n_initial = fr_data_ie$n_initial,                                           # data: number of prey provided initially, as integer
    complexity = fr_data_ie$complexity_level                                    # data: complexity levels
  ),
  control = list(reltol = 1e-12),                                               # fitting tolerance: the lower the better. typically around 1e-8 (see description of the optim function)
)

bbmle::summary(fit_ie_mod05r)                                                   # shows the summary table


################################################################################
## example for: mod05r_rrpe_scan.R

mod05r_rrpe_scan(
  n_eaten = fr_data_ie$n_eaten,                                                 # data: number of prey eaten, as integer
  n_initial = fr_data_ie$n_initial,                                             # data: number of prey provided initially, as integer
  complexity  = fr_data_ie$complexity_level,                                    # data: complexity levels
  p = 1,                                                                        # number of predators in the experiment: here 1 predator per vessel
  f_max_hab0_log10_range = log10(c(1, 1.5)),                                    # two values, the range, for maximum feeding rate (log10 scale) when habitat is absent
  f_max_hab1_log10_range = c(-0.05, 0.05),                                      # two values, the range, maximum feeding rate at complexity level 4 (log10 scale)
  n_half_log10_range = log10(c(1, 1.5)),                                        # two values, the range, half saturation density (log10 scale)
  t_end = 1,                                                                    # end time of the experiment, here 1 day
  no_lhs_samples = 100                                                          # number of latin hypercube samples that should be taken (i.e. 100 random values in the range of the above assigned ranges)
)

################################################################################
## example for: mod05r_rrpe_fit.R

mod05r_fit_ie <- mod05r_rrpe_fit(
  n_eaten = fr_data_ie$n_eaten,                                                 # data: number of prey eaten, as integer
  n_initial = fr_data_ie$n_initial,                                             # data: number of prey provided initially, as integer
  complexity = fr_data_ie$complexity_level                                      # data: complexity levels
)

bbmle::summary(mod05r_fit_ie)                                                   # shows the summary table of the best fit
bbmle::AIC(mod05r_fit_ie)                                                       # shows the AIC of the best fit
BIC(mod05r_fit_ie)                                                              # shows the BIC of the best fit

################################################################################
## example for: plot_mod05r.R

plot_mod05r(
  model_fit = mod05r_fit_ie,                                                    # the mle2 fit object
  include_habitat_pics = T,                                                     # include the habitat pictograms, default = True
  pic_x1 = c( 70.0,  95.0,  70.0,  95.0),                                       # lower (left) x values for the 4 habitat pictures, the vector has values for 4 pictograms
  pic_x2 = c( 95.0, 120.0,  95.0, 120.0),                                       # upper (right) x values for the 4 habitat pictures, the vector has values for 4 pictograms
  pic_y1 = c( 22.6,  22.6,  19.4,  19.4),                                       # lower (left) y values for the 4 habitat pictures, the vector has values for 4 pictograms
  pic_y2 = c( 25.6,  25.6,  22.4,  22.4),                                       # upper (right) y values for the 4 habitat pictures, the vector has values for 4 pictograms  
  ci_reps = 100,                                                                # number of samples for the confidence interval lines, default for publication is 10000
  ci_levels = c(0.025, 0.975),                                                  # lower and upper confidence limits
  x_res = 100,                                                                  # number of x values for regression line (200+ for a smooth shape, default = 1000)
)
