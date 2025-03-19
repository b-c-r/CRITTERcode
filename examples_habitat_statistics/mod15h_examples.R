################################################################################
#   Example code to fit mod15h to data                                         #
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
source(paste(gh_path, f_path, "mod15h_rrpe_nll.R", sep = ""))                   # calculates negative log likelihood
source(paste(gh_path, f_path, "mod15h_rrpe_scan.R", sep = ""))                  # calculates a set negative log likelihoods (nll) of random parameters and returns parameters from lowest nll
source(paste(gh_path, f_path, "mod15h_rrpe_fit.R", sep = ""))                   # fits functional response model to data
source(paste(gh_path, "functions_report/", "plot_mod15h.R", sep = ""))            # creates a nice plot

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

fit_ie_mod15h <- bbmle::mle2(                                                   # mle2 - the maximum likelihood estimator from the bbmle package
  minuslogl = mod15h_rrpe_nll,                                                  # mle2 requires the nll function for optimization
  start = list(
    t_h_0_log10 = -1,                                                           # start value for handling time at complexity level 0 (log10 scale)
    t_h_1_log10 = -1,                                                           # start value for handling time at complexity level 1 (log10 scale)
    t_h_2_log10 = -1,                                                           # start value for handling time at complexity level 2 (log10 scale)
    t_h_3_log10 = -1,                                                           # start value for handling time at complexity level 3 (log10 scale)
    t_h_4_log10 =  -1,                                                          # start value for handling time at complexity level 4 (log10 scale)
    a_intercept_log10 = -1,                                                     # start value for attack rate intercept (log10 scale, when habitat is absent)
    a_slope = 0                                                                 # start value for attack rate slope
  ),
  fixed = list(
    t_end = 1,                                                                  # end time of the experiment, here 1 day
    p = 1                                                                       # number of predators in the experiment: here 1 predator per vessel
  ),
  data = list(
    n_eaten = fr_data_ie$n_eaten,                                               # data: number of prey eaten, as integer
    n_initial = fr_data_ie$n_initial,                                           # data: number of prey provided initially, as integer
    n_rings = fr_data_ie$ring_count,                                            # data: number of habitat structural elements (rings)
    complexity = fr_data_ie$complexity_level                                    # data: complexity levels
  ),
  control = list(reltol = 1e-12),                                               # fitting tolerance: the lower the better. typically around 1e-8 (see description of the optim function)
)

bbmle::summary(fit_ie_mod15h)                                                   # shows the summary table


################################################################################
## example for: mod15h_rrpe_scan.R

mod15h_rrpe_scan(
  n_eaten = fr_data_ie$n_eaten,                                                 # data: number of prey eaten, as integer
  n_initial = fr_data_ie$n_initial,                                             # data: number of prey provided initially, as integer
  n_rings = fr_data_ie$ring_count,                                              # data: number of habitat structural elements (rings)
  complexity  = fr_data_ie$complexity_level,                                    # data: complexity levels
  p = 1,                                                                        # number of predators in the experiment: here 1 predator per vessel
  t_h_0_log10_range = log10(c(.1, 0.5)),                                        # two values, the range, handling time at complexity level 0 (log10 scale)
  t_h_1_log10_range = log10(c(.1,0.5)),                                         # two values, the range, handling time at complexity level 1 (log10 scale)
  t_h_2_log10_range = log10(c(.1,0.5)),                                         # two values, the range, handling time at complexity level 2 (log10 scale)
  t_h_3_log10_range  = log10(c(.1,0.5)),                                        # two values, the range, handling time at complexity level 3 (log10 scale)
  t_h_4_log10_range = log10(c(.1,0.5)),                                         # two values, the range, handling time at complexity level 4 (log10 scale)
  a_intercept_log10_range = log10(c(.01,0.5)),                                  # two values, the range, attack rate intercept (log10 scale, when habitat is absent)
  a_slope_range = c(-0.05,0.05),                                                # two values, the range, attack rate slope
  t_end = 1,                                                                    # end time of the experiment, here 1 day
  no_lhs_samples = 100                                                          # number of latin hypercube samples that should be taken (i.e. 100 random values in the range of the above assigned ranges)
)

################################################################################
## example for: mod15h_rrpe_fit.R

mod15h_fit_ie <- mod15h_rrpe_fit(
  n_eaten = fr_data_ie$n_eaten,                                                 # data: number of prey eaten, as integer
  n_initial = fr_data_ie$n_initial,                                             # data: number of prey provided initially, as integer
  n_rings = fr_data_ie$ring_count,                                              # data: number of habitat structural elements (rings)
  complexity = fr_data_ie$complexity_level                                      # data: complexity levels
)

bbmle::summary(mod15h_fit_ie)                                                   # shows the summary table of the best fit
bbmle::AIC(mod15h_fit_ie)                                                       # shows the AIC of the best fit
BIC(mod15h_fit_ie)                                                              # shows the BIC of the best fit

################################################################################
## example for: plot_mod15h.R

plot_mod15h(
  model_fit = mod15h_fit_ie,                                                    # the mle2 fit object
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
