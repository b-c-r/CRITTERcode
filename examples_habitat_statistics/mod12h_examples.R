################################################################################
#   Exmaple code to fit mod12h to data                                         #
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

################################################################################
## Setup

rm(list=ls())

library("foreach")
library("dplyr")

# import the functions from GitHub:
gh_path <- "https://raw.githubusercontent.com/b-c-r/CRITTERcode/refs/heads/main/"
f_path <- "functions_habitat_statistics"

source(paste(gh_path, f_path, "rrpe_sim.R"))                                    # simulates a feeding functional response based on parameters and initial prey density
source(paste(gh_path, f_path, "rrpe_nll_mod12h.R"))                             # calculates negative log likelihood
source(paste(gh_path, f_path, "rrpe_scan_mod12h.R"))                            # calculates a set negative log likelihoods (nll) of random parameters and returns parameters from lowest nll
source(paste(gh_path, f_path, "rrpe_fit_mod12h.R"))                             # fits functional response model to data

################################################################################
## Data

# import data from GitHub repository:
fr_data <- read.csv(
  "https://raw.githubusercontent.com/b-c-r/CRITTERdata/refs/heads/main/critter_data.csv"
)

# select a predator
fr_data_ie <- subset(fr_data, predator == "Ischnura elegans")

################################################################################
## example for: rrpe_nll_mod12h.R

fit_ie_mod12h <- bbmle::mle2(                                                   # mle2 - the maximum likelihood estimator from the bbmle package
 minuslogl = rrpe_nll_mod12h,                                                   # mle2 requires the nll function for optimization
 start = list(
   t_h_intercept_log10  = -1.5,                                                 # start value for handling time intercept (log10 scale)
   t_h_slope_log10      =    0,                                                 # start value for handling time (linear after log10 transformation)
   a_0_log10            =   -1,                                                 # start value for attack rate at complexity level 0 (log10 scale)
   a_1_log10            =   -1,                                                 # start value for attack rate at complexity level 1 (log10 scale)
   a_2_log10            =   -1,                                                 # start value for attack rate at complexity level 2 (log10 scale)
   a_3_log10            =   -1,                                                 # start value for attack rate at complexity level 3 (log10 scale)
   a_4_log10            =   -1                                                  # start value for attack rate at complexity level 4 (log10 scale)
 ),
 fixed = list(
   t_end = 1,                                                                   # end time of the experiment, here 1 day
   p = 1                                                                        # number of predators in the experiment: here 1 predator per vessel
 ),
 data = list(
   n_eaten = fr_data_ie$n_eaten,                                                # data: number of prey eaten, as integer
   n_initial = fr_data_ie$n_initial,                                            # data: number of prey provided initially, as integer
   n_rings = fr_data_ie$ring_count,                                             # data: number of habitat rings provided as structure
   complexity = fr_data_ie$complexity_level                                     # data: complexity levels
 ),
 control = list(reltol = 1e-12),                                                # fitting tolerance: the lower the better. typically around 1e-8 (see description of the optim function)
)

bbmle::summary(fit_ie_mod12h)                                                   # shows the summary table


################################################################################
## example for: rrpe_scan_mod12h.R

rrpe_scan_mod12h(
 n_eaten = fr_data_ie$n_eaten,                                                  # data: number of prey eaten, as integer
 n_initial = fr_data_ie$n_initial,                                              # data: number of prey provided initially, as integer
 n_rings = fr_data_ie$ring_count,                                               # data: number of habitat rings provided as structure
 complexity  = fr_data_ie$complexity_level,                                     # data: complexity levels
 p = 1,                                                                         # number of predators in the experiment: here 1 predator per vessel
 t_h_range_intercept_log10 = c(.05, .05),                                       # two values, the range, for handling time intercept (log10 scale)
 t_h_range_slope_log10 = c(-0.05, 0.05),                                        # two values, the range, handling time (linear after log10 transformation)
 a_range_0_log10 = c(.01,0.5),                                                  # two values, the range, attack rate at complexity level 0 (log10 scale)
 a_range_1_log10 = c(.01,0.5),                                                  # two values, the range, attack rate at complexity level 1 (log10 scale)
 a_range_2_log10 = c(.01,0.5),                                                  # two values, the range, attack rate at complexity level 2 (log10 scale)
 a_range_3_log10 = c(.01,0.5),                                                  # two values, the range, attack rate at complexity level 3 (log10 scale)
 a_range_4_log10 = c(.01,0.5),                                                  # two values, the range, attack rate at complexity level 4 (log10 scale)
 t_end = 1,                                                                     # end time of the experiment, here 1 day
 no_lhs_samples = 100                                                           # number of latin hypercube samples that should be taken (i.e. 100 random values in the range of the above assigned ranges)
)

################################################################################
## example for: rrpe_fit_mod12h.R

mod12h_fit_ie <- rrpe_fit_mod12h(
 n_eaten = fr_data_ie$n_eaten,                                                  # data: number of prey eaten, as integer
 n_initial = fr_data_ie$n_initial,                                              # data: number of prey provided initially, as integer
 n_rings = fr_data_ie$ring_count,                                               # data: number of habitat rings provided as structure
 complexity = fr_data_ie$complexity_level                                       # data: complexity levels
)

bbmle::summary(mod12h_fit_ie)                                                   # shows the summary table of the best fit
bbmle::AIC(mod12h_fit_ie)                                                       # shows the AIC of the best fit
BIC(mod12h_fit_ie)                                                              # shows the BIC of the best fit
