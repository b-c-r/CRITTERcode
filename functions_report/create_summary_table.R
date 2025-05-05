################################################################################
#    create_summary_table: creates a summary table for a pdf report            #
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
#' @description
#'     find the description including parameters here:
#'         https://github.com/b-c-r/CRITTERcode/blob/main/README.md
#'     
#'     find further details including the full statistics here:
#'         https://github.com/b-c-r/CRITTERstatistics/blob/main/statisticsReport.pdf
#'     
#'     if you are interested in the full scientific paper follow:
#'         https://doi.org/10.1101/2025.02.22.639633
#'     
#'     if you use this code, please cite:
#'         Rall et al. (2025): Habitat complexity reduces feeding strength of
#'         freshwater predators (CRITTER) - Code. Zenodo.
#'         https://doi.org/10.5281/zenodo.14894598
#'
#' @return Creates a nice table for a pdf report.
#'
#' @examples
#'
#' # find an example in our report:
#' # https://github.com/b-c-r/CRITTERstatistics/blob/main/statisticsReport.Rmd
#' 

create_summary_table <- function(
    model_fit,
    ci_reps = 10000,
    ci_levels = c(0.025, 0.975),
    dec_places = 3,
    par_names = c(
      "$T_{h}$",
      "$a$"
    ),
    unlog = c(TRUE, TRUE),
    caption_text = "Add caption text here."
){
  
  ci_samples <- MASS::mvrnorm(ci_reps, mu = model_fit@coef, Sigma = bbmle::vcov(model_fit))

  out_table <- foreach::foreach(i = 1: length(par_names), .combine = "rbind") %do% {
    if(unlog[i]){
      data.frame(
        par_names = par_names[i],
        par_value = round(10^bbmle::summary(model_fit)@coef[i], dec_places),
        par_value_low = round(10^as.numeric(quantile(ci_samples[,i], ci_levels[1])), dec_places),
        par_value_up = round(10^as.numeric(quantile(ci_samples[,i], ci_levels[2])), dec_places)
      )
    } else {
      data.frame(
        par_names = par_names[i],
        par_value = round(bbmle::summary(model_fit)@coef[i], dec_places),
        par_value_low = round(as.numeric(quantile(ci_samples[,i], ci_levels[1])), dec_places),
        par_value_up = round(as.numeric(quantile(ci_samples[,i], ci_levels[2])), dec_places)
      )
    }
  }
  
  out_table %>%
    knitr::kable(
      escape = FALSE,
      format = "latex",
      caption = caption_text,
      booktabs = TRUE,
      col.names = c(
        "Parameter Name",
        "Point Estimate",
        "Lower CI",
        "Upper CI"
      ),
      align = "c"
    ) %>%
    kableExtra::row_spec(
      row = 0,
      bold = TRUE
    ) %>%
    kableExtra::kable_styling(latex_options = "HOLD_position")
}



