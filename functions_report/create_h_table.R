################################################################################
#    create_h_table: create a nice table from from the hypotheses test results #
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
#'     if you prefer to download a pdf, including the full statistics, follow:
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
#' # see code in the statistical report (https://github.com/b-c-r/CRITTERstatistics)
#' 

create_h_table <- function(
  h_test_results,
  cut_after,
  caption_text = "Please add caption here as text string."
){
  
  model_names <- foreach::foreach(i = 1:(length(h_test_results)/2), .combine = "c") %do% {
    c(paste("Model ", i, "h", sep = ""), paste("Model ", i, "r", sep = ""))
  }
  
  aic_list <- foreach::foreach(
    i = 1:length(h_test_results),
    .combine = "rbind") %do% {
    data.frame(
      model_name_aic = model_names[i],
      df_aic = length(h_test_results[[i]]@coef),
      AIC = bbmle::AIC(h_test_results[[i]])
    )
  }
  
  aic_list <- aic_list %>%
    dplyr::arrange(AIC) %>%
    dplyr::mutate(dAIC = AIC-min(AIC)) %>%
    dplyr::mutate(
      AIC = round(AIC, 3),
      dAIC = round(dAIC, 3)
    )
  
  bic_list <- foreach::foreach(
    i = 1:length(h_test_results),
    .combine = "rbind") %do% {
      data.frame(
        model_name_bic = model_names[i],
        df_bic = length(h_test_results[[i]]@coef),
        BIC = BIC(h_test_results[[i]])
      )
    }
  
  bic_list <- bic_list %>%
    dplyr::arrange(BIC) %>%
    dplyr::mutate(dBIC = BIC-min(BIC)) %>%
    dplyr::mutate(
      BIC = round(BIC, 3),
      dBIC = round(dBIC, 3)
    )
  
  ic_list <- cbind(aic_list, bic_list)
  
  ic_list %>%
    dplyr::select(1,2,4,5,6,8) %>%
    dplyr::slice_head(n = cut_after) %>%
    kableExtra::kable(
      caption = caption_text,
      booktabs = TRUE,
      col.names = c(
        "Model (AIC)",
        "df",
        "dAIC",
        "Model (BIC)",
        "df",
        "dBIC"
      ),,
      align = "c"
    ) %>%
    kableExtra::row_spec(
      row = 0,
      bold = TRUE
    )
}
