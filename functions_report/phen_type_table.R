################################################################################
#    phen_type_table: create a nice table from phen_type_test results          #
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
#' @return Creates a nice table for a report.
#'
#' @examples
#'
#' # find an example in our report:
#' # https://github.com/b-c-r/CRITTERstatistics/blob/main/statisticsReport.Rmd
#' 

phen_type_table <- function(
    phen_test_results,
    add_caption = TRUE,
    caption_text = "Add caption text here."
    ){
  loops <- 1:length(phen_test_results)
  table <- foreach::foreach(i = loops, .combine = "rbind") %do% {

    ## get sign of linear term (linear GLM):
    sign_lin_lin <- sign(summary(phen_test_results[[i]][[2]]$modT2)[[12]][2])
    if(sign_lin_lin == 1){
      ssign_ll <- "+"
    }
    if(sign_lin_lin == 0){
      ssign_ll <- "0"
    }
    if(sign_lin_lin == -1){
      ssign_ll <- "-"
    }
    
    # get significance of linear term (linear GLM):
    sig_lin_lin <- summary(phen_test_results[[i]][[2]]$modT2)[[12]][8]
    if(sig_lin_lin >= 0.05){
      slev_ll <- "n.s."
    }
    if(sig_lin_lin < 0.05){
      slev_ll <- "*"
    }
    if(sig_lin_lin < 0.01){
      slev_ll <- "**"
    }
    if(sig_lin_lin < 0.001){
      slev_ll <- "***"
    }
    
    ## get sign of linear term (quadratic GLM):
    sign_lin <- sign(summary(phen_test_results[[i]][[2]]$modT3)[[12]][2])
    if(sign_lin == 1){
      ssign_l <- "+"
    }
    if(sign_lin == 0){
      ssign_l <- "0"
    }
    if(sign_lin == -1){
      ssign_l <- "-"
    }
    
    # get significance of linear term (quadratic GLM):
    sig_lin <- summary(phen_test_results[[i]][[2]]$modT3)[[12]][11]
    if(sig_lin >= 0.05){
      slev_l <- "n.s."
    }
    if(sig_lin < 0.05){
      slev_l <- "*"
    }
    if(sig_lin < 0.01){
      slev_l <- "**"
    }
    if(sig_lin < 0.001){
      slev_l <- "***"
    }
    
    ## get sign of quadratic term:
    sign_quad <- sign(summary(phen_test_results[[i]][[2]]$modT3)[[12]][3])
    if(sign_quad == 1){
      ssign_q <- "+"
    }
    if(sign_quad == 0){
      ssign_q <- "0"
    }
    if(sign_quad == -1){
      ssign_q <- "-"
    }
    
    # get significance of quadratic term
    sig_quad <- summary(phen_test_results[[i]][[2]]$modT3)[[12]][12]
    if(sig_quad >= 0.05){
      slev_q <- "n.s."
    }
    if(sig_quad < 0.05){
      slev_q <- "*"
    }
    if(sig_quad < 0.01){
      slev_q <- "**"
    }
    if(sig_quad < 0.001){
      slev_q <- "***"
    }
    
    # suggest type of functional response
    if(slev_q == "n.s."){
      stype <- "II"
    } else{
      if(ssign_q == "+"){
        stype <- "II"
      }
      if(ssign_q == "-"){
        stype <- "III"
      }
      if(ssign_q == "0"){
        stype <- "not clear"
      }
    }

    data.frame(
      predator = unique(phen_test_results[[i]][[3]][,2]),
      complexity = unique(phen_test_results[[i]][[3]][,3]),
      linear = paste(ssign_l, "(", slev_l, ")", sep =""),
      quadratic = paste(ssign_q, "(", slev_q, ")", sep =""),
      linear_linear = paste(ssign_ll, "(", slev_ll, ")", sep =""),
      suggested_type = stype)
  }
  
  if(add_caption){
    table %>%
      kableExtra::kable(
        booktabs = TRUE,
        caption = caption_text,
        col.names = c(
          "Predator",
          "Complexity",
          "Linear (Q)",
          "Quadratic (Q)",
          "Linear (L)",
          "Type"
        ),
        align = "c"
      ) %>%
      kableExtra::row_spec(
        row = 0,
        bold = TRUE
      ) %>%
      kableExtra::column_spec(
        column = 1,
        italic = TRUE
      ) %>%
      kableExtra::kable_styling(latex_options = "HOLD_position")
  } else{
    table %>%
      kableExtra::kable(
        booktabs = TRUE,
        col.names = c(
          "Predator",
          "Complexity",
          "Linear (Q)",
          "Quadratic (Q)",
          "Linear (L)",
          "Type"
        ),
        align = "c"
      ) %>%
      kableExtra::row_spec(
        row = 0,
        bold = TRUE
      ) %>%
      kableExtra::column_spec(
        column = 1,
        italic = TRUE
      ) %>%
      kableExtra::kable_styling(latex_options = "HOLD_position")
  }

  
}
