################################################################################
#   gen_fr_table: creates a good looking table for results from gen_fr_fit_all #
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
' 
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
#'
#' @return Creates a nice table for a report.
#'
#' @examples
#'
#' # find an example in our report:
#' # https://github.com/b-c-r/CRITTERstatistics/blob/main/statisticsReport.Rmd
#' 

gen_fr_table <- function(
    gen_fr_results,
    add_caption = TRUE,
    caption_text = "Add caption text here.",
    output_style = "github_document"
    ){
  loops <- 1:length(gen_fr_results)
  table <- foreach::foreach(i = loops, .combine = "rbind") %do% {
    
    # get significance level:
    sig_q <- bbmle::summary(gen_fr_results[[i]][[2]])@coef[[12]]
    if(sig_q  >= 0.05){
      slev_q <- "n.s."
    }
    if(sig_q  < 0.05){
      slev_q <- "*"
    }
    if(sig_q < 0.01){
      slev_q <- "**"
    }
    if(sig_q < 0.001){
      slev_q <- "***"
    }
    
    # suggest type of functional response
    if(slev_q == "n.s."){
      stype <- "II"
    } else{
      if(bbmle::summary(gen_fr_results[[i]][[2]])@coef[[3]] < 0){
        stype <- "II"
      } else {
        stype <- "III"
      }
    }
    
    if(output_style == "github_document"){
      data.frame(
        predator = paste("*", unique(gen_fr_results[[i]][[3]][,2]), "*", sep=""),
        complexity = unique(gen_fr_results[[i]][[3]][,3]),
        q = round(bbmle::summary(gen_fr_results[[i]][[2]])@coef[[3]],3),
        p_value = slev_q,
        suggested_type = stype)
    } else{
      if(output_style == "pdf"){
        data.frame(
          predator = unique(gen_fr_results[[i]][[3]][,2]),
          complexity = unique(gen_fr_results[[i]][[3]][,3]),
          q = round(bbmle::summary(gen_fr_results[[i]][[2]])@coef[[3]],3),
          p_value = slev_q,
          suggested_type = stype)
      } else{
        stop("No valid output style defined. Please select eihter github_document or pdf.")
      }
    }
  }

  if(output_style == "github_document"){
    
    if(add_caption){
      
      table %>%
        knitr::kable(
          caption = caption_text,
          col.names = c(
            "Predator Name",
            "Complexity",
            "*q*",
            "Significance",
            "Type-o-Response"
          )
        )
    } else {
      
      table %>%
        knitr::kable(
          col.names = c(
            "Predator Name",
            "Complexity",
            "*q*",
            "Significance",
            "Type-o-Response"
          )
        )
    }
    
  } else{
    if(output_style == "pdf"){
      
      if(add_caption){
        table %>%
          knitr::kable(
            escape = F,
            format = "latex",
            caption = caption_text,
            booktabs = TRUE,
            col.names = c(
              "Predator Name",
              "Complexity",
              "$q$",
              "Significance",
              "Type-o-Response"
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
          knitr::kable(
            escape = F,
            format = "latex",
            booktabs = TRUE,
            col.names = c(
              "Predator Name",
              "Complexity",
              "$q$",
              "Significance",
              "Type-o-Response"
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
      
    } else{
      stop("No valid output style defined. Please select eihter github_document or pdf.")
    }
  }
  

}
