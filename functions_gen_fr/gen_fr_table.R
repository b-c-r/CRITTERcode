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
#' 
#' @description “gen_fr_table” creates a nice-looking table from the output
#'     created by “gen_fr_fit_all”. Hard-coded and only useful in this project.
#'     
#'     Required packages and their dependencies to be installed:
#'       - `bbmle` (Bolker et al 2023)
#'       - `dplyr` (Wickham et al. 2023)
#'       - `foreach` (Microsoft and Weston 2022)
#'       - `kableExtra` (Zhu et al. 2024)
#'       - `knitr` (Xie 2024)
#'     Required packages to be attached:
#'       - `dplyr` (Wickham et al. 2023)
#'       - `foreach` (Microsoft and Weston 2022)
#' 
#' @references Bolker (2023) bbmle: tools for general Maximum Likelihood
#'     Estimation. Version 1.0.25.1. https://doi.org/10.32614/CRAN.package.bbmle
#' @references Microsoft and Weston (2022) foreach: provides foreach looping
#'     construct. Version 1.5.2. https://CRAN.R-project.org/package=foreach
#' @references Wickham et al. (2023) dplyr: a grammar of data manipulation.
#'     1.1.4. https://CRAN.R-project.org/package=dplyr
#' @references Xie (2024) knitr: a general-purpose package for dynamic report
#'     generation in R. Version 1.49. https://yihui.org/knitr/
#' @references Zhu et al. (2024) kableExtra: construct complex table with
#'     'kable' and pipe syntax. Version 1.4.0.
#'     https://doi.org/10.32614/CRAN.package.kableExtra
#'
#' 
#' @param gen_fr_results an object created by “gen_fr_fit_all”.
#' @param caption_text your caption test as a string.
#' @param output_style decide on the output format of the Rmd file. Can be
#'     "github_document" (default) or "pdf".
#' 

gen_fr_table <- function(
    gen_fr_results,
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
  } else{
    if(output_style == "pdf"){
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
      stop("No valid output style defined. Please select eihter github_document or pdf.")
    }
  }
  

}
