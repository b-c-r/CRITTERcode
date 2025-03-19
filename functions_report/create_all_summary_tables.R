################################################################################
#    create_all_summary_tables: creates a table including all habitat results  #
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

create_all_summary_tables <- function(
    fit_results,
    ci_samples = 10000,
    caption_text = "All 32 summary tables from \\textit{Ischnura elegans} model fits."
){
  
  out <- foreach::foreach(i = 1:32, .combine = "rbind") %do% {
    
    ci_samples <- MASS::mvrnorm(ci_samples, mu = fit_results[[i]]@coef, Sigma = bbmle::vcov(fit_results[[i]]))
    
    par_exp = 10^bbmle::summary(fit_results[[i]])@coef[,1]
    par_value_low <- foreach::foreach(j = 1:ncol(ci_samples), .combine = "c") %do% {10^as.numeric(quantile(ci_samples[,j], 0.025))}
    par_value_up  <- foreach::foreach(j = 1:ncol(ci_samples), .combine = "c") %do% {10^as.numeric(quantile(ci_samples[,j], 0.975))}
    
    data.frame(
      Parameter = rownames(bbmle::summary(fit_results[[i]])@coef),
      PointEstimate = round(bbmle::summary(fit_results[[i]])@coef[,1],3),
      SE = round(bbmle::summary(fit_results[[i]])@coef[,2],3),
      z = round(bbmle::summary(fit_results[[i]])@coef[,3],3),
      p = round(bbmle::summary(fit_results[[i]])@coef[,4],4),
      ParameterExp = round(par_exp,3),
      par_value_low = round(par_value_low,3),
      par_value_up = round(par_value_up,3),
      row.names = NULL
    )
  }
  
  non_exp <- c(13,16,46,50,66,69,72,76,80,82 ,84 , 86, 88, 95,133,140)
  
  for(i in 1:nrow(out)){
    if(any(non_exp==i)){
      out[i,6] <- round(log10(out[i,6]),3)
      out[i,7] <- round(log10(out[i,7]),3)
      out[i,8] <- round(log10(out[i,8]),3)
    }
  }
  
  # https://stackoverflow.com/questions/46251023/kableextra-continued-on-next-page-for-longtable
  out %>%
    knitr::kable(
      escape = TRUE,
      format = "latex",
      caption = caption_text,
      booktabs = TRUE,
      longtable = TRUE,
      col.names = c(
        "Name",
        "Orig. Est.",
        "SE",
        "z",
        "p",
        "Estimate",
        "CI low",
        "CI up"
      ),
      align = "c"
    ) %>%
    kableExtra::kable_styling() %>%
    kableExtra::pack_rows( "Model 1h",   1,   2) %>%
    kableExtra::pack_rows( "Model 1r",   3,   4) %>%
    kableExtra::pack_rows( "Model 2h",   5,   7) %>%
    kableExtra::pack_rows( "Model 2r",   8,  10) %>%
    kableExtra::pack_rows( "Model 3h",  11,  13) %>%
    kableExtra::pack_rows( "Model 3r",  14,  16) %>%
    kableExtra::pack_rows( "Model 4h",  17,  22) %>%
    kableExtra::pack_rows( "Model 4r",  23,  28) %>%
    kableExtra::pack_rows( "Model 5h",  29,  31) %>%
    kableExtra::pack_rows( "Model 5r",  32,  34) %>%
    kableExtra::pack_rows( "Model 6h",  35,  38) %>%
    kableExtra::pack_rows( "Model 6r",  39,  42) %>%
    kableExtra::pack_rows( "Model 7h",  43,  46) %>%
    kableExtra::pack_rows( "Model 7r",  47,  50) %>%
    kableExtra::pack_rows( "Model 8h",  51,  57) %>%
    kableExtra::pack_rows( "Model 8r",  58,  64) %>%
    kableExtra::pack_rows( "Model 9h",  65,  67) %>%
    kableExtra::pack_rows( "Model 9r",  68,  70) %>%
    kableExtra::pack_rows("Model 10h",  71,  74) %>%
    kableExtra::pack_rows("Model 10r",  75,  78) %>%
    kableExtra::pack_rows("Model 11h",  79,  82) %>%
    kableExtra::pack_rows("Model 11r",  83,  86) %>%
    kableExtra::pack_rows("Model 12h",  87,  93) %>%
    kableExtra::pack_rows("Model 12r",  94, 100) %>%
    kableExtra::pack_rows("Model 13h", 101, 106) %>%
    kableExtra::pack_rows("Model 13r", 107, 112) %>%
    kableExtra::pack_rows("Model 14h", 113, 119) %>%
    kableExtra::pack_rows("Model 14r", 120, 126) %>%
    kableExtra::pack_rows("Model 15h", 127, 133) %>%
    kableExtra::pack_rows("Model 15r", 134, 140) %>%
    kableExtra::pack_rows("Model 16h", 141, 150) %>%
    kableExtra::pack_rows("Model 16r", 151, 160) %>%
    kableExtra::kable_styling(
      latex_options = c("HOLD_position", "repeat_header"),
      repeat_header_continued = "\\textit{(continued on next page...)}"
    )
}
