################################################################################
#    create_package_info: creates a table with packages used                   #
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

create_package_info <- function(
    include_base = TRUE,
    dependencies = TRUE,
    caption_text = "All loaded R packages that we used in this report, including base packages and dependencies."
){
  dplyr::tibble(
    pname = sessioninfo::package_info(pkgs = "loaded", include_base = include_base, dependencies = dependencies)[,1],
    pver = sessioninfo::package_info(pkgs = "loaded", include_base = include_base, dependencies = dependencies)[,3],
    pattached = sessioninfo::package_info(pkgs = "loaded", include_base = include_base, dependencies = dependencies)[,6],
    pbase = sessioninfo::package_info(pkgs = "loaded", include_base = include_base, dependencies = dependencies)[,7]
  ) %>%
    dplyr::arrange(
      desc(pbase),
      desc(pattached)
    ) %>%
    knitr::kable(
      escape = TRUE,
      format = "latex",
      caption = caption_text,
      booktabs = TRUE,
      longtable = TRUE,
      col.names = c(
        "Package Name",
        "Version",
        "Is Package Attached?",
        "Is a Base Package?"
      ),
      align = "c"
    ) %>%
    kableExtra::kable_styling() %>%
    kableExtra::row_spec(
      row = 0,
      bold = TRUE
    ) %>%
    kableExtra::kable_styling(
      latex_options =  c("HOLD_position", "repeat_header"),
      repeat_header_continued = "\\textit{(continued on next page...)}"
    )
}
