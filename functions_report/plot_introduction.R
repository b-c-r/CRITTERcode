################################################################################
#    plot_introduction: creates a explanatory plot for the introduction        #
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

plot_introduction <- function(
  fmax = 20,                                                                    # maximum feeding rate
  nhalf = 50,                                                                   # half saturation density
  q = c(0, 0.2, 0.5, 1),                                                        # q values for gen FR
  maxN = 200,                                                                   # maximum prey density
  xres = 1000,                                                                  # number of points for lines
  ylim = c(0, 30),
  xlim = c(0, 200),
  journal_style = TRUE                                                          # should the style according to journal guidelines used? Default: FALSE
){
  
  ##############################################################################
  # simulate a type II functional response
  ##############################################################################
  
  N <- seq(0, maxN, length = 1000)
  t2 <- fmax * N / (nhalf + N)
  t2e1 <- 0.75*fmax * N / (nhalf + N)
  t2e2 <- 0.50*fmax * N / (nhalf + N)
  t2f1 <- fmax * N / (0.8*nhalf + N)
  t2f2 <- fmax * N / (0.4*nhalf + N)
  
  ##############################################################################
  # simulate a generalized functional response
  ##############################################################################
  
  tg <- foreach(i = 1:length(q)) %do% {
    fmax * N^(1+q[i]) / (nhalf^(1+q[i]) + N^(1+q[i]))
  }
  
  ##############################################################################
  # Plot Settings
  ##############################################################################
  
  par(
    mfrow=c(2,3),
    oma = c( .1,   .1,  .1,  .1),
    mar = c(4.25, 4.25, .25, .25),
    las = 1,
    pty = "s"
  )

  ##############################################################################
  # Plot (a)
  ##############################################################################
  
  plot(
    N,
    t2,
    xlab = expression(paste("prey density ", N, sep ="")),
    ylab = expression(paste("feeding rate ", F, sep ="")),
    type = "l",
    xlim = xlim,
    ylim = ylim,
    bty = ifelse(journal_style, "l", "o"),
    tck = ifelse(journal_style, -0.01, NA),
    lwd=3
  ) # create the plot
  
  # add saturation:
  lines(c(0, maxN), c(fmax, fmax))
  text(
    x = maxN*0.75,
    y = fmax+5,
    labels = expression(paste(F[max], " = ", frac(1,T[h]), sep =""))
  )
  
  # add half saturation:
  lines(c(nhalf, nhalf), c(0, fmax/2))
  #lines(c(0, nhalf), c(fmax/2, fmax/2))
  text(
    x = nhalf+5,
    y = fmax/3,
    labels = expression(paste(N[half], " = ", frac(1,aT[h]), sep ="")),
    pos = 4
  )
  
  # add attack rate:
  xa <- c(0, nhalf*0.8)
  ya <- fmax/nhalf * xa
  lines(xa, ya)
  text(
    x = xa[2],
    y = ya[2],
    labels = expression(paste(aN, sep ="")),
    pos = 4
  )

  # add plot identifier
  mtext(
    "(a)",
    adj = .05,
    line = -1.5
  )
  
  
  ##############################################################################
  # Plot (b)
  ##############################################################################
  
  plot(
    N,
    tg[[1]],
    xlab = expression(paste("prey density ", N, sep ="")),
    ylab = expression(paste("feeding rate ", F, sep ="")),
    type = "n",
    xlim = xlim,
    ylim = ylim,
    bty = ifelse(journal_style, "l", "o"),
    tck = ifelse(journal_style, -0.01, NA),
    lwd=3
  ) # create the plot
  
  
  for(i in 1:length(q)){
    lines(N, tg[[i]], col = rainbow(length(q))[i])
  }
  
  qnames <- foreach(
    i = 1:length(q),
    .combine = "c"
    ) %do% {
    paste("q = ", q[i], sep = "")
  }
  legend(
    "bottomright",
    legend = qnames,
    lty = 1,
    col = rainbow(length(q)),
    cex = 0.75
  )
  
  # add plot identifier
  mtext(
    "(b)",
    adj = .05,
    line = -1.5
  )
  
  
  ##############################################################################
  # Plot (c)
  ##############################################################################
  
  plot(
    N,
    tg[[1]]/N,
    xlab = expression(paste("prey density ", N, sep ="")),
    ylab = expression(paste("predation risk ", F, "/",N, sep ="")),
    type = "n",
    xlim = xlim,
    ylim = c(0,1),
    bty = ifelse(journal_style, "l", "o"),
    tck = ifelse(journal_style, -0.01, NA),
    lwd=3
  ) # create the plot
  
  
  for(i in 1:length(q)){
    lines(N, tg[[i]]/N, col = rainbow(length(q))[i])
  }
  
  qnames <- foreach(
    i = 1:length(q),
    .combine = "c"
  ) %do% {
    paste("q = ", q[i], sep = "")
  }
  legend(
    "right",
    legend = qnames,
    lty = 1,
    col = rainbow(length(q)),
    cex = 0.75
  )
  
  # add plot identifier
  mtext(
    "(c)",
    adj = .05,
    line = -1.5
  )

  ##############################################################################
  # Plot (d)
  ##############################################################################
  
  plot(
    c(0,1),
    c(0,1),
    xlab = "",
    ylab = "",
    xaxt = "n",
    yaxt = "n",
    type = "n",
    bty = "n",
    xaxs = "i"
  ) # create the plot
  
  gh_path <- "https://raw.githubusercontent.com/b-c-r/CRITTERdata/refs/heads/main/"
  f_path <- "pictures/"
  png_url <- paste(gh_path, f_path, "all_habitat_levels_numbered.png", sep ="")
  
  graphics::rasterImage(
    png::readPNG(RCurl::getURLContent(png_url)),
    0,
    0,
    1,
    1
  )
  
  # add plot identifier
  mtext(
    "(d)",
    adj = -.2,
    line = -1.5
  )
  
  ##############################################################################
  # Plot (e)
  ##############################################################################
  
  plot(
    N,
    t2,
    xlab = expression(paste("prey density ", N, sep ="")),
    ylab = expression(paste("feeding rate ", F, sep ="")),
    type = "l",
    xlim = xlim,
    ylim = ylim,
    bty = ifelse(journal_style, "l", "o"),
    tck = ifelse(journal_style, -0.01, NA),
    lwd=3
  ) # create the plot
  
  lines(N, t2e1, col = "darkgrey", lty = 2, lwd = 2)
  lines(N, t2e2, col = "lightgrey", lty = 3, lwd = 2)

  legend(
    "topright",
    legend = c(
      paste(expression(F[max]), " = ",     fmax, ", a = ",     fmax/nhalf, sep = ""),
      paste(expression(F[max]), " = ", .75*fmax, ", a = ", .75*fmax/nhalf, sep = ""),
      paste(expression(F[max]), " = ", .50*fmax, ", a = ", .50*fmax/nhalf, sep = ""),
      paste(expression(N[half]), " = ", nhalf, sep = "")
    ),
    lty = c(1,2,3, 1),
    col = c("black","darkgrey","lightgrey", "white"),
    cex = 0.5
  )
  
  # add attack rate:
  #xa  <- c(0, nhalf*0.8)
  #ya  <-      fmax/nhalf * xa
  #ya1 <- 0.75*fmax/nhalf * xa
  #ya2 <- 0.50*fmax/nhalf * xa
  #lines(xa, ya,  col = "black", lty = 1)
  #lines(xa, ya1, col = "darkgrey", lty = 2)
  #lines(xa, ya2, col = "lightgrey", lty = 3)
  
  # add plot identifier
  mtext(
    "(e)",
    adj = .05,
    line = -1.5
  )
  
  
  ##############################################################################
  # Plot (f)
  ##############################################################################
  
  plot(
    N,
    t2,
    xlab = expression(paste("prey density ", N, sep ="")),
    ylab = expression(paste("feeding rate ", F, sep ="")),
    type = "l",
    xlim = xlim,
    ylim = ylim,
    bty = ifelse(journal_style, "l", "o"),
    tck = ifelse(journal_style, -0.01, NA),
    lwd=3
  ) # create the plot
  
  lines(N, t2f1, col = "darkgrey", lty = 2, lwd = 2)
  lines(N, t2f2, col = "lightgrey", lty = 3, lwd = 2)
  
  legend(
    "topright",
    legend = c(
      paste(expression(N[half]), " = ",    nhalf, ", a = ", fmax/nhalf, sep = ""),
      paste(expression(N[half]), " = ", .8*nhalf, ", a = ", fmax/(.8*nhalf), sep = ""),
      paste(expression(N[half]), " = ", .4*nhalf, ", a = ", fmax/(.4*nhalf), sep = ""),
      paste(expression(F[max]), " = ", fmax, sep = "")
    ),
    lty = c(1,2,3, 1),
    col = c("black","darkgrey","lightgrey", "white"),
    cex = 0.5
  )
}
