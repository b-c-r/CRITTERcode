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
  
  fmax = 20
  nhalf = 50
  
  N <- seq(0, maxN, length = 1000)
  t2 <- fmax * N / (nhalf + N)
  t2d1 <- fmax * N / (0.8*nhalf + N)
  t2d2 <- fmax * N / (0.5*nhalf + N)
  t2e1 <- 0.75*fmax * N / (nhalf + N)
  t2e2 <- 0.50*fmax * N / (nhalf + N)
  
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
    mfrow=c(2,2),
    oma = c( .25,   .25,  .25,  .25),
    mar = c(4.25, 4.25, 1, 1),
    las = 1,
    pty = "m"
  )

  ##############################################################################
  # Plot (a)
  ##############################################################################
  
  plot(
    N,
    t2,
    xlab = "",
    ylab = "",
    type = "l",
    xlim = xlim,
    ylim = ylim,
    xaxs = "i",
    yaxs = "i",
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
  
  title(
    ylab = expression(paste("feeding rate ", F, sep ="")),
    line = 2.5,
  )
  title(
    xlab = expression(paste("prey density ", N, sep ="")),
    line = 2.5,
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
    xlab = "",
    ylab = "",
    type = "n",
    xlim = xlim,
    ylim = ylim,
    xaxs = "i",
    yaxs = "i",
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
    cex = 0.9
  )
  
  title(
    ylab = expression(paste("feeding rate ", F, sep ="")),
    line = 2.5,
  )
  title(
    xlab = expression(paste("prey density ", N, sep ="")),
    line = 2.5,
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
    t2,
    xlab = "",
    ylab = "",
    type = "l",
    xlim = xlim,
    ylim = ylim,
    xaxs = "i",
    yaxs = "i",
    bty = ifelse(journal_style, "l", "o"),
    tck = ifelse(journal_style, -0.01, NA),
    lwd=3
  ) # create the plot
  
  lines(N, t2d1, col = "darkgrey", lty = 2, lwd = 2)
  lines(N, t2d2, col = "lightgrey", lty = 3, lwd = 2)
  
  legend(
    "topright",
    legend = c(
      expression(paste(N[half], " = 50, a = 0.4", sep = "")),
      expression(paste(N[half], " = 40, a = 0.5", sep = "")),
      expression(paste(N[half], " = 25, a = 0.8", sep = "")),
      expression(paste(F[max], " = 20, ", T[h], " = 0.05", sep = ""))
    ),
    lty = c(1,2,3, 1),
    col = c("black","darkgrey","lightgrey", "white"),
    cex = 0.9
  )
  
  title(
    ylab = expression(paste("feeding rate ", F, sep ="")),
    line = 2.5,
  )
  title(
    xlab = expression(paste("prey density ", N, sep ="")),
    line = 2.5,
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
    N,
    t2,
    xlab = "",
    ylab = "",
    type = "l",
    xlim = xlim,
    ylim = ylim,
    xaxs = "i",
    yaxs = "i",
    bty = ifelse(journal_style, "l", "o"),
    tck = ifelse(journal_style, -0.01, NA),
    lwd=3
  ) # create the plot
  
  lines(N, t2e1, col = "darkgrey", lty = 2, lwd = 2)
  lines(N, t2e2, col = "lightgrey", lty = 3, lwd = 2)
  
  legend(
    "topright",
    legend = c(
      expression(paste(F[max], " = 20, ", T[h], " = 0.050, a = 0.4", sep = "")),
      expression(paste(F[max], " = 15, ", T[h], " = 0.066, a = 0.3", sep = "")),
      expression(paste(F[max], " = 10, ", T[h], " = 0.100, a = 0.2", sep = "")),
      expression(paste(N[half], " = 50", sep = ""))
    ),
    lty = c(1,2,3, 1),
    col = c("black","darkgrey","lightgrey", "white"),
    cex = 0.9
  )
  
  title(
    ylab = expression(paste("feeding rate ", F, sep ="")),
    line = 2.5,
  )
  title(
    xlab = expression(paste("prey density ", N, sep ="")),
    line = 2.5,
  )
  
  # add plot identifier
  mtext(
    "(d)",
    adj = .05,
    line = -1.5
  )
  
  
  
  
  ##############################################################################
  # Plot Inlay 
  ##############################################################################
  
  par(fig=c(0.7,0.95,0.83,0.99),
      new=T,
      oma = c(0,0,0,0),
      mar = c(2,2,2,2),
      pty = "m",
      cex.axis = 0.75)
  
  plot(
    N,
    tg[[1]]/N,
    xlab = expression(paste("prey density ", N, sep ="")),
    ylab = expression(paste("predation risk ", F, "/",N, sep ="")),
    type = "l",
    xlim = xlim,
    ylim = c(0,.41),
    bty = "l",
    tck = ifelse(journal_style, -0.01, NA),
    lwd=1,
    xaxt = "n",
    yaxt = "n"
  ) # create the plot
  
  
  for(i in 1:length(q)){
    lines(N, tg[[i]]/N, col = rainbow(length(q))[i])
  }
  
  axis(1, at = seq(0,200,50), labels = rep("",5), tck = -0.01)
  axis(2, at = seq(0,0.4,0.1), labels = rep("",5), tck = -0.01)
  
  text(
    seq(0,200,50),
    -0.075,
    adj = 0.5,
    labels = seq(0,200,50),
    xpd = T,
    cex = 0.75
  )
  text(
    -25,
    seq(0,0.4,0.1),
    adj = 0.5,
    labels = seq(0,0.4,0.1),
    xpd = T,
    cex = 0.75
  )
  title(
    ylab = "predation risk F/N",
    line = 1.2,
    cex.lab = 0.85
  )
  title(
    xlab = "prey density N",
    line = .5,
    cex.lab = 0.85
  )

  
}
