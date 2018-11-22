#' Makes simple slop plots of moderatedMediationSem object
#'
#' Description
#'
#' @param data data frame containg the variables of the model
#' @param xvar predictor variable name
#' @param yvar depedendent variable name
#' @param mod moderator name
#' @param mvars vector of mediators names
#' @param parEst parameter estimates from lavaan results
#' @param vorw character ("v" or "w") indicating for which path the simple slopes should be computed.
#' @param int estimate of interaction
#' @param vdichotomous indicates whether moderator is dichotomous (TRUE)
#' @param modLevels levels of dichotomous moderator
#' @param path which path is used
#' @param digits number of digits displayed
#' @import ggplot2
#' @return empty, directly plots all simple slopes and all indices of mediation
#' @export

simpleSlopes <- function(data,xvar,yvar,mod, mvars, parEst, vorw, int, vdichotomous,
                              modLevels, path = NULL, digits = 3) {
  
  xmin <- min(data[,xvar], na.rm = TRUE)
  xmax <- max(data[,xvar], na.rm = TRUE)
  miny <- min(data[,yvar], na.rm = TRUE)
  maxy <- max(data[,yvar], na.rm = TRUE)

  # compute simple slopes

    if (vdichotomous) {
      modmin <- 0;
      modmax <- 1
    } else {
      modmin <- mean(data[,mod]) - stats::sd(data[,mod]);
      modmax <- mean(data[,mod]) + stats::sd(data[,mod]);
    }

  a <- subset(parEst, grepl("a1", parEst$label))[,"est"]
  b <- subset(parEst, grepl("b", parEst$label))[,"est"]
  vw <- subset(parEst, grepl(vorw, parEst$label))[,"est"]
  int <- subset(parEst, grepl(int, parEst$label))[,c("ci.lower","est","ci.upper")]

  if (vorw == "v") vw <- rep(vw,length(mvars))

  title <- paste0("Simple slopes in ", path , " path for indirect effect ")

 # cat("\n","Simple slopes in ", path , " path(s) for indirect effect ", "\n")
 # cat(" ---------------------------------------------", "\n" );

  plotData <- data.frame(yIom=numeric(), moderator = numeric(), mediator = factor())
  moderator <- c(min(data[,mod]),max(data[,mod])) 
  
  plotDat2 <- data.frame(yv=numeric(), xv=numeric(), mov=numeric(),mev=factor())
  xv <- c(xmin,xmax,xmin,xmax)
  mov <- c(modmin, modmin, modmax, modmax)
  
  
  for (i in 1:length(mvars)) {

    incmin <- vw*modmin
    slopemin <- a[i]*b[i] + vw[i]*int[i,]*modmin
    incmax  <- vw*modmax
    slopemax <- a[i]*b[i] + vw[i]*int[i,]*modmax

    if (vdichotomous) {
      legendLabel <- modLevels
      maxLab <- max(stringr::str_length(modLevels[1]),stringr::str_length(modLevels[2]))
      minLabel <- stringr::str_pad(modLevels[1],maxLab, pad = " ")
      maxLabel <- stringr::str_pad(modLevels[2],maxLab, pad = " ")
    }
    else {
      legendLabel <- c("1 SD below mean", "1 SD above mean")
      minLabel <- c("for 1 sd below mean of moderator: ")
      maxLabel <- c("for 1 sd above mean of moderator: ")
      data$moderator <- data[,mod]

      yIom <- a[i]*b[i] + vw[i]*int[i,2]*moderator
      mediator <- c(mvars[i],mvars[i])
      plotDat0 <- data.frame(yIom,moderator,mediator);
      plotData <- rbind(plotData,plotDat0)

    }

    # tableRes <- matrix(c(round(sort(as.numeric(slopemin)), digits = 3),
    #                      round(sort(as.numeric(slopemax)), digits = 3)),
    #                      nrow = 2, byrow = TRUE)
    # colnames(tableRes) <- c("ci.lower", "est", "ci.upper")
    # rownames(tableRes) <- c(minLabel, maxLabel)
    # print(tableRes, digits = 3, quote = FALSE, row.names = TRUE)
    
    pred <- rep(0,4)
    pred[1] <- incmin  + slopemin[2] * xmin;
    pred[2] <- incmin + slopemin[2] * xmax;
    pred[3] <- incmax  + slopemax[2] * xmin;
    pred[4] <- incmax  + slopemax[2] * xmax
    pred <- unlist(pred)

    plotDat1 <- data.frame(cbind(yv = pred, xv = xv))
    
    plotDat1$mov <- as.factor(round(mov,1))                      
    plotDat1$mev <- as.factor(rep(mvars[i],4))
    
    plotDat2 <- rbind(plotDat2,plotDat1)
    
  

     }  # loop mvars
  
  if (!vdichotomous) {
    
    names(plotData) <- c('IMM', mod, "mediator")

    plot_indexOfmediation <- ggplot(plotData, aes_string(x=mod,y="IMM",colour = "mediator")) +
       geom_point(size=.5) + geom_line() +
       coord_cartesian(ylim=c(-0.5, 0.5)) +
       scale_y_continuous(breaks=seq(-0.5, 0.5, 0.1)) +
       ggtitle("Index of moderated mediation") +
       xlab(paste0("Moderator: ",mod))
  
    print(plot_indexOfmediation)
  }
  
    names(plotDat2) <- c(yvar, xvar, mod, "mediator")
  
    plot_simpleSlopes <- ggplot(plotDat2, aes_string(x=xvar,y=yvar,group= mod, colour=mod)) +
       geom_point() + geom_line() +
       #labs(x = xvar, y = yvar) +
       ylim(min(miny,min(plotDat2[,1])), max(maxy,max(plotDat2[,1]))) +
       theme(plot.title = ggplot2::element_text(lineheight=.8, face="bold")) +
       scale_colour_discrete(name  = mod, labels=legendLabel) 
       
      plot_simpleSlopes <- plot_simpleSlopes + facet_grid(rows=vars(mediator))
  
    print(plot_simpleSlopes)
  
  
  
  return()

} # end function





