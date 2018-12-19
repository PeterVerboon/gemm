#' Makes Index of Mediated Moderated plots 
#'
#' @param data data frame containg the variables of the model
#' @param xvar predictor variable name
#' @param yvar depedendent variable name
#' @param mod moderator name
#' @param mvars vector of mediators names
#' @param parEst parameter estimates from lavaan results
#' @param vdichotomous indicates whether moderator is dichotomous (TRUE)
#' @param modLevels levels of dichotomous moderator
#' @param path which path is used
#' @import ggplot2
#' @return empty, directly plots all simple slopes and all indices of mediation
#' @export

prepPlotIMM <- function(data,xvar,yvar,mod, mvars, parEst, vdichotomous,
                              modLevels, path = NULL) {
  
  xquant <- quantile(data[,xvar], c(.16,.84), na.rm = TRUE)
  yquant <- quantile(data[,yvar], c(.16,.84), na.rm = TRUE)
  
  # compute simple slopes

    if (vdichotomous) {
       modquant <- c(0,1)
    } else {
       modquant <- quantile(data[,mod], c(.16,.84), na.rm = TRUE)
    }

  if (path == "x-m") {
    vorw <- "w"
    inter <- "im"
    modmed <- "modmedx"
  } else {
    vorw <- "v"
    inter <- "iy"
    modmed <- "modmedm"
  }
  
  a <- subset(parEst, grepl("a", parEst$label))[,"est"]
  b <- subset(parEst, grepl("b", parEst$label))[,c("ci.lower","est","ci.upper")]
  ind <- subset(parEst, grepl("ind", parEst$label))[,c("ci.lower","est","ci.upper")]
  
  vw <- subset(parEst, grepl(vorw, parEst$label))[,c("ci.lower","est","ci.upper")]
  int <- subset(parEst, grepl(inter, parEst$label))[,c("ci.lower","est","ci.upper")]
  mm <- subset(parEst, grepl(modmed, parEst$label))[,c("ci.lower","est","ci.upper")]
  
  bw <- subset(parEst, grepl("bw", parEst$label))[,c("ci.lower","est","ci.upper")]
  gw <- subset(parEst, grepl("gw", parEst$label))[,c("ci.lower","est","ci.upper")]
  
  N <- dim(data)[1]
  
  if (vorw == "v") bw <- (matrix(as.numeric(vw), nrow=length(mvars), ncol = 3, byrow = TRUE ))
  
  if (vdichotomous) {
    legendLabel <- modLevels
    # maxLab <- max(stringr::str_length(modLevels[1]),stringr::str_length(modLevels[2]))
    # minLabel <- stringr::str_pad(modLevels[1],maxLab, pad = " ")
    # maxLabel <- stringr::str_pad(modLevels[2],maxLab, pad = " ")
  }
  else {
    legendLabel <- c("16th percentile", "84th percentile")
    # minLabel <- c("for 16th percentile of moderator: ")
    # maxLabel <- c("for 84th percentile of moderator: ")
  }

  # initialize data for index mediated moderation
 
  plotData <- data.frame(X1 = numeric(),X2 = numeric(),X3 = numeric(), 
                         moderator = numeric(), 
                         mediator = factor())
  moderator <- as.numeric(data[,mod])
  
  # initialize data for mediated simple slopes
  
  plotDat2 <- data.frame(yv=numeric(), xv=numeric(), mov=numeric(),mev=factor())
  xv <-  c(xquant[1],xquant[1],xquant[2],xquant[2])
  mov <- c(modquant,modquant)
 
  # loop over mediators 
  
  for (i in seq_along(mvars)) {

   # index of moderated mediation
      
        d1 <-  matrix(as.numeric(ind[i,]), nrow = N, ncol=3, byrow = TRUE)  
        d2 <- (moderator %o% as.numeric(mm[i,]))
        yIom <- d1+ d2
      mediator <- rep(mvars[i],nrow(data))
      plotDat0 <- data.frame(yIom,moderator,mediator);
      plotData <- rbind(plotData,plotDat0)

    # tableRes <- matrix(c(round(sort(as.numeric(slopemin)), digits = 3),
    #                      round(sort(as.numeric(slopemax)), digits = 3)),
    #                      nrow = 2, byrow = TRUE)
    # colnames(tableRes) <- c("ci.lower", "est", "ci.upper")
    # rownames(tableRes) <- c(minLabel, maxLabel)
    # print(tableRes, digits = 3, quote = FALSE, row.names = TRUE)
    
  # mediated simple slopes   
      
      d1 <-  matrix(as.numeric(ind[i,]), nrow = 2, ncol=3, byrow = TRUE)  
      d2 <-  (modquant %o% as.numeric(mm[i,]))
      yIom2 <- d1+ d2
     
    pred1 <-  as.numeric(modquant) %o%  as.numeric(bw[i,])  + yIom2*xquant[1]
    pred2 <-  as.numeric(modquant) %o%  as.numeric(bw[i,])  + yIom2*xquant[2]
    pred <- rbind(pred1, pred2)
    plotDat1 <- data.frame(cbind(pred, xv = xv))
    plotDat1$mov <- as.factor(round(mov,1))                      
    plotDat1$mev <- as.factor(rep(mvars[i],4))
    
    plotDat2 <- rbind(plotDat2,plotDat1)
    
     }  # loop mvars
  
  
  if (!vdichotomous) {
    
      names(plotData) <- c("IMM_lwr",'IMM',"IMM_upr", mod, "mediator")
      ymin <- min(plotData$IMM,plotData$IMM_lwr,plotData$IMM_upr, yquant, na.rm = TRUE)
      ymax <- max(plotData$IMM,plotData$IMM_lwr,plotData$IMM_upr, yquant, na.rm = TRUE)

      plot_indexOfmediation <- ggplot(plotData, aes_string(x=mod,y="IMM",colour = "mediator")) +
       geom_line(aes(colour = mediator, group = mediator)) +      
       coord_cartesian(ylim=c(ymin, ymax)) +
       ggtitle("Index of moderated mediation") +
       xlab(paste0("Moderator: ",mod))
      
      plot_indexOfmediation <- plot_indexOfmediation + 
        geom_ribbon(aes(ymin=IMM_lwr, ymax=IMM_upr), alpha=.3, linetype=0) 
    
      print(plot_indexOfmediation)
  }
  
  if (vdichotomous) {
    
    names(plotData) <- c("IMM_lwr",'IMM',"IMM_upr", mod, "mediator")
    ymin <- min(plotData$IMM,plotData$IMM_lwr,plotData$IMM_upr, yquant, na.rm = TRUE)
    ymax <- max(plotData$IMM,plotData$IMM_lwr,plotData$IMM_upr, yquant, na.rm = TRUE)
    pd <- position_dodge(0.1)
    
    plotData[,mod] <- as.factor(plotData[,mod])
    levels(plotData[,mod]) <- legendLabel
    
    plot_indexOfmediation <- ggplot(plotData, aes_string(x=mod,y="IMM",colour = "mediator")) +
      geom_point(aes(colour = mediator), position = pd, size=3) +  
      geom_errorbar(aes(ymin=IMM_lwr, ymax=IMM_upr), position = pd) +
      coord_cartesian(ylim=c(ymin, ymax)) +
      ggtitle("Index of moderated mediation") +
      xlab(paste0("Moderator: ",mod)) 

    
    print(plot_indexOfmediation)
  }
  
    names(plotDat2) <- c("lwr",yvar,"upr", xvar, mod, "mediator")
    ymin <- min(plotDat2$yvar, plotDat2$lwr,plotDat2$upr, yquant)
    ymax <- max(plotDat2$yvar, plotDat2$lwr,plotDat2$upr, yquant)
    
    plot_simpleSlopes <- ggplot(plotDat2, aes_string(x=xvar,y=yvar,group= mod, colour=mod)) +
       geom_point() + geom_line() +
       geom_ribbon(aes(ymin=lwr, ymax=upr),alpha=.3, linetype=0) +
       ylim(ymin,ymax) +
       theme(plot.title = ggplot2::element_text(lineheight=.4, face="italic")) +
       ggtitle(paste0("Simple slopes in ", path , " path for indirect effect ")) +
       scale_colour_discrete(name  = mod, labels=legendLabel) 
       
      plot_simpleSlopes <- plot_simpleSlopes + facet_grid(rows=vars(mediator))
  
    print(plot_simpleSlopes)
  
  
  
  return()

} # end function





