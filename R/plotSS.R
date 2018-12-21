

#' Makes simple slope plots of gemm object
#'
#' @param x   object moderatedMediationSem
#' @param ... optional
#'
#' @return simple slope plots for each mediator and simple slopes parameter estimates
#' @export
#'
plotSS <- function(x,...) {
  
  data <- x$intermediate$data
  xmmod <- x$input$xmmod
  mymod <- x$input$mymod
  xvar <- x$input$xvar
  yvar <- x$input$yvar
  mvars <- x$input$mvars
  parEst <- x$intermediate$parameterEstimates
  predLevels <- x$intermediate$predLevels
  
  if ((!length(xmmod)) & (!length(mymod))) {
    return(message("No plots can be given, because no moderators have been specified"))
  }
  
  xquant <- quantile(data[,xvar],c(.16,.84), na.rm = TRUE)
  if (xquant[1] == xquant[2])  {
    if (!is.null(predLevels))  {
      xquant <- c(0,1)
    } else {
      return(message("The 16th and 84th precentile of the predictor are the same. No plots are made."))
    }
  }
  
  ## test if moderator exists for x=m path and if it is dichotomous factor
  if (length(xmmod)) {
    xdichotomous <- FALSE
    if (is.factor(data[,xmmod])) {
      if (length(levels(data[,xmmod])) > 2) {
        stop("This function can not yet plot moderation with a moderator (x-m path) that is a factor with more than two levels.");
      }
      else {
        xmodLevels <- levels(data[,xmmod]);
        data[,xmmod] <- as.numeric(data[,xmmod]) - 1;
        xdichotomous <- TRUE;
      }
    }
    
    yquant <- quantile(data[,yvar], c(.16,.84), na.rm = TRUE)
    
    prepPlotSS(data=data, xvar=xvar, yvar = yvar, mod = xmmod, mvars = mvars, parEst = parEst, 
                 vdichotomous = xdichotomous, modLevels = xmodLevels, predLevels = predLevels, 
                 xquant = xquant, yquant = yquant, path = "x-m")
  }
  
  ## test if moderator exists for m=y path and if it is dichotomous factor
  
  if (length(mymod)) {
    ydichotomous <- FALSE
    if (is.factor(data[,mymod])) {
      if (length(levels(data[,mymod])) > 2) {
        stop("This function can not yet plot moderation with a moderator (x-y path) that is a factor with more than two levels.")}
      else {
        ymodLevels <- levels(data[,mymod]);
        data[,mymod] <- as.numeric(data[,mymod]) - 1;
        ydichotomous <- TRUE;
      }
    }
    prepPlotSS(data=data, xvar=xvar, yvar = yvar, mod = mymod, mvars = mvars, parEst = parEst, 
                 vdichotomous = ydichotomous, modLevels = ymodLevels,predLevels = predLevels, 
                 xquant = xquant, yquant = yquant, path = "m-y")
  }
  
  invisible()
  
}  # end function



