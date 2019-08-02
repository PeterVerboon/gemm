
#' print method of object of class gemm
#'
#' @param x object of class gemm
#' @param ... additional pars
#' @param digits  number of digits
#'
#' @export
#'
print.gemm <- function(x, ..., digits=2) {
  
  options(digits = digits)
  res <- list()
  
  cat("
The dependent variable is:     ", (x$input$yvar), "
The predictor variable is:     ", (x$input$xvar), "
The model contains", length(x$input$mvars),"mediators:",x$input$mvars,"\n")
  
  if(length(x$input$xmmod) > 0) {cat("The moderators for the x-m path(s):",x$input$xmmod,"\n")
  } else {cat("No moderators for the x-m path(s)","\n")}
  if(length(x$input$mymod) > 0) {cat("The moderators for the m-y path(s):",x$input$mymod,"\n")
  } else {cat("No moderators for the m-y path(s)","\n")}
  if(length(x$input$cmvars) > 0) {cat("The covariates for the mediators:",x$input$cmvars, "\n")
  } else {cat("No covariates for the mediators","\n")}
  if(length(x$input$cyvars) > 0) {cat("The covariates for the dependent:",x$input$cyvars, "\n")
  } else {cat("No covariates for the dependent variable","\n")}
  
  cat("\n\n")
  
  cat("Explained variance (R-square) of the mediators and dependent variable:\n");
  table1 <- data.frame(R2 = x$output$Rsq)
  row.names(table1) <- NULL
  table1$Dependent <- c(x$input$mvars, x$input$yvar)
  table1 <- table1[,c(2,1)]
  table1[,2] <- format(round(table1[,2], digits = 3), nsmall = 2)
  pander::pander(table1, justify = c("left", "right"))
  cat("\n")
  res$R2 <- table1
  
  cat("Estimate of total effect");
  lmout <- x$output$parameterEstimates.total
  a <- summary(lmout)$coefficients
  b <- confint(lmout)
  a <- as.data.frame(cbind(a,b))
  row.names(a) <- NULL
  terms <- c("intercept", paste0(x$input$xvar, " --> ", x$input$yvar))
  a[,"effect"] <- terms
  table2 <- a[,c(7,1:6)]
  table2[,c(2:7)] <- format(round(table2[,c(2:7)], digits = 3), nsmall = 2)
  colnames(table2) <- c("effect","est","se","t","pvalue", "ci.lower", "ci.upper")
  pander::pander(table2, justify = c("left", "right", "right", "right","right","right","right"))
  cat("\n\n")
  res$totalEff <- table2
  
  cat("Estimates of a-paths");
  table3 <- x$output$parameterEstimates.apath
  row.names(table3) <- NULL
  names(table3) <- c("path",names(table3)[-1])
  if (is.null(x$input$xmmod)) {
    terms <- paste0(x$input$xvar, " --> ", x$input$mvars) 
  } else {terms <- c(paste0(x$input$xvar, " --> ", x$input$mvars),
                     paste0(x$input$xmmod, " --> ", x$input$mvars),
                     paste0(x$input$xvar, " x ",x$input$xmmod, " --> ", x$input$mvars)) }
  table3[,1] <- terms
  table3[,c(2:7)] <- format(round(table3[,c(2:7)], digits = 3), nsmall = 2)
  pander::pander(table3, justify = c("left", "right", "right", "right","right","right","right"))
  cat("\n\n")
  res$aPaths <- table3
  
  cat("Estimates of b-paths");
  table4 <- x$output$parameterEstimates.bpath
  row.names(table4) <- NULL
  names(table4) <- c("path",names(table4)[-1])
  if (is.null(x$input$mymod)) {
    terms <- paste0(x$input$mvars, " --> ", x$input$yvar) 
  } else {terms <- c(paste0(x$input$mvars, " --> ", x$input$yvar),
                     paste0(x$input$mymod, " --> ", x$input$yvar),
                     paste0(x$input$mvars," x ", x$input$mymod, " --> ", x$input$yvar))}
  table4[,1] <- terms
  table4[,c(2:7)] <- format(round(table4[,c(2:7)], digits = 3), nsmall = 2)
  pander::pander(table4, justify = c("left", "right", "right", "right","right","right","right"))
  cat("\n")
  res$bPaths <- table4
  
  if (!is.null(x$input$cmvars) | !is.null(x$input$cyvars)) {
    cat("Estimates of covariates");
    table5 <- x$output$parameterEstimates.covs
    row.names(table5) <- terms1 <- terms2 <- NULL
    names(table5) <- c("path",names(table5)[-1])
    if (!is.null(x$input$cyvars)) {
      terms1 <- paste0(x$input$cyvars, " --> ", x$input$yvar) }
    if (!is.null(x$input$cmvars)) {
      terms2 <- paste0(x$input$cmvars, " --> ", x$input$mvars) }
    terms <- c(terms1,terms2)
    table5[,1] <- terms
    table5[,c(2:7)] <- format(round(table5[,c(2:7)], digits = 3), nsmall = 2)
    pander::pander(table5, justify = c("left", "right", "right", "right","right","right","right"))
    cat("\n")
    res$covs <- table5
  }
  
  
  cat("Direct effect (c') ");
  table6 <- x$output$parameterEstimates.direct
  row.names(table6) <- NULL
  table6[,1] <- paste0(x$input$xvar," --> ", x$input$yvar)
  names(table6) <- c("path",names(table6)[-1])
  table6[,c(2:7)] <- format(round(table6[,c(2:7)], digits = 3), nsmall = 2)
  pander::pander(table6, justify = c("left", "right", "right", "right","right","right","right"))
  cat("\n")
  res$directEff <- table6
  
  cat("Indirect effects (a*b) ");
  table7 <- x$output$parameterEstimates.indirect.raw
  row.names(table7) <- NULL
  names(table7) <- c("through",names(table7)[-1])
  table7[,1] <- c(x$input$mvars, "total")
  table7[,c(2:7)] <- format(round(table7[,c(2:7)], digits = 3), nsmall = 2)
  pander::pander(table7, justify = c("left", rep("right",6))) 
  cat("\n")
  res$indirectEff <- table7
  
  cat("Completely standardized effect sizes  ");
  table8 <- x$output$parameterEstimates.indirect.es_std
  row.names(table8) <- NULL
  table8[,7] <- c(x$input$mvars, "total")
  table8 <- table8[,c(7,1:6)]
  names(table8)[1] <- c("through")
  table8[,c(2:7)] <- format(round(table8[,c(2:7)], digits = 3), nsmall = 2)
  pander::pander(table8, justify = c("left", rep("right",6))) 
  cat("\n")
  res$csES <- table8
  
  cat("Ratio (indirect to total) based effect sizes  ");
  table9 <- x$output$parameterEstimates.indirect.es_rat
  row.names(table9) <- NULL
  table9[,7] <- c(x$input$mvars, "total")
  table9 <- table9[,c(7,1:6)]
  names(table9)[1] <- c("through")
  table9[,c(2:7)] <- format(round(table9[,c(2:7)], digits = 3), nsmall = 2)
  pander::pander(table9, justify = c("left", rep("right",6))) 
  cat("\n")
  res$ratES <- table9
  
  return(res)
}


