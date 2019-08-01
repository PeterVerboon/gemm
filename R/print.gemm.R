
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
  a <- data.frame(R2 = x$output$Rsq)
  row.names(a) <- NULL
  a$Dependent <- c(x$input$mvars, x$input$yvar)
  a <- a[,c(2,1)]
  a[,2] <- format(round(a[,2], digits = 3), nsmall = 1)
  pander::pander(a, justify = c("left", "right"))
  cat("\n")
  
  cat("Estimate of total effect");
  lmout <- x$output$parameterEstimates.total
  a <- summary(lmout)$coefficients
  b <- confint(lmout)
  a <- as.data.frame(cbind(a,b))
  row.names(a) <- NULL
  terms <- c("intercept", paste0(x$input$xvar, " --> ", x$input$yvar))
  a$effect <- terms
  a <- a[,c(7,1:6)]
  a[,c(2:7)] <- format(round(a[,c(2:7)], digits = 3), nsmall = 2)
  colnames(a) <- c("effect","est","se","t","pvalue", "ci.lower", "ci.upper")
  pander::pander(a, justify = c("left", "right", "right", "right","right","right","right"))
  cat("\n\n")
  
  
  cat("Estimates of a-paths");
  a <- x$output$parameterEstimates.apath
  row.names(a) <- NULL
  names(a) <- c("path",names(a)[-1])
  if (is.null(x$input$xmmod)) {
    terms <- paste0(x$input$xvar, " --> ", x$input$mvars) 
  } else {terms <- c(paste0(x$input$xvar, " --> ", x$input$mvars),
                     paste0(x$input$xmmod, " --> ", x$input$mvars),
                     paste0(x$input$xvar, " x ",x$input$xmmod, " --> ", x$input$mvars)) }
  a[,1] <- terms
  a[,c(2:7)] <- format(round(a[,c(2:7)], digits = 4), nsmall = 1)
  pander::pander(a, justify = c("left", "right", "right", "right","right","right","right"))
  cat("\n\n")
  
  cat("Estimates of b-paths");
  b <- x$output$parameterEstimates.bpath
  row.names(b) <- NULL
  names(b) <- c("path",names(b)[-1])
  if (is.null(x$input$mymod)) {
    terms <- paste0(x$input$mvars, " --> ", x$input$yvar) 
  } else {terms <- c(paste0(x$input$mvars, " --> ", x$input$yvar),
                     paste0(x$input$mymod, " --> ", x$input$yvar),
                     paste0(x$input$mvars," x ", x$input$mymod, " --> ", x$input$yvar))}
  b[,1] <- terms
  b[,c(2:7)] <- format(round(b[,c(2:7)], digits = 4), nsmall = 1)
  pander::pander(b, justify = c("left", "right", "right", "right","right","right","right"))
  cat("\n")
  
  
  if (!is.null(x$input$cmvars) | !is.null(x$input$cyvars)) {
    cat("Estimates of covariates");
    cov <- x$output$parameterEstimates.covs
    row.names(cov) <- terms1 <- terms2 <- NULL
    names(cov) <- c("path",names(cov)[-1])
    if (!is.null(x$input$cyvars)) {
      terms1 <- paste0(x$input$cyvars, " --> ", x$input$yvar) }
    if (!is.null(x$input$cmvars)) {
      terms2 <- paste0(x$input$cmvars, " --> ", x$input$mvars) }
    terms <- c(terms1,terms2)
    cov[,1] <- terms
    cov[,c(2:7)] <- format(round(cov[,c(2:7)], digits = 4), nsmall = 1)
    pander::pander(cov, justify = c("left", "right", "right", "right","right","right","right"))
    cat("\n")
  }
  
  cat("Direct effect (c') ");
  d <- x$output$parameterEstimates.direct
  row.names(d) <- NULL
  d[,1] <- paste0(x$input$xvar," --> ", x$input$yvar)
  d[,c(2:7)] <- format(round(d[,c(2:7)], digits = 4), nsmall = 1)
  pander::pander(d, justify = c("left", "right", "right", "right","right","right","right"))
  cat("\n")
  
  cat("Indirect effects (a*b) ");
  ind <- x$output$parameterEstimates.indirect.raw
  #ind$standardized <- x$output$parameterEstimates.indirect.standardized
  row.names(ind) <- NULL
  names(ind) <- c("through",names(ind)[-1])
  ind[,1] <- c(x$input$mvars, "total")
  ind[,c(2:7)] <- format(round(ind[,c(2:7)], digits = 4), nsmall = 1)
  pander::pander(ind, justify = c("left", rep("right",6))) 
  cat("\n")
  
  
  cat("Completely standardized effect sizes  ");
  es2 <- x$output$parameterEstimates.indirect.es_std
  row.names(es2) <- NULL
  es2[,7] <- c(x$input$mvars, "total")
  es2 <- es2[,c(7,1:6)]
  names(es2)[1] <- c("through")
  es2[,c(2:7)] <- format(round(es2[,c(2:7)], digits = 4), nsmall = 1)
  pander::pander(es2, justify = c("left", rep("right",6))) 
  cat("\n")
  
  cat("Ratio (indirect to total) based effect sizes  ");
  es1 <- x$output$parameterEstimates.indirect.es_rat
  row.names(es1) <- NULL
  es1[,7] <- c(x$input$mvars, "total")
  es1 <- es1[,c(7,1:6)]
  names(es1)[1] <- c("through")
  es1[,c(2:7)] <- format(round(es1[,c(2:7)], digits = 4), nsmall = 1)
  pander::pander(es1, justify = c("left", rep("right",6))) 
  cat("\n")
  
}


