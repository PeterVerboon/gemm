

#' Analyze moderated mediation model using SEM
#'
#' @param data data frame
#' @param xvar predictor variable, must be either numerical or dichotomous
#' @param mvars vector of names of mediator variables
#' @param yvar dependent variable, must be numerical
#' @param xmmod moderator of effect predictor on mediators, must be either numerical or dichotomous
#' @param mymod moderator of effect mediators on dependent variable, must be either numerical or dichotomous
#' @param cmvars covariates for mediators
#' @param cyvars covariates for dependent variable
#' @param estMethod estimation of standard errors method, bootstrap is default
#' @param nboot number of bootstrap samples
#'
#' @return gemm object
#' @export
#'
#' @examples
#' data("gemmDat")
#' res <- gemm(dat = gemmDat, xvar="x1", mvars= c("m1","m2","m3"),
#'        yvar = "y1", xmmod = "mod1", mymod= "bimod2",
#'        cmvars =c("c1","c2"), cyvars =c("c1","c2"), nboot=50)
#' print(res)
#' plot(res)

                 gemm <- function(data = NULL,
                                  xvar,
                                  mvars,
                                  yvar,
                                  xmmod = NULL,
                                  mymod = NULL,
                                  cmvars = NULL,
                                  cyvars = NULL,
                                  estMethod = "bootstrap",
                                  nboot = 1000) {

  res <- list(input = as.list(environment()),
              intermediate = list(),
              output = list());

  res$intermediate$dataName <- as.character(deparse(substitute(data)));

  res$intermediate$numberOfMediators <-
    nm <- length(mvars);
  
  ## check if predictor is dichotomous factor 
  if (is.factor(data[,xvar])) {
      if (nlevels(data[,xvar]) == 2) {
        res$intermediate$predLevels <- levels(data[,xvar])
        data[,xvar] <- as.numeric(data[,xvar])
      } else {
        return(message("Predictor is a factor with more than two levels"))
      }
  }
 
  ## initialize xmint and myint
  xmint <- NULL
  myint <- NULL
  
  
  ## check if there is a moderator for the x - m path
  if (!is.null(xmmod)) {
    if (length(xmmod) > 1) {
      return(message("This function can only handle one moderator for the x-m path."));
    }
    if (is.factor(data[,xmmod])) {
      if (nlevels(data[,xmmod]) > 2) {
        return(message("This function can not yet deal with categorical moderators 
                       with more than two levels."));
      } else {
        res$intermediate$xdichotomous <- xdichotomous <- TRUE;
        data[,"xmodOriginal"] <- data[,xmmod];
        data[,xmmod] <- as.numeric(data[,xmmod]) - 1;
      }
    } else {
      if (length(unique(data[,xmmod])) == 2) {
        res$intermediate$xdichotomous <-  xdichotomous <- TRUE
        levels(data[,xmmod]) <- unique(data[,xmmod])
        } else {
      res$intermediate$xdichotomous <- xdichotomous <- FALSE
        }
    }

    xmint <- paste0("xmInteraction",c(1:nm));
    xmterms <- paste0(paste0("data$",xmmod,"*","data$",mvars));

    for (i in 1:nm) {
      data[,xmint[i]] <- eval(parse(text = xmterms[i]));
    }

  }

  ### check if there is a moderator for the m - y path;
  if (!is.null(mymod)) {
    if (length(mymod) > 1) {
      return("This function can only handle one moderator for the m-y path.");
    }
    if (is.factor(data[,mymod])) {
      if (nlevels(data[,mymod]) > 2) {
        return("This function can not yet deal with categorical moderators with more than two levels.");
      } else {
        res$intermediate$ydichotomous <-  ydichotomous <- TRUE;
        data[,"ymodOriginal"] <- data[,mymod];
        data[,mymod] <- as.numeric(data[,mymod]) - 1;
      }
    } else {
      if (length(unique(data[,mymod])) == 2) {
        res$intermediate$ydichotomous <-  ydichotomous <- TRUE
        levels(data[,mymod]) <- unique(data[,mymod])
      } else {
      res$intermediate$ydichotomous <-  ydichotomous <- FALSE
      }
    }

    myint <- paste0("myInteraction",c(1:nm));
    myterms <- paste0(paste0("data$",mymod,"*","data$",mvars));

    for (i in 1:nm) {
      data[,myint[i]] <- eval(parse(text = myterms[i] ));
    }
  }

  res$intermediate$data <- data
  
  ### Build lavaan model
  res$intermediate$model <-
    buildModMedSemModel(xvar=xvar,
                        mvars= mvars,
                        yvar = yvar,
                        xmmod = xmmod,
                        mymod = mymod,
                        cmvars = cmvars,
                        cyvars = cyvars);

  ### Run SEM
  res$intermediate$result <- result <-
    lavaan::sem(res$intermediate$model,
        data=data,
        fixed.x = FALSE,
        std.lv = TRUE,
        se = estMethod,
        bootstrap=nboot);


  ### Extract R squared values
  res$output$Rsq <-
    lavaan::inspect(res$intermediate$result, "r2");
  
  ### Extract parameter estimates for a and b paths
  res$intermediate$parameterEstimates <- r1 <-
    lavaan::parameterestimates(result);
  res$output$parameterEstimates.apath <-
    r1[(r1[,"lhs"] %in% mvars & r1[,"rhs"] %in% c(xvar,xmmod,xmint)),-c(1:3)]
  res$output$parameterEstimates.bpath <-
    r1[(r1[,"lhs"] %in% yvar & r1[,"rhs"] %in% c(mvars,mymod,myint)),-c(1:3)]
  

  ### Extract parameter estimates for direct effects
  res$output$parameterEstimates.direct <-
      r1[(r1[,"lhs"] %in% yvar & r1[,"rhs"] %in% xvar),-c(1:3)]

  ### ... And for indirect effects
  a2 <- 1:nm;
  ind <- paste0("ind", a2 );
  indinter <- paste0("indinter", a2 );
  res$output$parameterEstimates.indirect.raw <- 
      r1[(r1[,"lhs"] %in% c(ind,indinter, "total")),-c(1:3)]

  ### ... And for the covariates
  res$output$parameterEstimates.covs <- 
    r1[(r1[,"rhs"] %in% c(cmvars,cyvars)) & (r1[,"label"] != ""),-c(1:3)]
  
  ### ... And the standardized indirect effects
   aa <- lavaan::lavInspect(result, "std")$beta[yvar,mvars] * lavaan::lavInspect(result, "std")$beta[mvars, xvar];
   names(aa) <- mvars
   res$output$parameterEstimates.indirect.standardized <- aa

  class(res) <- "gemm";

  return(res);

}



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
  
  cat("Direct effect (c) ");
  d <- x$output$parameterEstimates.direct
  row.names(d) <- NULL
  d[,1] <- paste0(x$input$xvar," --> ", x$input$yvar)
  d[,c(2:7)] <- format(round(d[,c(2:7)], digits = 4), nsmall = 1)
  pander::pander(d, justify = c("left", "right", "right", "right","right","right","right"))
  cat("\n")
  
  cat("Indirect effects (c') ");
  ind <- x$output$parameterEstimates.indirect.raw
  ind$standardized <- x$output$parameterEstimates.indirect.standardized
  row.names(ind) <- NULL
  names(ind) <- c("through",names(ind)[-1])
  ind[,1] <- x$input$mvars
  ind[,c(2:8)] <- format(round(ind[,c(2:8)], digits = 4), nsmall = 1)
  pander::pander(ind, justify = c("left", rep("right",7))) 
  cat("\n")
  

}









