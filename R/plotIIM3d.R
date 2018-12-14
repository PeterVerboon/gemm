
#' Makes plots of Index of moderated mediation of moderatedMediationSem object
#'
#' @param data data frame containg the moderators
#' @param xmmod moedraotr of x-m path  
#' @param mymod moderator of m-y path
#' @param mvars vector of mediators names
#' @param res results of moderatedMediationSem function
#' @return empty, directly plots all indices of mediation
#' @export
#' @example plotIIM3d(data=gemmDat, xmmod = "mod1", mymod = "mod2", mvars = mvars, res = res)

plotIIM3d <- function(data, xmmod, mymod, mvars = mvars, res = res) {
  
  Modxm <- quantile(as.numeric(data[,xmmod]), c(.10,.20,.40,.60,.80,.90))
  Modmy <- quantile(as.numeric(data[,mymod]), c(.10,.20,.40,.60,.80,.90))
  parEst <- lavaan::parameterestimates(res$intermediate$result)
  mm <- expand.grid(x = Modxm, y = Modmy)
  
  for (i in seq_along(mvars)){
    
    df <- IMM3d(mm$x, mm$y, parEst = parEst, i=i )
    z <- matrix(df[,2], nrow = length(Modxm), ncol = length(Modxm))
    upzlim <- max(max(z),.4) + .1
    lwzlim <- min(min(z),-.4) - .1
    
    persp(x=M1, y = M2, z = z, zlab = "IMM", xlab = "Mod1", ylab ="Mod2", 
          main = paste0("Index of Moderated Mediation for mediator: ", mvars[i]),
          theta = 30, phi = 30,axes = TRUE, scale = TRUE, zlim = c(lwzlim,upzlim),
          ticktype = "detailed",nticks = 4,
          cex.lab = .8, cex.axis = .5, col = "lightgrey", shade = 0.3, 
          expand = 1, d=2, r=5, border = NA)
  }
  invisible(z)
}


#' Computes Index of moderated mediation of moderatedMediationSem object
#'
#' @param data data frame containg the moderators
#' @param xmmod moedraotr of x-m path  
#' @param mymod moderator of m-y path
#' @param mvars vector of mediators names
#' @param parEst parameter estimates from lavaan results
#' @return vector of index of moderated mediation with 95%CI limits for a given mediator
#' @export

IMM3d <- function(M1, M2, parEst=parEst, i=1) {
  
  ind <- subset(parEst, grepl("ind", parEst$label))[,c("ci.lower","est","ci.upper")]
  mmx <- subset(parEst, grepl("modmedx", parEst$label))[,c("ci.lower","est","ci.upper")]
  mmy <- subset(parEst, grepl("modmedm", parEst$label))[,c("ci.lower","est","ci.upper")]
  im <- subset(parEst, grepl("im", parEst$label))[,c("ci.lower","est","ci.upper")]
  iy <- subset(parEst, grepl("iy", parEst$label))[,c("ci.lower","est","ci.upper")]
  
  N <- length(M1)
  p1 <- rep(1,N) %o% (unlist(ind[i,])) 
  p2 <- M1 %o% as.numeric(mmx[i,])
  p3 <- M2 %o% as.numeric(mmy[i,])
  p4 <- (M2 * M1) %o% (as.numeric(im[i,]) *  as.numeric(iy[i,2]))
  index <- p1 + p2 + p3 + p4

  return(index)
}













