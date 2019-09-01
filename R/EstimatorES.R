#' Wrapper function for computing the standard error of risk/performance measure estimators
#'
#' @param data xts object of the data
#' @param ... extra arugements to be passed to lower level functions, see details
#' @param estimator.fun a character string indicating which risk/performancem measure's standard error
#' should be computed. One of \code{"Mean"} (default), \code{"SD"}, \code{"VaR"}, \code{"ES"},
#' \code{"SR"}, \code{"SoR"}, \code{"STARR"}.
#' @param se.method a character string indicating which method should be used to compute
#' the standard error of the estimated standard deviation. One of \code{"none"} (default),
#' \code{"IFiid"}, \code{"IFcor"}, \code{"BOOTiid"}, \code{"BOOTcor"}. Currely, only \code{"IFiid"}
#' is implemented.
#'
#' @return the standard error of the specified risk/performance measure using the specified method
#'
#' @import PerformanceAnalytics
#' @import RPEIF
#'
#' @export
#' @author Xin Chen, \email{chenx26@uw.edu}
#'
#' @examples
#' data(edhec)
#' EstimatorSE(edhec,estimator.fun="SD",se.method="IFiid")
#' h2o.init()
#' EstimatorSE(edhec,estimator.fun="ES",se.method="IFcor")
#' # h2o.shutdown(prompt=FALSE)
#' EstimatorSE(edhec,estimator.fun="ES",se.method="BOOTiid", nsim=100)
#' EstimatorSE(edhec,estimator.fun="ES",se.method="BOOTcor", nsim=100,
#' sim = "fixed", l = round(nrow(edhec)/5))

EstimatorSE = function(data, ...,
                   estimator.fun = c("Mean","SD","VaR","ES","SR","SoR","STARR", "SoR", "LPM", "Omega", "SSD"),
                   se.method = c("none","IFiid","IFcor", "IFcor.h2o", "BOOTiid","BOOTcor"),
                   prewhiten = FALSE, adaptive = FALSE, a=0.3, b=0.7){

  if(adaptive)
    if(prewhiten)
      stop("Only one of prewhiten or adaptive can be TRUE.")

  estimator.fun = estimator.fun[1]
  se.method = se.method[1]
  myfun = switch(estimator.fun,
                 Mean = mean,
                 SD = sd,
                 VaR = VaR.hist,
                 ES = ES.hist,
                 SR = SR,
                 SoR = SoR,
                 STARR = STARR,
                 LPM = LPM,
                 OmegaRatio = OmegaRatio,
                 SSD = SSD,
                 RachevRatio = RachevRatio,
                 stop("The estimator.fun specified is not implemented yet, please contact Anthony Christidis (anthony.christidis@stat.ubc.ca) or submit an issue at the github repository")
  )
  myfun.IF = switch (estimator.fun,
                     Mean = IF.mean,
                     SD = IF.SD,
                     VaR = IF.VaR,
                     ES = IF.ES,
                     SR = IF.SR,
                     SoR = IF.SoR,
                     STARR = STARR,
                     LPM = IF.LPM,
                     OmegaRatio = IF.Omega,
                     SSD = IF.SSD,
                     RachevRatio = IF.RachR,
                     stop("The estimator.fun specified is not implemented yet, please contact Anthony Christidis (anthony.christidis@stat.ubc.ca) or submit an issue at the github repository")
  )

  res = switch (se.method,
    none = NULL,
    IFiid = SE.xts(data, SE.IF.iid, myfun, myfun.IF, ...),
    IFcor.h2o = SE.xts(data, SE.IF.cor.h2o, myfun, myfun.IF, ...),
    BOOTiid = SE.xts(data, SE.BOOT.iid, myfun, myfun.IF, ...),
    BOOTcor = SE.xts(data, SE.BOOT.cor, myfun, myfun.IF,...)
  )

  # Adaptive method
  if(se.method=="IFcor"){
    if(!adaptive){
      res = SE.xts(data, SE.IF.cor, myfun, myfun.IF, ...)
    } else{
        if(sum(ncol(data)>1)==0){
        ar1.param = arima(x=data, order=c(1,0,0), include.mean=TRUE)[[1]][1]
        IFcor = SE.xts(data, SE.IF.cor, myfun, myfun.IF, prewhiten=FALSE, ...)
        IFcorPW = SE.xts(data, SE.IF.cor, myfun, myfun.IF, prewhiten=TRUE, ...)
        if(0<=ar1.param & ar1.param<a)
          w = 0 else if(a<=ar1.param & ar1.param<=b)
            w = (ar1.param - a)/(b - a) else
              w = 1
        res = (1-w)*IFcor + w*IFcorPW
        names(res) = NULL
        } else{
          res = numeric()
          for(my.col in 1:ncol(data)){
            temp.data = data[, my.col]
            ar1.param = arima(x=temp.data, order=c(1,0,0), include.mean=TRUE)[[1]][1]
            IFcor = SE.xts(temp.data, SE.IF.cor, myfun, myfun.IF, prewhiten=FALSE, ...)
            IFcorPW = SE.xts(temp.data, SE.IF.cor, myfun, myfun.IF, prewhiten=TRUE, ...)
            if(0<=ar1.param & ar1.param<a)
              w = 0 else if(a<=ar1.param & ar1.param<=b)
                w = (ar1.param - a)/(b - a) else
                  w = 1
            res[my.col] = (1-w)*IFcor + w*IFcorPW
          }
          names(res) = colnames(data)
        }
      }
  }

  return(res)
}

#' compute the standard error for the xts object
#'
#' @param x the xts object of the data
#' @param se.fun the function used to compute the standard error
#' @param myfun the measure
#' @param myfun.IF the influence function of the measure
#' @param ... other parameters
#'
#' @return standard errors
#' @export
#' @author Xin Chen, \email{chenx26@uw.edu}
#'
#' @examples
#' data(edhec)
#' SE.xts(edhec, SE.IF.iid, sd, SD.IF)
SE.xts = function(x, se.fun, myfun, myfun.IF, prewhiten=FALSE, ...){
  if (is.vector(x) || is.null(ncol(x)) || ncol(x) == 1) {
    x <- as.numeric(x)
    #    if(na.rm) x <- na.omit(x)
    return(se.fun(x = x, myfun = myfun, myfun.IF = myfun.IF, prewhiten=FALSE, ...))
  }
  else {
    x <- coredata(x)
    #    if(na.rm) x <- na.omit(x)
    return(apply(x, 2, se.fun, myfun = myfun, myfun.IF = myfun.IF, prewhiten=FALSE, ... ))
  }
}

#' Compute the standard error using influence function approach for a vector
#'
#' @param x vector of data
#' @param myfun.IF influence function of the measure
#' @param ... other parameters used to compute the influence function
#'
#' @return standard error of the measure
#' @export
#' @author Xin Chen, \email{chenx26@uw.edu}
#'
#' @examples
#' SE.IF.iid(rnorm(100), SD.IF)
SE.IF.iid = function(x, myfun.IF, ...){
  N=length(x)
  x.IF = myfun.IF(x, ...)
  x.IF.2 = x.IF^2
  tmp = mean(x.IF.2)
  return(sqrt(tmp/N))
}

#' Compute the standard error using GLM-EN approach for serially correlated data using h2o package
#'
#' @param x the vector of data
#' @param myfun.IF the influene function of the measure
#' @param d maximum order of the polynomial
#' @param alpha.lasso weight for the elastic net
#' @param keep portion of frequencies to be used
#' @param ... other parameters
#'
#' @return the standard error of the measure
#' @export
#' @author Xin Chen, \email{chenx26@uw.edu}
#'

SE.IF.cor.h2o = function(x, myfun.IF, ..., d = 5, alpha.lasso = 0.5, keep = 1){
  data.IF = myfun.IF(x, ...)
  tmp = SE.GLM.LASSO(data.IF, d = d, alpha = alpha.lasso, keep = keep)
  return(sqrt(tmp))
}

#' Compute the standard error of the measure by iid bootstrapping
#'
#' @param x vector of data
#' @param myfun measure
#' @param nsim number of replicates
#' @param ... other parameters
#' @param myfun.IF not used
#'
#' @return standard error
#' @export
#' @author Xin Chen, \email{chenx26@uw.edu}
#'
#' @examples
#' SE.BOOT.iid(x = rnorm(100), myfun = mean, nsim = 100)
SE.BOOT.iid = function(x, myfun, myfun.IF, prewhiten=FALSE, ..., nsim = 100){
  res = boot(data = x, statistic = function(x,i,...) myfun(x[i],...), R = nsim, ... = ...)
  return(sd(res$t))
}

#' Compute the standard error of the measure by tsboot()
#'
#' @param x vector of data
#' @param myfun measure
#' @param nsim number of replicates
#' @param ... other parameters
#' @param myfun.IF not used
#' @param sim the type of simulation
#' @param l the length of the fixed block
#'
#' @return standard error
#' @export
#' @author Xin Chen, \email{chenx26@uw.edu}
#'
#' @examples
#' SE.BOOT.cor(x = rnorm(100), myfun = mean, nsim = 100)
SE.BOOT.cor = function(x, myfun, myfun.IF, prewhiten=FALSE, ..., nsim = 1000,
                       sim = "fixed", l = round(length(x)/5)){
  res = tsboot(tseries = x, statistic = function(x,...) myfun(x,...), R = nsim,
               sim = sim, l = l,...)
  return(sd(res$t))
}

#' Compute the standard error using GLM-EN approach for serially correlated data using glmnetRcpp
#'
#' @param x the vector of data
#' @param myfun.IF the influene function of the measure
#' @param d maximum order of the polynomial
#' @param alpha.lasso weight for the elastic net
#' @param keep portion of frequencies to be used
#' @param standardize whether to standardize data when computing standard error
#' @param ... other parameters
#'
#' @return the standard error of the measure
#' @export
#' @author Xin Chen, \email{chenx26@uw.edu}
#'

SE.IF.cor = function(x, myfun.IF, ..., return.coeffs = FALSE, d.GLM.EN = 5, alpha.EN = 0.5, keep = 1,
                     standardize = FALSE, prewhiten=FALSE){
  d = d.GLM.EN
  data.IF = myfun.IF(x, prewhiten=prewhiten, ...)
  tmp = SE.glmnet_exp(data.IF,  ..., standardize = standardize, return.coeffs = return.coeffs, d = d, alpha.EN = alpha.EN, keep = keep)
  if(return.coeffs){
    coeffs = tmp[[2]]
    tmp = tmp[[1]]
    return(list(sqrt(tmp), coeffs))
  }
  return(sqrt(tmp))
}









