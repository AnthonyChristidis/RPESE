#' Wrapper function that computes Rachev Ratio and the standard error of the estimate
#'
#' @param data data
#' @param \dots any other passthru parameters. This include two types of parameters.
#' The first type is parameters associated with the risk/performance measure, such as tail
#' probability for VaR and ES. The second type is the parameters associated with the metohd
#' used to compute the standard error. See \code{\link{SE.IF.iid}}, \code{\link{SE.IF.cor}},
#' \code{\link{SE.BOOT.iid}}, \code{\link{SE.BOOT.cor}} for details.
#' @param alpha Lower tail probability.
#' @param beta Upper tail probability.
#' @param rf Risk-free interest rate.
#' @param se.method a character string indicating which method should be used to compute
#' the standard error of the estimated standard deviation. One of \code{"none"} (default),
#' \code{"IFiid"}, \code{"IFcor"}, \code{"BOOTiid"}, \code{"BOOTcor"}.
#'
#' @return a vector or a list depending on se.method
#' @export
#'
#' @import PerformanceAnalytics
#' @import IFs
#'
#' @author Xin Chen, \email{chenx26@uw.edu}

Rachev.SE = function(data, ..., alpha=0.1, beta=0.1, rf=0, se.method = "none"){
  data = checkData(data)
  myRachevRatio = apply(data, 2, RachevRatio, const = const, k = k, ...)
  names(myRachevRatio) = colnames(data)
  if(se.method[1] == "none" & length(se.method)==1){
    return(myRachevRatio)
  } else {
    res=list(RachevRatio=myRachevRatio)
    # for each of the method specified in se.method, compute the standard error
    for(mymethod in se.method){
      res[[mymethod]]=EstimatorSE(data, estimator.fun = "RachevRatio",
                                  se.method = mymethod, alpha=0.1, beta=0.1, rf=0, ...)
    }
    return(res)
  }
}
