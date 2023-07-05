#' @export
#' @title getPortNTSParam
#' @description Portfolio return with capital allocation weight is \eqn{R_p=<w,r>},
#' which is a weighted sum of of elements in the N-dimensional NTS random vector.
#' \eqn{R_p} becomes an 1-dimensional NTS random variable.
#' \code{getPortNTSParam} find the parameters of \eqn{R_p}.
#'
#' @param strPMNTS Structure of parameters for the n-dimensional NTS distribution.
#'
#' \code{strPMNTS$ndim} : dimension
#'
#' \code{strPMNTS$mu} : \eqn{\mu} mean vector (column vector) of the input data.
#'
#' \code{strPMNTS$sigma} : \eqn{\sigma} standard deviation vector (column vector) of the input data.
#'
#' \code{strPMNTS$alpha} : \eqn{\alpha} of the std NTS distribution (X).
#'
#' \code{strPMNTS$theta} : \eqn{\theta} of the std NTS distribution (X).
#'
#' \code{strPMNTS$beta} : \eqn{\beta} vector (column vector) of the std NTS distribution (X).
#'
#' \code{res$Rho} : \eqn{\rho} matrix (Correlation) of the std NTS distribution (X).
#'
#' \code{res$Sigma} : Covariance \eqn{\Sigma} matrix of return data \eqn{r}.
#'
#' @param w Capital allocation weight vector.
#' @param stdform If \code{stdform} is \code{FALSE}, then the return parameter has the following representation
#'
#' \eqn{R_p = <w, r> = \mu + diag(\sigma) X},
#'
#' where
#'
#' \eqn{X} follows \eqn{stdNTS_1(\alpha, \theta, \beta, 1)}.
#'
#' If \code{stdform} is \code{TRUE}, then the return parameter has the following representation
#'
#' \eqn{R_p = <w, r>} follows \eqn{NTS_1(\alpha, \theta, \beta, \gamma, \mu, 1)}
#'
#' @usage
#' \code{res <- setPortfolioParam(strPMNTS,w)}
#' \code{res <- setPortfolioParam(strPMNTS,w, FALSE)}
#'
#' @return The weighted sum follows 1-dimensional NTS.
#'
#' \eqn{R_p = <w, r> = \mu + diag(\sigma) X},
#'
#' where
#'
#' \eqn{X} follows \eqn{stdNTS_1(\alpha, \theta, \beta, 1)}.
#'
#' Hence we obtain
#'
#' \code{res$mu} : \eqn{\mu} mean of \eqn{R_p}.
#'
#' \code{res$sigma} : \eqn{\sigma} standard deviation of \eqn{R_p}.
#'
#' \code{res$alpha} : \eqn{\alpha} of \eqn{X}.
#'
#' \code{res$theta} : \eqn{\theta} of \eqn{X}.
#'
#' \code{res$beta} : \eqn{\beta} of \eqn{X}.
#'
#' @references
#' Proposition 2.1 of
#' Kim, Y. S. (2020) Portfolio Optimization on the Dispersion Risk and the Asymmetric Tail Risk
#' \url{https://arxiv.org/pdf/2007.13972.pdf}
#'
#' @examples
#' library("temStaR")
#' strPMNTS <- list(ndim = 2,
#'                  mu = c( 9.876552e-05, 4.747343e-04 ),
#'                  sigma = c( 0.01620588, 0.02309643 ),
#'                  alpha = 0.1888129 ,
#'                  theta = 0.523042,
#'                  beta =  c( -0.04632938,  0.04063555 ),
#'                  Rho = matrix( data = c(1.0, 0.469883,
#'                                        0.469883, 1.0),
#'                                nrow = 2, ncol = 2)
#'                  CovMtx = matrix( data = c(0.0002626304, 0.0001740779,
#'                                          0.0001740779, 0.0005334452),
#'                                  nrow = 2, ncol = 2)
#'                  )
#' w <- c(0.3, 0.7)
#' res <- getPortNTSParam(strPMNTS,w)
#'
getPortNTSParam <- function(strPMNTS, w, stdform = TRUE){
  if(strPMNTS$ndim != length(w)){
    print("The dimension of the weight vector must be the same as strPMNTS$ndim")
    return(NULL)
  }
  w <- matrix(data = w, nrow = strPMNTS$ndim, ncol = 1)
  sigmaBar <- as.numeric(sqrt(t(w)%*%strPMNTS$CovMtx%*%w))
  mbar <- sum(w*strPMNTS$mu)
  betbar <- sum(w*strPMNTS$sigma*strPMNTS$beta)/sigmaBar
  stdntsparam <- c(as.numeric(strPMNTS$alpha), as.numeric(strPMNTS$theta), betbar)
  names(stdntsparam) <- c("alpha", "theta", "beta")
  param <- list(stdparam = stdntsparam, mu = mbar, sig = sigmaBar)
  if (stdform == FALSE){
    param <- change_stdntsparam2ntsparam(param$stdparam, mbar, sigmaBar)
    param <- param[1:5]
  }
  return(param)
}

#' @export
#' @title portfolioVaRmnts
#' @description
#' Calculate portfolio value at risk on the NTS market model
#'
#' @param strPMNTS Structure of parameters for the n-dimensional NTS distribution.
#'
#' \code{strPMNTS$ndim} : dimension
#'
#' \code{strPMNTS$mu} : \eqn{\mu} mean vector (column vector) of the input data.
#'
#' \code{strPMNTS$sigma} : \eqn{\sigma} standard deviation vector (column vector) of the input data.
#'
#' \code{strPMNTS$alpha} : \eqn{\alpha} of the std NTS distribution (X).
#'
#' \code{strPMNTS$theta} : \eqn{\theta} of the std NTS distribution (X).
#'
#' \code{strPMNTS$beta} : \eqn{\beta} vector (column vector) of the std NTS distribution (X).
#'
#' \code{res$Rho} : \eqn{\rho} matrix (Correlation) of the std NTS distribution (X).
#'
#' \code{res$Sigma} : Covariance \eqn{\Sigma} matrix of return data \eqn{r}.
#'
#' @param w Capital allocation weight vector.
#' @param eta significanlt level
#' @return portfolio value at risk on the NTS market model
#'
portfolioVaRmnts <- function(strPMNTS, w, eta){
  if(strPMNTS$ndim != length(w)){
    print("The dimension of the weight vector must be the same as strPMNTS$ndim")
    return(NULL)
  }
  ntsparam <- getPortNTSParam(strPMNTS, w)
  VaR <- (- ntsparam$mu - ntsparam$sig * qnts(eta, ntsparam$stdparam))
  return(VaR)
}


#' @export
#' @title portfolioCVaRmnts
#' @description
#' Calculate portfolio conditional value at risk (expected shortfall) on the NTS market model
#'
#' @param strPMNTS Structure of parameters for the n-dimensional NTS distribution.
#'
#' \code{strPMNTS$ndim} : dimension
#'
#' \code{strPMNTS$mu} : \eqn{\mu} mean vector (column vector) of the input data.
#'
#' \code{strPMNTS$sigma} : \eqn{\sigma} standard deviation vector (column vector) of the input data.
#'
#' \code{strPMNTS$alpha} : \eqn{\alpha} of the std NTS distribution (X).
#'
#' \code{strPMNTS$theta} : \eqn{\theta} of the std NTS distribution (X).
#'
#' \code{strPMNTS$beta} : \eqn{\beta} vector (column vector) of the std NTS distribution (X).
#'
#' \code{res$Rho} : \eqn{\rho} matrix (Correlation) of the std NTS distribution (X).
#'
#' \code{res$Sigma} : Covariance \eqn{\Sigma} matrix of return data \eqn{r}.
#'
#' @param w Capital allocation weight vector.
#' @param eta significanlt level
#' @return portfolio value at risk on the NTS market model
#'
portfolioCVaRmnts <- function(strPMNTS, w, eta){
  if(strPMNTS$ndim != length(w)){
    print("The dimension of the weight vector must be the same as strPMNTS$ndim")
    return(NULL)
  }
  ntsparam <- getPortNTSParam(strPMNTS, w)
  CVaR <- (- ntsparam$mu + ntsparam$sig * cvarnts(eta, ntsparam$stdparam))
  return(CVaR)
}


#' @export
#' @title mctStdDev
#' @description
#' Morginal contribution to Risk for Standard Deviation.
#'
#' @param n The targer stock to calculate the mctStdDev
#' @param w The capital allocation rate vector for the current portfolio
#' @param CovMtx Covariance matrix of return data.
#'
mctStdDev <- function(n, w, covMtx){
  sig <-  sqrt(w%*%covMtx%*%t(w))
  return( w%*%covMtx[,n]/sig )
}

#' @export
#' @title dBeta
#' @description
#' The first derivative of the beta.
#' Developer's version.
#'
dBeta <- function(n, w, betaArray, covMtx){
  barsig <-  sqrt(w%*%covMtx%*%t(w))
  barBeta <- sum(w*betaArray)
  mctsd <- mctStdDev(n, w, covMtx)
  sig <- sqrt(covMtx[n,n])
  #db <- w[n]*sig*betaArray[n]/(barsig^2)-(barBeta/barsig)*mctsd
  db <- sig*betaArray[n]/barsig-(barBeta/barsig)*mctsd
  return( db )
}

psi_stdNTS <- function( z, a, th, b ){
  sz <- (2-a)/(2*th)
  return(-1i*z+(1-1i*z*b/th+(1-b^2*sz)*z^2/(2*th))^(a/2-1)*(1i*z+b*sz*z^2))
}

int_phi_psi <- function(u, x, alpha, theta, beta, rho = 0.1){
  param <- c(alpha, theta, beta)
  res <- exp(-1i*u*x)*chf_stdNTS(u+rho*1i, param)/(rho-u*1i)
  res <- res*psi_stdNTS(u+1i*rho, alpha, theta, beta)
  return(res)
}

dFdBeta_eta <- function(x, eta, alpha, theta, beta, rho = 0.1, N = 20){
  fn <- pracma::integral(functional::Curry(int_phi_psi,
                                           x = x,
                                           alpha = alpha,
                                           theta = theta,
                                           beta = beta,
                                           rho = rho), 0, N)
  dFdB <- Re(fn)*exp(rho*x)/pi
  return(dFdB)
}

#' @export
#' @title dinvCdf_stdNTS_int
#' @description
#' The first derivative of inverse CDF for the beta parameter of the stdNTS.
#' Developer's version.
#'
dinvCdf_stdNTS_int <- function(eta, x = NULL, alpha, theta, beta){
  if( is.null(x) ) x <- ipnts(eta, c(alpha, theta, beta))

  denom <- dnts(x,  c(alpha, theta, beta))
  nom <- dFdBeta_eta(x, eta, alpha, theta, beta)

  return(-nom/denom)
}





#' @export
#' @title dCVaR_numint
#' @description
#' The first derivative of CVaR for the beta parameter of the stdNTS.
#' Developer's version.
#'
dCVaR_numint <- function( eta, alpha, theta, beta, N = 200, rho = 0.1 ){
  dFinv <- dinvCdf_stdNTS(eta, alpha, theta, beta)
  f <- Re(pracma::integral(functional::Curry(integrand_dCVaRstdNTS,
                                             eta = eta,
                                             alpha = alpha,
                                             theta = theta,
                                             beta = beta,
                                             rho = rho),
                           0, N))
  dcvar <- -dFinv-f/(pi*eta)
  return( dcvar )
}

#' @export
#' @title mctVaR_MNTS
#' @description Calculate the marginal contribution to VaR for the multivariate NTS market model:
#' the random vector \eqn{r} is
#'
#' \eqn{r = \mu + diag(\sigma) X}
#'
#' where
#'
#' \eqn{X} follows \eqn{stdNTS_N(\alpha, \theta, \beta, \Sigma)}
#'
#' @usage
#' \code{mctVaR_MNTS(n, eta, w, st)}
#'
#' @param n The target stock to calculate the mctVaR
#' @param eta Significant level of VaR.
#' @param w The capital allocation rate vector for the current portfolio
#' @param st Structure of parameters for the N-dimensional NTS distribution.
#' @param iCDFstd The inverst cdf of stdNTS at the significant level eta. If NULL, the function automatically find it. The default value is NULL
#'
#' \code{st$ndim} : Dimension of the model. Here \code{st$ndim=N}.
#'
#' \code{st$mu} : \eqn{\mu} mean vector (column vector) of the input data.
#'
#' \code{st$sigma} : \eqn{\sigma} standard deviation vector (column vector) of the input data.
#'
#' \code{st$alpha} : \eqn{\alpha} of the std NTS distribution (X).
#'
#' \code{st$theta} : \eqn{\theta} of the std NTS distribution (X).
#'
#' \code{st$beta} : \eqn{\beta} vector (column vector) of the std NTS distribution (X).
#'
#' \code{st$Rho} : \eqn{\rho} matrix of the std NTS distribution (X),
#'                     which is correlation matrix of epsilon.
#'
#' \code{st$CovMtx} : Covariance matrix of return data \eqn{r}.
#'
#' @references
#' Kim, Y. S. (2020) Portfolio Optimization on the Dispersion Risk and the Asymmetric Tail Risk
#' \url{https://arxiv.org/pdf/2007.13972.pdf}
#'
#' @examples
#' library(functional)
#' library(nloptr)
#' library(pracma)
#' library(spatstat)
#' library(Matrix)
#' library(quantmod)
#' library(mvtnorm)
#' library(temStaR)
#'
#' #Fix alpha and theta.
#' #Estimate alpha dna theta from DJIA and use those parameter for IBM, INTL parameter fit.
#' getSymbols("^DJI", src="yahoo", from = "2020-8-25", to = "2020-08-31")
#' prDJ <- as.numeric(DJI$DJI.Adjusted)
#' ret <- diff(log(prDJ))
#' ntsparam <-  fitnts(ret)
#' getSymbols("IBM", src="yahoo", from = "2016-1-1", to = "2020-08-31")
#' pr1 <- as.numeric(IBM$IBM.Adjusted)
#' getSymbols("INTL", src="yahoo", from = "2016-1-1", to = "2020-08-31")
#' pr2 <- as.numeric(INTL$INTL.Adjusted)
#'
#' returndata <- matrix(data = c(diff(log(pr1)),diff(log(pr2))),
#'                      ncol = 2, nrow = (length(pr1)-1))
#' st <- fitmnts( returndata = returndata,
#'                 n = 2,
#'                 alphaNtheta = c(ntsparam["alpha"], ntsparam["theta"])  )
#' w <- c(0.3, 0.7)
#' eta <- 0.01
#'
#' mctVaR_MNTS(1, eta, w, st) #MCT-VaR for IBM
#' mctVaR_MNTS(2, eta, w, st) #MCT-VaR for INTL
#'
mctVaR_MNTS <- function(n, eta, w, stmnts, iCDFstd = NULL ){
  barsig <-  sqrt(w%*%stmnts$CovMtx%*%t(w))
  mcts <- mctStdDev(n, w, stmnts$CovMtx)
  db <- dBeta(n, w, stmnts$beta, stmnts$CovMtx)
  barBeta <- sum(w*stmnts$beta)
  if (is.null(iCDFstd))
    iCDFstd <- ipnts(eta, ntsparam = c(stmnts$alpha, stmnts$theta, barBeta))
  dicdf <- dinvCdf_stdNTS_int(eta, iCDFstd, stmnts$alpha, stmnts$theta, barBeta)
  return(-stmnts$mu[n] - iCDFstd*mcts - barsig*db*dicdf )
}


integrand_dCVaRstdNTS <- function( u, eta, alpha, theta, beta, cv = NULL, v = NULL, rho = 0.0001 ){
  param <- c(alpha, theta, beta)
  Finv <- ipnts(eta, param)
  res <- (exp((1i*u+rho)*Finv)*chf_stdNTS(-u+rho*1i, param))
  res <- res*psi_stdNTS(-u+1i*rho, alpha, theta, beta)
  res <- res*(1/(rho*1i-u)^2+(cv+v)/(u*1i+rho))
  return(res)
}

#' @export
#' @title dCVaRstdNTS_numint
#' @description
#' Calculate the marginal contribution to CVaR for the multivariate stdNTS Model.
#' Developer's version.
#'
dCVaRstdNTS_numint <- function( eta, alpha, theta, beta, cv = NULL, v = NULL, N = 20, rho = 0.0001 ){

  if (is.null(cv)) cv <-  cvarnts(eps = eta, ntsparam = c(alpha, theta, beta))
  if (is.null(v)) v <- ipnts(u = eta, ntsparam = c(alpha, theta, beta))
  fn <- Re(pracma::integral(functional::Curry(integrand_dCVaRstdNTS,
                                              eta = eta,
                                              alpha = alpha,
                                              theta = theta,
                                              beta = beta,
                                              cv = cv,
                                              v = v,
                                              rho = rho),
                            0, N))
  dcvar <- -fn/(pi*eta)
  return( dcvar )
}

#' @export
#' @title mctCVaR_MNTS
#' @description Calculate the marginal contribution to CVaR for the multivariate NTS market model:
#' the random vector \eqn{r} is
#'
#' \eqn{r = \mu + diag(\sigma) X}
#'
#' where
#'
#' \eqn{X} follows \eqn{stdNTS_N(\alpha, \theta, \beta, \Sigma)}
#'
#' @usage
#' \code{mctCVaR_MNTS(eta, n, w, stmnts)}
#'
#' @param eta Significant level of CVaR.
#' @param n The targer stock to calculate the mctCVaR
#' @param w The capital allocation rate vector for the current portfolio
#' @param stmnts Structure of parameters for the N-dimensional NTS distribution.
#' @param CVaRstd CVaR Value of StdNTS residual. If NULL, the function automatically find it. The default value is NULL
#' @param dCVaRstd The first derivative of the stdNTS CVaR for beta. If NULL, the function automatically find it. The default value is NULL
#' @param iCDFstd The inverst cdf of stdNTS at the significant level eta. If NULL, the function automatically find it. The default value is NULL
#'
#' \code{st$ndim} : Dimension of the model. Here \code{st$ndim=N}.
#'
#' \code{st$mu} : \eqn{\mu} mean vector (column vector) of the input data.
#'
#' \code{st$sigma} : \eqn{\sigma} standard deviation vector (column vector) of the input data.
#'
#' \code{st$alpha} : \eqn{\alpha} of the std NTS distribution (X).
#'
#' \code{st$theta} : \eqn{\theta} of the std NTS distribution (X).
#'
#' \code{st$beta} : \eqn{\beta} vector (column vector) of the std NTS distribution (X).
#'
#' \code{st$Rho} : \eqn{\rho} matrix of the std NTS distribution (X),
#'                     which is correlation matrix of epsilon.
#'
#' \code{st$CovMtx} : Covariance matrix of return data \eqn{r}.
#'
#' @references
#' Kim, Y. S. (2020) Portfolio Optimization on the Dispersion Risk and the Asymmetric Tail Risk
#' \url{https://arxiv.org/pdf/2007.13972.pdf}
#'
#' @examples
#' library(functional)
#' library(nloptr)
#' library(pracma)
#' library(spatstat)
#' library(Matrix)
#' library(quantmod)
#' library(mvtnorm)
#' library("temStaR")
#'
#' #Fix alpha and theta.
#' #Estimate alpha dna theta from DJIA and use those parameter for IBM, INTL parameter fit.
#' getSymbols("^DJI", src="yahoo", from = "2020-8-25", to = "2020-08-31")
#' prDJ <- as.numeric(DJI$DJI.Adjusted)
#' ret <- diff(log(prDJ))
#' ntsparam <-  fitnts(ret)
#' getSymbols("IBM", src="yahoo", from = "2016-1-1", to = "2020-08-31")
#' pr1 <- as.numeric(IBM$IBM.Adjusted)
#' getSymbols("INTL", src="yahoo", from = "2016-1-1", to = "2020-08-31")
#' pr2 <- as.numeric(INTL$INTL.Adjusted)
#'
#' returndata <- matrix(data = c(diff(log(pr1)),diff(log(pr2))),
#'                      ncol = 2, nrow = (length(pr1)-1))
#' st <- fitmnts( returndata = returndata,
#'                 n = 2,
#'                 alphaNtheta = c(ntsparam["alpha"], ntsparam["theta"])  )
#' w <- c(0.3, 0.7)
#' eta <- 0.01
#'
#' mctCVaR_MNTS(1, eta, w, st) #MCT-CVaR for IBM
#' mctCVaR_MNTS(2, eta, w, st) #MCT-CVaR for INTL
#'
mctCVaR_MNTS <- function(n, eta, w, stmnts, CVaRstd=NULL, dCVaRstd=NULL, iCDFstd = NULL){
  barBeta <- sum(w*stmnts$beta)
  if (is.null(iCDFstd)) iCDFstd <- ipnts(eta, c(stmnts$alpha, stmnts$theta, barBeta))
  if (is.null(CVaRstd)) CVaRstd <- cvarnts(eps = eta, ntsparam = c(stmnts$alpha, stmnts$theta, barBeta))
  if (is.null(dCVaRstd))
    dCVaRstd <- dCVaRstdNTS_numint(eta, stmnts$alpha, stmnts$theta, barBeta, cv = CVaRstd, v = iCDFstd)

  barsig <-  sqrt(w%*%stmnts$CovMtx%*%t(w))
  mcts <- mctStdDev(n, w, stmnts$CovMtx)
  db <- dBeta(n, w, stmnts$beta, stmnts$CovMtx)
  return(-stmnts$mu[n] + CVaRstd*mcts + barsig*db*dCVaRstd )
}

#' @export
#' @title cvarGauss
#' @description
#' Calculate the CVaR for the normal distributed market model.
#' Developer's version.
#'
cvarGauss <- function( eta, mu = 0, sigma = 1){
  #var = -qnorm(eta, mu, sigma)
  stdK = qnorm(eta, 0, 1)
  cvar = sigma/(eta*sqrt(2*pi))*exp(-stdK^2/2)-mu
  return(cvar)
}

#' @title setPortfolioParam
#' @description Please use \code{getPortNTSParam} instead of \code{setPortfolioParam}.
#'
#' Portfolio return with capital allocation weight is \eqn{R_p=<w,r>},
#' which is a weighted sum of of elements in the N-dimensional NTS random vector.
#' \eqn{R_p} becomes an 1-dimensional NTS random variable.
#' \code{setPortfolioParam} find the parameters of \eqn{R_p}.
#'
#' @param strPMNTS Structure of parameters for the n-dimensional NTS distribution.
#'
#' \code{strPMNTS$ndim} : dimension
#'
#' \code{strPMNTS$mu} : \eqn{\mu} mean vector (column vector) of the input data.
#'
#' \code{strPMNTS$sigma} : \eqn{\sigma} standard deviation vector (column vector) of the input data.
#'
#' \code{strPMNTS$alpha} : \eqn{\alpha} of the std NTS distribution (X).
#'
#' \code{strPMNTS$theta} : \eqn{\theta} of the std NTS distribution (X).
#'
#' \code{strPMNTS$beta} : \eqn{\beta} vector (column vector) of the std NTS distribution (X).
#'
#' \code{strPMNTS$Rho} : \eqn{\Sigma} matrix of the std NTS distribution (X).
#'
#' @param w Capital allocation weight vector.
#'
#' @usage
#' \code{res <- setPortfolioParam(strPMNTS,w)}
#'
#' @return The weighted sum follows 1-dimensional NTS.
#'
#' \eqn{R_p = <w, r> = \mu + diag(\sigma) X},
#'
#' where
#'
#' \eqn{X} follows \eqn{stdNTS_1(\alpha, \theta, \beta, 1)}.
#'
#' Hence we obtain
#'
#' \code{res$mu} : \eqn{\mu} mean of \eqn{R_p}.
#'
#' \code{res$sigma} : \eqn{\sigma} standard deviation of \eqn{R_p}.
#'
#' \code{res$alpha} : \eqn{\alpha} of \eqn{X}.
#'
#' \code{res$theta} : \eqn{\theta} of \eqn{X}.
#'
#' \code{res$beta} : \eqn{\beta} \eqn{X}.
#'
#' @references
#' Kim, Y. S. (2020) Portfolio Optimization on the Dispersion Risk and the Asymmetric Tail Risk
#' \url{https://arxiv.org/pdf/2007.13972.pdf}
#'
#' @examples
#' library("temStaR")
#' strPMNTS <- list(ndim = 2,
#'                  mu = c( 9.876552e-05, 4.747343e-04 ),
#'                  sigma = c( 0.01620588, 0.02309643 ),
#'                  alpha = 0.1888129 ,
#'                  theta = 0.523042,
#'                  beta =  c( -0.04632938,  0.04063555 ),
#'                  Rho = matrix( data = c(1.0, 0.469883,
#'                                        0.469883, 1.0),
#'                                nrow = 2, ncol = 2)
#'                  CovMtx = matrix( data = c(0.0002626304, 0.0001740779,
#'                                          0.0001740779, 0.0005334452),
#'                                  nrow = 2, ncol = 2)
#' )
#' w <- c(0.3, 0.7)
#' res <- setPortfolioParam(strPMNTS,w)
#'
setPortfolioParam <- function(strPMNTS, w){
  if(strPMNTS$ndim != length(w)){
    print("The dimension of the weight vector must be the same as strPMNTS$ndim")
    return(NULL)
  }
  w <- matrix(data = w, nrow = strPMNTS$ndim, ncol = 1)
  m <- sum(w*strPMNTS$mu)
  b <- sum(w*strPMNTS$sigma*strPMNTS$beta)
  gammaMtx <- diag(as.numeric(sqrt(1-strPMNTS$beta^2*(2-strPMNTS$alpha)/(2*strPMNTS$theta))))
  sigMtx <- diag(as.numeric(strPMNTS$sigma))
  g <- as.numeric(sqrt(t(w)%*%sigMtx%*%gammaMtx%*%strPMNTS$Rho%*%gammaMtx%*%sigMtx%*%w))
  a <- as.numeric(strPMNTS$alpha)
  th <- as.numeric(strPMNTS$theta)
  param <- change_ntsparam2stdntsparam(c(a, th, b, g, m))
  return(param)
}



#' @export
#' @title portfolioVaRETmnts
#' @description
#' Calculate portfolio value at return on the NTS market model
#'
#' @param strPMNTS Structure of parameters for the n-dimensional NTS distribution.
#'
#' \code{strPMNTS$ndim} : dimension
#'
#' \code{strPMNTS$mu} : \eqn{\mu} mean vector (column vector) of the input data.
#'
#' \code{strPMNTS$sigma} : \eqn{\sigma} standard deviation vector (column vector) of the input data.
#'
#' \code{strPMNTS$alpha} : \eqn{\alpha} of the std NTS distribution (X).
#'
#' \code{strPMNTS$theta} : \eqn{\theta} of the std NTS distribution (X).
#'
#' \code{strPMNTS$beta} : \eqn{\beta} vector (column vector) of the std NTS distribution (X).
#'
#' \code{res$Rho} : \eqn{\rho} matrix (Correlation) of the std NTS distribution (X).
#'
#' \code{res$Sigma} : Covariance \eqn{\Sigma} matrix of return data \eqn{r}.
#'
#' @param w Capital allocation weight vector.
#' @param eta significanlt level
#' @return portfolio value at Return on the NTS market model
#'
portfolioVaRETmnts <- function(strPMNTS, w, eta){
  if(strPMNTS$ndim != length(w)){
    print("The dimension of the weight vector must be the same as strPMNTS$ndim")
    return(NULL)
  }
  ntsparam <- getPortNTSParam(strPMNTS, w)
  VaRet <- ( ntsparam$mu + ntsparam$sig * VaRetnts(eta, ntsparam$stdparam))
  return(VaRet)
}


#' @export
#' @title portfolioCVaRETmnts
#' @description
#' Calculate portfolio conditional value at Return on the NTS market model
#'
#' @param strPMNTS Structure of parameters for the n-dimensional NTS distribution.
#'
#' \code{strPMNTS$ndim} : dimension
#'
#' \code{strPMNTS$mu} : \eqn{\mu} mean vector (column vector) of the input data.
#'
#' \code{strPMNTS$sigma} : \eqn{\sigma} standard deviation vector (column vector) of the input data.
#'
#' \code{strPMNTS$alpha} : \eqn{\alpha} of the std NTS distribution (X).
#'
#' \code{strPMNTS$theta} : \eqn{\theta} of the std NTS distribution (X).
#'
#' \code{strPMNTS$beta} : \eqn{\beta} vector (column vector) of the std NTS distribution (X).
#'
#' \code{res$Rho} : \eqn{\rho} matrix (Correlation) of the std NTS distribution (X).
#'
#' \code{res$Sigma} : Covariance \eqn{\Sigma} matrix of return data \eqn{r}.
#'
#' @param w Capital allocation weight vector.
#' @param eta significanlt level
#' @return portfolio value at Return on the NTS market model
#'
portfolioCVaRETmnts <- function(strPMNTS, w, eta){
  if(strPMNTS$ndim != length(w)){
    print("The dimension of the weight vector must be the same as strPMNTS$ndim")
    return(NULL)
  }
  ntsparam <- getPortNTSParam(strPMNTS, w)
  CVaRet <- (ntsparam$mu + ntsparam$sig * cvaretnts(eta, ntsparam$stdparam))
  return(CVaRet)
}