#' @export
#' @title rmnts
#' @description \code{rmnts} generates random vector following the n dimensional NTS distribution using Cholesky decomposition.
#'
#' \eqn{r = \mu + diag(\sigma) X},
#'
#' where
#'
#' \eqn{X} follows \eqn{stdNTS_n(\alpha, \theta, \beta, \Sigma)}
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
#' \code{strPMNTS$Rho} : \eqn{\rho} matrix of the std NTS distribution (X).
#'
#' @param numofsample number of samples.
#'
#' @return Simulated NTS random vectors
#' @references
#' Kim, Y. S. (2020) Portfolio Optimization on the Dispersion Risk and the Asymmetric Tail Risk
#' \url{https://arxiv.org/pdf/2007.13972.pdf}
#'
#' @examples
#' strPMNTS <- list(ndim = 2,
#'                  mu = c( 0.00011, 0.00048 ),
#'                  sigma = c( 0.0162, 0.0231 ),
#'                  alpha = 1.23,
#'                  theta = 3.607,
#'                  beta =  c( -0.1209,  0.0905 ),
#'                  Rho = matrix( data = c(1.0, 0.55, 0.55, 1.0), nrow = 2, ncol = 2)
#' )
#' gensim <- rmnts( strPMNTS, 100 )
#' plot(gensim)
#'
rmnts <- function(strMnts, numofsample, u = NULL){
  N <- strMnts$ndim
  al <- strMnts$alpha
  th <- strMnts$theta
  beta <- strMnts$beta
  gamma <- sqrt(1-beta^2*(2-strMnts$alpha)/(2*strMnts$theta))
  L <- t(chol(strMnts$Rho))
  X <- matrix(nrow = numofsample, ncol = N)
  X[,1] = rnts(numofsample, c(al, th, beta[1]), u)
  for (n in 2:N){
    betatilde <- beta[n]-gamma[n]*sum(beta[1:(n-1)]*L[n,1:(n-1)]/gamma[1:(n-1)])
    gammatilde <- gamma[n]*L[n,n]
    Xtilde <- matrix(data=rnts(numofsample, c(al, th, betatilde, gammatilde, 0), u),
                     nrow = numofsample, ncol = 1)

    b <- matrix(data=gamma[n]*L[n,1:(n-1)]/gamma[1:(n-1)],
                nrow = n-1, ncol = 1)
    Xb <- X[,1:(n-1)]%*%b
    X[,n] <- Xtilde+Xb
  }
  Y <- matrix(nrow = numofsample, ncol = N)

  mu <- t(matrix(strMnts$mu, nrow = N, ncol = numofsample))
  sigma <- t(matrix(strMnts$sigma, nrow = N, ncol = numofsample))
  Y <- mu+sigma*X
  return(Y)
}

#' @export
#' @title rmnts_subord
#' @description \code{rmnts_subord} generates random vector following the n dimensional NTS distribution using subordination.
#'
#' \eqn{r = \mu + diag(\sigma) X},
#'
#' where
#'
#' \eqn{X} follows \eqn{stdNTS_n(\alpha, \theta, \beta, \Sigma)}
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
#' \code{strPMNTS$Rho} : \eqn{\rho} matrix of the std NTS distribution (X).
#'
#' @param numofsample number of samples.
#'
#' @return Simulated NTS random vectors
#' @references
#' Kim, Y. S. (2020) Portfolio Optimization on the Dispersion Risk and the Asymmetric Tail Risk
#' \url{https://arxiv.org/pdf/2007.13972.pdf}
#'
#' @examples
#' strPMNTS <- list(ndim = 2,
#'                  mu = c( 0.00011, 0.00048 ),
#'                  sigma = c( 0.0162, 0.0231 ),
#'                  alpha = 1.23,
#'                  theta = 3.607,
#'                  beta =  c( -0.1209,  0.0905 ),
#'                  Rho = matrix( data = c(1.0, 0.55, 0.55, 1.0), nrow = 2, ncol = 2)
#' )
#' gensim <- rmnts_subord( strPMNTS, 100 )
#' plot(gensim)
#'
rmnts_subord <- function( strPMNTS, numofsample, rW = NULL, rTau = NULL )
  # r = \mu + \sigma X
  # X = \beta-1+\sqrt(\tau) \gamma^T \epsilon
  # \mu = [\mu_1, \mu_2, \cdots, \mu_N]^T
  # \sigma = [\sigma_1,\sigma_2, \cdots, \sigma_N]^T
  # \beta = [\beta_1, \beta_2, \cdots, \beta_N]^T
  # \gamma = [gamma_1, gamma_2, \cdots, \gamma_N]^T
  # \gamma_n = \sqrt{1-\beta_n^2(2-\alpha)/(2\theta)}
  # \epsilon \sim N(0,\Rho)
  # \Rho Correlation metrix of N-dim standard normal distribution
{
    if( is.null(rW) ){
    rW <- pracma::randn(numofsample, strPMNTS$ndim)
  }

  if( is.null(rTau) ){
    rTau <- rsubTS(numofsample, c(strPMNTS$alpha, strPMNTS$theta) )
  }
  #print(rW)
  N <- strPMNTS$ndim
  beta <- matrix(data = strPMNTS$beta, nrow = 1, ncol = N)
  gamma <- sqrt(1-beta^2*(2-strPMNTS$alpha)/(2*strPMNTS$theta))
  gamma <- matrix(data = gamma, nrow = 1, ncol = N)
  reps <- t(t(chol(strPMNTS$Rho))%*%t(rW))
  rstdmnts <- (rTau-1)%*%beta+sqrt(rTau)*(reps%*%diag(as.numeric(gamma), N, N))
  #print(rstdmnts%*%diag(as.numeric(strPMNTS$sigma), N, N))
  res_rmnts <- rstdmnts%*%diag(as.numeric(strPMNTS$sigma), N, N) + t(matrix(data=as.numeric(strPMNTS$mu), nrow = N, ncol = numofsample))

  #print(res_rmnts)
  return( res_rmnts )
}

#' @export
#' @title fitmnts
#' @description \code{fitmnts} fit parameters
#' of the n-dimensional NTS distribution.
#'
#' \eqn{r = \mu + diag(\sigma) X}
#'
#' where
#'
#' \eqn{X} follows \eqn{stdNTS_n(\alpha, \theta, \beta, \Sigma)}
#'
#' @param returndata Raw data to fit the parameters.
#' The data must be given as a matrix form.
#' Each column of the matrix contains a sequence of asset returns.
#' The number of row of the matrix is the number of assets.
#'
#'@usage
#' \code{res <- fitmnts(returndata, n)}
#' \code{res <- fitmnts(returndata, n, alphaNtheta = c(alpha, theta))}
#' \code{res <- fitmnts(returndata, n, stdflag = TRUE ) }
#' \code{res <- fitmnts(returndata, n, alphaNtheta = c(alpha, theta), stdflag = TRUE)}
#'
#' @param n Dimension of the data. That is the number of assets.
#' @param alphaNtheta If \eqn{\alpha} and \eqn{\theta} are given,
#' then put those numbers in this parameter.
#' The function fixes those parameters and fits other remaining parameters.
#' If you set \code{alphaNtheta = NULL},
#' then the function fits all parameters including \eqn{\alpha} and \eqn{\theta}.
#'
#' @param stdflag If you want only standard NTS parameter fit, set this value be TRUE.
#'
#' @return Structure of parameters for the n-dimensional NTS distribution.
#'
#' \code{res$mu} : \eqn{\mu} mean vector of the input data.
#'
#' \code{res$sigma} : \eqn{\sigma} standard deviation vector of the input data.
#'
#' \code{res$alpha} : \eqn{\alpha} of the std NTS distribution (X).
#'
#' \code{res$theta} : \eqn{\theta} of the std NTS distribution (X).
#'
#' \code{res$beta} : \eqn{\beta} vector of the std NTS distribution (X).
#'
#' \code{res$Rho} : \eqn{\rho} matrix of the std NTS distribution (X),
#'                     which is correlation matrix of epsilon.
#'
#' \code{res$CovMtx} : Covariance matrix of return data \eqn{r}.
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
#' library(temStaR)
#'
#' getSymbols("^GSPC", src="yahoo", from = "2016-1-1", to = "2020-08-31")
#' pr1 <- as.numeric(GSPC$GSPC.Adjusted)
#' getSymbols("^DJI", src="yahoo", from = "2016-1-1", to = "2020-08-31")
#' pr2 <- as.numeric(DJI$DJI.Adjusted)
#'
#' returndata <- matrix(data = c(diff(log(pr1)),diff(log(pr2))),
#'                      ncol = 2, nrow = (length(pr1)-1))
#' res <- fitmnts( returndata = returndata, n=2 )
#'
#'
#' #Fix alpha and theta.
#' #Estimate alpha dna theta from DJIA and use those parameter for IBM, INTC parameter fit.
#' getSymbols("^DJI", src="yahoo", from = "2016-1-1", to = "2020-08-31")
#' prDJ <- as.numeric(DJI$DJI.Adjusted)
#' ret <- diff(log(prDJ))
#' ntsparam <-  fitnts(ret)
#' getSymbols("IBM", src="yahoo", from = "2016-1-1", to = "2020-08-31")
#' pr1 <- as.numeric(IBM$IBM.Adjusted)
#' getSymbols("INTC", src="yahoo", from = "2016-1-1", to = "2020-08-31")
#' pr2 <- as.numeric(INTC$INTC.Adjusted)
#'
#' returndata <- matrix(data = c(diff(log(pr1)),diff(log(pr2))),
#'                      ncol = 2, nrow = (length(pr1)-1))
#' res <- fitmnts( returndata = returndata,
#'                 n = 2,
#'                 alphaNtheta = c(ntsparam["alpha"], ntsparam["theta"])  )
#'
fitmnts <- function( returndata, n, alphaNtheta = NULL, stdflag = FALSE, PDflag = TRUE ){
  strPMNTS <- list(
    ndim = n,
    mu = matrix(data = 0, nrow = n, ncol = 1),
    sigma = matrix(data = 1, nrow = n, ncol = 1),
    alpha = 1,
    theta = 1,
    beta = matrix(data = 0, nrow = n, ncol = 1),
    Rho = matrix( nrow = n, ncol = n),
    CovMtx = matrix( nrow = n, ncol = n)
  )

  if( stdflag ){
    stdRetData <- returndata
    strPMNTS$CovMtx <- cor(returndata)
  }
  else{
    strPMNTS$CovMtx <- cov(returndata)
    strPMNTS$mu <- matrix(data = colMeans(returndata), nrow = n, ncol = 1)
    strPMNTS$sigma <- matrix(data = sqrt(diag(strPMNTS$CovMtx)), nrow = n, ncol = 1)
    muMtx <- t(matrix(data = strPMNTS$mu, nrow = n, ncol = length(returndata[,1])))
    sigmaMtx <- t(matrix(data = strPMNTS$sigma, nrow = n, ncol = length(returndata[,1])))
    stdRetData <- (returndata-muMtx)/sigmaMtx
  }

  athb <- matrix(nrow = n, ncol = 3)
  if (is.null(alphaNtheta)){
    for( k in 1:n ){
      stdntsparam <- fitstdnts(stdRetData[,k])
      athb[k,1] <- stdntsparam[1]
      athb[k,2] <- stdntsparam[2]
      athb[k,3] <- stdntsparam[3]
    }
    alphaNtheta = c(mean(athb[,1]), mean(athb[,2]))
  }
  strPMNTS$alpha <- alphaNtheta[1]
  strPMNTS$theta <- alphaNtheta[2]

  betaVec <- matrix(nrow = n, ncol = 1)
  for( k in 1:n ){
    betaVec[k,1] <- fitstdntsFixAlphaThata(stdRetData[,k], alpha = strPMNTS$alpha, theta = strPMNTS$theta)
  }
  strPMNTS$beta <- betaVec

  strPMNTS$Rho <- changeCovMtx2Rho(cov(stdRetData),
                                   strPMNTS$alpha,
                                   strPMNTS$theta,
                                   betaVec,
                                   PDflag)

  strPMNTS$mu <- as.numeric(strPMNTS$mu)
  strPMNTS$sigma <- as.numeric(strPMNTS$sigma)
  strPMNTS$beta <- as.numeric(strPMNTS$beta)
  return(strPMNTS)
}

#' @export
#' @title fitmnts
#' @description \code{fitmnts} fit parameters
#' of the n-dimensional NTS distribution. A parallel version of fitmnts()
#'
#' @examples
#' library(functional)
#' library(nloptr)
#' library(pracma)
#' library(spatstat)
#' library(Matrix)
#' library(foreach)
#' library(doParallel)
#' library(quantmod)
#' library(temStaR)
#'
#' getSymbols("^GSPC", src="yahoo", from = "2016-1-1", to = "2020-08-31")
#' pr1 <- as.numeric(GSPC$GSPC.Adjusted)
#' getSymbols("^DJI", src="yahoo", from = "2016-1-1", to = "2020-08-31")
#' pr2 <- as.numeric(DJI$DJI.Adjusted)
#'
#' returndata <- matrix(data = c(diff(log(pr1)),diff(log(pr2))),
#'                      ncol = 2, nrow = (length(pr1)-1))
#'
#' numofcluster  <- detectCores()
#' cl <- makePSOCKcluster(numofcluster)
#' registerDoParallel(cl)
#' res <- fitmnts_par( returndata = returndata, n=2, parallelSocketCluster = cl )
#' stopCluster(cl)
#'
fitmnts_par <- function( returndata, n, alphaNtheta = NULL, stdflag = FALSE, parallelSocketCluster = NULL, PDflag = TRUE ){
  flag_sock  <- FALSE
  if (is.null(parallelSocketCluster)){
    numofcluster  <- detectCores()
    parallelSocketCluster <- makePSOCKcluster(numofcluster)
    registerDoParallel(parallelSocketCluster)
    flag_sock <- TRUE
  }

  strPMNTS <- list(
    ndim = n,
    mu = matrix(data = 0, nrow = n, ncol = 1),
    sigma = matrix(data = 1, nrow = n, ncol = 1),
    alpha = 1,
    theta = 1,
    beta = matrix(data = 0, nrow = n, ncol = 1),
    Rho = matrix( nrow = n, ncol = n),
    CovMtx = matrix( nrow = n, ncol = n)
  )

  if( stdflag ){
    stdRetData <- returndata
    strPMNTS$CovMtx <- cor(returndata)
  }
  else{
    strPMNTS$CovMtx <- cov(returndata)
    strPMNTS$mu <- matrix(data = colMeans(returndata), nrow = n, ncol = 1)
    strPMNTS$sigma <- matrix(data = sqrt(diag(strPMNTS$CovMtx)), nrow = n, ncol = 1)
    muMtx <- t(matrix(data = strPMNTS$mu, nrow = n, ncol = length(returndata[,1])))
    sigmaMtx <- t(matrix(data = strPMNTS$sigma, nrow = n, ncol = length(returndata[,1])))
    stdRetData <- (returndata-muMtx)/sigmaMtx
  }

  athb <- matrix(nrow = n, ncol = 3)
  if (is.null(alphaNtheta)){

    res_par <- foreach(k = 1:n) %dopar% {
      library(temStaR)
      stdntsparam <- fitstdnts(stdRetData[,k])
      return(stdntsparam)
    }

    for( k in 1:n ){
      athb[k,1] <- res_par[[k]][1]
      athb[k,2] <- res_par[[k]][2]
      athb[k,3] <- res_par[[k]][3]
    }
    alphaNtheta = c(mean(athb[,1]), mean(athb[,2]))
  }
  strPMNTS$alpha <- alphaNtheta[1]
  strPMNTS$theta <- alphaNtheta[2]

  bet_par <- foreach(k = 1:n) %dopar% {
    library(temStaR)
    beta <- fitstdntsFixAlphaThata(stdRetData[,k], alpha = strPMNTS$alpha, theta = strPMNTS$theta)
    return(beta)
  }

  betaVec <- matrix(nrow = n, ncol = 1)
  for( k in 1:n ){
    betaVec[k,1] <- bet_par[[k]]
  }
  strPMNTS$beta <- betaVec

  strPMNTS$Rho <- changeCovMtx2Rho(cov(stdRetData),
                                   strPMNTS$alpha,
                                   strPMNTS$theta,
                                   betaVec,
                                   PDflag)

  strPMNTS$mu <- as.numeric(strPMNTS$mu)
  strPMNTS$sigma <- as.numeric(strPMNTS$sigma)
  strPMNTS$beta <- as.numeric(strPMNTS$beta)

  if(flag_sock){
    stopCluster(parallelSocketCluster)
  }
  return(strPMNTS)
}



#' @export
#' @title getGammaVec
#' @description beta to gamma in StdNTS
#'
getGammaVec <- function(alpha, theta, betaVec){
  n <- length(betaVec)
  gammaVec <- matrix(nrow = n, ncol = 1)
  for( k in 1:n ){
    gammaVec[k,1] <- sqrt(1-betaVec[k,1]^2*(2-alpha)/(2*theta))
  }
  return(gammaVec)
}

#' @export
#' @title changeCovMtx2Rho
#' @description Change covariance matrix to Rho matrix.
#'
changeCovMtx2Rho <- function(CovMtx, alpha, theta, betaVec, PDflag = TRUE){
  n <- length(betaVec)
  gammaVec <- getGammaVec(alpha, theta, betaVec)

  Rho <- (CovMtx-(2-alpha)/(2*theta)*(betaVec%*%t(betaVec)))
  igam <- diag(as.numeric(1/gammaVec))
  Rho <- igam%*%Rho%*%igam
  if( PDflag == TRUE ){
    pd <- nearPD(Rho, corr=TRUE)
    Rho <- matrix(data = as.numeric(pd$mat), ncol = n, nrow = n)
  }
  return(Rho)
}

#' @export
#' @title fitstdntsFixAlphaThata
#' @description Fit beta of stdNTS distribution with fixed alpha and theta.
#'
fitstdntsFixAlphaThata <- function( rawdat, alpha, theta, initialparam = NaN, maxeval = 100, ksdensityflag = 1){
  if (is.nan(sum(initialparam))){
    init <- 0
  } else {
    init = initialparam
  }

  if(ksdensityflag == 0){
    Femp = ecdf(rawdat)
    x = seq(from=min(rawdat), to = max(rawdat), by = 1000)
    y = Femp(x)
  } else{
    ks <- density(rawdat)
    cdfks <- spatstat.core::CDF(ks)
    x <- ks$x
    y <- cdfks(x)
  }

  ntsp <- nloptr::bobyqa(init,
                         functional::Curry(
                           llhfntsFixAlphaTheta,
                           alpha = alpha,
                           theta = theta,
                           x = x,
                           cemp = y,
                           dispF = 0),
                         lower = -0.9999,
                         upper = 0.9999,
                         control = nloptr::nl.opts(list(maxeval = maxeval)))

  beta  <-  ntsp$par*sz(alpha, theta)
  return( beta )
}

llhfntsFixAlphaTheta <- function(betaparam, alpha, theta, x, cemp, dispF = 0){
  betaparam <- betaparam*sz(alpha, theta)
  Fnts <- pnts(x, c(alpha, theta, betaparam) )
  MSE <- mean((cemp-Fnts)^2)
  #SSE = sqrt(sum((cemp-Fcts)^2))/length(x)

  if (dispF == 1){
    cat(ntsparam, sep = ',')
    cat(";")
    cat(MSE)
    cat("\n")
  }
  return( MSE )
}

#' @export
#' @title pmnts
#' @description \code{pmnts} calculates the cdf values of the multivariate NTS distribution:
#' \eqn{F(x_1, \cdots, x_n)=P(x_n<R_1, \cdots, x_n<R_n)}.
#' The multivariate NTS random vector \eqn{R = (R_1, \cdots, R_n)} is defined
#'
#' \eqn{R = \mu + diag(\sigma) X},
#'
#' where
#'
#' \eqn{X} follows \eqn{stdNTS_n(\alpha, \theta, \beta, \Sigma)}
#'
#' @param x array of the \eqn{(x_1, \cdots, x_n)}
#' @param st Structure of parameters for the n-dimensional NTS distribution.
#'
#' \code{st$ndim} : dimension
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
#' \code{st$Rho} : \eqn{\rho} matrix of the std NTS distribution (X).
#'
#' @param numofsample number of samples.
#'
#' @return Simulated NTS random vectors
#' @references
#' Kim, Y. S. (2020) Portfolio Optimization on the Dispersion Risk and the Asymmetric Tail Risk
#' \url{https://arxiv.org/pdf/2007.13972.pdf}
#'
#' @examples
#' library(mvtnorm)
#' library(temStaR)
#'
#' strPMNTS <- list(ndim = 2,
#'               mu = c( 0.5, -1.5 ),
#'               sigma = c( 2, 3 ),
#'               alpha = 0.1,
#'               theta = 3,
#'               beta =  c( 0.1, -0.3 ),
#'               Rho = matrix( data = c(1.0, 0.75, 0.75, 1.0),
#'                             nrow = 2, ncol = 2)
#' )
#' pmnts(c(0.6, -1.0), st = strPMNTS)
#'
#'
#' strPMNTS <- list(ndim = 2,
#'                  mu = c( 0, 0, 0 ),
#'                  sigma = c( 1, 1, 1 ),
#'                  alpha = 0.1,
#'                  theta = 3,
#'                  beta =  c( 0.1, -0.3, 0 ),
#'                  Rho = matrix(
#'                      data = c(1.0, 0.75, 0.1, 0.75, 1.0, 0.2, 0.1, 0.2, 1.0),
#'                      nrow = 3, ncol = 3)
#' )
#' pmnts(c(0,0,0), st = strPMNTS)
#' dmnts(c(0,0,0), st = strPMNTS)
#'
pmnts <- function( x, st, subTS = NULL ){
  #xadj <- (x-st$beta)/st$sigma
  xadj <- (x-st$mu)/st$sigma
  return(pMultiStdNTS(xadj, st, subTS))
}

#' @export
#' @title dmnts
#' @description \code{dmnts} calculates the density of the multivariate NTS distribution:
#' \eqn{f(x_1, \cdots, x_n)=\frac{d^n}{dx_1\cdots dx_n}P(x_n<R_1, \cdots, x_n<R_n)}.
#' The multivariate NTS random vector \eqn{R = (R_1, \cdots, R_n)} is defined
#'
#' \eqn{R = \mu + diag(\sigma) X},
#'
#' where
#'
#' \eqn{X} follows \eqn{stdNTS_n(\alpha, \theta, \beta, \Sigma)}
#'
#' @param x array of the \eqn{(x_1, \cdots, x_n)}
#' @param st Structure of parameters for the n-dimensional NTS distribution.
#'
#' \code{st$ndim} : dimension
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
#' \code{st$Rho} : \eqn{\rho} matrix of the std NTS distribution (X).
#'
#' @param numofsample number of samples.
#'
#' @return Simulated NTS random vectors
#' @references
#' Kim, Y. S. (2020) Portfolio Optimization on the Dispersion Risk and the Asymmetric Tail Risk
#' \url{https://arxiv.org/pdf/2007.13972.pdf}
#'
#' @examples
#' library("temStaR")
#' library(mvtnorm)
#' strPMNTS <- list(ndim = 2,
#'               mu = c( 0.5, -1.5 ),
#'               sigma = c( 2, 3 ),
#'               alpha = 0.1,
#'               theta = 3,
#'               beta =  c( 0.1, -0.3 ),
#'               Rho = matrix( data = c(1.0, 0.75, 0.75, 1.0),
#'                             nrow = 2, ncol = 2)
#' )
#' dmnts(c(0.6, -1.0), st = strPMNTS)
#'
#'
#' strPMNTS <- list(ndim = 2,
#'                  mu = c( 0, 0, 0 ),
#'                  sigma = c( 1, 1, 1 ),
#'                  alpha = 0.1,
#'                  theta = 3,
#'                  beta =  c( 0.1, -0.3, 0 ),
#'                  Rho = matrix(
#'                      data = c(1.0, 0.75, 0.1, 0.75, 1.0, 0.2, 0.1, 0.2, 1.0),
#'                      nrow = 3, ncol = 3)
#' )
#' pmnts(c(0,0,0), st = strPMNTS)
#' dmnts(c(0,0,0), st = strPMNTS)
#'
dmnts <- function( x, st, subTS = NULL ){
  #xadj <- (x-st$beta)/st$sigma
  xadj <- (x-st$mu)/st$sigma
  return(dMultiStdNTS(xadj, st, subTS)/prod(st$sigma))
}

#' @export
#' @title copulaStdNTS
#' @description \code{copulaStdNTS} calculates the stdNTS copula values
#' @references
#'	Y. S. Kim, D. Volkmann (2013), Normal Tempered Stable Copula, Applied Mathematics Letters,  26(7), 676-680
#'	\url{https://www.sciencedirect.com/science/article/pii/S0893965913000384}
#'
copulaStdNTS <- function(u, st, subTS = NULL){
  x <- matrix(nrow = length(u), ncol = 1)
  for (j in 1:length(u)){
    x[j] <- temStaR::qnts(u[j], c(st$alpha, st$theta, st$beta[j]))
  }
  return(pMultiStdNTS( as.numeric(x), st, subTS ))
}

#' @export
#' @title dcopulaStdNTS
#' @description \code{dcopulaStdNTS} calculates
#' density of the stdNTS copula.
#' @references
#'	Y. S. Kim, D. Volkmann (2013), Normal Tempered Stable Copula, Applied Mathematics Letters,  26(7), 676-680
#'	\url{https://www.sciencedirect.com/science/article/pii/S0893965913000384}
#'
dcopulaStdNTS <- function(u, st, subTS = NULL){
  x <- matrix(nrow = length(u), ncol = 1)
  y <- 1
  for (j in 1:length(u)){
    x[j] <- temStaR::qnts(u[j], c(st$alpha, st$theta, st$beta[j]))
    y <- y*temStaR::dnts(x[j], c(st$alpha, st$theta, st$beta[j]))
  }
  return(dMultiStdNTS( as.numeric(x), st, subTS )/y)
}

#' @export
#' @title importantSamplining
#' @description \code{importantSamplining} do the important sampling for the TS Subordinator.
#'
importantSamplining <- function( alpha, theta ){
  u  <- c(
    seq(from=0, to = 0.009, length.out = 10),
    seq(from=0.01, to = 0.09, length.out = 9),
    seq(from=0.1, to = 0.9, length.out = 9),
    0.95, 0.99, 0.999, 0.9999, 0.99999, 0.999999
  )
  ti <- temStaR::ipsubTS(u, c(alpha, theta), maxt = 20)
  subtsi <- temStaR::dsubTS(ti, c(alpha,theta))
  subTS <- list( ti = ti, subtsi = subtsi)
  return(subTS)
}

dMultiNorm_Subord <- function( tVec, x, beta, sig0 ){
  #gamma <- as.numeric(sqrt(1-(2-alpha)/(2*theta)*beta^2))
  #sig0 <- (cbind(gamma)%*%rbind(gamma))*rhoMtx
  re <- matrix(nrow = length(tVec), ncol = 1)
  for (i in 1:length(tVec) ){
    t <- tVec[i]
    mu <- beta*(t-1)
    Sig <- t*sig0

    re[i] <- dmvnorm(x, mean=mu, sigma = Sig )
  }
  return(re)
}

func_indegrand <- function(t, x, beta, sig0, ti, subtsi){
  fe <- dMultiNorm_Subord(t,
                          x = x,
                          beta = beta,
                          sig0 = sig0)
  ft <- pracma::pchip(ti, subtsi, t)
  return( fe*ft )
}

dMultiStdNTS <- function( x, st, subTS = NULL ){
  if ( is.null(subTS) ){
    subTS <- importantSamplining(st$alpha, st$theta)
  }
  ti <- subTS$ti
  subtsi <- subTS$subtsi

  gamma <- as.numeric(sqrt(1-(2-st$alpha)/(2*st$theta)*st$beta^2))
  sig0 <- (cbind(gamma)%*%rbind(gamma))*st$Rho

  d <- pracma::quad(
    functional::Curry(func_indegrand,
                      x = x,
                      beta = st$beta,
                      sig0 = sig0,
                      ti = ti,
                      subtsi = subtsi),
    xa = 0,
    xb = max(ti)
  )
  return(d)
}

pMultiNorm_Subord <- function( tVec, x, alpha, theta, beta, rhoMtx ){
  gamma <- as.numeric(sqrt(1-(2-alpha)/(2*theta)*beta^2))
  re <- matrix(nrow = length(tVec), ncol = 1)
  for (i in 1:length(tVec) ){
    t <- tVec[i]
    if (t == 0){
      re[i] <- 0
    }
    else{
      adjx <- (x-beta*(t-1))/(gamma*sqrt(t))
      re[i] <- pmvnorm(lower = rep(-Inf, length(x)), #c(-Inf, -Inf)
                       upper = adjx, mean = rep(0, length(x)), sigma = rhoMtx )
    }
  }
  return(re)
}

func_indegrand_cdf <- function(t, x, st, ti, subtsi){
  Gt <- pMultiNorm_Subord(t,
                          x = x,
                          alpha = st$alpha,
                          theta = st$theta,
                          beta = st$beta,
                          rhoMtx = st$Rho)
  ft <- pracma::pchip(ti, subtsi, t)
  #ft <- temStaR::dsubTS(t, c( st$alpha,  st$theta))
  return( Gt*ft )
}

pMultiStdNTS <- function( x, st, subTS = NULL ){
  if ( is.null(subTS) ){
    subTS <- importantSamplining(st$alpha, st$theta)
  }
  ti <- subTS$ti
  subtsi <- subTS$subtsi

  p <- pracma::quad(
    functional::Curry(func_indegrand_cdf,
                      x = x,
                      st = st,
                      ti = ti,
                      subtsi = subtsi),
    xa = 0,
    xb = max(ti)
  )
  p[p>1]=1
  return(p)
}

#' @export
#' @title pmarginalmnts
#' @description \code{pmarginalmnts} calculates
#' the marginal cdf of the \eqn{n}-th element
#' of the multivariate NTS distributed random variable.
#'
#' @param x the \eqn{x} such that \eqn{F(x) = P(X_n<x)}
#' @param n the \eqn{n}-th element to be calculated.
#' @param st Structure of parameters for the n-dimensional NTS distribution.
#'
#'
pmarginalmnts <- function(x, n, st){
  if (n>st$ndim){
    print("n must be less than or equal to st$ndim.")
    return(NULL)
  }
  if (n<1){
    print("n must be strictly positive integer.")
    return(NULL)
  }
  ntsparam  <- c(st$alpha, st$theta, st$beta[n])
  xdata <- (x-st$mu[n])/st$sigma[n]
  return(pnts(xdata, ntsparam))
}

#' @export
#' @title dmarginalmnts
#' @description \code{dmarginalmnts} calculates
#' the marginal density of the \eqn{n}-th element
#' of the multivariate NTS distributed random variable.
#'
#' @param x the \eqn{x} such that \eqn{f(x) = \frac{d}{dx}P(X_n<x)}
#' @param n the \eqn{n}-th element to be calculated.
#' @param st Structure of parameters for the n-dimensional NTS distribution.
#'
dmarginalmnts <- function(x, n, st){
  if (n>st$ndim){
    print("n must be less than or equal to st$ndim.")
    return(NULL)
  }
  if (n<1){
    print("n must be strictly positive integer.")
    return(NULL)
  }
  xdata <- (x-st$mu[n])/st$sigma[n]
  ntsparam  <- c(st$alpha, st$theta, st$beta[n])
  return(dnts(xdata, ntsparam)/st$sigma[n])
}

#' @export
#' @title qmarginalmnts
#' @description \code{qmarginalmnts} calculates
#' the quantile value of the \eqn{n}-th element
#' of the multivariate NTS distributed random variable.
#'
#' @param u vector of probabilities.
#' @param n the \eqn{n}-th element to be calculated.
#' @param st Structure of parameters for the n-dimensional NTS distribution.
#'
qmarginalmnts <- function(u, n, st){
  if (n>st$ndim){
    print("n must be less than or equal to st$ndim.")
    return(NULL)
  }
  if (n<1){
    print("n must be strictly positive integer.")
    return(NULL)
  }
  ntsparam  <- c(st$alpha, st$theta, st$beta[n])
  return(st$sigma[n]*qnts(u, ntsparam)+st$mu[n])
}

#' @export
#' @title cvarmarginalmnts
#' @description \code{cvarmarginalmnts} calculates
#' the CVaR of the \eqn{n}-th element
#' of the multivariate NTS distributed random variable.
#'
#' @param eta the significant level for CVaR. Real value between 0 and 1.
#' @param n the \eqn{n}-th element to be calculated.
#' @param st Structure of parameters for the n-dimensional NTS distribution.
#'
cvarmarginalmnts <- function(eta, n, st){
  if (n>st$ndim){
    print("n must be less than or equal to st$ndim.")
    return(NULL)
  }
  if (n<1){
    print("n must be strictly positive integer.")
    return(NULL)
  }
  ntsparam  <- c(st$alpha, st$theta, st$beta[n])
  return(st$sigma[n]*cvarnts(eta, ntsparam)-st$mu[n])
}
