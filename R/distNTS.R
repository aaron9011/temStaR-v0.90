
#' @export
#' @title dnts
#' @description \code{dnts} calculates pdf of the NTS distribution with parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu)}.
#' If only three parameters are given, it calculates pdf of the standard NTS distribution with parameter \eqn{(\alpha, \theta, \beta)}
#' If a time parameter value is given, it calculates pdf of the NTS profess
#' \eqn{f(x)dx=d(P((X(t+s)-X(s))<x))}, where X is the NTS process generated
#' by the NTS distribution with parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu)}.
#'
#' @param xdata An array of x
#' @param ntsparam A vector of the NTS parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu)}.
#' For the NTS process case it is a vector of parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu, t)}.
#' A vector of the standard NTS parameters \eqn{(\alpha, \theta, \beta)}.
#'
#' @return Density of NTS distribution
#'
#' @references
#' Kim, Y. S. (2020) Portfolio Optimization on the Dispersion Risk and the Asymmetric Tail Risk
#' \url{https://arxiv.org/pdf/2007.13972.pdf}
#'
#' @examples#
#' library(functional)
#' library(nloptr)
#' library(pracma)
#' library(spatstat)
#' library(Matrix)
#' library("temStaR")
#'
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' ntsparam <- c(alpha, theta, beta)
#' x <- seq(from = -6, to = 6, length.out = 101)
#' d <- dnts(x, ntsparam)
#' plot(x,d,type = 'l')
#'
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' gamma <- 0.3
#' mu <- 0.1
#' ntsparam <- c(alpha, theta, beta, gamma, mu)
#' x <- seq(from = -2, to = 2, by = 0.01)
#' d <- dnts(x, ntsparam)
#' plot(x,d,type = 'l')
#'
#'
#' #Annual based parameters
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' gamma <- 0.3
#' mu <- 0.1
#' #scaling annual parameters to one day
#' dt <- 1/250 #one day
#' ntsparam <- c(alpha, theta, beta, gamma, mu, dt)
#' x <- seq(from = -0.02, to = 0.02, length.out = 101)
#' d <- dnts(x, ntsparam)
#' plot(x,d,type = 'l')
dnts <- function( xdata, ntsparam ){
  newparam <- change_ntsparam2stdntsparam( ntsparam )
  pdf <- (1/newparam$sig)*pdf_FFT((xdata-newparam$mu)/newparam$sig,
                                 newparam$stdparam,  functional::Curry(chf_stdNTS))
  return( pdf )
}

#' @export
#' @title pnts
#' @description \code{pnts} calculates cdf of the NTS distribution with parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu)}.
#' If only three parameters are given, it calculates cdf of the standard NTS distribution with parameter \eqn{(\alpha, \theta, \beta)}
#' If a time parameter value is given, it calculates cdf of the  profess
#' \eqn{F(x)=P((X(t+s)-X(s))<x)}, where X is the NTS process generated
#' by the NTS distribution with parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu)}.
#'
#' @param xdata An array of x
#' @param ntsparam A vector of the NTS parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu)}.
#' For the NTS process case it is a vector of parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu, t)}.
#' A vector of the standard NTS parameters \eqn{(\alpha, \theta, \beta)}.
#'
#' @return Cumulative probability of the NTS distribution
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
#' library("temStaR")
#'
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' ntsparam <- c(alpha, theta, beta)
#' x <- seq(from = -6, to = 6, length.out = 101)
#' p <- pnts(x, ntsparam)
#' plot(x,p,type = 'l')
#'
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' gamma <- 0.3
#' mu <- 0.1
#' ntsparam <- c(alpha, theta, beta, gamma, mu)
#' x <- seq(from = -2, to = 2, by = 0.01)
#' p <- pnts(x, ntsparam)
#' plot(x,p,type = 'l')
#'
#'
#' #Annual based parameters
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' gamma <- 0.3
#' mu <- 0.1
#' #scaling annual parameters to one day
#' dt <- 1/250 #one day
#' ntsparam <- c(alpha, theta, beta, gamma, mu, dt)
#' x <- seq(from = -0.02, to = 0.02, length.out = 101)
#' p <- pnts(x, ntsparam)
#' plot(x,p,type = 'l')
pnts <- function( xdata, ntsparam,
                  dz = 2^-8, m = 2^12){
  newparam <- change_ntsparam2stdntsparam( ntsparam )
  cnts <- cdf_FFT_GilPelaez((xdata-newparam$mu)/newparam$sig,
                           newparam$stdparam,  functional::Curry(chf_stdNTS),
                           dz = 2^-8, m = 2^12)
  return( cnts )
}


#' @title qnts
#' @description \code{qnts} calculates quantile of the NTS distribution with parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu)}.
#' If only three parameters are given, it calculates quantile of the standard NTS distribution with parameter \eqn{(\alpha, \theta, \beta)}
#' If a time parameter value is given, it calculates quantile of NTS profess. That is it finds x such that
#' \eqn{u = P((X(t+s)-X(s))<x)}, where X is the NTS process generated
#' by the NTS distribution with parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu)}.
#'
#' @param u vector of probabilities.
#' @param ntsparam A vector of the NTS parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu)}.
#' For the NTS process case it is a vector of parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu, t)}.
#' A vector of standard NTS parameters \eqn{(\alpha, \theta, \beta)}.
#'
#' @return The quantile function of the NTS distribution
#' @export
#' @examples
#' library(functional)
#' library(nloptr)
#' library(pracma)
#' library(spatstat)
#' library(Matrix)
#' library("temStaR")
#'
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' ntsparam <- c(alpha, theta, beta)
#' u <- c(0.01,0.05,0.25,0.5, 0.75, 0.95, 0.99)
#' q <- qnts(u, ntsparam)
#'
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' gamma <- 0.3
#' mu <- 0.1
#' ntsparam <- c(alpha, theta, beta, gamma, mu)
#' u <- c(0.01,0.05,0.25,0.5, 0.75, 0.95, 0.99)
#' q <- qnts(u, ntsparam)
#'
#'
#' #Annual based parameters
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' gamma <- 0.3
#' mu <- 0.1
#' #scaling annual parameters to one day
#' dt <- 1/250 #one day
#' ntsparam <- c(alpha, theta, beta, gamma, mu, dt)
#' u <- c(0.01,0.05,0.25,0.5, 0.75, 0.95, 0.99)
#' q <- qnts(u, ntsparam)
qnts <- function(u, ntsparam){
  ipnts(u, ntsparam)
}

#' @export
#' @title rnts
#' @description \code{rnts} generates random numbers following NTS distribution with parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu)}.
#' If only three parameters are given, it generates random numbers of standard NTS distribution with parameter \eqn{(\alpha, \theta, \beta)}
#' If a time parameter value is given, it generates random numbers of increments of NTS profess for time interval t.
#'
#' @param n number of random numbers to be generated.
#' @param ntsparam A vector of NTS parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu)}.
#' For NTS process case it is a vector of parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu, t)}.
#' A vector of standard NTS parameters \eqn{(\alpha, \theta, \beta)}.
#'
#' @return NTS random numbers
#' @references
#' Kim, Y. S. (2020) Portfolio Optimization on the Dispersion Risk and the Asymmetric Tail Risk
#' \url{https://arxiv.org/pdf/2007.13972.pdf}
#' @examples
#' library(functional)
#' library(nloptr)
#' library(pracma)
#' library(spatstat)
#' library(Matrix)
#' library("temStaR")
#'
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' ntsparam <- c(alpha, theta, beta)
#' r <- rnts(100, ntsparam) #generate 100 NTS random numbers
#' plot(r)
#'
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' gamma <- 0.3
#' mu <- 0.1
#' ntsparam <- c(alpha, theta, beta, gamma, mu)
#' r <- rnts(100, ntsparam) #generate 100 NTS random numbers
#' plot(r)
#'
#'
#' #Annual based parameters
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' gamma <- 0.3
#' mu <- 0.1
#' #scaling annual parameters to one day
#' dt <- 1/250 #one day
#' ntsparam <- c(alpha, theta, beta, gamma, mu, dt)
#' r <- rnts(100, ntsparam) #generate 100 NTS random numbers
#' plot(r)
rnts <- function(n, ntsparam, u = NULL){
  if (is.null(u)){
    u <- pracma::rand(1,n)
  }
  r <- ipnts(u, ntsparam)
  return( c(r) )
}

#' @export
#' @title fitnts
#' @description \code{fitnts} fit parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu)} of the NTS distribution.
#' This function using the curvefit method between the empirical cdf and the NTS cdf.
#'
#' @usage
#' \code{fitnts(rawdat)}
#' \code{fitnts(rawdat), ksdensityflag = 1}
#' \code{fitnts(rawdat, initialparam = c(alpha, theta, beta, gamma, mu))}
#' \code{fitnts(rawdat, initialparam = c(alpha, theta, beta, gamma, mu)), ksdensityflag = 1}
#' \code{fitnts(rawdat, initialparam = c(alpha, theta, beta, gamma, mu)), maxeval = 100, ksdensityflag = 1}
#'
#'
#' @param rawdat Raw data to fit the parameters.
#' @param initialparam A vector of initial NTS parameters.
#' This function uses the \code{nloptr} package.
#' If it has a good initial parameter then estimation performs better.
#' If users do not know a good initial parameters, then just set it as \code{initialparam=NaN}, that is default.
#' The function \code{cffitnts()} may be helpful to find the initial parameters.
#' @param maxeval Maximum evaluation number for \code{nloptr}. The iteration stops on this many function evaluations.
#' @param ksdensityflag This function fit the parameters using the curvefit method between the empirical cdf and the NTS cdf.
#' If \code{ksdensityflag = 1} (default), then the empirical cdf is calculated by the kernel density estimation.
#' If \code{ksdensityflag = 0}, then the empirical cdf is calculated by the empirical cdf.
#'
#'
#' @return Estimated parameters
#' @references
#' Kim, Y. S. (2020) Portfolio Optimization on the Dispersion Risk and the Asymmetric Tail Risk
#' \url{https://arxiv.org/pdf/2007.13972.pdf}
#' @examples
#' library(evmix)
#' library(spatstat)
#' library("temStaR")
#' library("quantmod")
#' getSymbols("^GSPC", src="yahoo", from = "2013-1-1", to = "2023-12-31")
#' pr <- as.numeric(GSPC$GSPC.Adjusted)
#' ret <- diff(log(pr))
#' ntsparam <-  fitnts(ret)
#'
#' Femp = ecdf(ret)
#' x = seq(from=min(ret), to = max(ret), length.out = 100)
#' cemp = Femp(x)
#' ncdf = pnts(x, c(ntsparam))
#' plot(x,ncdf,type = 'l', col = "red")
#' points(x,cemp, type = 'l', col = "blue")
#' a = density(ret)
#' p = dnts(x,ntsparam)
#' plot(x,p,type = 'l', col = "red")
#' lines(a,type = 'l', col = "blue")
#'
fitnts <- function( rawdat, initialparam = NaN, maxeval = 100, ksdensityflag = 0){
  if (is.nan(sum(initialparam))){
    init = c(0.99, 2, 0)
  } else {
    sp = change_ntsparam2stdntsparam(initialparam)
    init = sp$stdparam
  }
  init[3] = init[3]/(sz(init[1],init[2]))

  mu = mean(rawdat)
  sig = sd(rawdat)
  obs = (rawdat-mu)/sig

  if(ksdensityflag == 1){
    ks <- density(obs)
    cdfks <- CDF(ks)
    x <- ks$x
    y <- cdfks(x)
    #y <- pkden(x,obs)
  } else{
    Femp = ecdf(obs)
    x = seq(from=min(obs), to = max(obs), length = 1000)
    y = Femp(x)
  }

  ntsp <- nloptr::bobyqa(init, functional::Curry(llhfnts, x = x, cemp = y, dispF = 0),
                         lower = c(0.0001, 0.0001, -0.9999),
                         upper = c(1.9999,  1000, 0.9999),
                         control = nloptr::nl.opts(list(maxeval = maxeval)))

  ntsp$par[3] <- ntsp$par[3]*sz(ntsp$par[1], ntsp$par[2])
  retparam <- change_stdntsparam2ntsparam(ntsp$par, mu, sig, 1)
  retparam  <- retparam[1:5]
  names(retparam) <- c("alpha", "theta", "beta", "gamma", "mu")
  return( retparam )
}

#' @export
#' @title fitstdnts
#' @description \code{fitstdnts} fit parameters \eqn{(\alpha, \theta, \beta)} of the standard NTS distribution.
#' This function using the curvefit method between the empirical cdf and the standard NTS cdf.
#'
#' @usage
#' \code{fitstdnts(rawdat)}
#' \code{fitstdnts(rawdat), ksdensityflag = 1}
#' \code{fitstdnts(rawdat, initialparam = c(alpha, theta, beta))}
#' \code{fitstdnts(rawdat, initialparam = c(alpha, theta, beta)), ksdensityflag = 1}
#' \code{fitstdnts(rawdat, initialparam = c(alpha, theta, beta)), maxeval = 100, ksdensityflag = 1}
#'
#'
#' @param rawdat Raw data to fit the parameters.
#' @param initialparam A vector of initial standard NTS parameters.
#' This function uses the \code{nloptr} package.
#' If it has a good initial parameter then estimation performs better.
#' If users do not know a good initial parameters, then just set it as \code{initialparam=NaN}, that is default.
#'
#' @param maxeval Maximum evaluation number for \code{nloptr}. The iteration stops on this many function evaluations.
#' @param ksdensityflag This function fit the parameters using the curvefit method between the empirical cdf and the standard NTS cdf.
#' If \code{ksdensityflag = 1} (default), then the empirical cdf is calculated by the kernel density estimation.
#' If \code{ksdensityflag = 0}, then the empirical cdf is calculated by the empirical cdf.
#'
#'
#' @return Estimated parameters
#' @references
#' Kim, Y. S. (2020) Portfolio Optimization on the Dispersion Risk and the Asymmetric Tail Risk
#' \url{https://arxiv.org/pdf/2007.13972.pdf}
#' @examples
#' library(functional)
#' library(nloptr)
#' library(pracma)
#' library(spatstat)
#' library(evmix)
#' library(Matrix)
#' library(quantmod)
#' library("temStaR")
#' library("quantmod")
#' getSymbols("^GSPC", src="yahoo", from = "2013-1-1", to = "2023-12-31")
#' pr <- as.numeric(GSPC$GSPC.Adjusted)
#' ret <- diff(log(pr))
#' stdret <- (ret-mean(ret))/sd(ret)
#' stdntsparam <-  fitstdnts(stdret)
#'
#' Femp = ecdf(stdret)
#' x = seq(from=min(stdret), to = max(stdret), length.out = 100)
#' cemp = Femp(x)
#' ncdf = pnts(x, c(stdntsparam))
#' plot(x,ncdf,type = 'l', col = "red")
#' lines(x,cemp, type = 'l', col = "blue")
#' a = density(stdret)
#' p = dnts(x,stdntsparam)
#' plot(x,p,type = 'l', col = "red", ylim = c(0, max(a$y, p)))
#' lines(a,type = 'l', col = "blue")
#'
fitstdnts <- function( rawdat, initialparam = NaN, maxeval = 100, ksdensityflag = 1){
  if (is.nan(sum(initialparam))){
    init = c(0.99, 2, 0)
  } else {
    init = initialparam
  }
  init[3] = init[3]/(sz(init[1],init[2]))

  if(ksdensityflag == 0){
    Femp = ecdf(rawdat)
    x = seq(from=min(rawdat), to = max(rawdat), by = 1000)
    y = Femp(x)
  } else{
    ks <- density(rawdat)
    cdfks <- CDF(ks)
    x <- ks$x
    y <- cdfks(x)
    #y <- pkden(x,rawdat)
  }

  ntsp <- nloptr::bobyqa(init, functional::Curry(llhfnts, x = x, cemp = y, dispF = 0),
                         lower = c(0.0001, 0.0001, -0.9999),
                         upper = c(1.9999,  1000, 0.9999),
                         control = nloptr::nl.opts(list(maxeval = maxeval)))

  ntsp$par[3]  <-  ntsp$par[3]*sz(ntsp$par[1], ntsp$par[2])
  stdntsparam <-  c(ntsp$par[1],ntsp$par[2],ntsp$par[3])
  names(stdntsparam) <- c("alpha", "theta", "beta")
  return( stdntsparam )
}

#' @export
#' @title chf_NTS
#' @description \code{chf_NTS} calculates Ch.F of the NTS distribution with parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu)}.
#' If a time parameter value is given, it calculates Ch.F of the NTS profess
#' \eqn{\phi(u)=E[\exp(i u (X(t+s)-X(s)) )]=\exp(t \log(E[\exp(i u X(1))]))}, where X is the NTS process generated
#' by the NTS distribution with parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu)}.
#'
#' @param u An array of u
#' @param ntsparam A vector of the NTS parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu)}.
#' For NTS process case it is a vector of parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu, t)}.
#'
#' @return Characteristic function of the NTS distribution
#' @examples
#' library(functional)
#' library(nloptr)
#' library(pracma)
#' library(spatstat)
#' library(Matrix)
#' library("temStaR")
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' gamma <- 0.3
#' mu <- 0.1
#' ntsparam <- c(alpha, theta, beta, gamma, mu)
#' u <- seq(from = -2*pi, to = 2*pi, length.out = 101)
#' phi <- chf_NTS(u, ntsparam)
#'
#'
#' #Annual based parameters
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' gamma <- 0.3
#' mu <- 0.1
#' #scaling annual parameters to one day
#' dt <- 1/250 #one day
#' ntsparam <- c(alpha, theta, beta, gamma, mu, dt)
#' u <- seq(from = -2*pi, to = 2*pi, length.out = 101)
#' phi <- chf_NTS(u, ntsparam)
#'
chf_NTS <- function(u, param){
  a = param[1]
  th = param[2]
  b = param[3]
  g = param[4]
  m = param[5]
  if ((length(param))>=6){
    dt = param[6]
  } else {
    dt = 1
  }

  y = exp( dt*(
            1i*(m-b)*u
            - 2*th^(1-a/2)/a
                *((th-1i*(b*u+1i*g^2*u^2/2))^(a/2)-th^(a/2) ) ) )
  return( y )
}

#' @export
#' @title chf_stdNTS
#' @description \code{chf_stdNTS} calculates Ch.F of the standard NTS distribution with parameters \eqn{(\alpha, \theta, \beta)}.
#' If a time parameter value is given, it calculates Ch.F of the standard NTS profess
#' \eqn{\phi(u)=E[\exp(i u (X(t+s)-X(s))]=\exp(t \log(E[\exp(i u X(1))]))}, where X is the standard NTS process generated
#' by the standard NTS distribution with parameters \eqn{(\alpha, \theta, \beta)}.
#'
#' @param u An array of u
#' @param ntsparam A vector of the standard NTS parameters \eqn{(\alpha, \theta, \beta)}.
#' For the standard NTS process case it is a vector of parameters \eqn{(\alpha, \theta, \beta, t)}.
#'
#' @return Characteristic function of the standatd NTS distribution
#' @examples
#' #' library(functional)
#' library(nloptr)
#' library(pracma)
#' library(spatstat)
#' library(Matrix)
#' library("temStaR")
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' ntsparam <- c(alpha, theta, beta)
#' u <- seq(from = -2*pi, to = 2*pi, length.out = 101)
#' phi <- chf_stdNTS(u, ntsparam)
#'
#'
#' #Annual based parameters
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' #scaling annual parameters to one day
#' dt <- 1/250 #one day
#' ntsparam <- c(alpha, theta, beta, gamma, mu, dt)
#' u <- seq(from = -2*pi, to = 2*pi, length.out = 101)
#' phi <- chf_stdNTS(u, ntsparam)
#'
chf_stdNTS = function(u, param){
  a = param[1]
  th = param[2]
  b = param[3]
  g = sqrt(abs(1-b^2*(2-a)/(2*th)));
  m  = 0;
  if (length(param)>=4){
    dt = param[4]
  } else {
    dt = 1
  }
  y = exp( dt*( 1i*(m-b)*u - 2*th^(1-a/2)/a*((th-1i*(b*u+1i*g^2*u^2/2))^(a/2)-th^(a/2) ) ) )
  return( y )
}


chf_NTSV2 <- function(u, param){
# param(1) = alpha
# param(2) = theta
# param(3) = b (not beta but b)
# param(4) = sigma
# param(5) = mu
# param(6) = dt (optional)
  a = param[1]
  th = param[2]
  bet  = param[3]*sqrt(2*th/(2-a));
  sig = param[4]
  m  = param[5]
  g = sqrt(1-param[2]^2)

  if (length(param)>=6){
    dt = param[6];
  } else {
    dt = 1
  }
  newparam = c(a, th, bet*sig, g*sig, m, dt)
  y = chf_NTS(u, newparam)
  return( y )
}

#' @export
#' @title moments_NTS
#' @description \code{moments_NTS} calculates mean, variance, skewness, and excess kurtosis
#' of the NTS distribution with parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu)}.
#'
#' @param param A vector of the NTS parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu)}.
#'
#' @return First 4 moments (Mean, Variance, Skewness, Excess Kurtosis) of NTS distribution.
#' The mean is always the same as the parameter \eqn{\mu}.
#'
#' @references
#' Kim, Y.S, K-H Roh, R. Douady (2020) Tempered Stable Processes with Time Varying Exponential Tails
#' \url{https://arxiv.org/pdf/2006.07669.pdf}
#'
#' @examples
#' library("temStaR")
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' gamma <- 0.3
#' mu <- 0.1
#' ntsparam <- c(alpha, theta, beta, gamma, mu)
#' moments_NTS(param = ntsparam)
moments_NTS <- function(param){
  # m4[1] mean
  # m4[2] variance
  # m4[3] skewness
  # m4[4] excess kurtosis

  m4 = c(0,0,0,0)
  cm = c(0,0,0,0)
  al = param[1];
  th = param[2];
  b = param[3];
  gm = param[4];
  m = param[5];

  #u = 0;
  #cm[1] = - b*1i + m*1i ...
  #    + th^(1 - al/2)*(- u*gm^2 + b*1i)*((gm^2*u^2)/2 - b*u*1i + th)^(al/2 - 1);
  #cm[2] = - gm^2*th^(1 - al/2)*((gm^2*u^2)/2 - b*u*1i + th)^(al/2 - 1)...
  #    - th^(1 - al/2)*(al/2 - 1)*(- u*gm^2 + b*1i)^2*((gm^2*u^2)/2 - b*u*1i + th)^(al/2 - 2);
  #cm[3] = th^(1 - al/2)*(al/2 - 1)*(al/2 - 2)*(- u*gm^2 + b*1i)^3*((gm^2*u^2)/2 - b*u*1i + th)^(al/2 - 3) ...
  #    + 3*gm^2*th^(1 - al/2)*(al/2 - 1)*(- u*gm^2 + b*1i)*((gm^2*u^2)/2 - b*u*1i + th)^(al/2 - 2);
  #cm[4] = - 3*gm^4*th^(1 - al/2)*(al/2 - 1)*((gm^2*u^2)/2 - b*u*1i + th)^(al/2 - 2)...
  #                    - th^(1 - al/2)*(al/2 - 1)*(al/2 - 2)*(al/2 - 3)*(- u*gm^2 + b*1i)^4 ...
  #                *((gm^2*u^2)/2 - b*u*1i + th)^(al/2 - 4) ...
  #        - 6*gm^2*th^(1 - al/2)*(al/2 - 1)*(al/2 - 2)*(- u*gm^2 + b*1i)^2 ...
  #                *((gm^2*u^2)/2 - b*u*1i + th)^(al/2 - 3);
  cm[1] = m*1i;
  cm[2] = b^2*th^(1 - al/2)*th^(al/2 - 2)*(al/2 - 1) - gm^2;
  cm[3] = - b^3*th^(1 - al/2)*th^(al/2 - 3)*(al/2 - 1)*(al/2 - 2)*1i + b*gm^2*th^(1 - al/2)*th^(al/2 - 2)*(al/2 - 1)*3i;

  cm[4] = 6*b^2*gm^2*th^(1 - al/2)*th^(al/2 - 3)*(al/2 - 1)*(al/2 - 2) - 3*gm^4*th^(1 - al/2)*th^(al/2 - 2)*(al/2 - 1) - b^4*th^(1 - al/2)*th^(al/2 - 4)*(al/2 - 1)*(al/2 - 2)*(al/2 - 3);

  cm = cm/c(1i, (1i)^2, (1i)^3, (1i)^4);

  m4[1] = cm[1];
  m4[2] = cm[2];
  m4[3] = cm[3]/((cm[2])^(3/2));
  m4[4] = cm[4]/((cm[2])^2);

  res <- Re(m4)
  names(res) <- c("mean", "variance", "skewness", "excess kurtosis")
  return(res)
}

#' @export
#' @title moments_stdNTS
#' @description \code{moments_stdNTS} calculates mean, variance, skewness, and excess kurtosis
#' of the standard NTS distribution with parameters \eqn{(\alpha, \theta, \beta)}.
#'
#' @param param A vector of the standard NTS parameters \eqn{(\alpha, \theta, \beta)}.
#'
#' @return First 4 moments (Mean, Variance, Skewness, Excess Kurtosis) of NTS distribution.
#' Of course, the mean and variance are always 0 and 1, respectively.
#'
#' @references
#' Kim, Y.S, K-H Roh, R. Douady (2020) Tempered Stable Processes with Time Varying Exponential Tails
#' \url{https://arxiv.org/pdf/2006.07669.pdf}
#'
#' @examples
#' library("temStaR")
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' ntsparam <- c(alpha, theta, beta)
#' moments_stdNTS(param = ntsparam)
moments_stdNTS <-function (param){
    a = param[1]
    th = param[2]
    bet = param[3] #beta
    gam = sqrt(abs(1-bet^2*(2-a)/(2*th)))
    m  = 0;
    m4 = moments_NTS(c(a,th,bet,gam,m))
    return( m4 )
}

# change stdntsparam  to ntsparam
#' @export
change_stdntsparam2ntsparam <-function(stdparam, mu, sig, dt = 1){
  a = stdparam[1]
  th = stdparam[2]
  bet = stdparam[3]*sig
  gam = sig*sqrt(1-stdparam[3]^2*(2-a)/(2*th))
  mu = mu
  ntsparam = c(a,th, bet, gam, mu, dt)
  names(ntsparam) <- c("alpha", "theta", "beta", "gamma", "mu", "dt")
  return( ntsparam )
}

# change ntsparam to stdntsparam
#' @export
change_ntsparam2stdntsparam <-function(ntsparam){
  if (length(ntsparam)==3){
    a = ntsparam[1] #alpha
    th_std = ntsparam[2]#theta
    bet_std = ntsparam[3] #beta
    sig_std = 1
    mu_std = 0
  } else {
    if(length(ntsparam)==4){
      a = ntsparam[1]
      th = ntsparam[2]*ntsparam[4]
      bet = ntsparam[3]*ntsparam[4]
      gam = sqrt(ntsparam[4]-bet^2*(2-a)/(2*th))
      m = 0
      dt = 1
    } else {
      a <- ntsparam[1] #alpha
      th <- ntsparam[2] #theta
      bet <- ntsparam[3] #beta
      gam <- ntsparam[4] #gamma
      m <- ntsparam[5] #mu
      if (length(ntsparam)>=6){
        dt = ntsparam[6]
      } else {
        dt = 1
      }
    }
    z  <- (2-a)/(2*th)
    th_std <- th*dt;
    bet_std <- bet*sqrt(dt/(gam^2+bet^2*z))
    sig_std <- sqrt(gam^2*dt+bet^2*dt*z)
    mu_std <- m*dt
  }
  stdntsparam <- c(a, th_std, bet_std)
  names(stdntsparam) <- c("alpha", "theta", "beta")
  param <- list(stdparam=stdntsparam, mu=mu_std, sig=sig_std)
  return( param )
}

#' @export
#' @title ipnts
#' @description \code{ipnts} calculates inverse cdf of the NTS distribution with parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu)}.
#' If only three parameters are given, it calculates inverse cdf of the standard NTS distribution with parameter \eqn{(\alpha, \theta, \beta)}
#'
#' @param u Real value between 0 and 1
#' @param ntsparam A vector of the NTS parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu)}.
#' A vector of the standard NTS parameters \eqn{(\alpha, \theta, \beta)}.
#'
#' @return Inverse cdf of the NTS distribution. It is the same as qnts function.
#' @references
#' Kim, Y. S. (2020) Portfolio Optimization on the Dispersion Risk and the Asymmetric Tail Risk
#' \url{https://arxiv.org/pdf/2007.13972.pdf}
#'
#' @examples
#' library("temStaR")
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' ntsparam <- c(alpha, theta, beta)
#' u <- seq(from = 0.01, to = 0.99, length.out = 99)
#' q <- ipnts(u, ntsparam)
#' plot(u,q,type = 'l')
#'
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' gamma <- 0.3
#' mu <- 0.1
#' ntsparam <- c(alpha, theta, beta, gamma, mu)
#' u <- seq(from = 0.01, to = 0.99, length.out = 99)
#' q <- ipnts(u, ntsparam)
#' plot(x,q,type = 'l')
#'
#'
#' #Annual based parameters
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' gamma <- 0.3
#' mu <- 0.1
#' #scaling annual parameters to one day
#' dt <- 1/250 #one day
#' ntsparam <- c(alpha, theta, beta, gamma, mu, dt)
#' u <- seq(from = 0.01, to = 0.99, length.out = 99)
#' q <- ipnts(u, ntsparam)
#' plot(x,q,type = 'l')
#'
ipnts <-function( u, ntsparam, maxmin = c(-10,10), du = 0.01){
  newparam = change_ntsparam2stdntsparam( ntsparam )
  x = invcdf_FFT_GilPelaez(u,newparam$stdparam, functional::Curry(chf_stdNTS), maxmin, du,
                           dz = 2^-8, m = 2^12)
  y = x*newparam$sig+newparam$mu
  names(y) <- u
  return( y )
}

#' @export
#' @title gensamplepathnts
#' @description \code{gensamplepathnts} generate sample paths
#' of the NTS process with parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu)}.
#' If only three parameters are given, it generate sample paths
#' of the standard NTS process with parameters \eqn{(\alpha, \theta, \beta)}.
#'
#' @param npath Number of sample paths
#' @param ntimestep number of time step
#' @param ntsparam A vector of the NTS parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu)}.
#' A vector of the standard NTS parameters \eqn{(\alpha, \theta, \beta)}.
#' @param dt the time length of one time step by the year fraction. "dt=1" means 1-year.
#'
#' @return Structure of the sample path.
#' Matrix of sample path. Column index is time.
#'
#' @examples
#' library("temStaR")
#' #standard NTS process sample path
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' ntsparam <- c(alpha, theta, beta)
#' npath <- 5
#' ntimestep <- 250
#' dt <- 1/250
#' simulation <- gensamplepathnts(npath, ntimestep, ntsparam, dt)
#' matplot(colnames(simulation), t(simulation), type = 'l')
#'
#' #NTS process sample path
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' gamma <- 0.3
#' mu <- 0.1
#' ntsparam <- c(alpha, theta, beta, gamma, mu)
#' npath <- 5
#' ntimestep <- 250
#' dt <- 1/250
#' simulation <- gensamplepathnts(npath, ntimestep, ntsparam, dt)
#' matplot(colnames(simulation), t(simulation), type = 'l')
#'
gensamplepathnts <- function(npath, ntimestep, ntsparam, dt){
  if (length(ntsparam)<=4){
    param <- ntsparam[1:3]
  } else if (length(ntsparam)>=5){
    param <- ntsparam[1:5]
  }
  u <- pracma::rand(npath,ntimestep)
  r <- ipnts(u, c(param,dt))
  x <- t(apply(r, 1, cumsum))
  colnames(x) <- seq(from = dt, to = dt*ntimestep, by = dt)
  return(x)
}


sz <- function(a,th){
  sqrt((2*th)/(2-a))
}

llhfnts <- function(ntsparam, x, cemp, dispF = 0){
  if (length(ntsparam) == 3){
    ntsparam[3] = ntsparam[3]*sz(ntsparam[1],ntsparam[2])
  }
  Fnts = pnts(x, ntsparam)
  MSE = mean((cemp-Fnts)^2)
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
#' @title cvarnts
#' @description \code{cvarnts} calculates Conditional Value at Risk (CVaR, or expected shortfall ES)
#' of the NTS market model with parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu)}.
#' If only three parameters are given, it calculates CVaR
#' of the standard NTS distribution with parameter \eqn{(\alpha, \theta, \beta)}
#'
#' @param eps the significant level for CVaR. Real value between 0 and 1.
#' @param ntsparam A vector of the NTS parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu)}.
#' A vector of the standard NTS parameters \eqn{(\alpha, \theta, \beta)}.
#'
#' @return CVaR of the NTS distribution.
#' @references
#' Y. S. Kim, S. T. Rachev, M. L. Bianchi, and F. J. Fabozzi (2010), Computing VaR and AVaR in infinitely divisible distributions, Probability and Mathematical Statistics, 30 (2), 223-245.
#'
#' S. T. Rachev, Y. S. Kim, M. L. Bianchi, and F. J. Fabozzi (2011), Financial Models with Levy Processes and Volatility Clustering, John Wiley & Sons
#'
#' @examples
#' library("temStaR")
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' ntsparam <- c(alpha, theta, beta)
#' u <- c(0.01,0.05)
#' q <- cvarnts(u, ntsparam)
#'
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' gamma <- 0.3
#' mu <- 0.1
#' ntsparam <- c(alpha, theta, beta, gamma, mu)
#' u <- c(0.01,0.05)
#' q <- cvarnts(u, ntsparam)
#'
#'
#' #Annual based parameters
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' gamma <- 0.3
#' mu <- 0.1
#' #scaling annual parameters to one day
#' dt <- 1/250 #one day
#' ntsparam <- c(alpha, theta, beta, gamma, mu, dt)
#' u <- c(0.01,0.05)
#' q <- cvarnts(u, ntsparam)
#'
cvarnts <- function( eps, ntsparam ){
  newparam = change_ntsparam2stdntsparam( ntsparam )
  varstd = -qnts(eps, newparam$stdparam)
  cvarstd = avar_numint(eps, varstd, newparam$stdparam, functional::Curry(chf_stdNTS))
  cvar = -newparam$mu + newparam$sig*cvarstd
  return(cvar)
}

#' @export
#' @title cvaretnts
#' @description \code{cvaretnts} calculates Conditional Value at Risk (CVaR, or expected shortfall ES)
#' of the NTS market model with parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu)}.
#' If only three parameters are given, it calculates CVaR
#' of the standard NTS distribution with parameter \eqn{(\alpha, \theta, \beta)}
#'
#' @param eps the significant level for CVaR. Real value between 0 and 1.
#' @param ntsparam A vector of the NTS parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu)}.
#' A vector of the standard NTS parameters \eqn{(\alpha, \theta, \beta)}.
#'
#' @return CVaReturn of the NTS distribution.
#' @references
#' Y. S. Kim, S. T. Rachev, M. L. Bianchi, and F. J. Fabozzi (2010), Computing VaR and AVaR in infinitely divisible distributions, Probability and Mathematical Statistics, 30 (2), 223-245.
#'
#' S. T. Rachev, Y. S. Kim, M. L. Bianchi, and F. J. Fabozzi (2011), Financial Models with Levy Processes and Volatility Clustering, John Wiley & Sons
#'
#' @examples
#' library("temStaR")
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' ntsparam <- c(alpha, theta, beta)
#' u <- c(0.01,0.05)
#' q <- cvaretnts(u, ntsparam)
#'
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' gamma <- 0.3
#' mu <- 0.1
#' ntsparam <- c(alpha, theta, beta, gamma, mu)
#' u <- c(0.01,0.05)
#' q <- cvaretnts(u, ntsparam)
#'
#'
#' #Annual based parameters
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' gamma <- 0.3
#' mu <- 0.1
#' #scaling annual parameters to one day
#' dt <- 1/250 #one day
#' ntsparam <- c(alpha, theta, beta, gamma, mu, dt)
#' u <- c(0.01,0.05)
#' q <- cvaretnts(u, ntsparam)
#'
cvaretnts <- function( eps, ntsparam ){
  newparam = change_ntsparam2stdntsparam( ntsparam )
  nparam <- c(newparam$stdparam[1], newparam$stdparam[2], -newparam$stdparam[3])
  varetstd = -qnts(eps, nparam)
  cvaretstd = avar_numint(eps, varetstd, nparam, functional::Curry(chf_stdNTS))
  cvaret = newparam$mu + newparam$sig*cvaretstd
  return(cvaret)
}

#' @export
#' @title VaRetnts
#' @description \code{VaRetnts} calculates Conditional Value at Risk (CVaR, or expected shortfall ES)
#' of the NTS market model with parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu)}.
#' If only three parameters are given, it calculates CVaR
#' of the standard NTS distribution with parameter \eqn{(\alpha, \theta, \beta)}
#'
#' @param eps the significant level for CVaR. Real value between 0 and 1.
#' @param ntsparam A vector of the NTS parameters \eqn{(\alpha, \theta, \beta, \gamma, \mu)}.
#' A vector of the standard NTS parameters \eqn{(\alpha, \theta, \beta)}.
#'
#' @return Value at Return of the NTS distribution.
#' @references
#' Y. S. Kim, S. T. Rachev, M. L. Bianchi, and F. J. Fabozzi (2010), Computing VaR and AVaR in infinitely divisible distributions, Probability and Mathematical Statistics, 30 (2), 223-245.
#'
#' S. T. Rachev, Y. S. Kim, M. L. Bianchi, and F. J. Fabozzi (2011), Financial Models with Levy Processes and Volatility Clustering, John Wiley & Sons
#'
#' @examples
#' library("temStaR")
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' ntsparam <- c(alpha, theta, beta)
#' u <- c(0.01,0.05)
#' q <- VaRetnts(u, ntsparam)
#'
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' gamma <- 0.3
#' mu <- 0.1
#' ntsparam <- c(alpha, theta, beta, gamma, mu)
#' u <- c(0.01,0.05)
#' q <- VaRetnts(u, ntsparam)
#'
#'
#' #Annual based parameters
#' alpha <- 1.2
#' theta <- 1
#' beta <- -0.2
#' gamma <- 0.3
#' mu <- 0.1
#' #scaling annual parameters to one day
#' dt <- 1/250 #one day
#' ntsparam <- c(alpha, theta, beta, gamma, mu, dt)
#' u <- c(0.01,0.05)
#' q <- VaRetnts(u, ntsparam)
#'
VaRetnts <- function( eps, ntsparam ){
  newparam = change_ntsparam2stdntsparam( ntsparam )
  nparam <- c(newparam$stdparam[1], newparam$stdparam[2], -newparam$stdparam[3])
  varetstd = -qnts(eps, nparam)
  varet = newparam$mu + newparam$sig*varetstd
  return(varet)
}



# NTS fit using Characteristic function
#' @export
cffitnts <- function(sample){

  # deal with the symmetric case
  sym <- c(sample, -sample)
  sigma0 <- sd(sym)
  sym <- sym/sigma0

  wk <- seq(from = 0.1, to = 1, by = 0.1)
  samplelogsymcf <- log(samplechf(wk, sym))

  initial <- c(1.4, 0.05)
  lower <- c(1.3, 0)

  out0 <- nloptr::bobyqa(x0 = initial, fn = function(x){return(errsymstcf(x, samplelogsymcf, wk))},
                         lower = lower)
  alpha <- out0$par[1]
  theta <- out0$par[2]

  # general case
  mu0 <- mean(sample)
  sigma0 <- sd(sample)
  sample <- (sample - mu0)/sigma0
  samplelogcf <- log(samplechf(wk, sample))
  initial <- c(theta, 1, 0)
  lower <- c(0, 0, -Inf)

  out <- nloptr::bobyqa(x0 = initial, fn = function(x){return(errstcf(x, samplelogcf, wk, alpha))},
                        lower = lower)
  theta <- out$par[1]
  beta <- out$par[3]*sigma0 + mu0
  gam <- sqrt(sigma0^2*out$par[2]^2 + (2 - alpha)*(sigma0^2*out$par[3]^2  - beta^2)/(2*theta))
  return(list(alpha = alpha, theta = theta, gamma = gam, beta = beta))
}

samplechf <- function(u, sample){
  f <- rep(0,length(u))
  for(k in 1:length(u)){
    f[k] <- mean(exp(1i*u[k]*sample))
  }
  return(f)
}

errsymstcf <- function(x, samplelogcf, wk){
  alpha <- x[1]
  theta <- x[2]

  theory <- -theta/(alpha/2)*(1 + wk^2/(2*theta))^(alpha/2) + theta/(alpha/2)
  f <- theory - Re(samplelogcf)
  f <- sum(f^2)
  return(f)
}

errstcf <- function(x, samplelogcf, wk, alpha){
  theta <- x[1]
  gamma <- x[2]
  beta <- x[3]

  theory <- -(2*theta^(1 - alpha/2)/alpha)*((theta - 1i*beta*wk + gamma^2*wk^2/2)^(alpha/2) - theta^(alpha/2))
  theory <- theory - 1i*beta*wk
  f <- c(Re(theory) - Re(samplelogcf), Im(theory) - Im(samplelogcf))
  f <- sum(f^2)
  return(f)
}

