library(functional)
library(nloptr)
library(pracma)
library(spatstat)
library(Matrix)

#' @export
pcts <-function( xdata, ctsparam,
                 dz = 2^-8, m = 2^12){
  newparam <- change_ctsparam2stdctsparam( ctsparam )
  ccts <- cdf_FFT_GilPelaez((xdata-newparam$mu)/newparam$sig,
                           newparam$stdparam,  functional::Curry(chf_stdCTS),
                           dz = 2^-8, m = 2^12)
  return( ccts )
}

#' @export
qcts <- function(u, ctsparam){
  ipcts(u, ctsparam)
}

#' @export
rcts <- function(n, ctsparam){
  u <- rand(1,n)
  r <- ipcts(u, ctsparam)
  return( c(r) )
}

#' @export
dcts <-function( xdata, ctsparam ){
  newparam = change_ctsparam2stdctsparam( ctsparam )
  pdf = (1/newparam$sig)*pdf_FFT((xdata-newparam$mu)/newparam$sig,
                                 newparam$stdparam, functional::Curry(chf_stdCTS))
  return( pdf )
}

#' @export
fitcts <- function( rawdat, initialparam = NaN, maxeval = 100, ksdensityflag = 1){
  if (is.nan(sum(initialparam))){
    init = c(0.99, 1, 1)
  } else {
    sp = change_ctsparam2stdctsparam(initialparam)
    init = sp$stdparam
  }

  mu = mean(rawdat)
  sig = sd(rawdat)
  obs = (rawdat-mu)/sig

  if(ksdensityflag == 0){
    Femp = ecdf(obs)
    x = seq(from=min(obs), to = max(obs), length = 1000)
    y = Femp(x)
  } else{
    ks <- density(obs)
    #cdfks <- CDF(ks)
    x <- ks$x
    #y <- cdfks(x)
    y <- pkden(x,obs)
  }

  ctsp <- nloptr::bobyqa(init, functional::Curry(llhfcts, x = x, cemp = y),
                 lower = c(0.0001, 0.0001, 0.0001),
                 upper = c(1.999,  1000, 1000),
                 control = nl.opts(list(maxeval = maxeval)))

  retparam = change_stdctsparam2ctsparam(ctsp$par, mu, sig, 1)
  retparam = retparam[1:5]
  return( retparam )
}

#' @export
chf_CTS <- function(u, param){
  alpha = param[1]
  C     = param[2]
  lmp   = param[3]
  lmm   = param[4]
  m     = param[5]-gamma(1-alpha)*C*(lmp^(alpha-1)-lmm^(alpha-1))
  if ((length(param))>=6){
    dt = param[6]
  } else {
    dt = 1
    }

  phi = exp(dt*( 1i*u*m + C*gamma(-alpha)*(
    (lmp - 1i*u )^alpha - lmp^alpha
    +(lmm + 1i*u )^alpha - lmm^alpha ) ) )
  return(phi)
}

normalizeCTSparam <-function( param ){
  c = 1/gamma(2-param[1])/(param[2]^(param[1]-2)+param[3]^(param[1]-2))
  normalizedparam = c(param[1], c, param[2], param[3], 0)
  return(normalizedparam)
}


chf_stdCTS = function(u, param){

  if (length(param)>=4){
    dt = param[4]
  } else {
    dt = 1
  }
  normalizedparam = normalizeCTSparam( param[1:3] )
  normalizedparam = c(normalizedparam, dt)
  phi = chf_CTS(u, normalizedparam)
  return(phi)
}

change_stdctsparam2ctsparam <-function(stdparam, mu, sig, dt){
  a = stdparam[1]
  C = sig^a/(gamma(2-a)*(stdparam[2]^(a-2)+stdparam[3]^(a-2)))
  lp = stdparam[2]/sig
  lm = stdparam[3]/sig
  ctsparam = c(a,C/dt,lp,lm, mu/dt,dt)
  return( ctsparam )
}

change_ctsparam2stdctsparam <-function(ctsparam){
  if (length(ctsparam)==3){
    a = ctsparam[1]
    lp_std = ctsparam[2]
    lm_std = ctsparam[3]
    sig_std = 1
    mu_std = 0
  } else {
    if (length(ctsparam)==4){
      dt = ctsparam[4]
      C = 1/(gamma(2-ctsparam[1])*(ctsparam[2]^(ctsparam[1]-2)+ctsparam[3]^(ctsparam[1]-2)))
      ctsparam = c(ctsparam[1],C,ctsparam[2],ctsparam[3],0,dt)
    } else if (length(ctsparam)==6){
      dt = ctsparam[6]
    } else {
      dt = 1
    }
    a = ctsparam[1] #alpha
    c = ctsparam[2]*dt #C
    lp = ctsparam[3] #lambda plus
    lm = ctsparam[4] #lambda minus
    m = ctsparam[5]*dt #m

    sig_std = sqrt(gamma(2-a)*(lp^(a-2)+lm^(a-2))*c)
    lm_std = lm*sig_std
    lp_std = lp*sig_std
    mu_std = m
    }
  param=list(stdparam=c(a,lp_std,lm_std), mu=mu_std, sig=sig_std)
  return( param )
}

ipcts <-function( u, ctsparam, maxmin = c(-10,10), du = 0.01){
  newparam = change_ctsparam2stdctsparam( ctsparam )
  x = invcdf_FFT_GilPelaez(u,newparam$stdparam,
                           functional::Curry(chf_stdCTS),
                           maxmin, du,
                           dz = 2^-8, m = 2^12)
  y = x*newparam$sig+newparam$mu
  return( y )
}

gensamplepathcts <- function(npath, ntimestep, ctsparam, dt){
  if (length(ctsparam)<=4){
    param = ctsparam[1:3]
  } else if (length(ctsparam)>=5){
    param = ctsparam[1:5]
  }
  u = rand(npath,ntimestep)
  r = ipcts(u, c(param,dt))
  pth = list(x = t(apply(r, 1, cumsum)), t = seq(from = dt, to = dt*ntimestep, by = dt))
  return( pth )
}


llhfcts <- function(ctsparam, x, cemp, dispF = 0){
  Fcts = pcts(x, ctsparam)
  MSE = mean((cemp-Fcts)^2)
  #SSE = sqrt(sum((cemp-Fcts)^2))/length(x)

  if (dispF == 1){
    cat(ctsparam, sep = ',')
    cat(";")
    cat(MSE)
    cat("\n")
  }
  return( MSE )
}

fitstdcts <- function( obs, initialparam = NaN, maxeval = 100, ksdensityflag = 1){
  if (is.nan(sum(initialparam))){
    init = c(0.99, 1, 1)
  } else {
    init = initialparam
  }

  if(ksdensityflag == 0){
    Femp = ecdf(obs)
    x = seq(from=min(obs), to = max(obs), length = 1000)
    y = Femp(x)
  } else{
    ks <- density(obs)
    #cdfks <- CDF(ks)
    x <- ks$x
    #y <- cdfks(x)
    y <- pkden(x,obs)
  }

  ctsp <- nloptr::bobyqa(init[1:3], functional::Curry(llhfcts, x = x, cemp = y),
                 lower = c(0.0001, 0.0001, 0.0001),
                 upper = c(1.999,  1000, 1000),
                 control = nl.opts(list(maxeval = maxeval)))
  retparam = ctsp$par
  return( retparam )
}

cvarcts <- function( eps, ctsparam ){
  newparam = change_ctsparam2stdctsparam( ctsparam )
  varstd = -qcts(eps, newparam$stdparam)
  cvarstd = avar_numint(eps, varstd, newparam$stdparam, functional::Curry(chf_stdCTS))
  cvar = -newparam$mu + newparam$sig*cvarstd
  return(cvar)
}
