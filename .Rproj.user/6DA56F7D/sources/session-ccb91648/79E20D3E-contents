
cleanupCdf <- function( cdfret, argout ){
  argout <- argout[cdfret>0]
  cdfret <- cdfret[cdfret>0]
  argout <- argout[cdfret<1]
  cdfret <- cdfret[cdfret<1]
  argout <- argout[diff(cdfret)!=0]
  cdfret <- cdfret[diff(cdfret)!=0]
  while( !pracma::isempty(cdfret[diff(cdfret)<=0]) ){
    argout <- argout(diff(cdfret)>0)
    cdfret <- cdfret(diff(cdfret)>0)
  }
  c$x <- argout
  c$y <- cdfret
  return( c )
}

pdf_FFT <- function(arg, param, chf, h = 2^-10,N = 2^17){
  s <- 1/(h*N)
  t1 <- 1:N
  t2 <- 2*pi*(t1 - 1 - N/2)*s

  cfvalues <- chf(u = t2, param = param)
  #chf_CTS( t2, param )
  x1 <- (-1)^(t1 - 1)*cfvalues

  pdfft <- Re(stats::fft(x1, inverse = FALSE))*(s*(-1)^(t1 - 1 - N/2))
  pdfft <- ifelse(pdfft>=0, pdfft, 0)

  x = (t1 - 1 - N/2)*h

  arg <- sort(arg)
  f <- pracma::pchip(x, pdfft, arg)
  return( f )
}

pdf_FFT2 <- function(argx, param, chf, dz = 2^-7, N = 2^17){
  h  <-  1/(dz*N)
  return( pdf_FFT(argx, param, chf, h, N ) )
  # z-grid for CF phi(z)
  #k <-  (0:(N-1))
  #K <-  dz*N

  #x <- (k-N/2)/K
  # z-grid for CDF
  #z <-  2*pi*k*dz

  # values of CF phi(z)
  #phi <-  (-1)^k*chf(z,param)

  #pdfft <- 2*Real(fft(phi)*dz)
  #pdfft <- ifelse(pdfft>=0, pdfft, 0)

  #f <- pracma::pchip(x, pdfft, argx)

}

cdf_FFT_GilPelaez <- function( arg, param, chf, dz = 2^-12, m = 2^17){

  # z-grid for CF phi(z)
  k = (0:(m-1))
  z = dz*(k+0.5)

  # z-grid for CDF
  x = pi*(2*k-m+1)/(m*dz)

  # values of CF phi(z)
  phi = chf(z,param)
  # point sequence for FFT
  seq = (phi/z)*exp(1i*pi*(m-1)/m*k)

  # CDF F(x) on x-grid
  F_fft = 0.5 - (1/pi)*Im(dz*exp(1i*pi*(m-1)/(2*m))*exp(-1i*pi/m*k)*stats::fft(seq, inverse = FALSE))

  # interpolation for requested arguments
  F = pracma::pchip(x, F_fft, arg)
  return( F )
}


cleanupCdf <- function( cdfret, argout ){
  argout <- argout[cdfret>0]
  cdfret <- cdfret[cdfret>0]
  argout <- argout[cdfret<1]
  cdfret <- cdfret[cdfret<1]
  argout <- argout[diff(cdfret)!=0]
  cdfret <- cdfret[diff(cdfret)!=0]
  while( !pracma::isempty(cdfret[diff(cdfret)<=0]) ){
    argout <- argout[diff(cdfret)>0]
    cdfret <- cdfret[diff(cdfret)>0]
  }
  x <- argout
  y <- cdfret
  c <- data.frame(x,y) # better to use an array for memory consumption reasons
  return( c )
}

invcdf_FFT_GilPelaez <- function(u, param, chf, maxmin, du, dz = 2^-12, m = 2^17){
  arg <- seq(from = maxmin[1], to = maxmin[2], by = du)
  c <- cdf_FFT_GilPelaez( arg, param, chf, dz, m )
  cc <- cleanupCdf(c, arg)

  y <- pracma::pchip(cc$y, cc$x, u)

  return( y )
}

integrandavar <- function( x, param, VaR, rho, chf){
  (exp(-1i*x*VaR)*chf(-x+rho*1i, param)/(rho*1i-x)^2)
}
avar_numint <- function( eps, VaR, param, chf, N = 10000, rho = 0.0001 ){
  f <- Re(pracma::integral(functional::Curry(integrandavar, param = param, VaR = VaR, rho = rho, chf = chf), 0, N))
  avar <- VaR-exp(-VaR*rho)*f/(pi*eps)
  return( avar )
}
