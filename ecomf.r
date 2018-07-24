#--- A set of tools to perform global and local multifracral analysis of community structure ---#

#' A function to compute moments
#' @param p abundance vector
#' @param q vector of orders
#' @keywords moment
#' @example mom(c(0.5, 0.3, 0.1, 0.1), 0:5)

mom <- function(p, q)
{
  p <- p[p > 0]
  p <- p/sum(p)
  m <- rep(0, length(q))
  for (ii in 1:length(q))
  {
    m[ii] <- sum(p^q[ii])
  }
  return(m)
}

#' A function to compute Shannon index
#' @param p abundance vector
#' @keywords Shannon
#' @example shannon(c(0.5, 0.3, 0.1, 0.1))

shannon <- function(p)
{
  p <- p[p > 0]
  p <- p/sum(p)
  return(-sum(p*log(p)))
}

#' A function to compute local derivative
#' @param vec vector
#' @param h interval size
#' @keywords derivative
#' @example derv(vec = seq(0,5, by = 0.1)^ 2, h = 0.1)

derv <- function(vec,h)
{
  n <- length(vec)
  f <- rep(0,n)
  f[2:(n-1)] <- (vec[3:n] - vec[1:(n-2)]) / (2*h)
  f[1] <- (-3*vec[1] + 4*vec[2] - vec[3]) / (2*h)
  f[n] <- (vec[n-2] - 4*vec[n-1] + 3*vec[n]) / (2*h)
  return(f)
}

#' A function to compute moments with a linear transect
#' @param ab matrix of abundances
#' @param q vector of orders

compute.moments.lin <- function(ab, q)
{
  ab <- as.matrix(ab)
  n <- dim(ab)[1]  # nr of samples
  nn <- sum(1:n)  # nr of cells
  
  m <- qD <- matrix(0, nrow = nn, ncol = length(q))
  pp <- matrix(0, nrow = nn, ncol = dim(ab)[2])
  a <- rep(0, nn)
  H <- rep(0, nn)
  pmin <- nmin <- rep(0, nn)
  
  counter <- 1
  for (ii in 1:(n-1))
  {
    for (jj in 1:(n-ii+1))
    {
      aa <- ab[jj:(jj+ii-1), ]
      
      if (is.matrix(aa))
      {
        aa <- colSums(aa)
      }
      if (sum(aa) > 0)
      {
        p <- aa/sum(aa)
        
        pp[counter,] <- p
        
        m[counter,] <- mom(p,q)
        a[counter] <- ii
        H[counter] <- shannon(p)
        pmin[counter] <- min(p[p > 0])
        nmin[counter] <- min(aa[aa > 0])
        counter <- counter + 1
      }
    }
  }
  
  aa <- colSums(ab)
  p <- aa/sum(aa)
  
  pp[counter, ] <- p
  
  m[counter,] <- mom(p,q)
  a[counter] <- n
  H[counter] <- shannon(p)
  pmin[counter] <- min(p[p > 0])
  nmin[counter] <- min(aa[aa > 0])
  
  qD <- m ^ (1/(1-matrix(q, nrow = dim(m)[1], ncol = dim(m)[2], byrow = T)))
  qD[, q == 1] <- exp(H)
  
  return(list(mom = m[1:counter, ], H = H[1:counter], qD = qD[1:counter, ], 
              a = a[1:counter], q = q, pmin = pmin[1:counter], nmin = nmin[1:counter],
              pp = pp[1:counter, ]))
}

#' A function to compute global multifractal spectrum
#' @param mom moments object returned by compute.moments.lin()

compute.spectrum <- function(mom)
{
  n <- length(mom$a)
  tau <- pp <- delta <- Dq <- rep(0,length(mom$q))
  xx <- log(mom$a)
  for (jj in 1:length(mom$q))
  {
    yy <- log(mom$qD[,jj])
    fin <- is.finite(yy)
    
    fit.lin <- lm(yy[fin] ~ xx[fin])
    fit.quad <- lm(yy[fin] ~ xx[fin] + I(xx[fin]^2))
    
    pp[jj] <- anova(fit.quad)[2,5]
    delta[jj] <- AIC(fit.quad) + (2*4*5/(n-4-1)) - AIC(fit.lin) - (2*3*4/(n-3-1))
    
    Dq[jj] <- fit.lin$coefficients[2]
  }
  
  tau <- Dq*(1-mom$q)
  
  alfa <- tau
  f <- tau
  alfa <- -derv(tau,mom$q[2]-mom$q[1])
  f <- mom$q*alfa + tau
  
  return(list(q = mom$q, tau = tau, Dq = Dq, alfa = alfa, f = f, pp = pp, delta = delta))
}

#' A function to compute local multifractal spectra
#' @param inp input object returned by compute.moments.lin()
#' @param sc vector of sample scales
#' @param smooth smoothing parameter

local.spectra <- function(inp, sc = 0, smooth = 1)
{
  if (exists("a", inp)) inp$A <- inp$a
  
  q <- inp$q
  
  if (sc[1] == 0)
  {
    fA <- log(max(inp$A))
    sA <- log(min(inp$A))
    dA <- fA - sA
    sc <- seq(sA + 0.05*dA, fA - 0.05*dA, length = 4)
  } else
  {
    sc <- log(sc)
  }
  
  tau <- rep(0,length(q)*length(sc))
  dim(tau) <- c(length(q),length(sc))
  alfa <- f <- Dq <- tau
  
  for (jj in 1:length(q))
  {
    xx <- log(inp$A)
    yy <- log(inp$qD[,jj])
    
    fin <- is.finite(yy)
    
    spl <- smooth.spline(x = xx[fin], y = yy[fin], spar = smooth, tol = .0001)
    
    der <- predict(spl, x = sc, deriv = 1)
    
    
    Dq[jj,] <- der$y
    tau[jj,] <- der$y*(1-q[jj])
  }
  
  for (ii in 1:length(sc))
  {
    alfa[,ii] <- -derv(tau[,ii],.1)
    f[,ii] <- q*alfa[,ii] + tau[,ii]
  }
  
  return(list(a = sc, q = q, tau = tau, Dq = Dq, alfa = alfa, f = f))
}

#' A function to truncate anomalous "horns" from multifractal spectra
#' @param mf spectra object returned by local.spectra()

trunc.spectra <- function(mf)
{
  n <- ncol(mf$f)
  len <- nrow(mf$f)
  start <- which(mf$q == 0)
  for (ii in 1:n)
  {
    pos <- start - 1
    while (mf$f[pos, ii] <= mf$f[pos+1, ii] & pos != 1) pos <- pos - 1
    if (pos != 1) 
    {
      mf$f[1:pos, ii] <- NA
      mf$alfa[1:pos, ii] <- NA
    }
    pos <- start + 1
    while (mf$f[pos, ii] <= mf$f[pos-1, ii] & pos != len) pos <- pos + 1
    if (pos != len) 
    {
      mf$f[pos:len, ii] <- NA
      mf$alfa[pos:len, ii] <- NA
    }
  }
  return(mf)
}