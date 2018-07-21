# Local multifractal analysis of a Fisher log-series model

# log-series parameter alpha
a <- 15

# scales to be analyzed
sc <- floor(10^seq(2,4, len = 4))

# abundances on 4 scales
fish1 <- a*(sc[1]/(sc[1]+a))^(1:sc[1])/(1:sc[1])
fish2 <- a*(sc[2]/(sc[2]+a))^(1:sc[2])/(1:sc[2])
fish3 <- a*(sc[3]/(sc[3]+a))^(1:sc[3])/(1:sc[3])
fish4 <- a*(sc[4]/(sc[4]+a))^(1:sc[4])/(1:sc[4])

# Figure 2a
si <- rbind(fish1[1:10], fish2[1:10], fish3[1:10], fish4[1:10])
sisc <- rep(1:4, each = 10)
barplot(si, beside = T, col = c("grey", "red", "blue", "black"), names = 1:10,
        xlab = "n", ylab = expression(s[n]))
legend("topright", legend = c("N = 100", "N = 464", "N = 2154", "N = 10000"),
       fill = c("grey", "red", "blue", "black"), bty = "n")

# Figure 2b
plot(fish4, log = "xy", type = "l", lwd = 3, col = "black",
     xlab = "n", ylab = expression(s[n]), bty = "L")
points(fish3, type = "l", lwd = 3, col = "blue")
points(fish2, type = "l", lwd = 3, col = "red")
points(fish1, type = "l", lwd = 3, col = "grey")
legend("topright", legend = c("N = 100", "N = 464", "N = 2154", "N = 10000"),
       col = c("grey", "red", "blue", "black"), lwd = 3, bty = "n")

#---#

# orders
q <- seq(-10, 10, by = .1)

# init vectors for total abundance N and Shannon index H
N <- H <- floor(10^(seq(1.5, 4.5, len = 20)))

# init matrices for moments and effective divercities
m <- qD <- matrix(0, nrow = length(N), ncol = length(q))

# compute moments and Shannon index
for (ii in 1:length(N))
{
  nn <- 1:N[ii]
  x <- N[ii]/(N[ii]+a)
  s <- a*x^nn/nn
  
  for (jj in 1:length(q))
  {
    m[ii,jj] <- sum(s*(nn/N[ii])^q[jj])
  }
  
  pp <- nn/N[ii]
  H[ii] <- -sum(pp*log(pp)*s)
}

# convert moments into effective diversities
qD <- m ^ (1/(1-matrix(q, nrow = 20, ncol = length(q), byrow = T)))
qD[, q == 1] <- exp(H)

#---#

# init matrices for moments and diversities on 4 local scales
m4 <- qD4 <- matrix(0, nrow = 4, ncol = length(q))
H4 <- H[1:4]

# compute moments and diversities on local scales
for (ii in 1:4)
{
  nn <- 1:sc[ii]
  x <- sc[ii]/(sc[ii]+a)
  s <- a*x^nn/nn
  
  for (jj in 1:length(q))
  {
    m4[ii,jj] <- sum(s*(nn/sc[ii])^q[jj])
  }
  
  pp <- nn/sc[ii]
  H4[ii] <- -sum(pp*log(pp)*s)
}

qD4 <- m4 ^ (1/(1-matrix(q, nrow = 4, ncol = length(q), byrow = T)))
qD4[, q == 1] <- exp(H4)

#---#

# Figure 3a
library(rgl)
persp3d(x = log(N), y = q, z = log(qD), col = "skyblue")

# Figure 3b
ind <- c(F,F,F,F,T)

plot(q[ind], qD4[4,ind], type = "o", log = "y", pch = 21, bg = "black",
     xlab = "q", ylab = expression("log "^q*D))
lines(q[ind], qD4[3,ind], type = "o", pch = 21, bg = "blue")
lines(q[ind], qD4[2,ind], type = "o", pch = 21, bg = "red")
lines(q[ind], qD4[1,ind], type = "o", pch = 21, bg = "grey")

legend("topright", legend = c("N = 100", "N = 464", "N = 2154", "N = 10000"),
       pt.bg = c("grey", "red", "blue", "black"), lwd = 1, pch = 21)

# Figure 3c
plot(N, qD[, q == -2], type = "o", log = "xy", ylim = c(5, max(qD[,q == -2])), 
     pch = 21, bg = "wheat", xlab = "N", ylab = expression("log "^q*D))
lines(N, qD[, q == 0], type = "o", pch = 22, bg = "darkolivegreen1")
lines(N, qD[, q == 2], type = "o", pch = 23, bg = "darkseagreen1")
lines(N, qD[, q == 5], type = "o", pch = 24, bg = "lavender")

abline(v = sc, col = c("grey", "red", "blue", "black"))
legend("topleft", legend = c("q = -2", "q = 0", "q = 2", "q = 5"),
       pt.bg = c("wheat", "darkolivegreen1", "darkseagreen1", "lavender"), 
       lwd = 1, pch = 21:24)

spar <- 0.5
spl <- smooth.spline(x = log(N), y = log(qD[, q == 0]), spar = spar, tol = .0001)
lin <- predict(spl)

lin1 <- predict(spl, x = log(sc[1]))
der1 <- predict(spl, deriv = 1, x = log(sc[1]))
xx <- log(10^c(1.5,2.8))
lines(exp(xx), exp(der1$y*xx + lin1$y-der1$y*lin1$x), lwd = 2, col = "grey", lty = 1)

lin4 <- predict(spl, x = log(sc[4]))
der4 <- predict(spl, deriv = 1, x = log(sc[4]))
xx <- log(10^c(3,4.5))
lines(exp(xx), exp(der4$y*xx + lin4$y-der4$y*lin4$x), lwd = 2, col = "black", lty = 1)

points(N, qD[, q == 0], pch = 22, bg = "darkolivegreen1")

#---#

# compute local multifractal spectra
source("ecomf.r")
mls <- local.spectra(list(q = q, a = N, mom = m, H = H, qD = qD), sc = sc, smooth = spar)

# Figure 4
plot(mls$alfa[,1], mls$f[,1], pch = 21, bg = "grey", type = "o", 
     xlim = c(0, 1), ylim = range(0,mls$f), xlab = "a", ylab = "f(a)")
points(mls$alfa[,2], mls$f[,2], pch = 21, bg = "red", type = "o")
points(mls$alfa[,3], mls$f[,3], pch = 21, bg = "blue", type = "o")
points(mls$alfa[,4], mls$f[,4], pch = 21, bg = "black", type = "o")

abline(h = 0, v = 0:1, col = "grey")

legend("topleft", legend = c("N = 100", "N = 464", "N = 2154", "N = 10000"),
       pt.bg = c("grey", "red", "blue", "black"), lwd = 1, pch = 21)