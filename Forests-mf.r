# Local multifractal analysis of forest layers

source("ecomf.r")

# loading plain forest data
load("D:/YandexDisk/Science/2016/China/R/ank-tree.rda")
load("D:/YandexDisk/Science/2016/China/R/ank-herbs.rda")
hc <- t(hc)
rm(sc, sh, tc, th)
ta <- ta[-51, ]
sa <- sa[-51, ]
hc <- hc[-51, ]

# compute moments with a linear scheme
m.at <- compute.moments.lin(ta, q = seq(-3, 3, .1))
m.as <- compute.moments.lin(sa, q = seq(-3, 3, .1))
m.ah <- compute.moments.lin(hc, q = seq(-3, 3, .1))

#--- diversity scaling plots in the plain forest ---#

# Figure 5a
ind <- m.at$q == 0

xx <- log(m.at$a*10)
yy <- log(m.at$qD[, ind])

plot(exp(xx), exp(yy), pch = 19, col = rgb(.7, .9, .7, .5), log = "xy",
     xlab = "L", ylab = expression(''^0*D))

spl <- smooth.spline(x = xx, y = yy, spar = 1.15, tol = .0001)
lines(exp(spl$x), exp(spl$y), lwd = 2, col = "red")

abline(v = 10*exp(1:4), col = c("grey", "red", "blue", "black"))

# Figure 5b
ind <- m.at$q == 3

xx <- log(m.at$a*10)
yy <- log(m.at$qD[, ind])

plot(exp(xx), exp(yy), pch = 19, col = rgb(.7, .9, .7, .5), log = "xy",
     xlab = "L", ylab = expression(''^3*D))

spl <- smooth.spline(x = xx, y = yy, spar = 1.15, tol = .0001)
lines(exp(spl$x), exp(spl$y), lwd = 2, col = "red")

abline(v = 10*exp(1:4), col = c("grey", "red", "blue", "black"))

# Figure 5c
ind <- m.ah$q == 0

xx <- log(m.ah$a*10)
yy <- log(m.ah$qD[, ind])

plot(exp(xx), exp(yy), pch = 19, col = rgb(.7, .9, .7, .5), log = "xy",
     xlab = "L", ylab = expression(''^0*D))

spl <- smooth.spline(x = xx, y = yy, spar = 1.15, tol = .0001)
lines(exp(spl$x), exp(spl$y), lwd = 2, col = "red")

abline(v = 10*exp(1:4), col = c("grey", "red", "blue", "black"))

# Figure 5d
ind <- m.ah$q == 3

xx <- log(m.ah$a*10)
yy <- log(m.ah$qD[, ind])

plot(exp(xx), exp(yy), pch = 19, col = rgb(.7, .9, .7, .5), log = "xy",
     xlab = "L", ylab = expression(''^3*D))

spl <- smooth.spline(x = xx, y = yy, spar = 1.15, tol = .0001)
lines(exp(spl$x), exp(spl$y), lwd = 2, col = "red")

abline(v = 10*exp(1:4), col = c("grey", "red", "blue", "black"))

#---#

# load mountain forest data
load("D:/YandexDisk/Science/2016/China/R/DLS.rda")

# compute moments
m.dt <- compute.moments.lin(ta, q = seq(-3,3,.1))
m.ds <- compute.moments.lin(sa, q = seq(-3,3,.1))
m.dh <- compute.moments.lin(hc, q = seq(-3,3,.1))

#---#

# compute and truncate local spectra
mfl.at <- trunc.spectra(local.spectra(m.at, sc = exp(1:4), smooth = 1.15))
mfl.as <- local.spectra(m.as, sc = exp(1:4), smooth = 1.15)
mfl.ah <- trunc.spectra(local.spectra(m.ah, sc = exp(1:4), smooth = 1.15))

mfl.dt <- local.spectra(m.dt, sc = exp(1:4), smooth = 1.15)
mfl.dt <- trunc.spectra(local.spectra(m.dt, sc = exp(1:4), smooth = 1.15))
mfl.ds <- local.spectra(m.ds, sc = exp(1:4), smooth = 1.15)
mfl.dh <- local.spectra(m.dh, sc = exp(1:4), smooth = 1.15)
mfl.dh <- trunc.spectra(local.spectra(m.dh, sc = exp(1:4), smooth = 1.15))

#---#

# Figure 6a
plot(mfl.at$alfa[,1], mfl.at$f[,1], type = "o", pch = 21, bg = "white", 
     ylim = c(-0.2, 0.55), xlim = range(-0.05, 1.5), xlab = "a", ylab = "f")
points(mfl.at$alfa[,2], mfl.at$f[,2], type = "o", pch = 22, bg = "red")
points(mfl.at$alfa[,3], mfl.at$f[,3], type = "o", pch = 24, bg = "blue")
points(mfl.at$alfa[,4], mfl.at$f[,4], type = "o", pch = 21, bg = "black")
abline(h = 0, v = 0)
abline(v = 1, lty = 2)
title(main = "Plain forest - tree layer", line = 0.25)

# Figure 6b
plot(mfl.dt$alfa[,1], mfl.dt$f[,1], type = "o", pch = 21, bg = "white", 
     ylim = c(-0.2, 0.55), xlim = range(-0.05, 1.5), xlab = "a", ylab = "f")
points(mfl.dt$alfa[,2], mfl.dt$f[,2], type = "o", pch = 22, bg = "red")
points(mfl.dt$alfa[,3], mfl.dt$f[,3], type = "o", pch = 24, bg = "blue")
points(mfl.dt$alfa[,4], mfl.dt$f[,4], type = "o", pch = 21, bg = "black")
abline(h = 0, v = 0)
abline(v = 1, lty = 2)
title(main = "Mountain forest - tree layer", line = 0.25)

# Figure 6c
plot(mfl.as$alfa[,1], mfl.as$f[,1], type = "o", pch = 21, bg = "white", 
     ylim = c(-0.45, 0.45), xlim = range(-0.05, 1.5), xlab = "a", ylab = "f")
points(mfl.as$alfa[,2], mfl.as$f[,2], type = "o", pch = 22, bg = "red")
points(mfl.as$alfa[,3], mfl.as$f[,3], type = "o", pch = 24, bg = "blue")
points(mfl.as$alfa[,4], mfl.as$f[,4], type = "o", pch = 21, bg = "black")
abline(h = 0, v = 0)
abline(v = 1, lty = 2)
title(main = "Plain forest - shrub layer", line = 0.25)

# Figure 6d
plot(mfl.ds$alfa[,1], mfl.ds$f[,1], type = "o", pch = 21, bg = "white", 
     ylim = c(-0.45, 0.45), xlim = range(-0.05, 1.5), xlab = "a", ylab = "f")
points(mfl.ds$alfa[,2], mfl.ds$f[,2], type = "o", pch = 22, bg = "red")
points(mfl.ds$alfa[,3], mfl.ds$f[,3], type = "o", pch = 24, bg = "blue")
points(mfl.ds$alfa[,4], mfl.ds$f[,4], type = "o", pch = 21, bg = "black")
abline(h = 0, v = 0)
abline(v = 1, lty = 2)
title(main = "Mountain forest - shrub layer", line = 0.25)

# Figure 6e
plot(mfl.ah$alfa[,1], mfl.ah$f[,1], type = "o", pch = 21, bg = "white", 
     ylim = c(-0.45, 0.55), xlim = range(-0.05, 1.5), xlab = "a", ylab = "f")
points(mfl.ah$alfa[,2], mfl.ah$f[,2], type = "o", pch = 22, bg = "red")
points(mfl.ah$alfa[,3], mfl.ah$f[,3], type = "o", pch = 24, bg = "blue")
points(mfl.ah$alfa[,4], mfl.ah$f[,4], type = "o", pch = 21, bg = "black")
abline(h = 0, v = 0)
abline(v = 1, lty = 2)
legend("topright", legend = paste0("L = ", round(10*exp(1:4))), 
       lwd = 1, pch = c(21,22,24,21), pt.bg = c("white", "red", "blue", "black"))
title(main = "Plain forest - herb layer", line = 0.25)

# Figure 6f
plot(mfl.dh$alfa[,1], mfl.dh$f[,1], type = "o", pch = 21, bg = "white", 
     ylim = c(-0.45, 0.55), xlim = range(-0.05, 1.5), xlab = "a", ylab = "f")
points(mfl.dh$alfa[,2], mfl.dh$f[,2], type = "o", pch = 22, bg = "red")
points(mfl.dh$alfa[,3], mfl.dh$f[,3], type = "o", pch = 24, bg = "blue")
points(mfl.dh$alfa[,4], mfl.dh$f[,4], type = "o", pch = 21, bg = "black")
abline(h = 0, v = 0)
abline(v = 1, lty = 2)
title(main = "Mountain forest - herb layer", line = 0.25)

#---#
