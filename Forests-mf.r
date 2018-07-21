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


