# Local multifractal analysis of forest layers

source("ecomf.r")

# loading plain forest data
library(readxl)
# r is for russian plain forest, tsh stands for tree, shrub and herb layers respectively
rt <- t(read_excel("Forest_data.xlsx", sheet = 1, range = "B2:CW12", col_names = F))
rs <- t(read_excel("Forest_data.xlsx", sheet = 2, range = "B2:CW11", col_names = F))
rh <- t(read_excel("Forest_data.xlsx", sheet = 3, range = "B2:CW43", col_names = F))

# compute moments with a linear scheme
m.rt <- compute.moments.lin(rt, q = seq(-3, 3, .1))
m.rs <- compute.moments.lin(rs, q = seq(-3, 3, .1))
m.rh <- compute.moments.lin(rh, q = seq(-3, 3, .1))

#--- diversity scaling plots in the plain forest ---#

# Figure 5a

# png("Fig5a.png", height = 800, width = 1200)
# opar <- par(cex = 2, mar = c(4, 4.2, 0.5, 0.5))

ind <- m.rt$q == 0

xx <- log(m.rt$a*10)
yy <- log(m.rt$qD[, ind])

plot(exp(xx), exp(yy), pch = 19, col = rgb(.7, .9, .7, .5), log = "xy",
     xlab = "L", ylab = expression(''^0*D))

spl <- smooth.spline(x = xx, y = yy, spar = 1.15, tol = .0001)
lines(exp(spl$x), exp(spl$y), lwd = 2, col = "red")

abline(v = 10*exp(1:4), col = c("grey", "red", "blue", "black"))

# par(opar)
# dev.off()

# Figure 5b

# png("Fig5b.png", height = 800, width = 1200)
# opar <- par(cex = 2, mar = c(4, 4.2, 0.5, 0.5))

ind <- m.rt$q == 3

xx <- log(m.rt$a*10)
yy <- log(m.rt$qD[, ind])

plot(exp(xx), exp(yy), pch = 19, col = rgb(.7, .9, .7, .5), log = "xy",
     xlab = "L", ylab = expression(''^3*D))

spl <- smooth.spline(x = xx, y = yy, spar = 1.15, tol = .0001)
lines(exp(spl$x), exp(spl$y), lwd = 2, col = "red")

abline(v = 10*exp(1:4), col = c("grey", "red", "blue", "black"))

# par(opar)
# dev.off()

# Figure 5c

# png("Fig5c.png", height = 800, width = 1200)
# opar <- par(cex = 2, mar = c(4, 4.2, 0.5, 0.5))

ind <- m.rh$q == 0

xx <- log(m.rh$a*10)
yy <- log(m.rh$qD[, ind])

plot(exp(xx), exp(yy), pch = 19, col = rgb(.7, .9, .7, .5), log = "xy",
     xlab = "L", ylab = expression(''^0*D))

spl <- smooth.spline(x = xx, y = yy, spar = 1.15, tol = .0001)
lines(exp(spl$x), exp(spl$y), lwd = 2, col = "red")

abline(v = 10*exp(1:4), col = c("grey", "red", "blue", "black"))

# par(opar)
# dev.off()

# Figure 5d

# png("Fig5d.png", height = 800, width = 1200)
# opar <- par(cex = 2, mar = c(4, 4.2, 0.5, 0.5))

ind <- m.rh$q == 3

xx <- log(m.rh$a*10)
yy <- log(m.rh$qD[, ind])

plot(exp(xx), exp(yy), pch = 19, col = rgb(.7, .9, .7, .5), log = "xy",
     xlab = "L", ylab = expression(''^3*D))

spl <- smooth.spline(x = xx, y = yy, spar = 1.15, tol = .0001)
lines(exp(spl$x), exp(spl$y), lwd = 2, col = "red")

abline(v = 10*exp(1:4), col = c("grey", "red", "blue", "black"))

# par(opar)
# dev.off()

#---#

# load mountain forest data
# c is for chinese mountain forest, tsh stands for tree, shrub and herb layers respectively
ct <- t(read_excel("Forest_data.xlsx", sheet = 4, range = "B2:CS21", col_names = F))
cs <- t(read_excel("Forest_data.xlsx", sheet = 5, range = "B2:CS42", col_names = F))
ch <- t(read_excel("Forest_data.xlsx", sheet = 6, range = "B2:CS194", col_names = F))

# compute moments
m.ct <- compute.moments.lin(ct, q = seq(-3,3,.1))
m.cs <- compute.moments.lin(cs, q = seq(-3,3,.1))
m.ch <- compute.moments.lin(ch, q = seq(-3,3,.1))

#---#

# compute and truncate local spectra
mfl.rt <- trunc.spectra(local.spectra(m.rt, sc = exp(1:4), smooth = 1.15))
mfl.rs <- local.spectra(m.rs, sc = exp(1:4), smooth = 1.15)
mfl.rh <- trunc.spectra(local.spectra(m.rh, sc = exp(1:4), smooth = 1.15))

mfl.ct <- trunc.spectra(local.spectra(m.ct, sc = exp(1:4), smooth = 1.15))
mfl.cs <- local.spectra(m.cs, sc = exp(1:4), smooth = 1.15)
mfl.ch <- trunc.spectra(local.spectra(m.ch, sc = exp(1:4), smooth = 1.15))

#--- local multifractal spectra plots ---#

# Figure 6a

# png("Fig6a.png", height = 800, width = 1200)
# opar <- par(cex = 2, mar = c(4, 4.2, 1, 0.5))

plot(mfl.rt$alfa[,1], mfl.rt$f[,1], type = "o", pch = 21, bg = "white", 
     ylim = c(-0.2, 0.55), xlim = range(-0.05, 1.5), xlab = "a", ylab = "f")
points(mfl.rt$alfa[,2], mfl.rt$f[,2], type = "o", pch = 22, bg = "red")
points(mfl.rt$alfa[,3], mfl.rt$f[,3], type = "o", pch = 24, bg = "blue")
points(mfl.rt$alfa[,4], mfl.rt$f[,4], type = "o", pch = 21, bg = "black")
abline(h = 0, v = 0)
abline(v = 1, lty = 2)
title(main = "Plain forest - tree layer", line = 0.25)

# par(opar)
# dev.off()

# Figure 6b

png("Fig6b.png", height = 800, width = 1200)
opar <- par(cex = 2, mar = c(4, 4.2, 1, 0.5))

plot(mfl.ct$alfa[,1], mfl.ct$f[,1], type = "o", pch = 21, bg = "white", 
     ylim = c(-0.2, 0.55), xlim = range(-0.05, 1.5), xlab = "a", ylab = "f")
points(mfl.ct$alfa[,2], mfl.ct$f[,2], type = "o", pch = 22, bg = "red")
points(mfl.ct$alfa[,3], mfl.ct$f[,3], type = "o", pch = 24, bg = "blue")
points(mfl.ct$alfa[,4], mfl.ct$f[,4], type = "o", pch = 21, bg = "black")
abline(h = 0, v = 0)
abline(v = 1, lty = 2)
title(main = "Mountain forest - tree layer", line = 0.25)

# par(opar)
# dev.off()

# Figure 6c

# png("Fig6c.png", height = 800, width = 1200)
# opar <- par(cex = 2, mar = c(4, 4.2, 1, 0.5))

plot(mfl.rs$alfa[,1], mfl.rs$f[,1], type = "o", pch = 21, bg = "white", 
     ylim = c(-0.45, 0.45), xlim = range(-0.05, 1.5), xlab = "a", ylab = "f")
points(mfl.rs$alfa[,2], mfl.rs$f[,2], type = "o", pch = 22, bg = "red")
points(mfl.rs$alfa[,3], mfl.rs$f[,3], type = "o", pch = 24, bg = "blue")
points(mfl.rs$alfa[,4], mfl.rs$f[,4], type = "o", pch = 21, bg = "black")
abline(h = 0, v = 0)
abline(v = 1, lty = 2)
title(main = "Plain forest - shrub layer", line = 0.25)

# par(opar)
# dev.off()

# Figure 6d

# png("Fig6d.png", height = 800, width = 1200)
# opar <- par(cex = 2, mar = c(4, 4.2, 1, 0.5))

plot(mfl.cs$alfa[,1], mfl.cs$f[,1], type = "o", pch = 21, bg = "white", 
     ylim = c(-0.45, 0.45), xlim = range(-0.05, 1.5), xlab = "a", ylab = "f")
points(mfl.cs$alfa[,2], mfl.cs$f[,2], type = "o", pch = 22, bg = "red")
points(mfl.cs$alfa[,3], mfl.cs$f[,3], type = "o", pch = 24, bg = "blue")
points(mfl.cs$alfa[,4], mfl.cs$f[,4], type = "o", pch = 21, bg = "black")
abline(h = 0, v = 0)
abline(v = 1, lty = 2)
title(main = "Mountain forest - shrub layer", line = 0.25)

# par(opar)
# dev.off()

# Figure 6e

# png("Fig6e.png", height = 800, width = 1200)
# opar <- par(cex = 2, mar = c(4, 4.2, 1, 0.5))

plot(mfl.rh$alfa[,1], mfl.rh$f[,1], type = "o", pch = 21, bg = "white", 
     ylim = c(-0.45, 0.55), xlim = range(-0.05, 1.5), xlab = "a", ylab = "f")
points(mfl.rh$alfa[,2], mfl.rh$f[,2], type = "o", pch = 22, bg = "red")
points(mfl.rh$alfa[,3], mfl.rh$f[,3], type = "o", pch = 24, bg = "blue")
points(mfl.rh$alfa[,4], mfl.rh$f[,4], type = "o", pch = 21, bg = "black")
abline(h = 0, v = 0)
abline(v = 1, lty = 2)
legend("topright", legend = paste0("L = ", round(10*exp(1:4))), 
       lwd = 1, pch = c(21,22,24,21), pt.bg = c("white", "red", "blue", "black"))
title(main = "Plain forest - herb layer", line = 0.25)

# par(opar)
# dev.off()

# Figure 6f

# png("Fig6f.png", height = 800, width = 1200)
# opar <- par(cex = 2, mar = c(4, 4.2, 1, 0.5))

plot(mfl.ch$alfa[,1], mfl.ch$f[,1], type = "o", pch = 21, bg = "white", 
     ylim = c(-0.45, 0.55), xlim = range(-0.05, 1.5), xlab = "a", ylab = "f")
points(mfl.ch$alfa[,2], mfl.ch$f[,2], type = "o", pch = 22, bg = "red")
points(mfl.ch$alfa[,3], mfl.ch$f[,3], type = "o", pch = 24, bg = "blue")
points(mfl.ch$alfa[,4], mfl.ch$f[,4], type = "o", pch = 21, bg = "black")
abline(h = 0, v = 0)
abline(v = 1, lty = 2)
title(main = "Mountain forest - herb layer", line = 0.25)

# par(opar)
# dev.off()

#--- analysis of minimal abundance scaling ---#

# Figure 7
spar <- 1.15

xx1 <- log(10*m.rt$a)
yy1 <- log(m.rt$nmin)
spl1 <- smooth.spline(x = xx1, y = yy1, spar = spar, tol = .0001)

xx2 <- log(10*m.rs$a)
yy2 <- log(m.rs$nmin)
spl2 <- smooth.spline(x = xx2, y = yy2, spar = spar, tol = .0001)

xx3 <- log(10*m.rh$a)
yy3 <- log(m.rh$nmin)
spl3 <- smooth.spline(x = xx3, y = yy3, spar = spar, tol = .0001)

# png("Fig7.png", height = 800, width = 1200)
# opar <- par(cex = 2, mar = c(4, 4.2, 0.5, 0.5))

plot(exp(spl1$x), exp(spl1$y), log = "xy", 
     type = "n", xlab = "L", ylab = expression(n[min]))

abline(v = 10*exp(1:4), col = c("grey", "red", "blue", "black"))

lines(exp(spl1$x), exp(spl1$y), lwd = 3, col = "darkred")
lines(exp(spl2$x), exp(spl2$y), lwd = 3, col = "darkgreen")
lines(exp(spl3$x), exp(spl3$y), lwd = 3, col = "darkblue")

legend("bottomleft", legend = c("tree layer", "shrub layer", "herb layer"), horiz = T,
       lwd = 3, col  = c("darkred", "darkgreen", "darkblue"))

# par(opar)
# dev.off()

#--- analysis of scaling of a dominant species ---#

# Figure 8
xx <- log(10*m.rt$a[m.rt$pp[, 1] > 0])
yy <- log(m.rt$pp[m.rt$pp[, 1] > 0, 1])

spl1 <- smooth.spline(x = xx, y = yy, spar = 1.15, tol = .0001)

# png("Fig8.png", height = 800, width = 1200)
# opar <- par(cex = 2, mar = c(4, 4.2, 0.5, 0.5))

plot(exp(xx), exp(yy), pch = 19, col = rgb(.5,.6,.1,.5), xlab = "L", 
     ylab = expression(p[Acer]), ylim = range(exp(yy)), log = "xy")
lines(exp(spl1$x), exp(spl1$y), lwd = 3, col = "darkgreen")

# par(opar)
# dev.off()