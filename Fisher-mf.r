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

# Fifure 2b
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

