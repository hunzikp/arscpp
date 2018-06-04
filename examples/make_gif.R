# Normal log pdf
f <- function(x, mu, sd) {
  -0.5*(sd^(-2))*(x-mu)^2
}

# Derivative of normal log pdf
f_prime <- function(x, mu, sd) {
  -1*(sd^(-2))*(x-mu)
}

# Set up
n_sample <- 100
mu <- 0
sd <- 1
sampler <- arscpp::ars(f, f_prime, c(-2, 5), mu = mu, sd = sd)
xr <- seq(-4, 4, length = 400)
xh <- exp(f(xr, mu, sd))
get_upper <- function() {exp(sapply(xr, function(x) sampler$get_hu(x)))}
get_lower <- function() {exp(sapply(xr, function(x) sampler$get_hl(x)))}

# Sample
if (!dir.exists("examples/gif")) {
  dir.create("examples/gif")
}
par(mar = c(2,2,1,1))
smp <- c()
for (i in 1:(n_sample+1)) {

  filename <- paste0('pc', formatC(i, width = 4, format = "d", flag = "0"), '.png')

  png(file.path("examples/gif", filename), width = 400, height = 400)
  up <- get_upper()
  lo <- get_lower()
  plot(xr, xh, type = 'l', xlab = '', ylab = '', ylim = c(0, 2), bty="n") # , yaxt = 'n'
  points(xr, up, type = 'l', col = 'red', lty = 2)
  points(xr, lo, type = 'l', col = 'blue', lty = 2)
  points(smp, y = rep(0, length(smp)), pch = "|", col = 'black')
  text(x = -2.5, y = 1.5, labels = paste0('Samples: ', i-1), cex = 1.25)
  dev.off()

  smp <- c(smp, sampler$sample(1))
}

# Make gif & clean up
system("convert examples/gif/*.png -delay 5 -loop 0 examples/gif/ars.gif")
unlink("examples/gif/*.png")
