gdfit = function (z, gd.ord, w.lambda = NULL, w.rho = 0.5) {
  if (is.null(w.lambda)) {
    w_prior = rep(0, gd.ord)
  } else {
    w_prior = w.lambda / sqrt(w.rho^(1 : gd.ord))
    w_prior[seq(1, gd.ord, by = 2)] = 0
  }
  hermite = Hermite(gd.ord)
  gd0.std = dnorm(z)
  matrix_lik_w = cbind(gd0.std)
  for (i in 1 : gd.ord) {
    gd.std = (-1)^i * hermite[[i]](z) * gd0.std / sqrt(factorial(i))
    matrix_lik_w = cbind(matrix_lik_w, gd.std)
  }
  w.fit = w.mosek(matrix_lik_w, w_prior, w.init = NULL)
  w.hat = c(1, w.fit$w)
  w.status = w.fit$status
  loglik.hat = sum(log(matrix_lik_w %*% w.hat))
  return(list(gd.ord = gd.ord, w = w.hat, loglik = loglik.hat, status = w.status))
}

gdfit.mom = function (z, gd.ord) {
  hermite.list = orthopolynom::hermite.he.polynomials(gd.ord)
  hermite.coef = orthopolynom::polynomial.coefficients(hermite.list)
  moments = c()
  for (i in 0 : gd.ord) {
    moments[i + 1] = mean(z^i)
  }
  w = c()
  for (i in 0 : gd.ord) {
    w[i + 1] = sum(moments[1 : (i + 1)] / sqrt(factorial(i)) * hermite.coef[[i + 1]]) * (-1)^i
  }
  return(list(gd.ord = gd.ord, w = w))
}

plot.gdfit = function (z, w, gd.ord, symm = TRUE, breaks = 100, std.norm = TRUE) {
  if (symm) {
    x.plot = seq(- max(abs(z)) - 2, max(abs(z)) + 2, length = 1000)
  } else {
    x.plot = seq(min(z) - 2, max(z) + 2, length = 1000)
  }
  hermite = Hermite(gd.ord)
  gd0.std = dnorm(x.plot)
  matrix_lik_plot = cbind(gd0.std)
  for (i in 1 : gd.ord) {
    gd.std = (-1)^i * hermite[[i]](x.plot) * gd0.std / sqrt(factorial(i))
    matrix_lik_plot = cbind(matrix_lik_plot, gd.std)
  }
  y.plot = matrix_lik_plot %*% w
  z.hist = hist(z, breaks, plot = FALSE)
  if (std.norm) {
    y.max = max(z.hist$density, y.plot, dnorm(0))
  } else {
    y.max = max(z.hist$density, y.plot)
  }
  hist(z, breaks, prob = TRUE, ylim = c(0, y.max))
  lines(x.plot, y.plot, col = "blue")
  legend("topright", lty = 1, col = "blue", "GD")
  if (std.norm) {
    lines(x.plot, dnorm(x.plot), col = "red")
    legend("topleft", lty = 1, col = "red", "N(0, 1)")
  }
}
