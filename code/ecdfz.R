require(EQL)
require(optimx)

ecdfz = function (z, ord = 5, w.init = NULL) {
  if (is.null(w.init)) {w.init = rep(0, ord)}
  if (ord == 0) {
    res = list(par = NA, value = 0, convergence = NA)
  } else {
    res = optim(w.init, loggaussderiv, ord = ord, z = z, control = list(fnscale = -1), method = "BFGS")
  }
  return(list(gaussianDerivOrder = ord,
              w = res$par,
              loglik.standardNormal = loglikn01z(z),
              loglik.gaussianDeriv = res$value,
              loglik.empiricalDist = res$value + loglikn01z(z),
              convergence = as.logical(!res$convergence)))
}

loggaussderiv = function (w, ord, z) {
  if (ord == 0) {
    return(0)
  } else {
    H = sapply(1 : ord, EQL::hermite, x = z)
    Hw = H %*% w + 1
    if (all(Hw > 0)) {
      return(sum(log(H %*% w + 1)))
    } else {
      return(-Inf)
    }
  }
}

loglikn01z = function (z) {
  return(sum(log(dnorm(z))))
}

loglikecdfz = function (w, ord, z) {
  return(loggaussderiv(w, ord, z) + loglikn01z(z))
}
