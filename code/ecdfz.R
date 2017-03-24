require(EQL)
require(optimx)

ecdfz = function (z, ord = 5, w.init = NULL, method = "BFGS") {
  if (is.null(w.init)) {w.init = rep(0, ord)}
  if (ord == 0) {
    res = list(par = NA, value = 0, convergence = NA)
  } else {
    res = optim(w.init, loggaussderiv, ord = ord, z = z, control = list(fnscale = -1), method = method)
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

require(EQL)
require(cvxr)

ecdfz.optimal = function (z, ord.max = 20, k = 2, alpha = 0.05, method = c("un", "con")) {
  z = as.numeric(z)
  H = log.lik.gd = c()
  log.lik.gd = c(0, log.lik.gd)
  res = list()
  ord = 1
  conv = "OPTIMAL"
  if (all(method == "con")) {w.cvxr = w.cvxr.cns} else {w.cvxr = w.cvxr.uncns}

  while ((ord <= k) & (conv == "OPTIMAL")) {
    H = cbind(H, hermite(x = z, n = ord))
    res[[ord]] = w.cvxr(H)
    conv = res[[ord]]$status
    log.lik.gd[1 + ord] = -res[[ord]]$optimal_value
    ord = ord + 1
  }

  if (conv == "OPTIMAL") {
    while (!w.stop(log.lik.gd, ord, k, alpha) & (ord <= (ord.max + k)) & (conv == "OPTIMAL")) {
      H = cbind(H, hermite(x = z, n = ord))
      res[[ord]] = w.cvxr(H)
      conv = res[[ord]]$status
      log.lik.gd[1 + ord] = -res[[ord]]$optimal_value
      ord = ord + 1
    }
  }

  ord.fitted = length(res)
  ord.optimal.found = w.stop(log.lik.gd, ord.fitted, k, alpha)
  if (ord.optimal.found) {
    ord.optimal = ord.fitted - k
    H.optimal = H[, (1 : ord.optimal)]
    w.optimal = as.numeric(res[[ord.optimal]]$primal_values[[1]])
    log.lik.gd.optimal = log.lik.gd[ord.optimal]
  } else {
    ord.optimal = NA
    H.optimal = NA
    w.optimal = NA
    log.lik.gd.optimal = NA
  }

  log.lik.n01 = loglikn01z(z)

  output.optimal <- list(ord.fitted = ord.fitted, ord.optimal.found = ord.optimal.found, ord.optimal = ord.optimal, w.optimal = w.optimal,
         log.lik.gd.optimal = log.lik.gd.optimal, log.lik.n01 = log.lik.n01, log.lik.ecdf.optimal = log.lik.gd.optimal + log.lik.n01,
         log.lik.gd = log.lik.gd)
  output <- list(optimal = output.optimal, H.optimal = H.optimal, H = H, res = res)

  class(output) <- "ecdfz.optimal"
  return(output)
}

summary.ecdfz.optimal <- function (object) {
  return(object$optimal)
}

print.ecdfz.optimal <- function (object) {
  print(summary(object))
}

w.cvxr.uncns = function (H) {
  p = ncol(H)
  w <- Variable(p)
  objective <- Maximize(SumEntries(Log(H %*% w + 1)))
  prob <- Problem(objective)
  capture.output(result <- solve(prob), file = "/dev/null")
  return(result)
}


w.cvxr.cns = function (H) {
  p = ncol(H)
  l = nrow(H)
  w <- Variable(p)
  y <- Variable(l)
  # Hp <- sapply(1 : p, EQL::hermite, x = c(min(z) - 1, max(z) + 1))
  objective <- Maximize(SumEntries(Log(y)))
  constraint <- list(H %*% w + 1 == y
                     #, Hp %*% w > -1
                     )
  prob <- Problem(objective, constraint)
  capture.output(result <- solve(prob), file = "/dev/null")
  return(result)
}

w.stop = function (log.lik.gd, ord, k, alpha) {
  all(2 * (log.lik.gd[ord - ((k - 1) : 0)] - log.lik.gd[ord - k]) <= qchisq(1 - alpha, 1:k))
}
