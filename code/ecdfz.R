require(EQL)
require(CVXR)

ecdfz = function (z, ord = 5, method = c("un", "con")) {
  z = as.numeric(z)
  ord.max = max(ord)
  H.max = sapply(1 : ord.max, EQL::hermite, x = z)
  res = output = list()
  status = log.lik.gd = c()
  ind = 1
  if (all(method == "con")) {w.cvxr = w.cvxr.cns} else {w.cvxr = w.cvxr.uncns}
  for (i in ord) {
    res[[ind]] = w.cvxr(cbind(H.max[, (1 : i)]))
    status[ind] = res[[ind]]$status
    log.lik.gd[ind] = -res[[ind]]$optimal_value
    output[[ind]] = list(status = res[[ind]]$status, log.lik.gd = -res[[ind]]$optimal_value, w = as.vector(res[[ind]]$primal_values[[1]]))
    ind = ind + 1
  }
  names(res) = paste("order =", ord)
  names(output) = paste("order =", ord)
  output = list(summary = output, status = status, log.lik.gd = log.lik.gd, res = res)
  class(output) <- "ecdfz"
  return(output)
}

summary.ecdfz = function (object) {
  print(object$summary)
}

print.ecdfz = function (object) {
  print(list(status = object$status, log.lik.gd = object$log.lik.gd))
}

ecdfz.optimal = function (z, ord.max = 20, k = 2, alpha = 0.05, method = c("un", "con"), firstk = FALSE) {
  z = as.numeric(z)
  H = log.lik.gd = c()
  log.lik.gd = c(0, log.lik.gd)
  res = list()
  ord = 1
  conv = "OPTIMAL"
  if (all(method == "con")) {w.cvxr = w.cvxr.cns} else {w.cvxr = w.cvxr.uncns}

  if (firstk) {
    while ((ord <= k)) {
      H = cbind(H, hermite(x = z, n = ord))
      res[[ord]] = w.cvxr(H)
      log.lik.gd[1 + ord] = -res[[ord]]$optimal_value
      ord = ord + 1
    }
    conv = "OPTIMAL"
  } else {
    while ((ord <= k) & (conv == "OPTIMAL")) {
      H = cbind(H, hermite(x = z, n = ord))
      res[[ord]] = w.cvxr(H)
      conv = res[[ord]]$status
      log.lik.gd[1 + ord] = -res[[ord]]$optimal_value
      ord = ord + 1
    }
  }

  while (!w.stop(log.lik.gd, ord, k, alpha) & (ord <= (ord.max + k)) & (conv == "OPTIMAL")) {
    H = cbind(H, hermite(x = z, n = ord))
    res[[ord]] = w.cvxr(H)
    conv = res[[ord]]$status
    log.lik.gd[1 + ord] = -res[[ord]]$optimal_value
    ord = ord + 1
  }

  ord.fitted = length(res)
  ord.optimal.found = w.stop(log.lik.gd, ord.fitted + 1, k, alpha) & (ord.fitted >= k)
  if (ord.fitted >=k & ord.optimal.found & conv == "OPTIMAL") {
    ord.optimal = ord.fitted - k
    H.optimal = H[, (1 : ord.optimal)]
    w.optimal = as.numeric(res[[ord.optimal]]$primal_values[[1]])
    log.lik.gd.optimal = log.lik.gd[1 + ord.optimal]
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


w.cvxr.uncns = function (H) {
  p = ncol(H)
  w <- Variable(p)
  objective <- Maximize(CVXR::sum_entries(Log(H %*% w + 1)))
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
  objective <- Maximize(CVXR::sum_entries(Log(y)))
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
