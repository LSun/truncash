require(EQL)
require(SQUAREM)
require(REBayes)
require(cvxr)
require(PolynomF)
require(Rmosek)

gdash = function (betahat, sebetahat, gd.ord = 10, primal = FALSE, gd.normalized = TRUE, w.lambda = 10, w.rho = 0.5, mixcompdist = "normal", method = "fdr", control = list(maxiter = 50)) {
  if (method == "fdr") {
    sd = c(0, autoselect.mixsd(betahat, sebetahat, mult = sqrt(2)))
    pi_prior = c(10, rep(1, length(sd) - 1))
  } else {
    sd = autoselect.mixsd(betahat, sebetahat, mult = sqrt(2))
    pi_prior = rep(1, length(sd))
  }
  array_F = array_f(betahat, sebetahat, sd, gd.ord, mixcompdist, gd.normalized)
  array_F = aperm(array_F, c(2, 3, 1))
  if (is.null(w.lambda)) {
    w_prior = rep(0, gd.ord)
  } else {
    w_prior = w.lambda / sqrt(w.rho^(1:gd.ord))
    w_prior[seq(1, gd.ord, by = 2)] = 0
  }
  res = biopt(array_F, pi_prior, w_prior, control, primal)
  pihat = res$pihat
  what = res$what
  fitted_g = normalmix(pi = pihat, mean = 0, sd = sd)
  return(list(fitted_g = fitted_g, w = what, niter = res$niter, converged = res$converged))
}

bifixpoint = function(pinw, array_F, pi_prior, w_prior, primal){
  Kpi = dim(array_F)[1]
  Lw = dim(array_F)[2]
  pi = pinw[1 : Kpi]
  w = pinw[-(1 : Kpi)]
  matrix_lik_w = apply(array_F * pi, 2, colSums)
  if (primal) {
    g_current = matrix_lik_w %*% w
    w_new = c(1, w.mosek.primal(matrix_lik_w, w_prior, w.init = c(g_current, w[-1]))$w)
  } else {
    w_new = c(1, w.mosek(matrix_lik_w, w_prior, w.init = w)$w)
  }
  matrix_lik_pi = apply(aperm(array_F, c(2,1,3)) * w_new, 2, colSums)
  pi_new = mixIP(matrix_lik_pi, pi_prior, pi)$pihat
  # w_new = c(1, w.cvxr.uncns(matrix_lik_w, w.init = w)$primal_values[[1]])
  return(c(pi_new, w_new))
}

w.cvxr.uncns = function (matrix_lik_w, w.init = NULL) {
  FF <- matrix_lik_w[, -1]
  f <- matrix_lik_w[, 1]
  p <- ncol(FF)
  w <- Variable(p)
  objective <- Maximize(SumEntries(Log(FF %*% w + f)))
  prob <- Problem(objective)
  if (is.null(w.init)) {
    capture.output(result <- solve(prob), file = "/dev/null")
  } else {
    capture.output(result <- solve(prob, warm_start = w.init[-1]), file = "/dev/null")
  }
  return(result)
}

w.mosek = function (matrix_lik_w, w_prior, w.init = NULL) {
  A = matrix_lik_w[, -1]
  a = matrix_lik_w[,1]
  m = ncol(A)
  n = nrow(A)
  P <- list(sense = "min")
  if (!is.null(w.init)) {
    g.init <- as.vector(matrix_lik_w %*% w.init)
    v.init <- 1 / g.init
    v.init.list <- list(xx = v.init)
    P$sol <- list(itr = v.init.list, bas = v.init.list)
  }
  P$c <- a
  P$A <- Matrix::Matrix(t(A), sparse = TRUE)
  if (is.null(w_prior) | all(w_prior == 0) | missing(w_prior)) {
    P$bc <- rbind(rep(0, m), rep(0, m))
  } else {
    P$bc <- rbind(-w_prior, w_prior)
  }
  P$bx <- rbind(rep(0, n), rep(Inf, n))
  opro <- matrix(list(), nrow = 5, ncol = n)
  rownames(opro) <- c("type", "j", "f", "g", "h")
  opro[1, ] <- as.list(rep("log", n))
  opro[2, ] <- as.list(1:n)
  opro[3, ] <- as.list(rep(-1, n))
  opro[4, ] <- as.list(rep(1, n))
  opro[5, ] <- as.list(rep(0, n))
  P$scopt <- list(opro = opro)
  z <- Rmosek::mosek(P, opts = list(verbose = 0, usesol = TRUE))
  status <- z$sol$itr$solsta
  w <- z$sol$itr$suc - z$sol$itr$slc
  list(w = w, status = status)
}


w.mosek.primal = function (matrix_lik_w, w_prior, w.init = NULL) {
  A = matrix_lik_w[, -1]
  a = matrix_lik_w[,1]
  m = ncol(A)
  n = nrow(A)
  P <- list(sense = "min")
  if (!is.null(w.init)) {
    w.init.list <- list(xx = w.init)
    P$sol <- list(bas = w.init.list)
  }
  P$c <- rep(0, n + m)
  P$A <- Matrix::Matrix(cbind(diag(n), -A), sparse = TRUE)
  P$bc <- rbind(a, a)
  P$bx <- rbind(c(rep(0, n), rep(-Inf, m)),
                c(rep(Inf, n), rep(Inf, m)))
  if (missing(w_prior) | is.null(w_prior) | all(w_prior == 0)) {
    opro <- matrix(list(), nrow = 5, ncol = n)
    rownames(opro) <- c("type", "j", "f", "g", "h")
    opro[1, ] <- as.list(rep("log", n))
    opro[2, ] <- as.list(1 : n)
    opro[3, ] <- as.list(rep(-1, n))
    opro[4, ] <- as.list(rep(1, n))
    opro[5, ] <- as.list(rep(0, n))
  } else {
    opro <- matrix(list(), nrow = 5, ncol = n + m)
    rownames(opro) <- c("type", "j", "f", "g", "h")
    opro[1, ] <- as.list(c(rep("log", n), rep("pow", m)))
    opro[2, ] <- as.list(1 : (n + m))
    opro[3, ] <- as.list(c(rep(-1, n), w_prior))
    opro[4, ] <- as.list(c(rep(1, n), rep(2, m)))
    opro[5, ] <- as.list(rep(0, n + m))
  }
  P$scopt <- list(opro = opro)
  z <- Rmosek::mosek(P, opts = list(verbose = 0, usesol = TRUE))
  status <- z$sol$itr$solsta
  w <- z$sol$itr$xx[-(1 : n)]
  list(w = w, status = status)
}


mixIP = function (matrix_lik, prior, pi_init = NULL, control = list()) {
  if(!requireNamespace("REBayes", quietly = TRUE)) {
    stop("mixIP requires installation of package REBayes")}
  control = set_control_mixIP(control)
  n = nrow(matrix_lik)
  k = ncol(matrix_lik)
  A = rbind(diag(length(prior)),matrix_lik) # add in observations corresponding to prior
  w = c(prior-1,rep(1,n))
  A = A[w!=0,] #remove zero weight entries, as these otherwise cause errors
  w = w[w!=0]
  #w = rep(1,n+k)
  res = REBayes::KWDual(A, rep(1,k), normalize(w), control=control)
  return(list(pihat = normalize(res$f), niter = NULL, converged=(res$status=="OPTIMAL"), control=control))
}

set_control_squarem=function(control,nobs){
  control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE)
  if (nobs > 50000) control.default$trace = TRUE
  control.default$tol = min(0.1/nobs,1.e-7) # set default convergence criteria to be more stringent for larger samples
  namc=names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  control=utils::modifyList(control.default, control)
  return(control)
}

set_control_mixIP=function(control){
  control.default=list(rtol=1e-6)
  namc=names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  control=utils::modifyList(control.default, control)
  return(control)
}

biopt = function (array_F, pi_prior, w_prior, control, primal) {
  control.default = list(K = 1, method = 3, square = TRUE,
                         step.min0 = 1, step.max0 = 1, mstep = 4, kr = 1, objfn.inc = 1,
                         tol = 1e-07, maxiter = 5000, trace = FALSE)
  namc = names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  controlinput = modifyList(control.default, control)
  Kpi = dim(array_F)[1]
  Lw = dim(array_F)[2]
  pinw_init = c(1, rep(0, Kpi - 1), 1, rep(0, Lw - 1))
  res = squarem(par = pinw_init, fixptfn = bifixpoint, objfn = binegpenloglik,
                array_F = array_F, pi_prior = pi_prior, w_prior = w_prior, primal = primal, control = controlinput)
  return(list(pihat = normalize(pmax(0, res$par[1 : Kpi])),
              what = res$par[-(1 : Kpi)],
              B = res$value.objfn,
              niter = res$fpevals,
              converged = res$convergence))
}

normalmix = function (pi, mean, sd) {
  structure(data.frame(pi, mean, sd), class = "normalmix")
}

binegpenloglik = function (pinw, array_F, pi_prior, w_prior, primal)
{
  return(-bipenloglik(pinw, array_F, pi_prior, w_prior))
}

bipenloglik = function (pinw, array_F, pi_prior, w_prior) {
  K = dim(array_F)[1]
  pi = pinw[1 : K]
  w = pinw[-(1 : K)]
  loglik = sum(log(pmax(0, colSums(t(apply(pi * array_F, 2, colSums)) * w))))
  subset = (pi_prior != 1)
  priordens = sum((pi_prior - 1)[subset] * log(pi[subset]))
  return(loglik + priordens - sum(abs(w[-1]) * w_prior))
}

normalize = function (x) {
  return(x/sum(x))
}


array_f = function (betahat, sebetahat, sd, gd.ord, mixcompdist, gd.normalized) {
  if (mixcompdist == "normal") {
    array_f = array_f.normal(betahat, sebetahat, sd, gd.ord, gd.normalized)
  } else {
    stop ("invalid prior mixture")
  }
  return(array_f)
}

## this function is vectorized for x
## more efficient if let it run for x at once
gauss.deriv = function(x, ord) {
  return(dnorm(x) * EQL::hermite(x, ord))
}

Hermite = function (gd.ord) {
  x <- polynom()
  H <- polylist(x, - 1 + x^2)
  if (gd.ord >= 3) {
    for(n in 2 : (gd.ord - 1))
      H[[n+1]] <- x * H[[n]] - n * H[[n-1]]
  }
  return(H)
}

array_f.normal = function (betahat, sebetahat, sd, gd.ord, gd.normalized) {
  sd.mat = sqrt(outer(sebetahat^2, sd^2, FUN = "+"))
  beta.std.mat = betahat / sd.mat
  temp2 = array(dim = c(dim(beta.std.mat), gd.ord + 1))
  temp2[, , 1] = dnorm(beta.std.mat)
  hermite = Hermite(gd.ord)
  if (gd.normalized) {
    for (i in 1 : gd.ord) {
      temp2[, , i + 1] = temp2[, , 1] * hermite[[i]](beta.std.mat) / sqrt(factorial(i))
    }
  } else {
    for (i in 1 : gd.ord) {
      temp2[, , i + 1] = temp2[, , 1] * hermite[[i]](beta.std.mat)
    }
  }
  # temp2.test = outer(beta.std.mat, 0:gd.ord, FUN = gauss.deriv)
  se.std.mat = sebetahat / sd.mat
  temp1 = exp(outer(log(se.std.mat), 0 : gd.ord + 1, FUN = "*"))
  array_f = temp1 * temp2 / sebetahat
  rm(temp1)
  rm(temp2)
  return(array_f)
}

autoselect.mixsd = function (betahat, sebetahat, mult)
{
  sebetahat = sebetahat[sebetahat != 0]
  sigmaamin = min(sebetahat)/10
  if (all(betahat^2 <= sebetahat^2)) {
    sigmaamax = 8 * sigmaamin
  }
  else {
    sigmaamax = 2 * sqrt(max(betahat^2 - sebetahat^2))
  }
  if (mult == 0) {
    return(c(0, sigmaamax/2))
  }
  else {
    npoint = ceiling(log2(sigmaamax/sigmaamin)/log2(mult))
    return(mult^((-npoint):0) * sigmaamax)
  }
}

plot.ghat = function (fitted.g, mixcompdist = "normal", xlim = c(-10, 10)) {
  pi = fitted.g$pi
  mean = fitted.g$mean
  sd = fitted.g$sd
  x = seq(xlim[1], xlim[2], 0.01)
  if (mixcompdist == "normal") {
    y = sum(pi * pnorm(x, mean, sd))
  }
}

ghat.cdf = function (x, fitted.g) {
  pi = fitted.g$pi
  mean = fitted.g$mean
  sd = fitted.g$sd
  return(sum(pi * pnorm(x, mean, sd)))
}
