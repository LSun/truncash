require(EQL)
require(SQUAREM)
require(REBayes)
require(cvxr)

gdash = function (betahat, sebetahat, gd.ord, mixcompdist = "normal", method = "fdr", control = list()) {
  if (method == "fdr") {
    sd = c(0, autoselect.mixsd(betahat, sebetahat, mult = sqrt(2)))
    pi_prior = c(10, rep(1, length(sd) - 1))
  } else {
    sd = autoselect.mixsd(betahat, sebetahat, mult = sqrt(2))
    pi_prior = rep(1, length(sd))
  }
  array_f = array_f(betahat, sebetahat, sd, gd.ord, mixcompdist)
  res = biopt(array_f, pi_prior, pi_init = NULL, w_prior = NULL, control)
  pihat = res$pihat
  what = res$what
  fitted.g = normalmix(pi = pihat, mean = 0, sd = sd)
  return(fitted.g = fitted.g)
}

bifixpoint = function(pinw, array_f, pi_prior){
  Kpi = dim(array_f)[1]
  Lw = dim(array_f)[2]
  pi = pinw[1:Kpi]
  w = pinw[-(1:Kpi)]
  matrix_lik_pi = apply(aperm(array_f, c(2,1,3)) * w, 2, colSums)
  pi_new = mixIP(matrix_lik_pi, pi_prior, pi)$pihat
  matrix_lik_w = apply(array_f * pi_new, 2, colSums)
  w_new = c(1, w.cvxr.uncns(matrix_lik_w)$primal_values[[1]])
  return(c(pi_new, w_new))
}

w.cvxr.uncns = function (matrix_lik_w) {
  p = ncol(matrix_lik_w) - 1
  w <- Variable(p)
  objective <- Maximize(SumEntries(Log(matrix_lik_w[, -1] %*% w + cbind(matrix_lik_w[, 1]))))
  prob <- Problem(objective)
  capture.output(result <- solve(prob), file = "/dev/null")
  return(result)
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

biopt = function (array_f, pi_prior, pi_init = NULL, w_prior = NULL, control) {
  control.default = list(K = 1, method = 3, square = TRUE,
                         step.min0 = 1, step.max0 = 1, mstep = 4, kr = 1, objfn.inc = 1,
                         tol = 1e-07, maxiter = 5000, trace = FALSE)
  namc = names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  controlinput = modifyList(control.default, control)
  Kpi = dim(array_f)[1]
  Lw = dim(array_f)[2]
  if (is.null(pi_init)) {
    pinw_init = c(rep(1/Kpi, Kpi), 1, rep(0, Lw - 1))
  } else {
    pinw_init = c(pi_init, 1, rep(0, Lw - 1))
  }
  res = squarem(par = pinw_init, fixptfn = bifixpoint, objfn = binegpenloglik,
                array_f = array_f, pi_prior = pi_prior, control = controlinput)
  return(list(pihat = normalize(pmax(0, res$par[1 : Kpi])),
              what = res$par[-(1 : Kpi)],
              B = res$value.objfn,
              niter = res$iter,
              converged = res$convergence))
}

normalmix = function (pi, mean, sd) {
  structure(data.frame(pi, mean, sd), class = "normalmix")
}

binegpenloglik = function (pinw, array_f, pi_prior)
{
  return(-bipenloglik(pinw, array_f, pi_prior))
}

bipenloglik = function (pinw, array_f, pi_prior) {
  K = dim(array_f)[1]
  pi = pinw[1 : K]
  w = pinw[-(1 : K)]
  loglik = sum(log(colSums(t(apply(pi * array_f, 2, colSums)) * w)))
  subset = (pi_prior != 1)
  priordens = sum((pi_prior - 1)[subset] * log(pi[subset]))
  return(loglik + priordens)
}

normalize = function (x) {
  return(x/sum(x))
}


array_f = function (betahat, sebetahat, sd, gd.ord, mixcompdist) {
  if (mixcompdist == "normal") {
    array_f = array_f.normal(betahat, sebetahat, sd, gd.ord)
  } else {
    stop ("invalid prior mixture")
  }
  return(array_f)
}

array_f.normal = function (betahat, sebetahat, sd, gd.ord) {
  n = length(betahat)
  K = length(sd)
  L = gd.ord + 1
  array_f = array(0, dim = c(K, L, n))
  for (j in 1:n) {
    for (k in 1:K) {
      for (l in 0:gd.ord) {
        array_f[k, (l + 1), j] =
          sebetahat[j]^l /
          sqrt(sd[k]^2 + sebetahat[j]^2)^(l + 1) *
          EQL::hermite(betahat[j] / sqrt(sd[k]^2 + sebetahat[j]^2), l) *
          dnorm(betahat[j] / sqrt(sd[k]^2 + sebetahat[j]^2))
      }
    }
  }
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
