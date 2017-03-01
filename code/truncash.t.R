library(ashr)
library(SQUAREM)

truncash = function (betahat, sebetahat, df, pval.thresh,
                     method = c("fdr", "shrink"),
                     mixcompdist = c("uniform", "halfuniform", "normal",
                                     "+uniform", "-uniform", "halfnormal"),
                     gridmult = sqrt(2), grange = c(-Inf,Inf),
                     nullweight = 10, pointmass = TRUE,
                     prior = c("nullbiased", "uniform", "unit")
                     ) {

  # check arguments for specifying priors
  if(!missing(pointmass) & !missing(method))
    stop("Specify either method or pointmass, not both")
  if(!missing(prior) & !missing(method))
    stop("Specify either method or prior, not both")
  if(!missing(method)){
    method = match.arg(method)
    if (method == "shrink"){pointmass = FALSE; prior = "uniform"}
    if (method == "fdr"){pointmass = TRUE; prior = "nullbiased"}
  }

  ## Check to see if is Inf, then switch to NULL.
  if (!is.null(df)) {
    if (df == Inf) {
      df <- NULL
    }
  }

  # First of all, separate the observations into 2 groups: moderate and extreme
  # given a p value threshold pval.thresh
  # Calculate marginal p values for each observation
  # Do normal or t calcultions according to df

  if (is.null(df)) {
    # calculate marginal p values for given data using normal distribution
    pval = (1 - pnorm(abs(betahat / sebetahat))) * 2
  } else {
    pval = (1 - pt(abs(betahat / sebetahat), df)) * 2
  }

  # break the observations into 2 groups: moderate and extreme
  I = (pval <= pval.thresh)
  betahat1 = betahat[I]
  sebetahat1 = sebetahat[I]
  m = length(betahat1)
  betahat2 = betahat[!I]
  sebetahat2 = sebetahat[!I]
  n = length(betahat2)

  # if no elements in the extreme group, prior is set to be a point mass at 0
  # or else the prior is estimated
  if (n == 0) {
    fitted.g = normalmix(pi = 1, mean = 0, sd = 0)
  } else {
    # the grid of sd estimated from both groups together
    mixsd = autoselect.mixsd(list(x = betahat, s = sebetahat), gridmult, mode = 0, grange, mixcompdist)
    if(pointmass){
      mixsd = c(0, mixsd)
    }
    k = length(sd)

    # the likelihood matrix for the moderate group
    sd.mat.1 = sqrt(outer(sebetahat1^2, sd^2, FUN = "+"))
    lik.mat.1 = apply(sd.mat.1, 2, FUN = function(x) {2 * pnorm(t * sebetahat1 / x) - 1})

    # the likelihood matrix for the extreme group
    sd.mat.2 = sqrt(outer(sebetahat2^2, sd^2, FUN = "+"))
    lik.mat.2 = apply(sd.mat.2, 2, FUN = function(x) {dnorm(betahat2, 0, x)})

    # combine the likelihood matrices for two groups
    lik.mat = rbind(lik.mat.1, lik.mat.2)

    # specify prior/penalty and initial value of pihat
    prior = c(10, rep(1, k - 1))
    pi_init = c(0.5, rep(0.5/(k - 1), k - 1))

    # the last of these conditions checks whether the gradient at the null is negative wrt pi0
    # to avoid running the optimization when the global null (pi0=1) is the optimal.
    if(max(prior[-1])>1 || min(gradient(matrix_lik = lik.mat)+prior[1]-1,na.rm=TRUE)<0){
      pihat.est = mixIP(matrix_lik = lik.mat, prior, pi_init)
      fitted.g = normalmix(pi = pihat.est$pihat, mean = 0, sd = sd)
    } else {
      fitted.g = normalmix(pi = 1, mean = 0, sd = 0)
    }

    # estimate pihat with IP
    # pihat.est = mixIP(matrix_lik = lik.mat, prior, pi_init)

    # fitted.g = normalmix(pi = pihat.est$pihat, mean = 0, sd = sd)
  }

  output = ash.workhorse(betahat, sebetahat, df = df, g = fitted.g, fixg = TRUE)
  return(output)
}

normalmix = function (pi, mean, sd) {
  structure(data.frame(pi, mean, sd), class = "normalmix")
}

# try to select a default range for the sigmaa values
# that should be used, based on the values of betahat and sebetahat
# mode is the location about which inference is going to be centered
# mult is the multiplier by which the sds differ across the grid
# grange is the user-specified range of mixsd
autoselect.mixsd = function(data,mult,mode,grange,mixcompdist){
  betahat = data$x - mode
  sebetahat = data$s
  exclude = get_exclusions(data)
  betahat = betahat[!exclude]
  sebetahat = sebetahat[!exclude]

  sigmaamin = min(sebetahat)/10 #so that the minimum is small compared with measurement precision
  if(all(betahat^2<=sebetahat^2)){
    sigmaamax = 8*sigmaamin #to deal with the occassional odd case where this could happen; 8 is arbitrary
  }else{
    sigmaamax = 2*sqrt(max(betahat^2-sebetahat^2)) #this computes a rough largest value you'd want to use, based on idea that sigmaamax^2 + sebetahat^2 should be at least betahat^2
  }

  if(mixcompdist=="halfuniform"){
    sigmaamax = min(max(abs(grange-mode)), sigmaamax)
  }else{
    sigmaamax = min(min(abs(grange-mode)), sigmaamax)
  }


  if(mult==0){
    return(c(0,sigmaamax/2))
  }else{
    npoint = ceiling(log2(sigmaamax/sigmaamin)/log2(mult))
    return(mult^((-npoint):0) * sigmaamax)
  }
}

get_exclusions = function (data) {
  return((data$s == 0 | data$s == Inf | is.na(data$x) | is.na(data$s)))
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

set_control_mixIP=function(control){
  control.default=list(rtol=1e-6)
  namc=names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  control=utils::modifyList(control.default, control)
  return(control)
}

normalize = function(x){return(x/sum(x))}

mixEM = function(matrix_lik,prior,pi_init=NULL,control=list()){
  control = set_control_squarem(control,nrow(matrix_lik))
  k=dim(matrix_lik)[2]
  if(is.null(pi_init)){
    pi_init = rep(1/k,k)# Use as starting point for pi
  }
  res = squarem(par=pi_init,fixptfn=fixpoint, objfn=negpenloglik,matrix_lik=matrix_lik, prior=prior, control=control)
  return(list(pihat = normalize(pmax(0,res$par)), B=res$value.objfn, niter = res$iter, converged=res$convergence, control=control))
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

fixpoint = function(pi, matrix_lik, prior){
  pi = normalize(pmax(0,pi)) #avoid occasional problems with negative pis due to rounding
  m = t(pi * t(matrix_lik)) # matrix_lik is n by k; so this is also n by k
  m.rowsum = rowSums(m)
  classprob = m/m.rowsum #an n by k matrix
  pinew = normalize(colSums(classprob) + prior - 1)
  return(pinew)
}

negpenloglik = function(pi,matrix_lik,prior){return(-penloglik(pi,matrix_lik,prior))}

penloglik = function(pi, matrix_lik, prior){
  pi = normalize(pmax(0,pi))
  m = t(pi * t(matrix_lik)) # matrix_lik is n by k; so this is also n by k
  m.rowsum = rowSums(m)
  loglik = sum(log(m.rowsum))
  subset = (prior != 1.0)
  priordens = sum((prior-1)[subset]*log(pi[subset]))
  return(loglik+priordens)
}

gradient = function(matrix_lik){
  n = nrow(matrix_lik)
  grad = n - ColsumModified(matrix_lik)
  return(grad)
}

ColsumModified = function(matrix_l){
  small = abs(matrix_l) < 10e-100
  matrix_l[small] = matrix_l[small]+10e-100
  colSums(matrix_l/matrix_l[,1])
}
