library(ashr)
library(SQUAREM)

truncash = function(betahat, sebetahat, t) {
	# break the observations into 2 groups: moderate and extreme
	I = (abs(betahat/sebetahat) <= t)
	betahat1 = betahat[I]
	sebetahat1 = sebetahat[I]
	betahat2 = betahat[!I]
	sebetahat2 = sebetahat[!I]
	m = length(betahat1)
	n = length(betahat2)
	
	# if no elements in the extreme group, prior is set to be a point mass at 0
	# or else the prior is estimated
	if (n == 0) {
		fitted.g = normalmix(pi = 1, mean = 0, sd = 0)
	} else {
		# the grid of sd estimated from the extreme group only
		sd = c(0, autoselect.mixsd(betahat2, sebetahat2, mult = sqrt(2)))
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
		
		# estimate pihat with IP
		pihat.est = mixIP(matrix_lik = lik.mat, prior, pi_init)
		
		fitted.g = normalmix(pi = pihat.est$pihat, mean = 0, sd = sd)
	}
	output = ash.workhorse(betahat, sebetahat, g = fitted.g, fixg = TRUE)
	return(output)
}

normalmix = function (pi, mean, sd) {
	structure(data.frame(pi, mean, sd), class = "normalmix")
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