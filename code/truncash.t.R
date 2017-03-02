library(ashr)
library(SQUAREM)

truncash.t = function (betahat, sebetahat, df = NULL, pval.thresh,
                     method = c("fdr", "shrink"),
                     mixcompdist = c("uniform", "halfuniform", "normal",
                                     "+uniform", "-uniform", "halfnormal"),
                     gridmult = sqrt(2), grange = c(-Inf,Inf),
                     nullweight = 10, pointmass = TRUE, mode = 0,
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
  # group1 is moderate group; group2 is extreme group
  I = (pval <= pval.thresh)
  betahat1 = betahat[I]
  sebetahat1 = sebetahat[I]
  m = length(betahat1)
  betahat2 = betahat[!I]
  sebetahat2 = sebetahat[!I]
  n = length(betahat2)
  if (!is.null(df)) {df1 = df[I]; df2 = df[!I]}

  # if no elements in the extreme group, prior is set to be a point mass at 0
  # or else the prior is estimated
  if (n == 0) {
    fitted.g = normalmix(pi = 1, mean = 0, sd = 0)
  } else {

    ## Generating mixture distribution g

    # the grid of sd estimated from both groups together
    mixsd = autoselect.mixsd(list(x = betahat, s = sebetahat), gridmult, mode, grange, mixcompdist)
    if(pointmass){
      mixsd = c(0, mixsd)
    }
    null.comp = which.min(mixsd) #which component is the "null"

    k = length(mixsd)
    prior = setprior(prior, k, nullweight, null.comp)
    pi = initpi(k, length(betahat), null.comp)

    if(!is.element(mixcompdist,c("normal","uniform","halfuniform","+uniform","-uniform","halfnormal")))
      stop("Error: invalid type of mixcompdist")
    if(mixcompdist == "normal") g=normalmix(pi,rep(mode,k),mixsd)
    if(mixcompdist == "uniform") g=unimix(pi,mode - mixsd,mode + mixsd)
    if(mixcompdist == "+uniform") g = unimix(pi,rep(mode,k),mode+mixsd)
    if(mixcompdist == "-uniform") g = unimix(pi,mode-mixsd,rep(mode,k))
    if(mixcompdist == "halfuniform"){

      if(min(mixsd)>0){ #simply reflect the components
        pi = c(pi[mode-mixsd>=min(grange)],pi[mode+mixsd<=max(grange)])
        pi = pi/sum(pi)
        g = unimix(pi,c((mode-mixsd)[mode-mixsd>=min(grange)],rep(mode,sum(mode+mixsd<=max(grange)))),
                   c(rep(mode,sum(mode-mixsd>=min(grange))),(mode+mixsd)[mode+mixsd<=max(grange)]))
        prior = c(prior[mode-mixsd>=min(grange)], prior[mode+mixsd<=max(grange)])
      } else { #define two sets of components, but don't duplicate null component
        null.comp=which.min(mixsd)
        tmp = (mode+mixsd)[-null.comp]
        pi = c(pi[mode-mixsd>=min(grange)],(pi[-null.comp])[tmp<=max(grange)])
        pi = pi/sum(pi)
        g = unimix(pi,
                   c((mode-mixsd)[mode-mixsd>=min(grange)],rep(mode,sum(tmp<=max(grange)))),
                   c(rep(mode,sum(mode-mixsd>=min(grange))),tmp[tmp<=max(grange)]))
        prior = c(prior[mode-mixsd>=min(grange)],(prior[-null.comp])[tmp<=max(grange)])
        #pi = c(pi,pi[-null.comp])
      }
    }

    #check that all prior are >=1 (as otherwise have problems with infinite penalty)
    if(!all(prior>=1)){
      stop("Error: prior must all be >=1")}

    ##3. Fitting the mixture
    pi_init = g$pi
    k=ncomp(g)

    # the likelihood matrix for the extreme group
    if (is.null(df)) {
      # normal likelihood
      if (mixcompdist == "normal") {
        # normal mixture prior
        sd.mat.2 = sqrt(outer(sebetahat2^2, g$sd^2, FUN = "+"))
        matrix_lik2 = apply(sd.mat.2, 2, FUN = function(x) {dnorm(betahat2, 0, x)})
      } else {
        # uniform mixture prior
        matrix_lik2 = matrix(nrow = n, ncol = k)
        for (i in 1:n) {
          beta.hat = betahat1[i]
          sebeta.hat = sebetahat1[i]
          for (j in 1:k) {
            a = g$a[j]
            b = g$b[j]
            matrix_lik2[i, j] = (pnorm(b, beta.hat, sebeta.hat) - pnorm(a, beta.hat, sebeta.hat)) / (b - a)
          }
        }
      }
    } else {
      # t likelihood
      t = qt(1 - pval.thresh / 2, df1)
      if (mixcompdist == "normal") {
        # normal mixture prior
        stop("Error: t likelihood and normal mixture prior not implemented")
      } else {
        # uniform mixture prior
        for (i in 1:m) {
          beta.hat = betahat1[i]
          sebeta.hat = sebetahat1[i]
          pbeta = function (beta) {
            pt(t[i] - beta / sebeta.hat, df[i]) - pt(-t[i] - beta / sebeta.hat, df[i])
          }
          for (j in 1:k) {
            a = g$a[j]
            b = g$b[j]
            matrix_lik1[i, j] = integrate(pbeta, a, b)$value / (b - a)
          }
        }
      }
    }

    # the likelihood matrix for the moderate group

    if (is.null(df)) {
      t = qnorm(1 - pval.thresh / 2)
      # normal likelihood
      if (mixcompdist == "normal") {
        # normal mixture prior
        sd.mat.1 = sqrt(outer(sebetahat1^2, g$sd^2, FUN = "+"))
        matrix_lik1 = apply(sd.mat.1, 2, FUN = function(x) {2 * pnorm(t * sebetahat1 / x) - 1})
      } else {
        # uniform mixture prior
        matrix_lik1 = matrix(nrow = m, ncol = k)
        for (i in 1:m) {
          beta.hat = betahat1[i]
          sebeta.hat = sebetahat1[i]
          for (j in 1:k) {
            a = g$a[j]
            b = g$b[j]
            pbetahat = function(beta.hat) {
              pbeta = (pnorm(b, beta.hat, sebeta.hat) - pnorm(a, beta.hat, sebeta.hat)) / (b - a)
              return(pbeta)
            }
            matrix_lik1[i, j] = integrate(pbetahat, -t * sebeta.hat, t * sebeta.hat)$value
          }
        }
      }
    } else {
      # t likelihood
      t = qt(1 - pval.thresh / 2, df1)
      if (mixcompdist == "normal") {
        # normal mixture prior
        stop("Error: t likelihood and normal mixture prior not implemented")
      } else {
        # uniform mixture prior
        for (i in 1:m) {
          beta.hat = betahat1[i]
          sebeta.hat = sebetahat1[i]
          pbeta = function (beta) {
            pt(t[i] - beta / sebeta.hat, df[i]) - pt(-t[i] - beta / sebeta.hat, df[i])
          }
          for (j in 1:k) {
            a = g$a[j]
            b = g$b[j]
            matrix_lik1[i, j] = integrate(pbeta, a, b)$value / (b - a)
          }
        }
      }
    }

    # combine the likelihood matrices for two groups
    matrix_lik = rbind(matrix_lik1, matrix_lik2)

    # the last of these conditions checks whether the gradient at the null is negative wrt pi0
    # to avoid running the optimization when the global null (pi0=1) is the optimal.
    if(max(prior[-1])>1 || min(gradient(matrix_lik)+prior[1]-1,na.rm=TRUE)<0){
      pihat.est = mixIP(matrix_lik, prior, pi_init)
    } else {
      pihat.est = c(1, rep(0, k-1))
    }
    g$pi=pihat.est
    # estimate pihat with IP
    # pihat.est = mixIP(matrix_lik = lik.mat, prior, pi_init)

    # fitted.g = normalmix(pi = pihat.est$pihat, mean = 0, sd = sd)
  }

  output = ash.workhorse(betahat, sebetahat, df = df, g = g, fixg = TRUE)
  return(output)
}

############################### METHODS FOR normalmix class ###########################

#' @title Constructor for normalmix class
#'
#' @description Creates an object of class normalmix (finite mixture
#'     of univariate normals)
#'
#' @details None
#'
#' @param pi vector of mixture proportions
#' @param mean vector of means
#' @param sd vector of standard deviations
#'
#' @return an object of class normalmix
#'
#' @export
#'
#' @examples normalmix(c(0.5,0.5),c(0,0),c(1,2))
#'

normalmix = function (pi, mean, sd) {
  structure(data.frame(pi, mean, sd), class = "normalmix")
}

unimix = function (pi, a, b) {
  structure(data.frame(pi, a, b), class = "unimix")
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

setprior = function (prior, k, nullweight, null.comp) {
  if (!is.numeric(prior)) {
    if (prior == "nullbiased") {
      prior = rep(1, k)
      prior[null.comp] = nullweight
    }
    else if (prior == "uniform") {
      prior = rep(1, k)
    }
    else if (prior == "unit") {
      prior = rep(1/k, k)
    }
  }
  if (length(prior) != k | !is.numeric(prior)) {
    stop("invalid prior specification")
  }
  return(prior)
}

initpi = function (k, n, null.comp, randomstart = FALSE) {
  if (randomstart) {
    pi = stats::rgamma(k, 1, 1)
  }
  else {
    if (k < n) {
      pi = rep(1, k)/n
      pi[null.comp] = (n - k + 1)/n
    }
    else {
      pi = rep(1, k)/k
    }
  }
  pi = normalize(pi)
  return(pi)
}

normal_lik = function () {
  list(name = "normal", const = TRUE, lcdfFUN = function(x) {
    stats::pnorm(x, log = TRUE)
  }, lpdfFUN = function(x) {
    stats::dnorm(x, log = TRUE)
  }, etruncFUN = function(a, b) {
    my_etruncnorm(a, b)
  }, e2truncFUN = function(a, b) {
    my_e2truncnorm(a, b)
  })
}

t_lik = function (df) {
  list(name = "t", const = (length(unique(df)) == 1), lcdfFUN = function(x) {
    stats::pt(x, df = df, log = TRUE)
  }, lpdfFUN = function(x) {
    stats::dt(x, df = df, log = TRUE)
  }, etruncFUN = function(a, b) {
    etrunct::e_trunct(a, b, df = df, r = 1)
  }, e2truncFUN = function(a, b) {
    etrunct::e_trunct(a, b, df = df, r = 2)
  })
}

set_data = function (betahat, sebetahat, lik = NULL, alpha = 0) {
  if (length(sebetahat) == 1L) {
    sebetahat = rep(sebetahat, length(betahat))
  }
  data = list(x = betahat/(sebetahat^alpha), s = sebetahat^(1 -
                                                              alpha), alpha = alpha, s_orig = sebetahat)
  if (is.null(lik)) {
    lik = normal_lik()
  }
  data$lik = lik
  return(data)
}

#' @title log_comp_dens_conv.normalmix
#' @description returns log-density of convolution of each component
#'     of a normal mixture with N(0,s^2) or s*t(v) at x. Note that
#'     convolution of two normals is normal, so it works that way
#' @inheritParams comp_dens_conv.normalmix
#' @return a k by n matrix
log_comp_dens_conv.normalmix = function(m,data){
  if(!is_normal(data$lik)){
    stop("Error: normal mixture for non-normal likelihood is not yet implemented")
  }
  sdmat = sqrt(outer(data$s^2,m$sd^2,"+")) #n by k matrix of standard deviations of convolutions
  return(t(stats::dnorm(outer(data$x,m$mean,FUN="-")/sdmat,log=TRUE) - log(sdmat)))
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
