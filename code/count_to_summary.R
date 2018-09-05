## counts: ngene * nsample RNA-seq count matrix
## design: design matrix for group identities of samples

count_to_summary = function (counts, design) {
  dgecounts = edgeR::calcNormFactors(edgeR::DGEList(counts = counts, group = design[, 2]))
  v = limma::voom(dgecounts, design, plot = FALSE)
  lim = limma::lmFit(v)
  r.ebayes = limma::eBayes(lim)
  p = r.ebayes$p.value[, 2]
  t = r.ebayes$t[, 2]
  z = -sign(t) * qnorm(p/2)
  betahat = lim$coefficients[,2]
  sebetahat = betahat / z
  return (list(betahat = betahat, sebetahat = sebetahat, z = z))
}
