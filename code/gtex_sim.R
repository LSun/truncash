require(limma)
require(edgeR)
require(qvalue)
require(ashr)

# Load in the gtex liver data

r = read.csv("../data/liver.csv")
r = r[, -(1 : 2)] # remove gene name and description

#extract top g genes from G by n matrix X of expression

top_genes_index = function (g, X)
{return(order(rowSums(X), decreasing = TRUE)[1 : g])
}

lcpm = function (r) {
  R = colSums(r)
  t(log2(((t(r) + 0.5) / (R + 1)) * 10^6))
}

Y = lcpm(r)
subset = top_genes_index(10000, Y)
Y = Y[subset,]
r = r[subset,]

# transform counts to z scores
# these z scores are marginally N(0, 1) under null

counts_to_z = function (counts, condition) {
  design = model.matrix(~condition)
  dgecounts = calcNormFactors(DGEList(counts = counts, group = condition))
  v = voom(dgecounts, design, plot = FALSE)
  lim = lmFit(v)
  r.ebayes = eBayes(lim)
  p = r.ebayes$p.value[, 2]
  t = r.ebayes$t[, 2]
  z = sign(t) * qnorm(1 - p/2)
  return (z)
}


# Define voom transform (using code from Mengyin Lu)

voom_transform = function(counts, condition){
  dgecounts = calcNormFactors(DGEList(counts = counts, group = condition))
  design = model.matrix(~condition)
  v = voom(dgecounts, design, plot = FALSE)
  lim = lmFit(v)
  betahat.voom = lim$coefficients[, 2]
  sebetahat.voom = lim$stdev.unscaled[,2] * lim$sigma
  df.voom = length(condition) - 2
  return(list(v = v, lim = lim, betahat = betahat.voom, sebetahat = sebetahat.voom, df = df.voom, v = v))
}

# Make 2 groups of size n, and repeat this random sampling for m times
# to generate a m * 10k matrix of p values, each row is a random sampling of 10k null genes.

set.seed(777)
n = 5 # number in each group
m = 1000
betahat = sebetahat = p = tscore = z = matrix(nrow = m, ncol = 10000)

for(i in 1:m){
  counts = r[,sample(1:ncol(r),2*n)]
  condition = c(rep(0,n),rep(1,n))
  r.voom = voom_transform(counts,condition)
  r.ebayes = eBayes(r.voom$lim)
  p[i, ] = r.ebayes$p.value[,2]
  tscore[i, ] = r.ebayes$t[,2]
  z[i, ] = sign(tscore[i, ]) * qnorm(1 - p[i, ]/2)
  betahat[i, ] = r.voom$betahat
  sebetahat[i, ] = r.voom$sebetahat
}

write.table(p, file = "../output/p_null_liver_777.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(tscore, file = "../output/t_null_liver_777.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(z, file = "../output/z_null_liver_777.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(betahat, file = "../output/betahat_null_liver_777.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(sebetahat, file = "../output/sebetahat_null_liver_777.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)


dat = list()
z = read.table("../output/z_null_liver_777.txt")
z1 = as.numeric(z)

m = dim(z)[1]
n = dim(z)[2]

num = matrix(0, nrow = m, ncol = n)
for (i in 1:m) {
  for (j in 1:n) {
    num[i, j] = is.numeric(z[i, j])
  }
}

res_ash = res_truncash = list()

for (i in 1:m) {
  betahat = as.numeric(z[i, ])
  sebetahat = rep(1, n)
  res_truncash[[i]] = truncash(betahat, sebetahat, t = qnorm(0.975))
  res_ash[[i]] = ash(betahat, sebetahat, mixcompdist = "normal", method = "fdr")
}

pihat0_ash = pihat0_truncash = c()
for (i in 1:m) {
  pihat0_ash[i] = get_fitted_g(res_ash[[i]])$pi[1]
  pihat0_truncash[i] = get_fitted_g(res_truncash[[i]])$pi[1]
}
