## selection top expressed genes
top_genes_index = function (ngene, data) {return(order(rowSums(data), decreasing = TRUE)[1 : ngene])}

lcpm = function(data) {R = colSums(data); t(log2(((t(data) + 0.5) / (R + 1)) * 10^6))}

top_gene_selection = function (ngene, data) {
  Y = lcpm(data)
  subset = top_genes_index(ngene, data)
  data = data[subset, ]
  return(list(data = data, ngene = ngene))
}

# this script generates a 1k * 10k matrix with p-values
# of correlatednull data simulated from GTEx/Liver rna-seq
# each row is a random simulation run of n genes
#
# then record the number of extreme observations
# that is, observations with p-values <= T (pre-specified)
#


# Load in the gtex liver data

library(limma)
library(edgeR)
library(qvalue)
library(ashr)
r = read.csv("../data/liver.csv")
r = r[,-(1:2)] # remove gene name and description
#extract top g genes from G by n matrix X of expression
top_genes_index=function(g,X){return(order(rowSums(X),decreasing =TRUE)[1:g])}
lcpm = function(r){R = colSums(r); t(log2(((t(r)+0.5)/(R+1))* 10^6))}
Y=lcpm(r)
subset = top_genes_index(10000,Y)
Y = Y[subset,]
r = r[subset,]

# Define voom transform (using code from Mengyin Lu)

voom_transform = function(counts, condition, W=NULL){
  dgecounts = calcNormFactors(DGEList(counts=counts,group=condition))
  #dgecounts = DGEList(counts=counts,group=condition)
  if (is.null(W)){
    design = model.matrix(~condition)
  }else{
    design = model.matrix(~condition+W)
  }

  v = voom(dgecounts,design,plot=FALSE)
  lim = lmFit(v)
  betahat.voom = lim$coefficients[,2]
  sebetahat.voom = lim$stdev.unscaled[,2]*lim$sigma
  df.voom = length(condition)-2-!is.null(W)

  return(list(v=v,lim=lim,betahat=betahat.voom, sebetahat=sebetahat.voom, df=df.voom, v=v))
}

# Make 2 groups of size n, and repeat this random sampling for m times
# to generate a m * 10k matrix of p values, each row is a random sampling of 10k null genes.


set.seed(101)
n = 5 # number in each group
m = 1000
p = matrix(nrow = m, ncol = 10000)

for(i in 1:m){
  counts = r[,sample(1:ncol(r),2*n)]
  condition = c(rep(0,n),rep(1,n))
  r.voom = voom_transform(counts,condition)
  r.ebayes = eBayes(r.voom$lim)
  p[i, ] = r.ebayes$p.value[,2]
}

write.table(p, file = "../output/p_null_liver.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)



# record the number of extreme observations
# with 22 thresholds

p = read.table("../output/p_null_liver.txt")
p = read.table(gzcon("../output/p_null_liver.tar.gz"))


extreme_p = matrix(nrow = dim(p)[1], ncol = 22)
for (i in 1:dim(p)[1]) {
  extreme_p[i, 1] = sum(p[i, ] <= 5e-4)
  extreme_p[i, 2] = sum(p[i, ] <= 0.001)
  extreme_p[i, 3] = sum(p[i, ] <= 0.005)
  extreme_p[i, 4] = sum(p[i, ] <= 0.01)
  extreme_p[i, 5] = sum(p[i, ] <= 0.02)
  extreme_p[i, 6] = sum(p[i, ] <= 0.03)
  extreme_p[i, 7] = sum(p[i, ] <= 0.04)
  extreme_p[i, 8] = sum(p[i, ] <= 0.05)
  extreme_p[i, 9] = sum(p[i, ] <= 0.06)
  extreme_p[i, 10] = sum(p[i, ] <= 0.07)
  extreme_p[i, 11] = sum(p[i, ] <= 0.08)
  extreme_p[i, 12] = sum(p[i, ] <= 0.09)
  extreme_p[i, 13] = sum(p[i, ] <= 0.10)
  extreme_p[i, 14] = sum(p[i, ] <= 0.15)
  extreme_p[i, 15] = sum(p[i, ] <= 0.20)
  extreme_p[i, 16] = sum(p[i, ] <= 0.25)
  extreme_p[i, 17] = sum(p[i, ] <= 0.30)
  extreme_p[i, 18] = sum(p[i, ] <= 0.35)
  extreme_p[i, 19] = sum(p[i, ] <= 0.40)
  extreme_p[i, 20] = sum(p[i, ] <= 0.45)
  extreme_p[i, 21] = sum(p[i, ] <= 0.50)
  extreme_p[i, 22] = sum(p[i, ] <= 0.75)
}

names(extreme_p) = c(5e-4, 0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.75)
write.table(extreme_p, file = "../output/p_null_liver_extreme.txt", quote = FALSE, col.names = TRUE, row.names = FALSE)
