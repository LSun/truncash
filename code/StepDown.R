# This script applies two step-down multiple comparison procedures to
# simulated, correlated global null data
# Holm-Bonferroni & Sidak-Holm


Sidak <- function(vecP)
  #
  # This function corrects a vector of probabilities for multiple testing
  # using the Bonferroni (1935) and Sidak (1967) corrections.
  #
  # References: Bonferroni (1935), Sidak (1967), Wright (1992).
  #
  # Bonferroni, C. E. 1935. Il calcolo delle assicurazioni su gruppi di teste.
  # Pp. 13-60 in: Studi in onore del Professore Salvatore Ortu Carboni. Roma.
  #
  # Sidak, Z. 1967. Rectangular confidence regions for the means of multivariate
  # normal distributions. Journal of the American Statistical Association 62:626-633.
#
# Wright, S. P. 1992. Adjusted P-values for simultaneous inference.
# Biometrics 48: 1005-1013.
#
#                  Pierre Legendre, May 2007
{
  k = length(vecP)

  vecPB = 0
  vecPS = 0

  for(i in 1:k) {
    bonf = vecP[i]*k
    if(bonf > 1) bonf=1
    vecPB = c(vecPB, bonf)
    vecPS = c(vecPS, (1-(1-vecP[i])^k))
  }
  #
  return(list(OriginalP=vecP, BonfP=vecPB[-1], SidakP=vecPS[-1]))
}


p = read.table("../output/p_null_liver.txt")

holm = sidak = matrix(nrow = nrow(p), ncol = ncol(p))
for (i in 1:dim(p)[1]) {
  holm[i, ] = p.adjust(p[i, ], method = "holm")
  sidak[i, ] = Sidak(p[i, ])$SidakP
}

write.table(holm, file = "../output/p_null_liver_holm.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(sidak, file = "../output/p_null_liver_sidak.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
