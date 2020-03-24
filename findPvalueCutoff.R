# This function takes a list of P-values, and calculates what the first P-value
# above a certain FDR threshold would have had to be in order to have an FDR
# less than said threshold. Uses an iterative approach to approximate, but can
# be made arbitraritly precise.
# Designed to be used when generating volcano-plots

findPvalueCutoff <- function(x, cutoff = 0.05)
{
  if (any(x < 0)) stop("Negative P-values not allowed")
  pvals <- sort(x)
  fdr <- p.adjust(pvals, method = "fdr")
  
  if (!any(fdr < cutoff)){
    return(cutoff/(length(pvals)))
  }
  
  firstFailing <- max(which(fdr < cutoff)) + 1
  
  diffFn <- function(x)
  {
    p.adjust(c(x, pvals[-(firstFailing)]), method = "fdr")[[1]] - cutoff
  }
  optimize(diffFn, pvals[c(firstFailing - 1, firstFailing)])$minimum
}
