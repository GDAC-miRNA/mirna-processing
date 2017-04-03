library(bigrquery)
library(NMF)

#' This function clusters samples based on their miRNA signature
#'
#' @export
mirnaCluster = function(df)
{
  #dataframe must contain miRNAs along rows, samples along columns
  #TODO: how to best choose rank
  nmfFit = nmf(df[1:20,1:30],rank=10,method = "brunet",seed=123456)
  showRNG(nmfFit)
  basismap(nmfFit)
}
