# Use pspectrum function to compute long run variance
SE.pspectrum = function(x, ..., verbose = FALSE, plot = FALSE){
  n = length(x)
  res.pspectrum = pspectrum(x, verbose = verbose, plot = plot)
  return(res.pspectrum$spec[1]/n/2)
}
