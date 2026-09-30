withProgress <- function(message=NULL,value=0,expr,...) eval(substitute(expr),parent.frame())
incProgress <- function(...) invisible(NULL)
suppressPackageStartupMessages({
  library(magicaxis)
  for(f in c('functions.R','periodograms.R','periodoframe.R','mcmc_func.R','sofa.R','orbit.R')) source(f)
})
set.seed(61995)
renew <- TRUE; Nsamp <- 1
outdir <- 'results/HD61995_20260919'
tab <- read.table(file.path(outdir,'HD61995_KECK.vels'),header=TRUE)[,1:3]
# Conditional eccentric refinement of the previously identified long-period peak.
pp <- list(ns='RV',ofac=1,frange=c(1/600,1/350),per.type='BFP',per.target='HD61995_KECK',
  sequence=FALSE,Nsig.max=1,per.type.seq='BFP',Niter=0,SigType='kepler',Nh=2,
  Nmas=1,Nars=0,Inds=list(0),lnBF.min=5)
x <- calc.1Dper(1,'RV',pp,list(HD61995_KECK=tab),Ncores=2)
saveRDS(x,file.path(outdir,'kepler_refinement.rds'))
write.csv(x$par.list$RV,file.path(outdir,'kepler_parameters.csv'))
write.csv(x$phase.list$RV,file.path(outdir,'kepler_fit_data.csv'),row.names=FALSE)
for(kind in c('phase','fit','residual')) {
  png(file.path(outdir,paste0('kepler_',kind,'.png')),width=1400,height=850,res=160)
  plot1D.single(x,kind,'RV',1,SigType='kepler')
  dev.off()
}
print(x$par.list$RV)
cat('DONE\n')
