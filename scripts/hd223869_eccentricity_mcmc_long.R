withProgress <- function(message=NULL,value=0,expr,...) eval(substitute(expr),parent.frame())
incProgress <- function(...) invisible(NULL)
suppressPackageStartupMessages({
  library(magicaxis)
  for(f in c('functions.R','periodograms.R','periodoframe.R','mcmc_func.R','sofa.R','orbit.R')) source(f)
})
outdir <- 'results/HD223869_20260919/eccentricity_mcmc_long'
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
tab <- read.table('results/HD223869_20260919/HD223869_KECK.vels',header=TRUE)[,1:3]
x <- readRDS('results/HD223869_20260919/kepler_refinement.rds')
per <- list(par.opt=x$par.list$RV,Popt=unname(x$par.list$RV['P1']))
fits <- parallel::mclapply(1:2,function(i){
  set.seed(223969+i)
  sink(file.path(outdir,paste0('chain',i,'.log')))
  on.exit(sink())
  registerDoMC(1)
  fit <- mcfit(per,tab,tsim=seq(0,diff(range(tab$BJD)),length.out=1000),
    Niter=40000,SigType='kepler',Ncores=1,mcmc.verbose=TRUE,
    e.prior=list(type='beta',a=0.867,b=3.03))
  saveRDS(fit,file.path(outdir,paste0('chain',i,'.rds')))
  fit
},mc.cores=2)
stopifnot(all(vapply(fits,function(f)is.list(f)&&!is.null(f$mc),logical(1))))
mc <- do.call(rbind,lapply(fits,`[[`,'mc'))
write.csv(mc,file.path(outdir,'posterior.csv'),row.names=FALSE)
q <- quantile(mc[,'e1'],c(.025,.16,.5,.84,.95,.975))
print(q)
write.csv(data.frame(probability=names(q),eccentricity=as.numeric(q)),file.path(outdir,'eccentricity_quantiles.csv'),row.names=FALSE)
saveRDS(fits,file.path(outdir,'fits.rds'))
png(file.path(outdir,'eccentricity_posterior.png'),width=1400,height=800,res=160)
par(mfrow=c(1,2),mar=c(4,4,2,1))
hist(mc[,'e1'],breaks=50,freq=FALSE,xlab='Eccentricity',main='576-day signal: posterior',col='lightblue',border='white')
abline(v=q[c(2,3,4)],lty=c(2,1,2),col='firebrick')
plot(fits[[1]]$mc[,'e1'],type='l',col='#0072B280',xlab='Post-burn-in draw',ylab='Eccentricity',main='Independent cold chains',ylim=c(0,1))
lines(fits[[2]]$mc[,'e1'],col='#D55E0080')
dev.off()
cat('DONE\n')
