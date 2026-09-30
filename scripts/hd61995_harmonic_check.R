withProgress <- function(message=NULL,value=0,expr,...) eval(substitute(expr),parent.frame())
incProgress <- function(...) invisible(NULL)
suppressPackageStartupMessages({
  library(magicaxis)
  for(f in c('functions.R','periodograms.R','periodoframe.R','mcmc_func.R','sofa.R','orbit.R')) source(f)
})
set.seed(619952)
registerDoMC(2)
outdir <- 'results/HD61995_20260919/harmonic_check'
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
tab <- read.table('results/HD61995_20260919/HD61995_KECK.vels',header=TRUE)[,1:3]
old <- readRDS('results/HD61995_20260919/diagnosis.rds')
span <- diff(range(tab$BJD)); runs <- list(); rows <- list(); fits <- list()
for(m in c('W','MA','AR')) {
  cat(format(Sys.time()),'two-harmonic BFP',m,'\n')
  z <- BFP(tab$BJD,tab$RV,tab$eRV,Nma=as.integer(m=='MA'),Nar=as.integer(m=='AR'),
    Nh=2,ofac=2,fmin=1/span,fmax=1/1.1,progress=FALSE,renew=TRUE,Nsamp=1)
  runs[[m]] <- z
  saveRDS(z,file.path(outdir,paste0(m,'_Nh2.rds')))
  write.csv(data.frame(P=z$P,lnBF=z$power),file.path(outdir,paste0(m,'_Nh2.csv')),row.names=FALSE)
  for(region in c('200d','469d','global')) {
    ii <- switch(region,'200d'=which(z$P>=150&z$P<=250),'469d'=which(z$P>=350&z$P<=600),'global'=seq_along(z$P))
    i <- ii[which.max(z$power[ii])]
    rows[[length(rows)+1]] <- data.frame(model=m,region=region,P=z$P[i],lnBF=z$power[i])
    if(region!='global') {
      p <- z; p$Popt <- z$P[i]; p$par.opt <- z$pars[i,]
      k <- KeplerFit(p,basis='natural')
      fits[[paste(m,region,sep='_')]] <- k
      write.csv(as.data.frame(as.list(k$ParKep)),file.path(outdir,paste0(m,'_',region,'_kepler.csv')),row.names=FALSE)
      cat(m,region,'Kepler P=',k$ParKep$P1,'e=',k$ParKep$e1,'logL=',k$LogLike,'\n')
    }
  }
  write.csv(do.call(rbind,rows),file.path(outdir,'peak_summary.csv'),row.names=FALSE)
  saveRDS(fits,file.path(outdir,'kepler_fits.rds'))
}
saveRDS(runs,file.path(outdir,'harmonic_runs.rds'))
# Targeted residual test: subtract only the 469-day Keplerian signal, then
# re-fit offset, trend, jitter and the chosen AR/MA noise at each trial period.
resrows <- list()
for(m in names(runs)) {
  k <- fits[[paste(m,'469d',sep='_')]]
  sig <- KeplerSig(k$ParKep,k$df,basis='natural')$ysig
  cat(format(Sys.time()),'200d residual search',m,'\n')
  z <- BFP(tab$BJD,tab$RV-sig,tab$eRV,Nma=as.integer(m=='MA'),Nar=as.integer(m=='AR'),
    Nh=2,ofac=.5,fmin=1/250,fmax=1/150,progress=FALSE,renew=TRUE,Nsamp=1)
  saveRDS(z,file.path(outdir,paste0(m,'_200d_residual.rds')))
  i <- which.max(z$power)
  resrows[[m]] <- data.frame(model=m,P=z$P[i],lnBF=z$power[i])
}
write.csv(do.call(rbind,resrows),file.path(outdir,'residual_summary.csv'),row.names=FALSE)
png(file.path(outdir,'harmonic_comparison.png'),width=1800,height=1500,res=160)
par(mfrow=c(3,2),mar=c(4,4,2,1))
for(m in names(runs)) {
  z <- runs[[m]]; circ <- old$combined[[m]]$per.list$RV
  for(xr in list(c(100,1200),c(150,250))) {
    i<-order(z$P); j<-order(circ[,1]); a<-z$P>=xr[1]&z$P<=xr[2]; b<-circ[,1]>=xr[1]&circ[,1]<=xr[2]
    plot(z$P[i],z$power[i],type='l',log='x',xlim=xr,ylim=range(0,5,z$power[a],circ[b,2]),
      xlab='Period [days]',ylab='ln(BF), BIC approximation',main=paste(m,'noise'),col='#D55E00')
    lines(circ[j,1],circ[j,2],col='#0072B2')
    abline(h=5,lty=2,col='grey50'); abline(v=c(200,469.4),lty=3,col='grey40')
    legend('topright',c('Two harmonics','One harmonic'),col=c('#D55E00','#0072B2'),lty=1,bty='n',cex=.8)
  }
}
dev.off()
print(do.call(rbind,rows)); print(do.call(rbind,resrows)); cat('DONE\n')
