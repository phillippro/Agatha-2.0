withProgress <- function(message=NULL,value=0,expr,...) eval(substitute(expr),parent.frame())
incProgress <- function(...) invisible(NULL)
suppressPackageStartupMessages({
  library(magicaxis)
  for(f in c('functions.R','periodograms.R','periodoframe.R','mcmc_func.R','sofa.R','orbit.R')) source(f)
})
set.seed(61995)
renew <- TRUE
Nsamp <- 1
outdir <- 'results/HD61995_20260919'
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
src <- '../data/combined/HD61995/HD61995_KECK.vels'
tab <- read.table(src,header=TRUE)
stopifnot(all(is.finite(as.matrix(tab[,1:3]))),all(tab$eRV>0),!anyDuplicated(tab$BJD))
tab <- tab[order(tab$BJD),]
file.copy(src,file.path(outdir,basename(src)),overwrite=TRUE)
capture.output(sessionInfo(),file=file.path(outdir,'sessionInfo.txt'))
writeLines(paste(names(tools::md5sum(c(src,'functions.R','periodoframe.R','mcmc_func.R'))),
                 tools::md5sum(c(src,'functions.R','periodoframe.R','mcmc_func.R'))),file.path(outdir,'input_code_md5.txt'))
span <- diff(range(tab$BJD)); fr <- c(1/span,1/1.1)
cat('N=',nrow(tab),'span=',span,'RV rms=',sd(tab$RV),'median error=',median(tab$eRV),'\n')
data <- list(HD61995_KECK=tab[,1:3])
d <- NULL
for(model in c('MA','W','AR')) {
  one <- diagnose.signals(data,names(data),noise.models=model,Nsig.max=3,
    ofac=2,frange=fr,Nh=1,Ncores=2,Nwin=5,proxies=FALSE,moving=FALSE,
    progress=function(msg) cat(format(Sys.time()),msg,'\n'))
  if(is.null(d)) d <- one else d$combined[[model]] <- one$combined[[model]]
  d$noise.models <- names(d$combined)
  saveRDS(d,file.path(outdir,'checkpoint.rds'))
}

saveRDS(d,file.path(outdir,'diagnosis_rv.rds'))
# Missing activity measurements use -1; never treat these as observations.
for(idx in names(tab)[4:7]) {
  keep <- is.finite(tab[[idx]]) & tab[[idx]]>0
  z <- tab[keep,c('BJD','RV','eRV',idx)]
  cat('Proxy',idx,'N=',nrow(z),'\n')
  pp <- list(ns=c('RV',idx),ofac=2,frange=fr,per.type='BFP',per.target='HD61995_KECK',
    sequence=FALSE,Nsig.max=1,per.type.seq='BFP',Niter=0,SigType='circular',Nh=1,
    Nmas=0,Nars=0,Inds=list(0),lnBF.min=5)
  d$proxies[[paste0('HD61995_',idx)]] <- calc.1Dper(1,idx,pp,list(HD61995_KECK=z),Ncores=2)
}
# Moving windows use the MA reference fit after subtracting other signals.
if(length(d$periods)>0) {
  mp <- list(ns=c('RV','Window Function'),ofac=2,frange=fr,per.type='BFP',
    per.target=names(data),sequence=FALSE,Nmas=0L,Nars=0L,Inds=list(0),
    Nsig.max=1L,per.type.seq='BFP',Niter=0,SigType='circular',Nh=1,
    lnBF.min=5,noise.model='ARMA',gp.Prot=NA,gp.tau=NA,files=NULL,
    Dt=signif(span/2,3),Nbin=5L,alpha=5,scale=TRUE,pmin.zoom=1/fr[2],
    pmax.zoom=1/fr[1],show.signal=TRUE,use.fit=TRUE,adaptive=TRUE)
  d$moving.par <- mp; d$moving <- list()
  fit <- d$combined$MA; ph <- fit$phase.list$RV
  for(k in seq_along(d$periods)) {
    cat(format(Sys.time()),'moving signal',k,'\n')
    fitk <- fit
    cols <- intersect(paste0('ysig_sig',setdiff(seq_along(d$periods),k)),colnames(ph))
    if(length(cols)) fitk$phase.list$RV[,'y_all'] <- ph[,'y_all']-rowSums(ph[,cols,drop=FALSE])
    d$moving[[k]] <- per2D.data('RV',mp,data,fit1D=fitk)
    saveRDS(d,file.path(outdir,'checkpoint.rds'))
  }
}
saveRDS(d,file.path(outdir,'diagnosis.rds'))
write.csv(diagnosis.table(d),file.path(outdir,'diagnosis.csv'),row.names=FALSE)
for(m in names(d$combined)) {
  x <- d$combined[[m]]
  write.csv(x$model.comp$RV,file.path(outdir,paste0(m,'_model_comparison.csv')),row.names=FALSE)
  write.csv(x$per.list$RV,file.path(outdir,paste0(m,'_periodograms.csv')),row.names=FALSE)
  write.csv(x$par.list$RV,file.path(outdir,paste0(m,'_parameters.csv')))
  write.csv(x$phase.list$RV,file.path(outdir,paste0(m,'_fit_data.csv')),row.names=FALSE)
  for(kind in c('phase','fit','residual')) {
    png(file.path(outdir,paste0(m,'_',kind,'.png')),width=1400,height=850,res=160)
    plot1D.single(x,kind,'RV',1)
    dev.off()
  }
}
png(file.path(outdir,'diagnostics.png'),width=2100,height=1900,res=160)
diagnosis.plot(d,ncol=3)
dev.off()
if(length(d$moving)) for(k in seq_along(d$moving)) {
  if(is.null(d$moving[[k]])) next
  png(file.path(outdir,paste0('moving_signal',k,'.png')),width=1600,height=1400,res=160)
  plotMP(d$moving[[k]],d$moving.par)
  dev.off()
}
print(diagnosis.table(d))
cat('DONE\n')
