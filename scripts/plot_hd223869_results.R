suppressPackageStartupMessages({library(magicaxis); source('functions.R')})
outdir <- 'results/HD223869_20260919'
d <- readRDS(file.path(outdir,'diagnosis.rds'))
for(m in names(d$combined)) {
  x <- d$combined[[m]]
  for(kind in c('phase','fit','residual')) {
    png(file.path(outdir,paste0(m,'_',kind,'.png')),width=1400,height=850,res=160)
    plot1D.single(x,kind,'RV',1)
    dev.off()
  }
}
for(k in seq_along(d$moving)) if(!is.null(d$moving[[k]])) {
  png(file.path(outdir,paste0('moving_signal',k,'.png')),width=1600,height=1400,res=160)
  plotMP(d$moving[[k]],d$moving.par)
  dev.off()
}
png(file.path(outdir,'diagnostics.png'),width=2100,height=1900,res=160)
diagnosis.plot(d,ncol=3)
dev.off()
