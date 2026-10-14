# Run from the repository root: Rscript docs/examples/quickstart.R
suppressPackageStartupMessages({
  library(magicaxis)
  source('periodograms.R')
  source('periodoframe.R')
  source('functions.R')
})
renew <- TRUE
Nsamp <- 1
set.seed(20260930)
tab <- read.table('data/HD210193_PFS.vels', header=TRUE)[,1:3]
tab <- tab[order(tab[,1]),]
stopifnot(all(is.finite(as.matrix(tab))), all(tab[,3]>0))
data <- list(HD210193_PFS=tab)
pp <- list(ns='RV', ofac=0.5, frange=c(1/3000,1/2),
  per.type='BFP', per.target=names(data),
  sequence=FALSE, Nsig.max=1, per.type.seq='BFP',
  Niter=0, SigType='circular', Nh=1,
  Nmas=1, Nars=0, Inds=list(0), lnBF.min=5)
out <- calc.1Dper(1, 'RV', pp, data, Ncores=1)
dest <- 'results/documentation_quickstart'
dir.create(dest, recursive=TRUE, showWarnings=FALSE)
saveRDS(out, file.path(dest,'result.rds'))
write.csv(out$per.list$RV, file.path(dest,'periodogram.csv'),
          row.names=FALSE)
write.csv(par.table(out$par.list$RV),
          file.path(dest,'parameters.csv'), row.names=FALSE)
pdf(file.path(dest,'periodogram.pdf'), width=7, height=4.5)
plot1D.single(out, 'periodogram', 'RV', 1)
dev.off()
capture.output(sessionInfo(), file=file.path(dest,'sessionInfo.txt'))
print(out$model.comp$RV)
cat('Quickstart completed:', dest, '\n')
