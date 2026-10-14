# Training exercise only. Run from the repository root.
# Usage: Rscript docs/examples/validation_smoke.R
suppressPackageStartupMessages({
  library(magicaxis)
  source('periodograms.R')
  source('periodoframe.R')
  source('functions.R')
})
renew <- TRUE
Nsamp <- 1
Ncores <- 1
dest <- file.path('results', 'student_validation',
                  paste0('smoke_', format(Sys.time(), '%Y%m%d_%H%M%S')))
if (dir.exists(dest)) stop('Output directory already exists: ', dest)
dir.create(dest, recursive=TRUE)
capture.output(sessionInfo(), file=file.path(dest, 'sessionInfo.txt'))
writeLines(system('git rev-parse HEAD', intern=TRUE),
           file.path(dest, 'commit.txt'))

set.seed(1730)
t <- sort(runif(80, 0, 180))
dy <- rep(1, length(t))
noise <- rnorm(length(t), sd=dy)
Ptrue <- 17.3
signal <- 5*cos(2*pi*t/Ptrue + 0.7)
series <- list(null=noise, planet=noise+signal)
truth <- data.frame(case=c('null','planet'), nplanet=c(0,1),
                    P=c(NA,Ptrue), K=c(0,5), e=c(NA,0))
write.csv(truth, file.path(dest, 'training_truth.csv'), row.names=FALSE)

pp <- list(ns='RV', ofac=2, frange=c(1/60,1/2),
  per.type='BFP', per.target='training', sequence=FALSE,
  Nsig.max=1, per.type.seq='BFP', Niter=0,
  SigType='circular', Nh=1, Nmas=0, Nars=0,
  Inds=list(0), lnBF.min=5, noise.model='ARMA')
saveRDS(pp, file.path(dest, 'settings.rds'))

rows <- lapply(names(series), function(id) {
  d <- data.frame(Time=t, RV=series[[id]], eRV=dy)
  write.table(d, file.path(dest, paste0(id, '.dat')),
              row.names=FALSE, quote=FALSE)
  start <- proc.time()[['elapsed']]
  warnings <- character()
  fit <- tryCatch(withCallingHandlers(
    calc.1Dper(1, 'RV', pp, list(training=d), Ncores=Ncores),
    warning=function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart('muffleWarning')
    }), error=function(e) e)
  elapsed <- proc.time()[['elapsed']]-start
  writeLines(warnings, file.path(dest, paste0(id, '_warnings.txt')))
  if (inherits(fit, 'error')) {
    return(data.frame(case=id, status='error', accepted=NA,
      P=NA, lnBF=NA, seconds=elapsed, message=conditionMessage(fit)))
  }
  saveRDS(fit, file.path(dest, paste0(id, '_result.rds')))
  cmp <- fit$model.comp$RV
  write.csv(cmp, file.path(dest, paste0(id, '_comparison.csv')),
            row.names=FALSE)
  pdf(file.path(dest, paste0(id, '_periodogram.pdf')), 7, 4.5)
  plot1D.single(fit, 'periodogram', 'RV', 1)
  dev.off()
  data.frame(case=id, status='ok', accepted=cmp$accepted[1],
    P=cmp$P[1], lnBF=cmp$lnBF[1], seconds=elapsed,
    message=paste(warnings, collapse=' | '))
})
summary <- do.call(rbind, rows)
write.csv(summary, file.path(dest, 'summary.csv'), row.names=FALSE)
inputs <- file.path(dest, paste0(names(series), '.dat'))
write.csv(data.frame(file=basename(inputs), md5=unname(tools::md5sum(inputs))),
          file.path(dest, 'input_hashes.csv'), row.names=FALSE)
print(summary)
cat('Training outputs:', dest, '\n')
if (any(summary$status!='ok')) stop('Training run failed; inspect summary.csv')
hit <- subset(summary, case=='planet')
stopifnot(hit$accepted, abs(1/hit$P-1/Ptrue)<1/(2*diff(range(t))))
cat('Strong-signal recovery check passed. One null run is not FPR calibration.\n')
