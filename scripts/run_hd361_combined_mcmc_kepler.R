withProgress <- function(message=NULL, value=0, expr, ...) {
  eval(substitute(expr), parent.frame())
}
incProgress <- function(amount=0, detail=NULL, ...) invisible(NULL)

suppressPackageStartupMessages({
  library(magicaxis)
  source("functions.R")
  source("periodoframe.R")
  source("periodograms.R")
  source("mcmc_func.R")
  source("sofa.R")
  source("orbit.R")
})

input_dir <- "/Users/ffeng/Documents/projects/dwarfs/agatha2/data/HD361"
out_dir <- file.path("results", "HD361_combined_mcmc_kepler")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

files <- sort(Sys.glob(file.path(input_dir, "*")))
dataset <- tools::file_path_sans_ext(basename(files))
tabs <- lapply(files, read.agatha.table)
names(tabs) <- dataset

rows <- do.call(rbind, Map(function(tab, id) {
  data.frame(t=tab[[1]], y=tab[[2]], dy=tab[[3]], set.id=id)
}, tabs, dataset))
rows <- rows[order(rows$t),]
rows$set.id <- factor(rows$set.id, levels=dataset)

ofac <- 1
fmin <- 1 / 1e5
fmax <- 1 / 1.1
Niter <- 1000
Nma <- 1
Nar <- 0

log_con <- file(file.path(out_dir, "run.log"), open="wt")
sink(log_con)
sink(log_con, type="message")

cat("Combined HD361 Keplerian Agatha run\n")
cat("datasets:", paste(dataset, collapse=", "), "\n")
cat("Ndata:", nrow(rows), "Niter:", Niter, "Nma:", Nma, "Nar:", Nar, "\n")

per <- BFP.multiset(t=rows$t, y=rows$y, dy=rows$dy, set.id=rows$set.id,
                    Nma=Nma, Nar=Nar, ofac=ofac, fmin=fmin, fmax=fmax)
per_tab <- data.frame(Period.day=per$P, logBF=per$power)
per_tab <- per_tab[order(per_tab$Period.day),]
write.table(per_tab, file.path(out_dir, "HD361_combined_BFP_MA1_periodogram.txt"),
            quote=FALSE, row.names=FALSE)

pdf(file.path(out_dir, "HD361_combined_BFP_MA1_periodogram.pdf"), width=7, height=5)
plot(per_tab$Period.day, per_tab$logBF, type="l", log="x", xaxt="n",
     xlab="Period [day]", ylab="ln(BF)",
     main="HD361 combined BFP, MA(1), shared signal")
magaxis(side=1, tcl=-0.5)
abline(v=per$Popt, col="red", lwd=2)
abline(h=5, lty=2, col="grey50")
legend("topright", legend=paste0("P=", signif(per$Popt, 6), " d"),
       bty="n", text.col="red")
dev.off()

## Set up the lower-level Agatha MCMC globals for a multi-instrument Keplerian fit.
target <- "HD361"
prior.type <- "mt"
period.par <- "logP"
time.unit <- 365.25
offset <- TRUE
ins <- dataset
tmin <- min(rows$t)
trv.all <- rows$t
out <- list()
out$trv.all <- trv.all
out$ins <- ins
out$prior.type <- prior.type

start_index <- 1
for (k in seq_along(dataset)) {
  tab <- tabs[[dataset[k]]]
  n <- nrow(tab)
  out[[dataset[k]]] <- list()
  out[[dataset[k]]]$RV <- as.matrix(tab[, 1:3])
  out[[dataset[k]]]$index <- start_index:(start_index + n - 1)
  out[[dataset[k]]]$noise <- list(nqp=c(0, Nma, Nar))
  start_index <- start_index + n
}
## Match trv.all to the concatenated per-instrument indexes used above.
trv.all <- unlist(lapply(tabs, function(tab) tab[[1]]), use.names=FALSE)
out$trv.all <- trv.all
tmin <- min(trv.all)

A <- unname(per$par.opt["A"])
B <- unname(per$par.opt["B"])
K0 <- max(sqrt(A^2 + B^2), 1e-3)
phase0 <- as.numeric(xy2phi(A, B)) %% (2 * pi)
trend0 <- unname(per$par.opt["beta"]) * time.unit

startvalue <- c(log(per$Popt), K0, 0.1, 0, phase0, trend0)
names(startvalue) <- c("per1", "K1", "e1", "omega1", "Mo1", "a11")
for (k in seq_along(dataset)) {
  gamma_name <- paste0("gamma_", make.names(dataset[k]))
  ma_name <- paste0("ma_", make.names(dataset[k]), "_lag1")
  b0 <- unname(per$par.opt[gamma_name])
  s0 <- max(0.1, 0.5 * median(tabs[[dataset[k]]][[3]], na.rm=TRUE))
  w0 <- max(-0.9, min(0.9, unname(per$par.opt[ma_name])))
  beta0 <- log(max(1, median(diff(sort(tabs[[dataset[k]]][[1]])), na.rm=TRUE)))
  tmp <- c(b0, s0, w0, beta0)
  names(tmp) <- c(paste0("b", k), paste0("s", k), paste0("w", k, "1"), paste0("beta", k))
  startvalue <- c(startvalue, tmp)
}

span <- max(trv.all) - min(trv.all)
Dt <- span / time.unit
Esd <- 0.1
phi.min <- wmin <- -1
phi.max <- wmax <- 1
alpha.min <- beta.min <- log(1 / 24)
alpha.max <- beta.max <- log(span)
tol <- 1e-16
tol1 <- 1e-12
Sd <- 2.4^2 / length(startvalue)
Npar <- length(startvalue)

par.min <- startvalue - abs(startvalue)
par.max <- startvalue + abs(startvalue)
par.min["per1"] <- log(max(1.1, per$Popt * 0.5))
par.max["per1"] <- log(min(1e5, per$Popt * 2))
par.min["K1"] <- 0
par.max["K1"] <- max(2 * K0, K0 + 20)
par.min["e1"] <- 0
par.max["e1"] <- 0.95
par.min[c("omega1", "Mo1")] <- 0
par.max[c("omega1", "Mo1")] <- 2 * pi
par.min["a11"] <- -10 * sd(rows$y) / Dt
par.max["a11"] <- 10 * sd(rows$y) / Dt
for (k in seq_along(dataset)) {
  rv_sd <- sd(tabs[[dataset[k]]][[2]], na.rm=TRUE)
  par.min[paste0("b", k)] <- min(tabs[[dataset[k]]][[2]], na.rm=TRUE) - 5 * rv_sd
  par.max[paste0("b", k)] <- max(tabs[[dataset[k]]][[2]], na.rm=TRUE) + 5 * rv_sd
  par.min[paste0("s", k)] <- 0
  par.max[paste0("s", k)] <- max(1, 2 * rv_sd)
  par.min[paste0("w", k, "1")] <- wmin
  par.max[paste0("w", k, "1")] <- wmax
  par.min[paste0("beta", k)] <- beta.min
  par.max[paste0("beta", k)] <- beta.max
}

cat("BFP Popt:", per$Popt, "logBF:", max(per$power, na.rm=TRUE), "\n")
cat("MCMC start:\n")
print(startvalue)

opt <- optim(startvalue, function(p) -loglikelihood(p, bases=rep("natural", 10)),
             method="L-BFGS-B", lower=par.min, upper=par.max,
             control=list(maxit=500))
startvalue <- opt$par
cat("Optimized start negative logL:", opt$value, "\n")
print(startvalue)

cov.start <- diag(length(startvalue)) * 1e-5
diag(cov.start)[names(startvalue) == "per1"] <- 1e-6
diag(cov.start)[names(startvalue) == "K1"] <- 1e-3
diag(cov.start)[names(startvalue) == "e1"] <- 1e-5
diag(cov.start)[grepl("^omega|^Mo", names(startvalue))] <- 1e-3

set.seed(361)
mc_out <- run.metropolis.MCMC(startvalue=startvalue, cov.start=cov.start,
                              iterations=Niter, tem=1, bases=rep("natural", 10))
mc <- mc_out$out
write.table(mc, file.path(out_dir, "HD361_combined_kepler_MA1_MCposterior.txt"),
            quote=FALSE, row.names=FALSE)

ll <- mc[, "loglike"]
par_cols <- setdiff(colnames(mc), c("logpost", "loglike"))
mc_par <- mc[, par_cols, drop=FALSE]
mc_export <- mc_par
mc_export[, "per1"] <- exp(mc_export[, "per1"])
colnames(mc_export)[colnames(mc_export) == "per1"] <- "P1"
par_stat <- t(apply(mc_export, 2, function(x) {
  c(xopt=x[which.max(ll)], x1per=quantile(x, 0.01), x99per=quantile(x, 0.99),
    x10per=quantile(x, 0.10), x90per=quantile(x, 0.90),
    med=median(x), mean=mean(x), sd=sd(x))
}))
write.csv(par_stat, file.path(out_dir, "HD361_combined_kepler_MA1_OptPar.csv"),
          quote=FALSE)

best <- mc[which.max(ll), par_cols]
rv_all <- RV.kepler(best, bases=rep("natural", 10))
rv_sig <- RV.kepler(best, kep.only=TRUE, bases=rep("natural", 10))
fit_rows <- do.call(rbind, lapply(seq_along(dataset), function(k) {
  id <- dataset[k]
  tab <- tabs[[id]]
  model <- as.numeric(rv_all[[id]])
  signal <- as.numeric(rv_sig[[id]])
  nqp <- out[[id]]$noise$nqp
  ma <- if (nqp[2] > 0) {
    arma(t=tab[[1]], ymodel=model, ydata=tab[[2]], pars=best,
         ind.set=k, p=nqp[3], q=nqp[2])$ma
  } else {
    rep(0, nrow(tab))
  }
  data.frame(dataset=id, BJD=tab[[1]], RV=tab[[2]], eRV=tab[[3]],
             kepler_signal=signal, model=model + ma, ma=ma,
             residual=tab[[2]] - model - ma)
}))
write.table(fit_rows, file.path(out_dir, "HD361_combined_kepler_MA1_fit_residuals.txt"),
            quote=FALSE, row.names=FALSE)

P_best <- exp(best["per1"])
t_grid <- seq(min(trv.all), max(trv.all), length.out=2000)
sim_signal <- RV.kepler(best, tt=t_grid, kep.only=TRUE, bases=rep("natural", 10))
sim <- data.frame(BJD=t_grid, phase=(t_grid - tmin) %% P_best, kepler_signal=as.numeric(sim_signal))
write.table(sim, file.path(out_dir, "HD361_combined_kepler_MA1_sim_signal.txt"),
            quote=FALSE, row.names=FALSE)

pdf(file.path(out_dir, "HD361_combined_kepler_MA1_phase_fit.pdf"), width=7, height=7)
layout(matrix(c(1, 2), nrow=2), heights=c(2, 1))
phase <- (fit_rows$BJD - tmin) %% P_best
sim_phase <- sim$phase
ord <- order(sim_phase)
plot(phase, fit_rows$RV - (fit_rows$model - fit_rows$kepler_signal),
     xlab="", ylab="RV minus offsets/noise [m/s]", xaxt="n",
     main=paste0("HD361 combined Keplerian fit, P=", signif(P_best, 6), " d"))
axis(side=1, labels=FALSE)
arrows(phase, fit_rows$RV - (fit_rows$model - fit_rows$kepler_signal) - fit_rows$eRV,
       phase, fit_rows$RV - (fit_rows$model - fit_rows$kepler_signal) + fit_rows$eRV,
       length=0.03, angle=90, code=3, col="grey60")
points(phase, fit_rows$RV - (fit_rows$model - fit_rows$kepler_signal),
       pch=19, cex=0.65, col=as.integer(factor(fit_rows$dataset)))
lines(sim_phase[ord], sim$kepler_signal[ord], col="red", lwd=2)
legend("topright", legend=dataset, pch=19, col=seq_along(dataset), bty="n", cex=0.8)
plot(phase, fit_rows$residual, xlab="Orbital phase [day]", ylab="O-C [m/s]",
     pch=19, cex=0.65, col=as.integer(factor(fit_rows$dataset)))
abline(h=0, lty=2, col="grey50")
dev.off()

summary <- data.frame(
  dataset="HD361_combined",
  status="ok",
  n_data=nrow(rows),
  n_sets=length(dataset),
  bfp_peak_days=per$Popt,
  bfp_logBF=max(per$power, na.rm=TRUE),
  mcmc_best_days=exp(best["per1"]),
  mcmc_median_days=median(mc_export[, "P1"]),
  max_loglike=max(ll, na.rm=TRUE),
  acceptance_percent=mc_out$acc,
  stringsAsFactors=FALSE
)
write.csv(summary, file.path(out_dir, "summary.csv"), quote=FALSE, row.names=FALSE)
print(summary)

sink(type="message")
sink()
close(log_con)
