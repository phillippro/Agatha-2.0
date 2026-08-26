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
out_dir <- file.path("results", "HD361_combined_2signal_mcmc_kepler")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

files <- sort(Sys.glob(file.path(input_dir, "*")))
dataset <- tools::file_path_sans_ext(basename(files))
tabs <- lapply(files, read.agatha.table)
names(tabs) <- dataset

rows_from_tabs <- function(tabs, values=NULL) {
  out_rows <- do.call(rbind, Map(function(tab, id) {
    y <- if (is.null(values)) tab[[2]] else values[[id]]
    data.frame(t=tab[[1]], y=y, dy=tab[[3]], set.id=id)
  }, tabs, names(tabs)))
  out_rows <- out_rows[order(out_rows$t),]
  out_rows$set.id <- factor(out_rows$set.id, levels=names(tabs))
  out_rows
}

periodogram_plot <- function(per, path, title) {
  tab <- data.frame(Period.day=per$P, logBF=per$power)
  tab <- tab[order(tab$Period.day),]
  pdf(path, width=7, height=5)
  plot(tab$Period.day, tab$logBF, type="l", log="x", xaxt="n",
       xlab="Period [day]", ylab="ln(BF)", main=title)
  magaxis(side=1, tcl=-0.5)
  abline(v=per$Popt, col="red", lwd=2)
  abline(h=5, lty=2, col="grey50")
  legend("topright", legend=paste0("P=", signif(per$Popt, 6), " d"),
         bty="n", text.col="red")
  dev.off()
  tab
}

make_start <- function(periods, amps, phases, gammas, mas, jitters, beta0, trend0=0) {
  start <- c()
  for (i in seq_along(periods)) {
    tmp <- c(log(periods[i]), amps[i], 0.1, 0, phases[i] %% (2 * pi))
    names(tmp) <- paste0(c("per", "K", "e", "omega", "Mo"), i)
    start <- c(start, tmp)
  }
  start <- c(start, a11=trend0)
  for (k in seq_along(dataset)) {
    tmp <- c(gammas[k], jitters[k], mas[k], beta0[k])
    names(tmp) <- c(paste0("b", k), paste0("s", k), paste0("w", k, "1"), paste0("beta", k))
    start <- c(start, tmp)
  }
  start
}

configure_bounds <- function(startvalue, periods, amps) {
  span <- max(trv.all) - min(trv.all)
  Dt <- span / time.unit
  par.min <<- startvalue - abs(startvalue)
  par.max <<- startvalue + abs(startvalue)
  for (i in seq_along(periods)) {
    par.min[paste0("per", i)] <<- log(max(1.1, periods[i] * 0.4))
    par.max[paste0("per", i)] <<- log(min(1e5, periods[i] * 2.5))
    par.min[paste0("K", i)] <<- 0
    par.max[paste0("K", i)] <<- max(2 * amps[i], amps[i] + 20)
    par.min[paste0("e", i)] <<- 0
    par.max[paste0("e", i)] <<- 0.95
    par.min[paste0("omega", i)] <<- 0
    par.max[paste0("omega", i)] <<- 2 * pi
    par.min[paste0("Mo", i)] <<- 0
    par.max[paste0("Mo", i)] <<- 2 * pi
  }
  par.min["a11"] <<- -10 * sd(raw_rows$y) / Dt
  par.max["a11"] <<- 10 * sd(raw_rows$y) / Dt
  for (k in seq_along(dataset)) {
    rv_sd <- sd(tabs[[dataset[k]]][[2]], na.rm=TRUE)
    par.min[paste0("b", k)] <<- min(tabs[[dataset[k]]][[2]], na.rm=TRUE) - 5 * rv_sd
    par.max[paste0("b", k)] <<- max(tabs[[dataset[k]]][[2]], na.rm=TRUE) + 5 * rv_sd
    par.min[paste0("s", k)] <<- 0
    par.max[paste0("s", k)] <<- max(1, 2 * rv_sd)
    par.min[paste0("w", k, "1")] <<- wmin
    par.max[paste0("w", k, "1")] <<- wmax
    par.min[paste0("beta", k)] <<- beta.min
    par.max[paste0("beta", k)] <<- beta.max
  }
  Npar <<- length(startvalue)
  Sd <<- 2.4^2 / length(startvalue)
}

fit_mcmc <- function(startvalue, label) {
  configure_bounds(startvalue, exp(startvalue[grepl("^per", names(startvalue))]),
                   startvalue[grepl("^K", names(startvalue))])
  opt <- optim(startvalue, function(p) -loglikelihood(p, bases=rep("natural", 10)),
               method="L-BFGS-B", lower=par.min, upper=par.max,
               control=list(maxit=500))
  startvalue <- opt$par
  cat(label, "optimized negative logL:", opt$value, "\n")
  print(startvalue)

  cov.start <- diag(length(startvalue)) * 1e-5
  diag(cov.start)[grepl("^per", names(startvalue))] <- 1e-6
  diag(cov.start)[grepl("^K", names(startvalue))] <- 1e-3
  diag(cov.start)[grepl("^e", names(startvalue))] <- 1e-5
  diag(cov.start)[grepl("^omega|^Mo", names(startvalue))] <- 1e-3
  run.metropolis.MCMC(startvalue=startvalue, cov.start=cov.start,
                      iterations=Niter, tem=1, bases=rep("natural", 10))
}

export_fit <- function(mc_out, label) {
  mc <- mc_out$out
  write.table(mc, file.path(out_dir, paste0("HD361_combined_", label, "_MCposterior.txt")),
              quote=FALSE, row.names=FALSE)
  ll <- mc[, "loglike"]
  par_cols <- setdiff(colnames(mc), c("logpost", "loglike"))
  mc_par <- mc[, par_cols, drop=FALSE]
  mc_export <- mc_par
  per_cols <- grep("^per", colnames(mc_export))
  mc_export[, per_cols] <- exp(mc_export[, per_cols])
  colnames(mc_export)[per_cols] <- gsub("^per", "P", colnames(mc_export)[per_cols])
  par_stat <- t(apply(mc_export, 2, function(x) {
    c(xopt=x[which.max(ll)], x1per=quantile(x, 0.01), x99per=quantile(x, 0.99),
      x10per=quantile(x, 0.10), x90per=quantile(x, 0.90),
      med=median(x), mean=mean(x), sd=sd(x))
  }))
  write.csv(par_stat, file.path(out_dir, paste0("HD361_combined_", label, "_OptPar.csv")),
            quote=FALSE)
  best <- mc[which.max(ll), par_cols]
  list(mc=mc, best=best, par_stat=par_stat, llmax=max(ll), acc=mc_out$acc)
}

fit_rows_for <- function(best, label) {
  rv_all <- RV.kepler(best, bases=rep("natural", 10))
  rv_sig <- RV.kepler(best, kep.only=TRUE, bases=rep("natural", 10))
  out_rows <- do.call(rbind, lapply(seq_along(dataset), function(k) {
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
  write.table(out_rows, file.path(out_dir, paste0("HD361_combined_", label, "_fit_residuals.txt")),
              quote=FALSE, row.names=FALSE)
  out_rows
}

phase_plot <- function(best, fit_rows, label, path) {
  P_best <- exp(best["per1"])
  t_grid <- seq(min(trv.all), max(trv.all), length.out=2000)
  sim_signal <- RV.kepler(best, tt=t_grid, kep.only=TRUE, bases=rep("natural", 10))
  sim <- data.frame(BJD=t_grid, phase=(t_grid - tmin) %% P_best,
                    kepler_signal=as.numeric(sim_signal))
  write.table(sim, file.path(out_dir, paste0("HD361_combined_", label, "_sim_signal.txt")),
              quote=FALSE, row.names=FALSE)
  pdf(path, width=7, height=7)
  layout(matrix(c(1, 2), nrow=2), heights=c(2, 1))
  phase <- (fit_rows$BJD - tmin) %% P_best
  sim_phase <- sim$phase
  ord <- order(sim_phase)
  y_signal_space <- fit_rows$RV - (fit_rows$model - fit_rows$kepler_signal)
  plot(phase, y_signal_space, xlab="", ylab="RV minus offsets/noise [m/s]", xaxt="n",
       main=paste0("HD361 ", label, ", phase to P1=", signif(P_best, 6), " d"))
  axis(side=1, labels=FALSE)
  arrows(phase, y_signal_space - fit_rows$eRV, phase, y_signal_space + fit_rows$eRV,
         length=0.03, angle=90, code=3, col="grey60")
  points(phase, y_signal_space, pch=19, cex=0.65, col=as.integer(factor(fit_rows$dataset)))
  lines(sim_phase[ord], sim$kepler_signal[ord], col="red", lwd=2)
  legend("topright", legend=dataset, pch=19, col=seq_along(dataset), bty="n", cex=0.8)
  plot(phase, fit_rows$residual, xlab="Orbital phase [day]", ylab="O-C [m/s]",
       pch=19, cex=0.65, col=as.integer(factor(fit_rows$dataset)))
  abline(h=0, lty=2, col="grey50")
  dev.off()
}

ofac <- 1
fmin <- 1 / 1e5
fmax <- 1 / 1.1
Niter <- 1000
Nma <- 1
Nar <- 0

log_con <- file(file.path(out_dir, "run.log"), open="wt")
sink(log_con)
sink(log_con, type="message")

raw_rows <- rows_from_tabs(tabs)

target <- "HD361"
prior.type <- "mt"
period.par <- "logP"
time.unit <- 365.25
offset <- TRUE
ins <- dataset
trv.all <- unlist(lapply(tabs, function(tab) tab[[1]]), use.names=FALSE)
tmin <- min(trv.all)
out <- list(trv.all=trv.all, ins=ins, prior.type=prior.type)
start_index <- 1
for (k in seq_along(dataset)) {
  tab <- tabs[[dataset[k]]]
  n <- nrow(tab)
  out[[dataset[k]]] <- list(RV=as.matrix(tab[, 1:3]),
                            index=start_index:(start_index + n - 1),
                            noise=list(nqp=c(0, Nma, Nar)))
  start_index <- start_index + n
}

span <- max(trv.all) - min(trv.all)
Esd <- 0.1
phi.min <- wmin <- -1
phi.max <- wmax <- 1
alpha.min <- beta.min <- log(1 / 24)
alpha.max <- beta.max <- log(span)
tol <- 1e-16
tol1 <- 1e-12

cat("Combined HD361 two-signal Keplerian Agatha run\n")
cat("datasets:", paste(dataset, collapse=", "), "\n")
cat("Ndata:", nrow(raw_rows), "Niter:", Niter, "Nma:", Nma, "Nar:", Nar, "\n")

per1 <- BFP.multiset(t=raw_rows$t, y=raw_rows$y, dy=raw_rows$dy, set.id=raw_rows$set.id,
                     Nma=Nma, Nar=Nar, ofac=ofac, fmin=fmin, fmax=fmax)
per1_tab <- periodogram_plot(per1, file.path(out_dir, "HD361_combined_periodogram_signal1.pdf"),
                             "HD361 combined BFP signal 1")
write.table(per1_tab, file.path(out_dir, "HD361_combined_periodogram_signal1.txt"),
            quote=FALSE, row.names=FALSE)

beta0 <- sapply(seq_along(dataset), function(k) log(max(1, median(diff(sort(tabs[[dataset[k]]][[1]])), na.rm=TRUE))))
jitters <- sapply(tabs, function(tab) max(0.1, 0.5 * median(tab[[3]], na.rm=TRUE)))
gammas1 <- sapply(dataset, function(id) unname(per1$par.opt[paste0("gamma_", make.names(id))]))
mas1 <- sapply(dataset, function(id) max(-0.9, min(0.9, unname(per1$par.opt[paste0("ma_", make.names(id), "_lag1")]))))
A1 <- unname(per1$par.opt["A"])
B1 <- unname(per1$par.opt["B"])
start1 <- make_start(periods=per1$Popt, amps=max(sqrt(A1^2 + B1^2), 1e-3),
                     phases=as.numeric(xy2phi(A1, B1)), gammas=gammas1,
                     mas=mas1, jitters=jitters, beta0=beta0,
                     trend0=unname(per1$par.opt["beta"]) * time.unit)
set.seed(361)
fit1 <- export_fit(fit_mcmc(start1, "one_signal"), "one_signal_kepler_MA1")
fit_rows1 <- fit_rows_for(fit1$best, "one_signal_kepler_MA1")
phase_plot(fit1$best, fit_rows1, "one_signal_kepler_MA1",
           file.path(out_dir, "HD361_combined_one_signal_kepler_MA1_phase_fit.pdf"))

residual1 <- split(fit_rows1$residual, fit_rows1$dataset)
per2_rows <- rows_from_tabs(tabs, residual1)
per2 <- BFP.multiset(t=per2_rows$t, y=per2_rows$y, dy=per2_rows$dy, set.id=per2_rows$set.id,
                     Nma=Nma, Nar=Nar, ofac=ofac, fmin=fmin, fmax=fmax)
per2_tab <- periodogram_plot(per2, file.path(out_dir, "HD361_combined_periodogram_signal2.pdf"),
                             "HD361 combined BFP after signal 1")
write.table(per2_tab, file.path(out_dir, "HD361_combined_periodogram_signal2.txt"),
            quote=FALSE, row.names=FALSE)

A2 <- unname(per2$par.opt["A"])
B2 <- unname(per2$par.opt["B"])
planet1 <- fit1$best[c("per1", "K1", "e1", "omega1", "Mo1")]
planet2 <- c(log(per2$Popt), max(sqrt(A2^2 + B2^2), 1e-3), 0.1, 0, as.numeric(xy2phi(A2, B2)) %% (2 * pi))
names(planet2) <- c("per2", "K2", "e2", "omega2", "Mo2")
noise2 <- c(a11=unname(fit1$best["a11"]))
for (k in seq_along(dataset)) {
  tmp <- fit1$best[c(paste0("b", k), paste0("s", k), paste0("w", k, "1"), paste0("beta", k))]
  names(tmp) <- c(paste0("b", k), paste0("s", k), paste0("w", k, "1"), paste0("beta", k))
  noise2 <- c(noise2, tmp)
}
start2 <- c(planet1, planet2, noise2)
set.seed(362)
fit2 <- export_fit(fit_mcmc(start2, "two_signal"), "two_signal_kepler_MA1")
fit_rows2 <- fit_rows_for(fit2$best, "two_signal_kepler_MA1")
phase_plot(fit2$best, fit_rows2, "two_signal_kepler_MA1",
           file.path(out_dir, "HD361_combined_two_signal_kepler_MA1_phase_fit.pdf"))

residual2 <- split(fit_rows2$residual, fit_rows2$dataset)
per3_rows <- rows_from_tabs(tabs, residual2)
per3 <- BFP.multiset(t=per3_rows$t, y=per3_rows$y, dy=per3_rows$dy, set.id=per3_rows$set.id,
                     Nma=Nma, Nar=Nar, ofac=ofac, fmin=fmin, fmax=fmax)
per3_tab <- periodogram_plot(per3, file.path(out_dir, "HD361_combined_periodogram_residual_after_2signals.pdf"),
                             "HD361 combined BFP residual after 2 signals")
write.table(per3_tab, file.path(out_dir, "HD361_combined_periodogram_residual_after_2signals.txt"),
            quote=FALSE, row.names=FALSE)

pdf(file.path(out_dir, "HD361_combined_three_periodograms.pdf"), width=7, height=9)
par(mfrow=c(3, 1), mar=c(4, 4, 2, 1))
for (obj in list(list(per=per1, title="Signal 1"),
                 list(per=per2, title="Signal 2 after subtracting signal 1"),
                 list(per=per3, title="Residual after two signals"))) {
  tab <- data.frame(P=obj$per$P, logBF=obj$per$power)
  tab <- tab[order(tab$P),]
  plot(tab$P, tab$logBF, type="l", log="x", xaxt="n",
       xlab="Period [day]", ylab="ln(BF)", main=obj$title)
  magaxis(side=1, tcl=-0.5)
  abline(v=obj$per$Popt, col="red", lwd=2)
  abline(h=5, lty=2, col="grey50")
}
dev.off()

summary <- data.frame(
  step=c("signal1_raw", "signal2_after_signal1", "residual_after_2signals"),
  periodogram_peak_days=c(per1$Popt, per2$Popt, per3$Popt),
  max_logBF=c(max(per1$power, na.rm=TRUE), max(per2$power, na.rm=TRUE), max(per3$power, na.rm=TRUE)),
  mcmc_best_P1_days=c(exp(fit1$best["per1"]), exp(fit2$best["per1"]), NA),
  mcmc_best_P2_days=c(NA, exp(fit2$best["per2"]), NA),
  mcmc_median_P1_days=c(fit1$par_stat["P1", "med"], fit2$par_stat["P1", "med"], NA),
  mcmc_median_P2_days=c(NA, fit2$par_stat["P2", "med"], NA),
  max_loglike=c(fit1$llmax, fit2$llmax, NA),
  acceptance_percent=c(fit1$acc, fit2$acc, NA),
  n_data=nrow(raw_rows),
  n_sets=length(dataset),
  stringsAsFactors=FALSE
)
write.csv(summary, file.path(out_dir, "summary.csv"), quote=FALSE, row.names=FALSE)
print(summary)

sink(type="message")
sink()
close(log_con)
