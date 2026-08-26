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
base_dir <- "results/HD361_combined_2signal_nightly_ma_modelselect_mcmc_kepler"
out_dir <- file.path("results", "HD361_1000d_candidate_checks")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

nightly_bin <- function(tab) {
  night <- floor(tab[[1]] - 0.5)
  binned <- do.call(rbind, lapply(split(seq_len(nrow(tab)), night), function(ii) {
    w <- 1 / tab[[3]][ii]^2
    c(BJD=sum(tab[[1]][ii] * w) / sum(w),
      RV=sum(tab[[2]][ii] * w) / sum(w),
      eRV=sqrt(1 / sum(w)))
  }))
  binned <- as.data.frame(binned)
  rownames(binned) <- NULL
  binned
}

files <- sort(Sys.glob(file.path(input_dir, "*")))
dataset <- tools::file_path_sans_ext(basename(files))
tabs <- lapply(files, function(f) nightly_bin(read.agatha.table(f)))
names(tabs) <- dataset

res1 <- read.table(file.path(base_dir, "HD361_combined_one_signal_kepler_ARMA02_fit_residuals.txt"),
                   header=TRUE, check.names=FALSE)
res1 <- res1[order(res1$BJD),]
res1$set.id <- factor(res1$dataset, levels=dataset)

peak_summary <- data.frame(method=character(), global_period=numeric(), global_power=numeric(),
                           peak_500_1500=numeric(), power_500_1500=numeric(),
                           peak_300_3000=numeric(), power_300_3000=numeric(),
                           power_at_1000=numeric(), stringsAsFactors=FALSE)

extract_peaks <- function(P, power, method) {
  finite <- is.finite(P) & is.finite(power)
  P <- P[finite]
  power <- power[finite]
  ord <- order(P)
  P <- P[ord]
  power <- power[ord]
  band1 <- P >= 500 & P <= 1500
  band2 <- P >= 300 & P <= 3000
  p1000 <- power[which.min(abs(P - 1000))]
  data.frame(method=method,
             global_period=P[which.max(power)],
             global_power=max(power),
             peak_500_1500=if (any(band1)) P[band1][which.max(power[band1])] else NA_real_,
             power_500_1500=if (any(band1)) max(power[band1]) else NA_real_,
             peak_300_3000=if (any(band2)) P[band2][which.max(power[band2])] else NA_real_,
             power_300_3000=if (any(band2)) max(power[band2]) else NA_real_,
             power_at_1000=p1000,
             stringsAsFactors=FALSE)
}

write_periodogram <- function(tab, name, ylab) {
  write.table(tab, file.path(out_dir, paste0(name, ".txt")), quote=FALSE, row.names=FALSE)
  pdf(file.path(out_dir, paste0(name, ".pdf")), width=7, height=5)
  plot(tab$Period.day, tab$power, type="l", log="x", xaxt="n",
       xlab="Period [day]", ylab=ylab, main=name)
  magaxis(side=1, tcl=-0.5)
  abline(v=c(500, 1000, 1500), col=c("grey70", "red", "grey70"), lty=c(3, 2, 3))
  dev.off()
}

## Selected red-noise model: MA(2), combined multi-set residual BFP/MLP.
bfp <- BFP.multiset(t=res1$BJD, y=res1$residual, dy=res1$eRV, set.id=res1$set.id,
                    Nma=2, Nar=0, ofac=1, fmin=1 / 10000, fmax=1 / 10)
bfp_tab <- data.frame(Period.day=bfp$P, power=bfp$power)
write_periodogram(bfp_tab[order(bfp_tab$Period.day),], "BFP_MA2_one_signal_residuals", "ln(BF)")
peak_summary <- rbind(peak_summary, extract_peaks(bfp$P, bfp$power, "BFP_MA2"))

mlp <- MLP.multiset(t=res1$BJD, y=res1$residual, dy=res1$eRV, set.id=res1$set.id,
                    Nma=2, Nar=0, ofac=1, fmin=1 / 10000, fmax=1 / 10)
mlp_tab <- data.frame(Period.day=mlp$P, power=mlp$power)
write_periodogram(mlp_tab[order(mlp_tab$Period.day),], "MLP_MA2_one_signal_residuals", "relative log(ML)")
peak_summary <- rbind(peak_summary, extract_peaks(mlp$P, mlp$power, "MLP_MA2"))

## Offset-removed residual checks with simpler periodogram methods.
res_simple <- res1
res_simple$y <- res_simple$residual - ave(res_simple$residual, res_simple$dataset, FUN=function(z) weighted.mean(z, 1 / res_simple$eRV[match(z, res_simple$residual)]^2))
## The ave/weighted match above is fragile when values repeat; compute robustly by split.
for (id in dataset) {
  ii <- which(res_simple$dataset == id)
  res_simple$y[ii] <- res_simple$residual[ii] - weighted.mean(res_simple$residual[ii], 1 / res_simple$eRV[ii]^2)
}
t0 <- res_simple$BJD - min(res_simple$BJD)
simple_methods <- list(
  GLST=glst(t=res_simple$BJD, y=res_simple$y, err=res_simple$eRV, ofac=1, fmin=1 / 10000, fmax=1 / 10),
  GLS=gls(t=t0, y=res_simple$y, err=res_simple$eRV, ofac=1, fmin=1 / 10000, fmax=1 / 10),
  BGLS=bgls(t=t0, y=res_simple$y, err=res_simple$eRV, ofac=1, fmin=1 / 10000, fmax=1 / 10),
  LS=lsp(times=t0, x=res_simple$y, ofac=1, from=1 / 10000, to=1 / 10)
)
for (nm in names(simple_methods)) {
  obj <- simple_methods[[nm]]
  power <- if (!is.null(obj$power)) obj$power else obj$y
  tab <- data.frame(Period.day=obj$P, power=power)
  write_periodogram(tab[order(tab$Period.day),], paste0(nm, "_one_signal_residuals"), "power")
  peak_summary <- rbind(peak_summary, extract_peaks(obj$P, power, nm))
}

write.csv(peak_summary, file.path(out_dir, "periodogram_peak_summary.csv"),
          quote=FALSE, row.names=FALSE)

## Direct likelihood comparison: two-Keplerian MA(2) fit initialized at 15 d vs long-period candidates.
setup_globals <- function() {
  assign("target", "HD361", envir=.GlobalEnv)
  assign("prior.type", "mt", envir=.GlobalEnv)
  assign("period.par", "logP", envir=.GlobalEnv)
  assign("time.unit", 365.25, envir=.GlobalEnv)
  assign("offset", TRUE, envir=.GlobalEnv)
  assign("ins", dataset, envir=.GlobalEnv)
  assign("trv.all", unlist(lapply(tabs, function(tab) tab[[1]]), use.names=FALSE), envir=.GlobalEnv)
  assign("tmin", min(trv.all), envir=.GlobalEnv)
  tmp_out <- list(trv.all=trv.all, ins=ins, prior.type=prior.type)
  start_index <- 1
  for (k in seq_along(dataset)) {
    tab <- tabs[[dataset[k]]]
    n <- nrow(tab)
    tmp_out[[dataset[k]]] <- list(RV=as.matrix(tab[, 1:3]),
                                  index=start_index:(start_index + n - 1),
                                  noise=list(nqp=c(0, 2, 0)))
    start_index <- start_index + n
  }
  assign("out", tmp_out, envir=.GlobalEnv)
  Esd <<- 0.1
  phi.min <<- wmin <<- -1
  phi.max <<- wmax <<- 1
  span <- max(trv.all) - min(trv.all)
  alpha.min <<- beta.min <<- log(1 / 24)
  alpha.max <<- beta.max <<- log(span)
  tol <<- 1e-16
  tol1 <<- 1e-12
}

setup_bounds <- function(startvalue, p2_center) {
  par.min <<- startvalue - abs(startvalue)
  par.max <<- startvalue + abs(startvalue)
  for (i in 1:2) {
    p0 <- exp(startvalue[paste0("per", i)])
    par.min[paste0("per", i)] <<- log(max(10, p0 * 0.5))
    par.max[paste0("per", i)] <<- log(min(10000, p0 * 2))
    par.min[paste0("K", i)] <<- 0
    par.max[paste0("K", i)] <<- max(80, 3 * startvalue[paste0("K", i)])
    par.min[paste0("e", i)] <<- 0
    par.max[paste0("e", i)] <<- 0.95
    par.min[paste0("omega", i)] <<- 0
    par.max[paste0("omega", i)] <<- 2 * pi
    par.min[paste0("Mo", i)] <<- 0
    par.max[paste0("Mo", i)] <<- 2 * pi
  }
  par.min["per2"] <<- log(max(10, p2_center * 0.75))
  par.max["per2"] <<- log(min(10000, p2_center * 1.25))
  Dt <- (max(trv.all) - min(trv.all)) / time.unit
  all_rv <- unlist(lapply(tabs, function(tab) tab[[2]]), use.names=FALSE)
  par.min["a11"] <<- -10 * sd(all_rv) / Dt
  par.max["a11"] <<- 10 * sd(all_rv) / Dt
  for (k in seq_along(dataset)) {
    rv_sd <- sd(tabs[[dataset[k]]][[2]], na.rm=TRUE)
    par.min[paste0("b", k)] <<- min(tabs[[dataset[k]]][[2]], na.rm=TRUE) - 5 * rv_sd
    par.max[paste0("b", k)] <<- max(tabs[[dataset[k]]][[2]], na.rm=TRUE) + 5 * rv_sd
    par.min[paste0("s", k)] <<- 0
    par.max[paste0("s", k)] <<- max(1, 2 * rv_sd)
    par.min[paste0("w", k, 1:2)] <<- -1
    par.max[paste0("w", k, 1:2)] <<- 1
    par.min[paste0("beta", k)] <<- beta.min
    par.max[paste0("beta", k)] <<- beta.max
  }
}

build_start_two <- function(p2) {
  one <- read.csv(file.path(base_dir, "HD361_combined_one_signal_kepler_ARMA02_OptPar.csv"),
                  row.names=1, check.names=FALSE)
  par1 <- c(log(one["P1", "xopt"]), one["K1", "xopt"], one["e1", "xopt"],
            one["omega1", "xopt"], one["Mo1", "xopt"])
  names(par1) <- c("per1", "K1", "e1", "omega1", "Mo1")
  p2_phase <- 0
  p2_amp <- 1
  band <- bfp$P >= p2 * 0.75 & bfp$P <= p2 * 1.25
  if (any(band)) {
    ind <- which(band)[which.max(bfp$power[band])]
    p2 <- bfp$P[band][which.max(bfp$power[band])]
    opt <- bfp$par.opt
    p2_amp <- max(sqrt(unname(opt["A"])^2 + unname(opt["B"])^2), 1e-3)
    p2_phase <- as.numeric(xy2phi(unname(opt["A"]), unname(opt["B"])))
  }
  par2 <- c(log(p2), p2_amp, 0.1, 0, p2_phase %% (2 * pi))
  names(par2) <- c("per2", "K2", "e2", "omega2", "Mo2")
  noise_names <- c("a11")
  for (k in seq_along(dataset)) {
    noise_names <- c(noise_names, paste0("b", k), paste0("s", k),
                     paste0("w", k, 1:2), paste0("beta", k))
  }
  noise <- one[noise_names, "xopt"]
  names(noise) <- noise_names
  c(par1, par2, noise)
}

setup_globals()
candidate_periods <- c(15.2, peak_summary$peak_500_1500[peak_summary$method == "BFP_MA2"],
                       1000, peak_summary$peak_300_3000[peak_summary$method == "BFP_MA2"])
candidate_periods <- unique(round(candidate_periods[is.finite(candidate_periods)], 6))
opt_params <- list()
likelihood_checks <- do.call(rbind, lapply(candidate_periods, function(p2) {
  start <- build_start_two(p2)
  setup_bounds(start, p2)
  opt <- optim(start, function(p) -loglikelihood(p, bases=rep("natural", 10)),
               method="L-BFGS-B", lower=par.min, upper=par.max,
               control=list(maxit=3000))
  opt_params[[as.character(p2)]] <<- opt$par
  data.frame(seed_P2=p2, opt_P1=exp(opt$par["per1"]), opt_P2=exp(opt$par["per2"]),
             K2=opt$par["K2"], neg_loglike=opt$value, loglike=-opt$value,
             convergence=opt$convergence, stringsAsFactors=FALSE)
}))
likelihood_checks <- likelihood_checks[order(likelihood_checks$neg_loglike),]
write.csv(likelihood_checks, file.path(out_dir, "two_signal_likelihood_seed_comparison.csv"),
          quote=FALSE, row.names=FALSE)

best_seed <- as.character(likelihood_checks$seed_P2[1])
best_start <- opt_params[[best_seed]]
setup_bounds(best_start, likelihood_checks$opt_P2[1])
Npar <- length(best_start)
Sd <- 2.4^2 / Npar
cov.start <- diag(length(best_start)) * 1e-5
diag(cov.start)[grepl("^per", names(best_start))] <- 1e-6
diag(cov.start)[grepl("^K", names(best_start))] <- 1e-3
diag(cov.start)[grepl("^e", names(best_start))] <- 1e-5
diag(cov.start)[grepl("^omega|^Mo", names(best_start))] <- 1e-3
set.seed(1037)
long_mc <- run.metropolis.MCMC(startvalue=best_start, cov.start=cov.start,
                               iterations=1000, tem=1, bases=rep("natural", 10))
write.table(long_mc$out, file.path(out_dir, "long_candidate_1000d_MCMCposterior.txt"),
            quote=FALSE, row.names=FALSE)
ll <- long_mc$out[, "loglike"]
par_cols <- setdiff(colnames(long_mc$out), c("logpost", "loglike"))
mc_export <- long_mc$out[, par_cols, drop=FALSE]
per_cols <- grep("^per", colnames(mc_export))
mc_export[, per_cols] <- exp(mc_export[, per_cols])
colnames(mc_export)[per_cols] <- gsub("^per", "P", colnames(mc_export)[per_cols])
par_stat <- t(apply(mc_export, 2, function(x) {
  c(xopt=x[which.max(ll)], x1per=quantile(x, 0.01), x99per=quantile(x, 0.99),
    x10per=quantile(x, 0.10), x90per=quantile(x, 0.90),
    med=median(x), mean=mean(x), sd=sd(x))
}))
write.csv(par_stat, file.path(out_dir, "long_candidate_1000d_OptPar.csv"),
          quote=FALSE)

pdf(file.path(out_dir, "periodogram_method_comparison.pdf"), width=8, height=10)
par(mfrow=c(3, 2), mar=c(4, 4, 2, 1))
for (nm in c("BFP_MA2", "MLP_MA2", "GLST", "GLS", "BGLS", "LS")) {
  path <- file.path(out_dir, paste0(nm, "_one_signal_residuals.txt"))
  if (file.exists(path)) {
    tab <- read.table(path, header=TRUE)
    plot(tab$Period.day, tab$power, type="l", log="x", xaxt="n",
         xlab="Period [day]", ylab="power", main=nm)
    magaxis(side=1, tcl=-0.5)
    abline(v=c(500, 1000, 1500), col=c("grey70", "red", "grey70"), lty=c(3, 2, 3))
  }
}
dev.off()

print(peak_summary)
print(likelihood_checks)
