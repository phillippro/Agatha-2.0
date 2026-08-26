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
check_dir <- "results/HD361_1000d_candidate_checks"
out_pdf <- file.path(check_dir, "HD361_1000d_candidate_summary.pdf")

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

mc <- read.table(file.path(check_dir, "long_candidate_1000d_MCMCposterior.txt"),
                 header=TRUE, check.names=FALSE)
best <- unlist(mc[which.max(mc$loglike), setdiff(colnames(mc), c("logpost", "loglike"))])
best <- as.numeric(best)
names(best) <- setdiff(colnames(mc), c("logpost", "loglike"))
p1 <- exp(best["per1"])
p2 <- exp(best["per2"])

pars_p1 <- best
pars_p1["K2"] <- 0
pars_p2 <- best
pars_p2["K1"] <- 0

fit_rows <- do.call(rbind, lapply(seq_along(dataset), function(k) {
  id <- dataset[k]
  tab <- tabs[[id]]
  model_all <- as.numeric(RV.kepler(best, bases=rep("natural", 10))[[id]])
  signal2 <- as.numeric(RV.kepler(pars_p2, kep.only=TRUE, bases=rep("natural", 10))[[id]])
  red <- arma(t=tab[[1]], ymodel=model_all, ydata=tab[[2]], pars=best,
              ind.set=k, p=0, q=2)$arma
  data.frame(dataset=id, BJD=tab[[1]], RV=tab[[2]], eRV=tab[[3]],
             model=model_all + red, signal2=signal2,
             residual=tab[[2]] - model_all - red)
}))

phase <- ((fit_rows$BJD - tmin) %% p2) / p2
y_phase <- fit_rows$RV - (fit_rows$model - fit_rows$signal2)

t_grid <- seq(min(trv.all), max(trv.all), length.out=4000)
sim2 <- as.numeric(RV.kepler(pars_p2, tt=t_grid, kep.only=TRUE, bases=rep("natural", 10)))
sim_phase <- ((t_grid - tmin) %% p2) / p2
ord <- order(sim_phase)

bfp <- read.table(file.path(check_dir, "BFP_MA2_one_signal_residuals.txt"), header=TRUE)
bgls <- read.table(file.path(check_dir, "BGLS_one_signal_residuals.txt"), header=TRUE)
ls <- read.table(file.path(check_dir, "LS_one_signal_residuals.txt"), header=TRUE)

pdf(out_pdf, width=8, height=10)
layout(matrix(c(1, 2, 3), ncol=1), heights=c(1, 1.3, 0.8))

plot(bfp$Period.day, bfp$power, type="l", log="x", xaxt="n",
     xlim=c(300, 3000), xlab="Period [day]", ylab="ln(BF)",
     main="HD361 one-signal residuals: long-period candidate")
magaxis(side=1, tcl=-0.5)
abline(v=c(722.53, 1000, p2, 1094.05), col=c("grey60", "grey50", "red", "grey60"),
       lty=c(3, 2, 1, 3), lwd=c(1, 1, 2, 1))
legend("topright", bty="n",
       legend=c("BFP MA(2)", paste0("MCMC P2=", round(p2, 2), " d")),
       col=c("black", "red"), lwd=c(1, 2))

cols <- seq_along(dataset)
plot(phase, y_phase, pch=19, cex=0.7, col=cols[match(fit_rows$dataset, dataset)],
     xlab=paste0("Phase at P2 = ", round(p2, 2), " d"),
     ylab="RV minus P1/noise/offset [m/s]",
     main="Two-Keplerian MA(2) fit, folded on candidate P2")
arrows(phase, y_phase - fit_rows$eRV, phase, y_phase + fit_rows$eRV,
       length=0.03, angle=90, code=3, col="grey70")
lines(sim_phase[ord], sim2[ord], col="red", lwd=2)
legend("topright", legend=dataset, col=cols, pch=19, bty="n", cex=0.8)

plot(phase, fit_rows$residual, pch=19, cex=0.7, col=cols[match(fit_rows$dataset, dataset)],
     xlab=paste0("Phase at P2 = ", round(p2, 2), " d"),
     ylab="O-C [m/s]", main=paste0("Residuals; RMS = ", round(sd(fit_rows$residual), 2), " m/s"))
abline(h=0, lty=2, col="grey50")

dev.off()
cat(out_pdf, "\n")
