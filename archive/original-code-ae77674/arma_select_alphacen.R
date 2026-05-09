## ARMA(p,q) model selection for alpha Cen residuals using Agatha's null-model
## likelihood. Grid p=0..4, q=0..4 (25 combos). Rank by BIC and AICc.
##
## Usage: Rscript arma_select_alphacen.R <input.vels> <outdir>
## Produces: <outdir>/<target>_ARMA_select.txt with the full grid table.

suppressPackageStartupMessages({
    library(magicaxis); library(foreach); library(doMC); library(parallel); library(minpack.lm)
})
source('periodograms.R'); source('periodoframe.R')
source('mcmc_func.R'); source('sofa.R'); source('orbit.R')
options(warn = 1)

args <- commandArgs(trailingOnly = TRUE)
f <- args[1]
outdir <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
target <- sub("\\.vels$", "", basename(f))

tab <- read.table(f, header = TRUE)
tmin <- min(tab[,1]); if (tmin < 2400000) tmin <- tmin + 2400000
t <- tab[,1] - tmin; y <- tab[,2]; ey <- tab[,3]
N <- length(y)
cat(sprintf("%s  N=%d  RMS=%.3f  baseline=%.1f d\n", target, N, sqrt(mean(y^2)), diff(range(t))))

fmin <- 1/max(diff(range(t)), 5000); fmax <- 1/1.1

grid <- expand.grid(p = 0:4, q = 0:4)
res <- data.frame()
for (i in seq_len(nrow(grid))) {
    p <- grid$p[i]; q <- grid$q[i]
    label <- sprintf("ARMA(%d,%d)", p, q)
    t0 <- Sys.time()
    # Ultra-coarse frequency grid so BFP just computes null model quickly.
    per <- tryCatch(
        BFP(t, y, ey, Nma = q, Nar = p, Indices = NULL, ofac = 0.005,
            model.type = "man", fmin = fmin, fmax = fmax, quantify = FALSE,
            progress = FALSE, GP = FALSE, gp.par = rep(NA, 3),
            noise.only = FALSE, Nsamp = 1, sampling = "combined",
            par.opt = NULL, renew = TRUE),
        error = function(e) { cat("  ", label, "FAILED:", conditionMessage(e), "\n"); NULL }
    )
    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (is.null(per)) next
    ll0 <- per$LogLike0
    llmax <- per$llmax
    k <- p + q + 2          # ARMA coefficients + 2 mean params (gamma, beta)
    BIC <- k * log(N) - 2 * ll0
    AIC <- 2 * k - 2 * ll0
    AICc <- AIC + 2 * k * (k + 1) / max(N - k - 1, 1)
    cat(sprintf("  %-10s  k=%2d  ll0=%10.2f  BIC=%10.2f  AICc=%10.2f  (%.1fs)\n",
                label, k, ll0, BIC, AICc, dt))
    res <- rbind(res, data.frame(p = p, q = q, k = k, ll0 = ll0,
                                  BIC = BIC, AIC = AIC, AICc = AICc,
                                  llmax = llmax, time_s = dt))
}

res$dBIC <- res$BIC - min(res$BIC, na.rm = TRUE)
res$dAICc <- res$AICc - min(res$AICc, na.rm = TRUE)
res <- res[order(res$BIC), ]

cat("\n=== Ranked by BIC ===\n")
print(res, row.names = FALSE, digits = 5)

best <- res[1, ]
cat(sprintf("\nBest ARMA(%d,%d)  BIC=%.2f  ll0=%.2f\n", best$p, best$q, best$BIC, best$ll0))

fout <- file.path(outdir, paste0(target, "_ARMA_select.txt"))
sink(fout)
cat("# ARMA(p,q) model selection for", target, "\n")
cat("# N=", N, "  RMS=", round(sqrt(mean(y^2)), 3), " m/s\n")
cat("# Columns: p q k ll0 BIC AIC AICc llmax dBIC dAICc time_s\n")
write.table(res, file = "", row.names = FALSE, quote = FALSE)
sink()
cat("Wrote", fout, "\n")
