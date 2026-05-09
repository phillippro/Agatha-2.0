## Run Agatha BFP with white-noise model only (Nma=0, Nar=0).
## Usage: Rscript run_bfp_white.R <input.vels> <outdir> <ofac>

suppressPackageStartupMessages({
    library(magicaxis); library(foreach); library(doMC); library(parallel); library(minpack.lm)
})
source('periodograms.R'); source('periodoframe.R')
source('mcmc_func.R'); source('sofa.R'); source('orbit.R')
options(warn = 1)

args <- commandArgs(trailingOnly = TRUE)
f <- args[1]; outdir <- args[2]
ofac <- as.numeric(args[3]); if (is.na(ofac)) ofac <- 1
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
target <- sub("\\.vels$", "", basename(f))

tab <- read.table(f, header = TRUE)
tmin <- min(tab[,1]); if (tmin < 2400000) tmin <- tmin + 2400000
t <- tab[,1] - tmin; y <- tab[,2]; ey <- tab[,3]
N <- length(y)
cat(sprintf("%s  N=%d  RMS=%.3f  baseline=%.1f d  ofac=%.2f\n",
            target, N, sqrt(mean(y^2)), diff(range(t)), ofac))

fmin <- 1/max(diff(range(t)), 5000); fmax <- 1/1.1

t0 <- Sys.time()
per <- BFP(t, y, ey, Nma = 0, Nar = 0, Indices = NULL, ofac = ofac,
           model.type = "man", fmin = fmin, fmax = fmax, quantify = TRUE,
           progress = FALSE, GP = FALSE, gp.par = rep(NA, 3),
           noise.only = FALSE, Nsamp = 1, sampling = "combined",
           par.opt = NULL, renew = TRUE)
dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

P <- per$P
lnbf <- if (!is.null(per$lnBFs)) per$lnBFs else per$power
isort <- order(P)

is_peak <- c(FALSE, (lnbf[2:(length(lnbf)-1)] > lnbf[1:(length(lnbf)-2)]) &
                    (lnbf[2:(length(lnbf)-1)] > lnbf[3:length(lnbf)]), FALSE)
ord <- order(lnbf, decreasing = TRUE)
cand_p <- c(); cand_b <- c()
for (i in ord) {
    if (!is_peak[i]) next
    if (length(cand_p) >= 20) break
    if (length(cand_p) == 0 || all(abs(log10(P[i] / cand_p)) > 0.015)) {
        cand_p <- c(cand_p, P[i]); cand_b <- c(cand_b, lnbf[i])
    }
}
cat(sprintf("  BFP+white: time=%.1fs  max(lnBF)=%.2f at P=%.4fd  llmax=%.2f  ll0=%.2f  dll=%.2f\n",
            dt, max(lnbf), P[which.max(lnbf)], per$llmax, per$LogLike0, per$llmax - per$LogLike0))
cat("  Top peaks (lnBF > 5):\n")
for (i in seq_along(cand_p)) {
    flag <- if (cand_b[i] > 5) "*" else " "
    cat(sprintf("  %s %2d  P=%12.4f d  lnBF=%8.3f\n", flag, i, cand_p[i], cand_b[i]))
}

f_asc <- file.path(outdir, paste0(target, "_BFP_white_periodogram.txt"))
write.table(cbind(Period.day = P[isort], lnBF = lnbf[isort]),
            file = f_asc, row.names = FALSE, quote = FALSE)

f_sum <- file.path(outdir, paste0(target, "_BFP_white_summary.txt"))
sink(f_sum)
cat("# Agatha BFP + white-noise model for", target, "\n")
cat("# N=", N, "  RMS=", round(sqrt(mean(y^2)), 3), "m/s\n")
cat("# ofac=", ofac, "  time=", round(dt, 1), "s\n")
cat("# max(lnBF)=", round(max(lnbf), 3), " at P=", round(P[which.max(lnbf)], 4), "d\n")
cat("# ll0=", round(per$LogLike0, 2), "  llmax=", round(per$llmax, 2),
    "  delta_ll=", round(per$llmax - per$LogLike0, 2), "\n")
cat("#\n# Top 20 peaks (period_d, lnBF):\n")
for (i in seq_along(cand_p))
    cat(sprintf("  %2d  %12.4f  %8.3f\n", i, cand_p[i], cand_b[i]))
sink()

fpdf <- file.path(outdir, paste0(target, "_BFP_white_periodogram.pdf"))
pdf(fpdf, 11, 8.5)
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1), mgp = c(2.2, 0.7, 0))
plot(P[isort], lnbf[isort], type = "l", log = "x", lwd = 0.5,
     xlab = "Period (d)", ylab = "ln(BF)",
     main = sprintf("%s   BFP + white noise   N=%d   max lnBF=%.2f at P=%.3f d",
                    target, N, max(lnbf), P[which.max(lnbf)]))
abline(h = 5, lty = 2, col = "red", lwd = 1.2)
abline(h = c(0, 10), lty = 3, col = "grey70")
for (win in c(365.25, 365.25/2, 365.25/3, 1.0027, 0.997)) abline(v = win, lty = 3, col = "grey60")
magaxis(side = c(1, 2), tcl = -0.3)

msel <- P > 50 & P < 300
plot(P[isort][msel[isort]], lnbf[isort][msel[isort]], type = "l", lwd = 0.7,
     xlab = "Period (d)", ylab = "ln(BF)",
     main = "Zoom: 50-300 d")
abline(h = 5, lty = 2, col = "red")
abline(h = 0, lty = 3, col = "grey70")
abline(v = 170, lty = 3, col = "orange")
magaxis(side = c(1, 2), tcl = -0.3)
dev.off()

cat("\nWrote", fpdf, "\n       ", f_sum, "\n       ", f_asc, "\n")
