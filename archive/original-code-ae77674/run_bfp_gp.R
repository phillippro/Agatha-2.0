## Run Agatha BFP with GP (SHO celerite) noise model to absorb stellar rotation.
## The SHO kernel models quasi-periodic rotation at ~38 d for alpha Cen B.
## Usage: Rscript run_bfp_gp.R <input.vels> <outdir> <ofac> [Prot]
##   Prot: rotation period to fix in GP (default 38.4 d); use "free" to let GP find it.

suppressPackageStartupMessages({
    library(magicaxis); library(foreach); library(doMC); library(parallel); library(minpack.lm)
})
source('periodograms.R'); source('periodoframe.R')
source('mcmc_func.R'); source('sofa.R'); source('orbit.R')
options(warn = 1)

args <- commandArgs(trailingOnly = TRUE)
f <- args[1]; outdir <- args[2]
ofac <- as.numeric(args[3]); if (is.na(ofac)) ofac <- 0.5
Prot_arg <- if (length(args) >= 4) args[4] else "38.4"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
target <- sub("\\.vels$", "", basename(f))

tab <- read.table(f, header = TRUE)
tmin <- min(tab[,1]); if (tmin < 2400000) tmin <- tmin + 2400000
t <- tab[,1] - tmin; y <- tab[,2]; ey <- tab[,3]
N <- length(y)

# GP parameters: c(sigmaGP, logProt, logtauGP)
# NA = free parameter; fixed value = locked
if (Prot_arg == "free") {
    gp.par <- rep(NA, 3)  # all free
    gp_label <- "GP_free"
    cat(sprintf("%s  N=%d  GP: all free  ofac=%.2f\n", target, N, ofac))
} else {
    Prot <- as.numeric(Prot_arg)
    gp.par <- c(NA, log(Prot), NA)  # fix rotation period, free amplitude & damping
    gp_label <- sprintf("GP_Prot%.0f", Prot)
    cat(sprintf("%s  N=%d  GP: Prot=%.1fd (fixed)  ofac=%.2f\n", target, N, Prot, ofac))
}

fmin <- 1/max(diff(range(t)), 5000); fmax <- 1/1.1

t0 <- Sys.time()
per <- tryCatch(
    BFP(t, y, ey, Nma = 0, Nar = 0, Indices = NULL, ofac = ofac,
        model.type = "man", fmin = fmin, fmax = fmax, quantify = TRUE,
        progress = FALSE, GP = TRUE, gp.par = gp.par,
        noise.only = FALSE, Nsamp = 1, sampling = "combined",
        par.opt = NULL, renew = TRUE),
    error = function(e) { cat("  BFP+GP FAILED:", conditionMessage(e), "\n"); NULL }
)
dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
if (is.null(per)) quit(status = 1)

P <- per$P
lnbf <- if (!is.null(per$lnBFs)) per$lnBFs else per$power
isort <- order(P)

# Find peaks
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

cat(sprintf("  BFP+%s: time=%.1fs  max(lnBF)=%.2f at P=%.4fd\n",
            gp_label, dt, max(lnbf), P[which.max(lnbf)]))
cat(sprintf("  ll0=%.2f  llmax=%.2f  dll=%.2f\n", per$LogLike0, per$llmax, per$llmax - per$LogLike0))

# Check specific periods
for (Pcheck in c(38.4, 170, 365.25)) {
    m <- which(P > Pcheck*0.9 & P < Pcheck*1.1)
    if (length(m) > 0) {
        cat(sprintf("  at ~%.0fd: lnBF=%.3f at P=%.2f\n", Pcheck, max(lnbf[m]), P[m][which.max(lnbf[m])]))
    }
}

cat("  Top 20 peaks:\n")
for (i in seq_along(cand_p)) {
    flag <- if (cand_b[i] > 5) "*" else " "
    cat(sprintf("  %s %2d  P=%12.4f d  lnBF=%8.3f\n", flag, i, cand_p[i], cand_b[i]))
}

# Save periodogram ASCII
tag <- paste0(gp_label, "_ofac", ofac)
f_asc <- file.path(outdir, paste0(target, '_BFP_', tag, '_periodogram.txt'))
write.table(cbind(Period.day = P[isort], lnBF = lnbf[isort]),
            file = f_asc, row.names = FALSE, quote = FALSE)

# Save summary
f_sum <- file.path(outdir, paste0(target, '_BFP_', tag, '_summary.txt'))
sink(f_sum)
cat("# Agatha BFP +", gp_label, "for", target, "\n")
cat("# N=", N, "  RMS=", round(sqrt(mean(y^2)), 3), "m/s  ofac=", ofac, "\n")
cat("# max(lnBF)=", round(max(lnbf), 3), " at P=", round(P[which.max(lnbf)], 4), "d\n")
cat("# ll0=", round(per$LogLike0, 2), "  llmax=", round(per$llmax, 2),
    "  delta_ll=", round(per$llmax - per$LogLike0, 2), "\n")
cat("#\n")
for (i in seq_along(cand_p))
    cat(sprintf("  %2d  %12.4f  %8.3f\n", i, cand_p[i], cand_b[i]))
sink()

# PDF
fpdf <- file.path(outdir, paste0(target, '_BFP_', tag, '_periodogram.pdf'))
pdf(fpdf, 11, 8.5)
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1), mgp = c(2.2, 0.7, 0))
plot(P[isort], lnbf[isort], type = "l", log = "x", lwd = 0.5,
     xlab = "Period (d)", ylab = "ln(BF)",
     main = sprintf("%s  BFP + %s  max lnBF=%.2f at P=%.2f d",
                    target, gp_label, max(lnbf), P[which.max(lnbf)]))
abline(h = 5, lty = 2, col = "red", lwd = 1.2)
abline(v = c(38.4, 170, 365.25), lty = 3, col = c("brown", "orange", "grey60"), lwd = c(1.5, 2, 1))
magaxis(side = c(1, 2), tcl = -0.3)
legend("topright", legend = c("ln(BF)=5", "38d rot", "170d", "1yr"),
       col = c("red", "brown", "orange", "grey60"), lty = c(2, 3, 3, 3), bty = "n", cex = 0.8)
# Short-period zoom
msel <- P > 1.5 & P < 200
plot(P[isort][msel[isort]], lnbf[isort][msel[isort]], type = "l", lwd = 0.7,
     xlab = "Period (d)", ylab = "ln(BF)",
     main = "Short-period zoom (1.5–200 d)")
abline(h = 5, lty = 2, col = "red")
abline(v = c(38.4, 170), lty = 3, col = c("brown", "orange"), lwd = c(1.5, 2))
magaxis(side = c(1, 2), tcl = -0.3)
dev.off()
cat("\nWrote", fpdf, "\n       ", f_sum, "\n       ", f_asc, "\n")
