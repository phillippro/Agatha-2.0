## Custom Agatha BFP runner for alpha Cen residuals.
## Computes BFP with white / MA(1) / AR(1) noise models, compares lnBF,
## saves periodograms and top-peak tables. Avoids Circ2kep follow-up fit
## (which has a recycling bug on data with no proxies).
##
## Usage: Rscript run_bfp_alphacen.R <input.vels> <outdir> <ofac>

suppressPackageStartupMessages({
    library(magicaxis); library(foreach); library(doMC); library(parallel)
})
source('periodograms.R'); source('periodoframe.R')
source('mcmc_func.R'); source('sofa.R'); source('orbit.R')
options(warn = 1)

args <- commandArgs(trailingOnly = TRUE)
f <- args[1]
outdir <- args[2]
ofac <- as.numeric(args[3])
if (is.na(ofac)) ofac <- 0.3

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
target <- sub("\\.vels$", "", basename(f))
cat("target:", target, "  ofac:", ofac, "\n")

tab <- read.table(f, header = TRUE)
tmin <- min(tab[, 1]); if (tmin < 2400000) tmin <- tmin + 2400000
t <- tab[, 1] - tmin
y <- tab[, 2]; ey <- tab[, 3]
cat(sprintf("  N=%d  RMS=%.3f m/s  baseline=%.1f d\n",
            length(y), sqrt(mean(y^2)), diff(range(t))))

fmin <- 1 / max(diff(range(t)), 5000)
fmax <- 1 / 1.1
Nsamp <- 1

results <- list()
for (noise in c("white", "MA", "AR")) {
    cat(sprintf("\n=== Noise model: %s ===\n", noise))
    q <- 0; p <- 0
    if (noise == "MA") q <- 1
    if (noise == "AR") p <- 1
    t0 <- Sys.time()
    per <- tryCatch({
        BFP(t, y, ey, Nma = q, Nar = p, Indices = NULL, ofac = ofac,
            model.type = "man", fmin = fmin, fmax = fmax, quantify = TRUE,
            progress = FALSE, GP = FALSE, gp.par = rep(NA, 3),
            noise.only = FALSE, Nsamp = Nsamp, sampling = "combined",
            par.opt = NULL, renew = TRUE)
    }, error = function(e) {
        cat("  BFP failed:", conditionMessage(e), "\n"); NULL
    })
    if (is.null(per)) next
    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    P <- per$P
    lnbf <- if (!is.null(per$lnBFs)) per$lnBFs else per$power
    isort <- order(P)
    # Top peaks (log BF > 1, local maxima, min period sep)
    is_peak <- c(FALSE, (lnbf[2:(length(lnbf) - 1)] > lnbf[1:(length(lnbf) - 2)]) &
                        (lnbf[2:(length(lnbf) - 1)] > lnbf[3:length(lnbf)]), FALSE)
    ord <- order(lnbf, decreasing = TRUE)
    cand_p <- c(); cand_b <- c()
    for (i in ord) {
        if (!is_peak[i]) next
        if (length(cand_p) >= 15) break
        if (length(cand_p) == 0 || all(abs(log10(P[i] / cand_p)) > 0.015)) {
            cand_p <- c(cand_p, P[i]); cand_b <- c(cand_b, lnbf[i])
        }
    }
    # Global maximum log-likelihood of the full model (per$LogLike0 = null model)
    llmax <- if (!is.null(per$llmax)) per$llmax else NA
    ll0   <- if (!is.null(per$LogLike0)) per$LogLike0 else NA
    cat(sprintf("  time=%.1fs  max(lnBF)=%.2f at P=%.3f d  llmax=%.2f  ll0=%.2f\n",
                dt, max(lnbf, na.rm = TRUE), P[which.max(lnbf)], llmax, ll0))
    cat("  top peaks (P day, lnBF):\n")
    for (i in seq_along(cand_p))
        cat(sprintf("    %2d  %10.4f  %8.3f\n", i, cand_p[i], cand_b[i]))

    results[[noise]] <- list(P = P, lnbf = lnbf, cand_p = cand_p, cand_b = cand_b,
                              llmax = llmax, ll0 = ll0, time = dt)
    # Save ASCII
    f_asc <- file.path(outdir, paste0(target, "_BFP_", noise, "_periodogram.txt"))
    write.table(cbind(Period.day = P[isort], lnBF = lnbf[isort]),
                file = f_asc, row.names = FALSE, quote = FALSE)
    cat("  wrote", f_asc, "\n")
}

if (length(results) == 0) { cat("No successful runs.\n"); quit() }

# ---- Comparison PDF ----
fpdf <- file.path(outdir, paste0(target, "_BFP_noise_compare.pdf"))
pdf(fpdf, 10, 10)
par(mfrow = c(4, 1), mar = c(3, 4, 2, 1), mgp = c(2, 0.6, 0))
cols <- c(white = "black", MA = "blue", AR = "red")
ymax <- max(sapply(results, function(r) max(r$lnbf, na.rm = TRUE)), na.rm = TRUE)
for (noise in names(results)) {
    r <- results[[noise]]
    plot(r$P, r$lnbf, type = "l", log = "x", col = cols[noise], lwd = 0.6,
         xlab = "Period (d)", ylab = "ln BF",
         main = sprintf("%s BFP (%s noise)  N=%d  max lnBF=%.2f at P=%.3fd",
                        target, noise, length(y), max(r$lnbf), r$P[which.max(r$lnbf)]),
         ylim = c(min(-2, min(r$lnbf)), max(10, ymax) * 1.05))
    abline(h = 5, lty = 2, col = "grey50")
    abline(h = 0, lty = 3, col = "grey70")
    for (i in seq_along(r$cand_p))
        abline(v = r$cand_p[i], col = adjustcolor(cols[noise], 0.3), lwd = 0.5)
    magaxis(side = c(1, 2), tcl = -0.3)
    legend("topleft",
           legend = c(sprintf("llmax=%.1f", r$llmax),
                      sprintf("ll0=%.1f", r$ll0),
                      sprintf("Δll=%.1f", r$llmax - r$ll0)),
           bty = "n", cex = 0.85)
}
# Overlay panel
plot(NA, xlim = range(results[[1]]$P), ylim = c(0, ymax * 1.05),
     log = "x", xlab = "Period (d)", ylab = "ln BF",
     main = paste0(target, ": overlay of noise models"))
abline(h = 5, lty = 2, col = "grey50")
for (noise in names(results)) {
    lines(results[[noise]]$P, pmax(results[[noise]]$lnbf, 0),
          col = adjustcolor(cols[noise], 0.8), lwd = 0.6)
}
legend("topright", legend = names(results), col = cols[names(results)], lwd = 2, bty = "n")
magaxis(side = c(1, 2), tcl = -0.3)
dev.off()
cat("\nWrote", fpdf, "\n")

# ---- Noise-model summary table ----
f_sum <- file.path(outdir, paste0(target, "_BFP_noise_compare_summary.txt"))
sink(f_sum)
cat("# Agatha BFP noise-model comparison for", target, "\n")
cat("# Columns: noise  N_freq  max_lnBF  P_at_max(d)  llmax  ll0  delta_ll  time_s\n")
for (noise in names(results)) {
    r <- results[[noise]]
    cat(sprintf("%-5s  %7d  %8.3f  %12.4f  %10.2f  %10.2f  %9.2f  %7.1f\n",
                noise, length(r$P), max(r$lnbf), r$P[which.max(r$lnbf)],
                r$llmax, r$ll0, r$llmax - r$ll0, r$time))
}
cat("\n# Top 10 peaks per noise model:\n")
for (noise in names(results)) {
    r <- results[[noise]]
    cat("#", noise, "\n")
    for (i in seq_len(min(10, length(r$cand_p))))
        cat(sprintf("  %-5s  %12.4f  %8.3f\n", noise, r$cand_p[i], r$cand_b[i]))
}
sink()
cat("Wrote", f_sum, "\n")
