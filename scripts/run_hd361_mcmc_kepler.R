withProgress <- function(message=NULL, value=0, expr, ...) {
  eval(substitute(expr), parent.frame())
}
incProgress <- function(amount=0, detail=NULL, ...) invisible(NULL)

suppressPackageStartupMessages({
  library(magicaxis)
  source("functions.R")
  source("periodograms.R")
  source("periodoframe.R")
  source("mcmc_func.R")
  source("sofa.R")
  source("orbit.R")
})

input_dir <- "/Users/ffeng/Documents/projects/dwarfs/agatha2/data/HD361"
out_dir <- file.path("results", "HD361_mcmc_kepler")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

files <- sort(list.files(input_dir, pattern="\\.(dat|vels)$", full.names=TRUE))
run_summary <- data.frame(
  dataset=character(),
  status=character(),
  n_data=integer(),
  p_bfp_days=numeric(),
  p_mcmc_days=numeric(),
  logbf_max=numeric(),
  llmax=numeric(),
  stringsAsFactors=FALSE
)

write_one <- function(x, path, csv=FALSE) {
  if (csv) {
    write.csv(x, path, quote=FALSE)
  } else {
    write.table(x, path, quote=FALSE, row.names=FALSE)
  }
}

for (path in files) {
  dataset <- tools::file_path_sans_ext(basename(path))
  log_path <- file.path(out_dir, paste0(dataset, "_run.log"))
  status <- "ok"
  n_data <- NA_integer_
  p_bfp <- NA_real_
  p_mcmc <- NA_real_
  logbf_max <- NA_real_
  llmax <- NA_real_

  log_con <- file(log_path, open="wt")
  sink(log_con)
  sink(log_con, type="message")
  result <- tryCatch({
    tab <- read.agatha.table(path)
    n_data <- nrow(tab)
    data <- list(tab)
    names(data) <- dataset
    per_par <- list(
      ns=colnames(tab),
      ofac=1,
      frange=c(1 / 1e5, 1 / 1.1),
      per.type="BFP",
      per.target=dataset,
      SigType="kepler",
      sequence=FALSE,
      Nsig.max=1,
      per.type.seq="BFP",
      Niter=1000,
      Nmas=list(1),
      Nars=list(0),
      Inds=list(0),
      Nsamp=1,
      renew=TRUE
    )

    Nmax.plots <<- 1
    out <- calc.1Dper(1, "RV", per_par, data, Ncores=1)
    prefix <- file.path(out_dir, paste0(dataset, "_", out$fname))

    write_one(out$per.list$RV, paste0(prefix, "_Periodogram.txt"))
    write_one(out$phase.list$RV, paste0(prefix, "_PhaseData.txt"))
    write_one(out$sim.list$RV, paste0(prefix, "_SimFit.txt"))
    write_one(out$mc.list$RV, paste0(prefix, "_MCposterior.txt"))
    write_one(out$par.list$RV, paste0(prefix, "_OptPar.csv"), csv=TRUE)
    saveRDS(out, paste0(prefix, ".rds"))

    pdf(paste0(prefix, "_Periodogram.pdf"), width=7, height=5)
    per1D.plot(out$per.list, out$tits, out$pers, out$levels,
               ylabs=out$ylabs, SigType="kepler")
    dev.off()

    pdf(paste0(prefix, "_PhaseFit.pdf"), width=7, height=8)
    phase1D.plot(out$phase.list, out$sim.list, out$tits, repar=TRUE)
    dev.off()

    per_tab <- out$per.list$RV
    power_col <- ncol(per_tab)
    p_bfp <- per_tab[which.max(per_tab[, power_col]), 1]
    logbf_max <- max(per_tab[, power_col], na.rm=TRUE)
    if (!is.null(out$mc) && length(out$mc) > 0 && "loglike" %in% colnames(out$mc)) {
      llmax <- max(out$mc[, "loglike"], na.rm=TRUE)
    }
    par_tab <- out$par.list$RV
    par_row <- if ("med" %in% rownames(par_tab)) "med" else if ("median" %in% rownames(par_tab)) "median" else NA
    if (!is.na(par_row) && !is.null(dim(par_tab)) && "P1" %in% colnames(par_tab)) {
      p_mcmc <- par_tab[par_row, "P1"]
    } else if (!is.na(par_row) && !is.null(dim(par_tab)) && "P" %in% colnames(par_tab)) {
      p_mcmc <- par_tab[par_row, "P"]
    }
    TRUE
  }, error=function(e) {
    status <<- paste("error:", conditionMessage(e))
    FALSE
  })
  sink(type="message")
  sink()
  close(log_con)

  run_summary <- rbind(run_summary, data.frame(
    dataset=dataset,
    status=status,
    n_data=n_data,
    p_bfp_days=p_bfp,
    p_mcmc_days=p_mcmc,
    logbf_max=logbf_max,
    llmax=llmax,
    stringsAsFactors=FALSE
  ))
}

write.csv(run_summary, file.path(out_dir, "summary.csv"), quote=FALSE, row.names=FALSE)
print(run_summary)
