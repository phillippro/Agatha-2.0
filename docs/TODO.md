# TODO

Completed in the multi-set implementation:

1. Implemented simultaneous multi-dataset RV support in the 1D periodogram path.
1. Extended selected data sets to fit shared circular signal parameters with one offset per data set.
1. Updated `BFP` and `MLP` handling so multi-dataset inputs use the combined model instead of residual concatenation.
1. Updated `server.R` so the 1D UI exposes multi-set `BFP`/`MLP` and keeps per-data-set AR/MA order controls available.
1. Added per-data-set AR/MA lag nuisance terms to the multi-set `BFP`/`MLP` fit.
1. Added `tests/test_multiset_periodogram.R` for period recovery, offset recovery, and `calc.1Dper()` integration.
1. Refreshed `README.md` with the app workflow and command-line limitation.

Completed in the Keplerian, stochastic and PT-MCMC implementation:

1. Added a shared Keplerian model for multi-set fits (`KeplerFit.multiset`), fitted with the
   linear-in-`K cos w`/`K sin w` parameterization so only `P`, `e` and `Mo` are optimized.
1. Added a purely stochastic multi-set model (`multiset_noise_periodogram`) that scans the
   AR/MA kernel time scale instead of a signal period, against a white-noise baseline.
1. Gave `multiset_lag_terms` an exponential `exp(-|dt|/tau)` kernel, defaulting to the previous
   plain lags so the signal periodogram is unchanged.
1. Replaced the multi-set branch of `sigfit()` with `sigfit.multiset()`, which dispatches on
   circular, Keplerian and stochastic signal types.
1. Exposed all three signal types for multi-set fits in `server.R`.
1. Replaced the sequential adaptive tempering of `hot_chain.R` with parallel tempering
   (`run.ptmcmc`) as the default MCMC, keeping the old path under `mcmc.method`.
1. Added `tests/test_ptmcmc.R` and extended `tests/test_multiset_periodogram.R`.
1. Made the PT ladder automatic (Vousden et al. 2016): `tem.min` from pilot prior draws,
   rung count from the dimension, rungs adapted during burn-in, convergence-driven extension.

Completed in the Fourier/harmonic implementation (Delisle et al. 2016):

1. `BFP`, `MLP` and the multi-set periodograms take `Nh` harmonics; `CircularSig` and
   `par.integral` now build their normal equations / marginal likelihood from a general
   design matrix (`marginal_logL`), which reproduces the previous closed forms for `Nh=1`.
1. Analytical Keplerian elements from `V1`, `V2` (`fourier_to_kepler`, Hansen coefficients by
   eq. A.6) seed `KeplerFit`, `KeplerFit.multiset` and the MCMC start.
1. P/2P disambiguation at the peak (`harmonic_period_check`).
1. Added `tests/test_fourier_kepler.R`.

Completed in the GP implementation:

1. `noise.model='GP'` (shared SHO kernel) for single- and multi-set BFP/MLP, the stochastic
   periodogram (scan of `Prot`) and the Keplerian refits, on a dense Cholesky solver
   (`gp_sho_cov`, `gp_gls`, `gp_res`, `gp_fit_hyper`, `multiset_gp_periodogram`).
1. The single-set GP likelihood no longer uses the celerite R port, which disagrees with a
   direct solve by ~15 in lnL on a 60-point test and whose `sho.term`/`celerite` were being
   overwritten by `mcmc_func.R` in the app anyway.
1. P/2P disambiguation of multi-harmonic peaks now compares single-sinusoid likelihoods at
   P and P/2 instead of harmonic amplitudes (amplitudes misled on sparse data).
1. Multi-set BFP fits the GP hyperparameters jointly with the signal per frequency
   (`gp.fit='joint'`); the fixed-hyperparameter scan remains for MLP / as an option.
1. GP time scales can be fixed from the panel (`gp.Prot`, `gp.tau`) or left free.
1. Added `tests/test_gp_periodogram.R`.

Completed for figures:

1. Publication-style single-panel plotting layer (`panel.*`, `plot1D.single`, `save.single.plot`)
   used by the screen, the bundled PDF and per-figure downloads in PDF/PNG/JPG with size and dpi.
1. Added `tests/test_plots.R`.

Completed for the moving (2D) periodogram, after testing it across six data sets and all
periodogram types:

1. `LS` read a field `lsp()` does not return (empty powers, plot crash) - fixed.
1. Windows with too few points for the model are skipped (blank column) instead of aborting;
   every window's periodogram is wrapped so one failure cannot kill the run; a clear error is
   raised only if no window has enough data.
1. Singular normal equations (small window + MA + proxies) no longer leak a `try-error` into
   arithmetic (`solve.try` returns zeros -> finite poor likelihood).
1. `combine.data` (two or more data sets in the 2D tab) called `global.notation`/`sopt` with an
   obsolete argument order - fixed; proxies are scaled per set as in the 1D path.
1. `MP_plot.R` tolerates NA columns and a single detected peak; the 2D oversampling slider can
   no longer be 0.
1. The 2D tab has the ARMA/GP noise-model choice with fixed or free GP time scales.
1. Added `tests/test_moving_periodogram.R`.

Remaining follow-up work:

1. MCMC with the GP noise model (the `mcfit` data wrapper is ARMA-only).

1. Extend the multi-set model to proxy terms.
1. Extend MCMC refinement to multi-set fits; `sigfit.multiset()` still warns and returns the
   weighted linear fit when MCMC is requested.
1. Add a public bundled-data demo with documented expected output.
1. Verify against the unavailable draft model in `paper/abfp/bkp/astro_periodogram2.pdf` if that paper tree is restored.

Backup and safety:

- Keep `archive/original-code-ae77674/` unchanged as the rollback snapshot.
- Create a fresh backup before any large second-stage rewrite if the model changes further.
