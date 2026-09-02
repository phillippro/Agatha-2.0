# Session Handoff

Current state:

- Repository is organized and public as `phillippro/Agatha-2.0`.
- DACE RV test data are already included under `data/`.
- A backup of the pre-multi-dataset code is stored in `archive/original-code-ae77674/`.

Implemented in this session:

- Added simultaneous 1D multi-dataset periodograms for `BFP` and `MLP`.
- The implemented model is a shared circular periodic signal with one RV offset per data set and a shared linear trend.
- Per-data-set AR/MA order controls are available in the 1D UI and enter the multi-set fit as per-set lag nuisance terms.
- The previous residual-concatenation path is no longer used for 1D multi-set `BFP`/`MLP`.
- Added `tests/test_multiset_periodogram.R`.

Implemented in the following session:

- Shared Keplerian fitting for multi-set data (`KeplerFit.multiset`, `multiset_kepler_*`).
- Purely stochastic multi-set fitting (`multiset_noise_periodogram`), selected by passing
  `noise.only=TRUE` to `BFP.multiset`.
- `sigfit.multiset()` dispatches circular/Keplerian/stochastic for multi-set periodograms.
- Parallel-tempering MCMC (`run.ptmcmc` in `mcmc_func.R`) is now the default sampler used to
  constrain periodogram signals; `hot_chain.R` is kept behind `mcmc.method`.
- The PT ladder is automatic: `tem.min` from pilot draws over the prior, rung count from the
  dimension, Vousden et al. 2016 adaptation during burn-in, extension until split `Rhat<1.1`.
- Two latent bugs fixed on the way: `multiset_design()` used `model.matrix()` and failed for a
  single data set, and `KeplerFit.multiset` has to read `per$df$data` because `sigfit()`
  overwrites `per$data` with the unshifted input table.

Implemented after that (Fourier/harmonic periodograms):

- `Nh` harmonics in BFP/MLP/multi-set periodograms; the linear solves in `CircularSig` and
  `par.integral` are now general-matrix (`marginal_logL`) and regression-tested against the
  previous hand-written forms for `Nh=1`.
- Analytical Keplerian seed from the fundamental and first harmonic (Delisle et al. 2016),
  used by `KeplerFit`, `KeplerFit.multiset` and `par.a2m`.
- Known limit: for `e >~ 0.9` the harmonic series converges slowly and `V1`, `V2` from sparse
  data are biased; the seed is then flagged (`ok=FALSE`) and the numerical fit multi-starts.

Implemented after that (GP red noise):

- `noise.model` ('ARMA'|'GP') through the UI, `calc.1Dper`, BFP/MLP and the multi-set path.
- Dense SHO-kernel GLS core in `periodoframe.R`; the old celerite R port is wrong (kept but
  unused). Known cost: O(N^3) per evaluation.
- Multi-harmonic P/2P check now likelihood-based (`harmonic_period_check(P, logL1)`).
- HIP 88962 PFS (47 points): with white noise the BFP finds ~1950-2000 d, the harmonic seed and
  LM fit give e~0.6, K~60 m/s; with MA(1) the noise model absorbs part of the signal and the
  result is unstable between runs (multimodal `sopt` restarts) - worth a larger `Nsamp`.

Remaining engineering targets:

- MCMC with the GP noise model.

- Extend multi-set support to proxy terms.
- Multi-set MCMC implemented (`mcfit.multiset`); sequential multi-signal MCMC for multi-set
  fits and single-set GP MCMC remain open.
- Verify against section 2.6 of `paper/abfp/bkp/astro_periodogram2.pdf` if that paper tree is restored.

Touched implementation points:

- `mcmc_func.R`
- `functions.R`
- `periodoframe.R`
- `server.R`
- `README.md`
- `ui.R`
- `tests/test_multiset_periodogram.R`
- `tests/test_ptmcmc.R`
- `tests/test_fourier_kepler.R`
- `tests/test_gp_periodogram.R`
- `tests/test_plots.R`
- `tests/test_moving_periodogram.R`
- `tests/test_multiset_mcmc.R`
- `MP_plot.R`
- `additional_signals.R`

Tracking list:

- `docs/TODO.md`

Notes for the next session:

- Keep the original-code backup untouched unless a new backup is explicitly requested.
- Treat the current backup directory as the rollback point for the multi-dataset rewrite.
- Update the user-facing docs only after the behavior changes land.
