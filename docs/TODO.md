# TODO

Completed in the multi-set implementation:

1. Implemented simultaneous multi-dataset RV support in the 1D periodogram path.
1. Extended selected data sets to fit shared circular signal parameters with one offset per data set.
1. Updated `BFP` and `MLP` handling so multi-dataset inputs use the combined model instead of residual concatenation.
1. Updated `server.R` so the 1D UI exposes multi-set `BFP`/`MLP` and locks unsupported multi-set options to safe values.
1. Added `tests/test_multiset_periodogram.R` for period recovery, offset recovery, and `calc.1Dper()` integration.
1. Refreshed `README.md` with the app workflow and command-line limitation.

Remaining follow-up work:

1. Extend the multi-set model beyond white noise to per-data-set AR/MA and proxy terms.
1. Add Keplerian and MCMC refinement for multi-set fits.
1. Add a public bundled-data demo with documented expected output.
1. Verify against the unavailable draft model in `paper/abfp/bkp/astro_periodogram2.pdf` if that paper tree is restored.

Backup and safety:

- Keep `archive/original-code-ae77674/` unchanged as the rollback snapshot.
- Create a fresh backup before any large second-stage rewrite if the model changes further.
