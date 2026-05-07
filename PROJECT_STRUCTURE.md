# Project Structure

Agatha keeps the Shiny app entry points at the repository root because Shiny automatically discovers `ui.R` and `server.R` there.

## Main Files

- `ui.R` - Shiny user interface.
- `server.R` - Shiny server logic.
- `agatha2.R` - command-line Agatha workflow.
- `functions.R`, `periodoframe.R`, `periodograms.R`, `mcmc_func.R`, `orbit.R`, `sofa.R` - core R implementation used by the app and scripts.

## Directories

- `data/` - bundled RV and activity-indicator input data. The DACE test data are documented in `data/DACE_TEST_DATA.md`.
- `results/` - default runtime output directory. Generated files are ignored by Git.
- `examples/results/` - historical example outputs kept for reference.
- `scripts/` - utility scripts, including DACE data fetching.
- `docs/` - developer notes and exploratory documents.
- `archive/backups/` - legacy editor backup files retained for traceability.
- `rsconnect/` - ShinyApps deployment metadata.

## Runtime Notes

Run the Shiny app from the repository root:

```sh
R -e "shiny::runApp('.', launch.browser = TRUE)"
```

Run the CLI workflow from the repository root:

```sh
Rscript agatha2.R BFP kepler 2 0.1 MA data HD210193_PFS.vels HD103949_PFS.vels
```
