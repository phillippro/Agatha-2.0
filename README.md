# Agatha v2.0

## Advanced version of Agatha software for periodic signal diagnostics

### Compared with Agatha v1.0 ([Shiny App](https://phillippro.shinyapps.io/Agatha/) or [GitHub source code](https://github.com/phillippro/agatha)), v2.0 provides the following new features:

> 1. The autoregressive and Gaussian process noise model can be used to calculate BFP

> 2. The Keplerian fit using LM algorithm is able to constrain the circular signal identified through BFP

> 3. The phase curve for Keplerian fit and residuals visualize the fit

> 4. Metrix ln(BF) is calculated in GLST to provide better signal detection threshold (i.e. ln(BF)>5)

> 5. add transit periodogram

## 1. Installation

Install most recent R and R packages from CRAN, for example, [https://cran.r-project.org/bin/macosx/](https://cran.r-project.org/bin/macosx/).

Then install R packages in the R console using

```r
install.packages(c(
  "shiny",
  "minpack.lm",
  "magicaxis",
  "foreach",
  "doMC",
  "ramify"
))
```

`parallel` is included with base R and normally does not need to be installed separately.

## 2. Run the Agatha Shiny app

From a terminal, enter the repository directory and start the app:

```sh
cd /path/to/Agatha-2.0
R -e "shiny::runApp('.', launch.browser = TRUE)"
```

If you are already in the repository directory, this shorter command is enough:

```sh
R -e "shiny::runApp('.', launch.browser = TRUE)"
```

The app opens in a local browser window. If it prints a URL such as `http://127.0.0.1:xxxx`, open that URL manually.

### Input data format

Agatha expects a plain-text table. The first three columns are required:

1. Observation time
2. Radial velocity or other observable
3. Measurement uncertainty

Additional columns are treated as noise proxies or activity indicators. Column names are recommended. Uploaded filenames should follow the pattern `star_instrument.ext`, for example `HD210193_PFS.vels`.

Bundled example data are available in the `data/` directory and can be selected directly in the app.

The repository also includes public DACE-derived RV test sets:

- `HD40307_DACE_HARPS03`
- `GJ436_DACE_HARPS03`
- `HD189567_DACE_HARPS03`

These are HARPS RV time series fetched with `dace-query` and saved in Agatha's table format. See `data/DACE_TEST_DATA.md` for provenance and `scripts/fetch_dace_test_data.py` for the reproducible download script.

### Basic app workflow

1. Open the `Choose File` tab.
2. Keep `Upload Type` set to `Select from the list` for bundled demo data, or choose `Upload files` for your own table.
3. Select one or more data sets.
4. Click `upload and show data`.
5. Use the `Scatter Plot` tab to inspect columns and obvious outliers.
6. Use the `Model Comparison` tab to compare AR, MA, and proxy noise models.
7. Use the `1D Periodogram` tab to calculate BFP, MLP, or related periodograms.
8. Use the `2D Periodogram` tab to calculate moving periodograms and check whether a signal is stable in time.
9. Download plots or data tables from the download buttons shown in each tab.

### Demo: bundled RV data

This demo uses data already included in the repository.

1. Start the app:

```sh
cd /path/to/Agatha-2.0
R -e "shiny::runApp('.', launch.browser = TRUE)"
```

2. In `Choose File`, set `Upload Type` to `Select from the list`.
3. Select `HD210193` and click `upload and show data`.
4. Open `Scatter Plot`.
5. Select `HD210193` as the target, choose time for the x-axis and RV for the y-axis, then click `show scatter plot`.
6. Open `1D Periodogram`.
7. Select `HD210193` as the data set.
8. Select `BFP` as the periodogram type.
9. Use `Circular` as the signal type for a quick first pass, then `Keplerian` to refine the orbit or `Stochastic` to fit correlated noise with no signal.
10. Set `Number of AR components` to `0`.
11. Set `Number of MA components` to `1`.
12. Keep the default period range for a first run.
13. Set `Oversampling factor` to `0.5` for a faster demo, or `1` for a denser period grid.
14. Select `RVs` as the observable.
15. Click `Calculate periodograms`.

The app will render the periodogram in the main panel. Peaks with larger Bayes-factor power are candidate periodic signals. The downloaded data can be used to remake publication-style figures outside the app.

To demo with DACE-derived data instead, choose `HD40307_DACE_HARPS03` in `Choose File` and use the same 1D periodogram settings.

### Demo: multi-set 1D periodogram in the app

The Shiny 1D periodogram can fit multiple selected RV data sets simultaneously. For multiple data sets, `BFP` and `MLP` fit a shared signal with one fitted RV offset per data set and a shared linear trend. The AR/MA controls remain available per data set and are included as per-set lag nuisance terms in the combined fit.

Three signal types are available for multi-set fits:

- `Circular shared signal` fits a shared sinusoid, and is the fastest first pass.
- `Keplerian shared signal` refits the periodogram peak as a shared eccentric orbit (`P`, `K`, `e`, `omega`, `Mo`) while keeping the per-set offsets, the trend and the AR/MA terms. The Keplerian is written so that the model is linear in `K*cos(omega)` and `K*sin(omega)`, which leaves only `P`, `e` and `Mo` to be optimized and makes the multi-start fit robust.
- `Stochastic` (BFP only) fits a purely stochastic model with no periodic signal at all. The periodogram then scans the time scale of the exponential AR/MA kernel rather than the period of a signal, and the Bayes factor is measured against the white-noise multi-set model. At least one AR or MA component must be selected; otherwise the periodic model is fitted instead.

Proxy terms and MCMC refinement remain single-data-set features in this workflow.

1. Start the app.
2. In `Choose File`, select two or more RV data sets and click `upload and show data`.
3. Open `1D Periodogram`.
4. Select the same data sets under `Data sets`.
5. Choose `BFP` or `MLP`.
6. Choose `Signal type`: `Circular shared signal`, `Keplerian shared signal`, or `Stochastic` (BFP with at least one AR/MA component).
7. Click `Calculate periodograms`.

### Demo: command-line run

The repository also includes a command-line workflow. This example runs a BFP search with a Keplerian signal model, two signals, oversampling factor `0.1`, and an MA noise model:

```sh
Rscript agatha2.R BFP kepler 2 0.1 MA data HD210193_PFS.vels HD103949_PFS.vels
```

The command writes periodograms, phase plots, fitted parameters, residuals, and tabular plotting data to `results/`.

## 3. Command-line usage

After downloading the source code and entering the repository directory, you can run the following command line in your terminal:

`Rscript agatha2.R BFP kepler 2 0.1 MA data HD210193_PFS.vels HD103949_PFS.vels`

>The first argument `BFP` specifies the type of periodograms used, it could be BFP or GLST and more types will be implemented soon. Note that the GLST is the most efficient periodogram but it does not account for white noise jitter nor red noise. BFP is typically recommended because it models white jitter, red nosie as well as a floating linear trend. In principle it could include a linear function of noise proxies as done in Agatha v1.0. However, this linear function may also introduce spurious signals due to nonlinear correlation between RVs and noise proxies. So the overlap between RV and activity signals might be a better way to diagnose RV signals. 

>The second argument `kepler` is the type of signal to be constrained. There are two options, `kepler` or `circular`. 

>The third argument `2` is the number of signals you want to find.

>The fourth argument `0.1` is the oversampling parameter (ofac), ofac>1 is recommended although lower sampling rate is efficient for test.

>The fifth argument `MA` specifies the noise model used to calculate BFPs. There are three options in this version, `white`, `MA` and `AR`. The Gaussian process version would be tested soon although GP is found to remove red noise as well as signals (Feng et al. 2016 and Ribas et al. 2018). 

>The sixth argument `data` is the relative path where the data files are put. If the directory is where the agatha2.R file is located, use `.` instead.

>The last arguments provide the data files to be analyzed. The command-line runner still treats multiple data files independently. Use the Shiny 1D Periodogram tab for simultaneous multi-set BFP/MLP fitting with shared circular, Keplerian or purely stochastic model parameters and one offset per data set.

By running the above commandline, the output would be

```
target: HD210193 

 results/HD210193_BFP_MA_periodogram_sig1.pdf 
results/HD210193_BFP_MA_periodogram_sig1.txt 

 results/HD210193_BFP_MA_phase_sig1.pdf 
results/HD210193_BFP_MA_phase_sig1_OptPar.txt 
results/HD210193_BFP_MA_phase_sig1_DataPhase.txt 
results/HD210193_BFP_MA_phase_sig1_SimPhase.txt 
results/HD210193_BFP_MA_phase_sig1_RawRes.txt 
results/HD210193_BFP_MA_phase_sig1_bin.txt 

 results/HD210193_BFP_MA_periodogram_sig2.pdf 
results/HD210193_BFP_MA_periodogram_sig2.txt 

 results/HD210193_BFP_MA_phase_sig2.pdf 
results/HD210193_BFP_MA_phase_sig2_OptPar.txt 
results/HD210193_BFP_MA_phase_sig2_DataPhase.txt 
results/HD210193_BFP_MA_phase_sig2_SimPhase.txt 
results/HD210193_BFP_MA_phase_sig2_RawRes.txt 
results/HD210193_BFP_MA_phase_sig2_bin.txt 

 results/HD210193_BFP_fit_allsig.pdf 
results/HD210193_BFP_fit_allsig_data.txt 
results/HD210193_BFP_fit_allsig_fit.txt 
results/HD210193_BFP_fit_allsig_RawRes.txt 
results/HD210193_BFP_fit_allsig_BinRes.txt 

 results/HD210193_BFP_MA_periodogram_res.pdf 
results/HD210193_BFP_MA_periodogram_res.txt 

 results/HD210193_BFP_MA_periodogram_Sindex.pdf 
results/HD210193_BFP_MA_periodogram_Sindex.txt 

 results/HD210193_BFP_MA_periodogram_Halpha.pdf 
results/HD210193_BFP_MA_periodogram_Halpha.txt 

 results/HD210193_BFP_MA_periodogram_PhotonCount.pdf 
results/HD210193_BFP_MA_periodogram_PhotonCount.txt 

 results/HD210193_BFP_MA_periodogram_ObsTime.pdf 
results/HD210193_BFP_MA_periodogram_ObsTime.txt 

 results/HD210193_BFP_MA_periodogram_window.pdf 
results/HD210193_BFP_MA_periodogram_window.txt 

target: HD103949 
...
```

These files provide you plots and relevant data which store the x, y and probably ey values for each plot. The meaning of file names are as follows:

>results/HD210193_BFP_MA_periodogram_sig1 - plot and data for periodogram calculated using BFP+MA(1) for the first signal

>results/HD210193_BFP_MA_periodogram_sig2 - plot and data for periodogram for the second signal

>results/HD210193_BFP_MA_phase_sig1 - phase plot and data for the first signal

>results/HD210193_BFP_MA_phase_sig2 - phase plot and data for the second signal

>results/HD210193_BFP_MA_phase_allsig - phase plot and data for all signals

>results/HD210193_BFP_MA_periodogram_res - residual periodogram

>results/HD210193_BFP_MA_periodogram_xxx - periodograms for columns 4,5,... in the data file. These columns store activity indices in the case of RV set. 

## 4. Eccentric orbits in BFP and MLP: harmonic (Fourier) fitting

BFP and MLP originally fitted a single sinusoid at each trial period. An eccentric Keplerian
moves power out of the fundamental into its harmonics, so a sinusoidal periodogram loses
sensitivity as the eccentricity grows (the power in the fundamental drops to zero as
`e -> 1`; Delisle & Segransan 2022, Fig. 2). Following Delisle, Segransan, Buchschacher &
Alesina (2016, A&A 590, A134), the periodograms can now fit `Nh` harmonics
`cos(k n t), sin(k n t)`, `k = 1..Nh`, at each trial period:

- `Number of harmonics of the signal` in the 1D UI (default 2); `Nh` argument of `BFP()`,
  `MLP()`, `BFP.multiset()` and `MLP.multiset()` (default 1, the previous behaviour).
- The harmonic amplitudes enter the same weighted least-squares solve as the offset, trend,
  proxies and AR/MA terms, so the red-noise treatment is unchanged. The BIC penalty of BFP
  counts `2*Nh` amplitudes.
- With `Nh >= 2`, the fitted fundamental and first harmonic (`V1`, `V2`) give `K`, `e`,
  `omega` and `M0` analytically (Paper I, eqs. 20-34, with the Newton-Raphson refinement of
  eqs. 38-40). That solution seeds the Keplerian fit (`KeplerFit`, `KeplerFit.multiset`) and
  the MCMC start, replacing the previous blind random starts. When the first harmonic is
  stronger than any Keplerian allows (Paper I, Sect. 4), the seed is flagged and the numerical
  fit multi-starts over eccentricity instead.
- A multi-harmonic fit is ambiguous between `P` and `2P` (a sinusoid at `P` fits the `2P`
  model through its first harmonic). At the peak, if the first harmonic dominates the
  fundamental the period is halved and refitted.

On a simulated `e = 0.8` orbit the two-harmonic BFP gains about 40 per cent in
`Delta chi^2` at the true period compared with the sinusoidal one, and nothing at `e = 0`.
The analytical estimate is exact on exact Fourier coefficients; on real data its accuracy is
limited by how well `V1`, `V2` can be measured from sparse sampling, which is why it is used
as a seed rather than as the final answer.

## 5. Red noise models: ARMA or Gaussian process

`Red noise model` in the 1D UI (`noise.model` in `calc.1Dper()`) selects how correlated noise
is treated in BFP and MLP:

- `ARMA` (default): the moving-average / auto-regressive terms of the original Agatha, with
  per-data-set orders.
- `GP`: a Gaussian process with the SHO (stochastically driven damped harmonic oscillator)
  kernel of Foreman-Mackey et al. (2017), parametrized by `sigmaGP` (S0), `logProt` and
  `logtauGP` (`Q = tauGP*pi/Prot`). For several data sets one GP is shared, because the
  activity belongs to the star, with one offset per data set and a shared trend.

How the GP enters each computation:

- Single data set, BFP: the hyperparameters are free at every trial period and fitted together
  with the signal by `sopt()`, exactly as the ARMA parameters are. MLP: the GP is fitted on
  the signal-free model, its conditional mean is subtracted, and the marginalized periodogram
  is computed on the residual (`MLP.type='sub'`).
- Several data sets (`multiset_gp_periodogram`): for BFP the hyperparameters are refitted
  together with the harmonics (`gp.fit='joint'`), so that the noise and the signal compete at
  each frequency exactly as in the single-set BFP and in the ARMA case. Because the
  hyperparameters vary smoothly with frequency, the refit is done on a coarse grid of about
  5 per cent of the trial frequencies (warm-started along the grid) and interpolated in
  between, where only the generalized least squares is solved; the ten highest peaks and the
  reported peak are then refitted exactly. On a 130-point test this takes 29 s at `ofac=2`
  against 218 s for a refit at every frequency, with identical peak values (lnBF 40.0 at
  the true period, against 12.4 for the fixed-hyperparameter scan of the same data). Hyperparameters fitted on the signal-free model
  absorb signal power and suppress the peak, which is why the fixed-hyperparameter
  (GP-whitened) scan is kept only as an option (`gp.fit='fixed'`) and as the MLP default,
  where the noise is fixed before marginalizing by construction. The Keplerian refit uses the
  covariance found at the peak.
- `Stochastic` signal type with `GP`: no periodic signal; the rotation period of the kernel is
  scanned and `sigma`, `tauGP` are refitted at each trial period against the white-noise
  baseline, so the periodogram shows the evidence for quasi-periodic red noise.
- Phase plots and residuals have the GP conditional mean removed.
- When `GP` is selected the panel offers `GP rotation period` and `GP coherence time scale`.
  Leave a field empty to fit that hyperparameter; enter a value (e.g. a photometric rotation
  period) to hold it fixed. The amplitude is always fitted. In the API this is
  `gp.par=c(sigmaGP, logProt, logtauGP)` with `NA` for free, on `BFP()`, `MLP()`,
  `BFP.multiset()` and `MLP.multiset()`, and `gp.Prot`/`gp.tau` (in days) in `calc.1Dper()`.
  With the `Stochastic` signal type the rotation period is what is scanned, so a fixed value
  is ignored there (with a warning).

The GP likelihood is computed with a dense Cholesky solve (`gp_sho_cov`, `gp_gls`,
`gp_res`), validated against an independent kernel implementation and a direct solve. The
celerite R port that an earlier, unreachable GP path relied on disagrees with a direct solve
and is no longer used. The dense solve costs `O(N^3)` per likelihood evaluation, which is
immediate for a few hundred points; for several thousand points prefer `ARMA` or a coarser
`ofac`.

MCMC sampling is not yet available with the GP noise model; the maximum-likelihood fit is
returned with a warning.

## 6. MCMC sampling

When `MCMC sample size` is set to a non-zero value, the signal found by the periodogram is
constrained with a parallel-tempering (replica-exchange) MCMC. A ladder of chains samples
`loglike * tem + logprior` for tempering parameters `tem` running geometrically from 1 (the
cold chain) down to `tem.min`. Neighbouring chains periodically propose to exchange states,
with an acceptance ratio that depends only on the likelihoods and the tempering parameters
because the priors cancel. The hot chains see a flattened likelihood and therefore roam over
period aliases and over eccentricity, and replica exchange carries whatever they find down to
the cold chain. Only the cold chain is returned, so the reported posterior is the untempered
one. Each replica adapts both its proposal covariance and a step-size scale towards a 23 per
cent acceptance rate, so no per-parameter tuning is required.

The ladder is set automatically, following Vousden, Farr & Mandel (2016, MNRAS 455, 1919):

- the hottest temperature `tem.min` is chosen from pilot draws over the prior box so that the
  tempered log-likelihood varies by order unity across the prior (the hottest replica then
  samples close to the prior and can reach any alias or eccentricity);
- the number of rungs follows from the dimension of the problem and `tem.min`, using the
  geometric spacing of their Table 1;
- during burn-in the interior rungs move so that neighbouring swap rates become equal (their
  eq. 12, with a rate decaying as `t0/(t+t0)`), and the ladder is frozen for the sampling
  phase so the cold chain is a valid stationary sample;
- if the cold chain has not converged after the requested length (split Gelman-Rubin
  `Rhat > 1.1`), sampling continues in blocks of the same length, up to three extra blocks.

The behaviour is controlled by the arguments of `mcfit()` and `sigfit()`:

- `mcmc.method` - `'PT'` (default) for parallel tempering, or anything else to fall back to
  the older sequential adaptive-tempering scheme in `hot_chain.R`.
- `Ntem`, `tem.min` - `NULL` (default) for the automatic choices above, or fixed values.
- `swap.interval` - iterations between replica-exchange sweeps (default 10).
- `mcmc.verbose` - print the chosen and adapted ladder, per-replica acceptance rates, swap
  rates, step scales and the final `Rhat`.

`run.ptmcmc()` also accepts `adapt.ladder`, `adapt.window`, `adapt.nu`, `adapt.t0`,
`max.extend` and `Rhat.max` for finer control.

## 7. Roadmap

>Develop a R markdown code to visualize the results

>Develope a Shiny app for graphic application of Agatha 2.0

>Enable Agatha to analyze multiple data sets

>Combine Agatha with MCMC to give parameter uncertainty and posterior# agatha2
