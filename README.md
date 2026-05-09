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
9. Use `Circular` as the signal type for a quick first pass.
10. Set `Number of AR components` to `0`.
11. Set `Number of MA components` to `1`.
12. Keep the default period range for a first run.
13. Set `Oversampling factor` to `0.5` for a faster demo, or `1` for a denser period grid.
14. Select `RVs` as the observable.
15. Click `Calculate periodograms`.

The app will render the periodogram in the main panel. Peaks with larger Bayes-factor power are candidate periodic signals. The downloaded data can be used to remake publication-style figures outside the app.

To demo with DACE-derived data instead, choose `HD40307_DACE_HARPS03` in `Choose File` and use the same 1D periodogram settings.

### Demo: multi-set 1D periodogram in the app

The Shiny 1D periodogram can fit multiple selected RV data sets simultaneously. For multiple data sets, `BFP` and `MLP` use a shared circular signal with one fitted RV offset per data set. AR/MA noise, proxy terms, Keplerian refinement, and MCMC refinement remain single-data-set features in this workflow.

1. Start the app.
2. In `Choose File`, select two or more RV data sets and click `upload and show data`.
3. Open `1D Periodogram`.
4. Select the same data sets under `Data sets`.
5. Choose `BFP` or `MLP`.
6. Keep `Signal type` as `Circular shared signal`.
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

>The last arguments provide the data files to be analyzed. The command-line runner still treats multiple data files independently. Use the Shiny 1D Periodogram tab for simultaneous multi-set BFP/MLP fitting with shared circular signal parameters and one offset per data set.

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

## 4. Roadmap

>Develop a R markdown code to visualize the results

>Develope a Shiny app for graphic application of Agatha 2.0

>Enable Agatha to analyze multiple data sets

>Combine Agatha with MCMC to give parameter uncertainty and posterior# agatha2
