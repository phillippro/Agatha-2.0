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

`install.packages(c('minpack.lm','magicaxis','foreach','doMC','parallel'))`

## 2. Usage

After downloading the source code and entering the `agatha2/` directory, you can run the following commandline in your terminal,

`Rscript agatha2.R BFP kepler 2 0.1 MA data HD210193_PFS.vels HD103949_PFS.vels`

>The first argument `BFP` specifies the type of periodograms used, it could be BFP or GLST and more types will be implemented soon. Note that the GLST is the most efficient periodogram but it does not account for white noise jitter nor red noise. BFP is typically recommended because it models white jitter, red nosie as well as a floating linear trend. In principle it could include a linear function of noise proxies as done in Agatha v1.0. However, this linear function may also introduce spurious signals due to nonlinear correlation between RVs and noise proxies. So the overlap between RV and activity signals might be a better way to diagnose RV signals. 

>The second argument `kepler` is the type of signal to be constrained. There are two options, `kepler` or `circular`. 

>The third argument `2` is the number of signals you want to find.

>The fourth argument `0.1` is the oversampling parameter (ofac), ofac>1 is recommended although lower sampling rate is efficient for test.

>The fifth argument `MA` specifies the noise model used to calculate BFPs. There are three options in this version, `white`, `MA` and `AR`. The Gaussian process version would be tested soon although GP is found to remove red noise as well as signals (Feng et al. 2016 and Ribas et al. 2018). 

>The sixth argument `data` is the relative path where the data files are put. If the directory is where the agatha2.R file is located, use `.` instead.

>The last arguments provide the data files to be analyzed. Agatha 2.0 only support analysis of one set per target and thus multiple data files would be treated independently. 

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

## 3. Future developement

>Develop a R markdown code to visualize the results

>Develope a Shiny app for graphic application of Agatha 2.0

>Enable Agatha to analyze multiple data sets

>Combine Agatha with MCMC to give parameter uncertainty and posterior# agatha2
