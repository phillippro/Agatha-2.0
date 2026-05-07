ui.R
IO interface

server.R
call functions but also do rendering

functions.R, periodograms.R, periodoframe.R
functions to be called by server.R

demo:
all peroidogram data will be stored in output$per, the relation between different objects and functions are as follows:
output$per -> per1D.plot() -> per1D.data() -> calc.1Dper() -> BFP() or MLP() or ...

First, get parameters from the unser input on the app using
periodogram.var: user choose which columns or time series data to be analyzed, including RVs, and various noise proxies
This function will return vars which will be used by calc.1Dper()

Second, get data from the input file or archived data file, the eventReactive() function will handle this. The data object will be used by calc.1Dper().

Third, get periodogram-related parameters from the user input. This will be save to the "per.par" object and will be used by calc.1Dper(). These parameters inlucde the periodogram type and orders, e.g. per.type=MA and Nma=3 for MA(q=3).
Task: change the input to enable per.type=AR and Nar from 0 to 10. 

Fourth, get the maximum number of signals to be found, i.e. Nmax.plots, which will be called by cala.1Dper().

After getting all parameters, calc.1Dper() will calculate peroidograms and then return data to per1D.data() which will be used by per1D.plot() to generate plots and other data. All these results will be stored in output$per. 

