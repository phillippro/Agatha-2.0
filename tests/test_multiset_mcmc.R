source('periodograms.R')
source('periodoframe.R')
source('functions.R')
source('mcmc_func.R')
withProgress <- function(message,value,expr,...) force(expr)
incProgress <- function(...) invisible(NULL)

expect_true <- function(value,message){
    if(!isTRUE(value)){
        stop(message,call.=FALSE)
    }
}

set.seed(4)
P <- 19.4; K <- 6; e <- 0.3; omega <- 1.2; Mo <- 0.7
kepler_rv <- function(t){
    M <- (Mo+2*pi*t/P)%%(2*pi)
    E <- kep.mt2(M,e)
    nu <- 2*atan(sqrt((1+e)/(1-e))*tan(E/2))
    K*(cos(omega+nu)+e*cos(omega))
}

####################################################
## Two data sets, Keplerian signal, MCMC refinement through calc.1Dper
####################################################
t1 <- sort(runif(60,0,300)); t2 <- sort(runif(55,30,330))
mk <- function(tt,off) data.frame(Time=tt,RV=kepler_rv(tt)+off+rnorm(length(tt),0,1),eRV=rep(1,length(tt)))
d2 <- list(s1=mk(t1,12),s2=mk(t2,-7))
renew <- TRUE; Nsamp <- 1
pp <- list(ns=c('RV','Window Function'),ofac=2,frange=c(1/40,1/8),per.type='BFP',
           per.target=c('s1','s2'),sequence=FALSE,Nmas=c(0,0),Nars=c(0,0),Inds=list(0,0),
           Nsig.max=1,per.type.seq='BFP',Niter=2000,SigType='kepler',Nh=2,noise.model='ARMA')
out <- calc.1Dper(Nmax.plots=2,vars='RV',per.par=pp,data=d2,Ncores=1)
ps <- out$par.list$RV
expect_true(is.matrix(ps),'multi-set MCMC did not return a parameter table')
expect_true(all(c('P1','K1','e1','omega1','Mo1','gamma_s1','gamma_s2','beta','sj_s1','sj_s2')%in%colnames(ps)),
            paste('multi-set MCMC parameter table misses columns; has:',paste(colnames(ps),collapse=',')))
expect_true(abs(ps['mode','P1']-P)<0.5,paste('multi-set MCMC gave P=',ps['mode','P1'],'instead of',P))
expect_true(abs(ps['mode','e1']-e)<0.2,paste('multi-set MCMC gave e=',ps['mode','e1'],'instead of',e))
expect_true(abs(ps['mode','K1']-K)<1.5,paste('multi-set MCMC gave K=',ps['mode','K1'],'instead of',K))
expect_true(abs((ps['mode','gamma_s1']-ps['mode','gamma_s2'])-19)<3,
            'multi-set MCMC did not recover the offset separation')
expect_true(nrow(out$mc.list$RV)>500,'the multi-set MCMC returned too few samples')
####uncertainties must be present and sane
expect_true(ps['sd','P1']>0 && ps['sd','P1']<1,'the period uncertainty is missing or absurd')

####################################################
## Five data sets, one of them tiny, circular signal
####################################################
offs <- c(10,-5,3,-8,0)
d5 <- list()
set.seed(6)
for(j in 1:5){
    n <- c(40,35,30,3,25)[j]
    tt <- sort(runif(n,0,300))
    d5[[paste0('set',j)]] <- data.frame(Time=tt,RV=5*cos(2*pi*tt/P)+offs[j]+rnorm(n,0,1),eRV=rep(1,n))
}
pp5 <- pp
pp5$per.target <- names(d5); pp5$Nmas <- rep(0,5); pp5$Nars <- rep(0,5); pp5$Inds <- rep(list(0),5)
pp5$SigType <- 'circular'; pp5$Nh <- 1; pp5$Niter <- 1000
out5 <- calc.1Dper(Nmax.plots=2,vars='RV',per.par=pp5,data=d5,Ncores=1)
ps5 <- out5$par.list$RV
expect_true(all(c('P1','A1','B1',paste0('gamma_set',1:5),paste0('sj_set',1:5))%in%colnames(ps5)),
            'five-set circular MCMC misses parameters')
expect_true(abs(ps5['mode','P1']-P)<0.5,paste('five-set circular MCMC gave P=',ps5['mode','P1']))
expect_true(max(abs(ps5['mode',paste0('gamma_set',1:5)]-offs))<3,
            'five-set circular MCMC did not recover the offsets')

####################################################
## GP noise model: MCMC on the whitened data, no jitter parameters
####################################################
ppg <- pp; ppg$noise.model <- 'GP'; ppg$Niter <- 1000; ppg$Nh <- 2
outg <- calc.1Dper(Nmax.plots=2,vars='RV',per.par=ppg,data=d2,Ncores=1)
psg <- outg$par.list$RV
expect_true(all(c('P1','K1','e1')%in%colnames(psg)) && !any(grepl('^sj_',colnames(psg))),
            'multi-set GP MCMC should have no per-set jitters (the GP covariance is fixed)')
expect_true(abs(psg['mode','P1']-P)<0.5,paste('multi-set GP MCMC gave P=',psg['mode','P1']))

####################################################
## Sequential signals with several data sets and MCMC: each signal is refined
## by its own multi-set MCMC and the combined model is assembled from the
## per-signal fits (this exact combination used to die on colnames(NULL))
####################################################
pps <- pp
pps$sequence <- TRUE; pps$Nsig.max <- 2; pps$per.type.seq <- 'BFP'; pps$Niter <- 1000
outs <- calc.1Dper(Nmax.plots=50,vars='RV',per.par=pps,data=d2,Ncores=1)
phs <- outs$phase.list$RV
expect_true(all(c('ysig_sig1','ysig_sig2','y_all','res_all')%in%colnames(phs)),
            'sequential multi-set MCMC did not return the per-signal and combined columns')
expect_true(all(is.finite(phs[,'res_all'])),'sequential multi-set MCMC produced non-finite residuals')
expect_true(ncol(outs$per.list$RV)>=3,'the second signal periodogram is missing')

cat('multi-set MCMC tests passed\n')
