source('periodograms.R')
source('periodoframe.R')
source('functions.R')

set.seed(42)

make_data <- function(){
    P <- 12.34
    A <- 2.2
    B <- -1.4
    t1 <- sort(runif(45,0,120))
    t2 <- sort(runif(50,10,130))
    dy1 <- rep(0.35,length(t1))
    dy2 <- rep(0.40,length(t2))
    y1 <- 10+A*cos(2*pi*t1/P)+B*sin(2*pi*t1/P)+rnorm(length(t1),0,dy1)
    y2 <- -5+A*cos(2*pi*t2/P)+B*sin(2*pi*t2/P)+rnorm(length(t2),0,dy2)
    list(P=P,A=A,B=B,
         d1=data.frame(time=t1,RV=y1,eRV=dy1),
         d2=data.frame(time=t2,RV=y2,eRV=dy2))
}

expect_true <- function(value,message){
    if(!isTRUE(value)){
        stop(message,call.=FALSE)
    }
}

sim <- make_data()
t <- c(sim$d1$time,sim$d2$time)
y <- c(sim$d1$RV,sim$d2$RV)
dy <- c(sim$d1$eRV,sim$d2$eRV)
set.id <- c(rep('set1',nrow(sim$d1)),rep('set2',nrow(sim$d2)))

fit <- BFP.multiset(t=t,y=y,dy=dy,set.id=set.id,ofac=6,fmin=1/20,fmax=1/8)
expect_true(abs(fit$Popt-sim$P)<0.25,
            paste('BFP.multiset recovered period',fit$Popt,'instead of',sim$P))
expect_true(all(c('A','B','gamma_set1','gamma_set2')%in%names(fit$par.opt)),
            'BFP.multiset did not return shared amplitudes plus per-set offsets')
expect_true(abs((fit$par.opt['gamma_set1']-fit$par.opt['gamma_set2'])-15)<1,
            'BFP.multiset did not recover the expected offset separation')
expect_true(sd(fit$res)<0.7,'BFP.multiset residual scatter is too large')

fit.arma <- BFP.multiset(t=t,y=y,dy=dy,set.id=set.id,Nma=c(1,0),Nar=c(0,1),
                         ofac=6,fmin=1/20,fmax=1/8)
expect_true(abs(fit.arma$Popt-sim$P)<0.25,
            paste('BFP.multiset with ARMA orders recovered period',fit.arma$Popt,'instead of',sim$P))
expect_true(identical(as.integer(fit.arma$Nma),c(1L,0L)) && identical(as.integer(fit.arma$Nar),c(0L,1L)),
            'BFP.multiset did not preserve selected per-set ARMA orders')
expect_true(any(grepl('ma_set1_lag1|ar_set2_lag1',colnames(fit.arma$pars))),
            'BFP.multiset did not include selected ARMA lag terms')

mlp <- MLP.multiset(t=t,y=y,dy=dy,set.id=set.id,ofac=6,fmin=1/20,fmax=1/8)
expect_true(abs(mlp$Popt-sim$P)<0.25,
            paste('MLP.multiset recovered period',mlp$Popt,'instead of',sim$P))

renew <- TRUE
Nsamp <- 1
per.par <- list(ns=colnames(sim$d1),ofac=6,frange=c(1/20,1/8),
                per.type='BFP',per.target=c('set1','set2'),
                SigType='circular',sequence=FALSE,Nmas=c(1,0),Nars=c(0,1),
                Inds=list(0,0),Nsig.max=1,per.type.seq='BFP',Niter=0)
out <- calc.1Dper(Nmax.plots=2,vars='RV',per.par=per.par,
                  data=list(set1=sim$d1,set2=sim$d2),Ncores=1)
expect_true(length(out$per.list)==1,'calc.1Dper did not return one observable')
expect_true(ncol(out$per.list$RV)==2,'calc.1Dper did not return BFP powers')
expect_true(any(grepl('combined',out$tits)),'calc.1Dper did not label the multi-set fit as combined')
expect_true(all(c('ma_set1_lag1','ar_set2_lag1')%in%names(out$par.list$RV)),
            'calc.1Dper did not return ARMA nuisance parameters for multi-set BFP')

####################################################
## Keplerian fitting of a multi-data-set time series
####################################################
make_kepler_data <- function(P,K,e,omega,Mo){
    t1 <- sort(runif(60,0,150))
    t2 <- sort(runif(65,20,170))
    tt <- c(t1,t2)
    sid <- c(rep('set1',length(t1)),rep('set2',length(t2)))
    ord <- order(tt)
    tt <- tt[ord]
    sid <- sid[ord]
    M <- (Mo+2*pi*(tt-min(tt))/P)%%(2*pi)
    E <- kep.mt2(M,e)
    nu <- 2*atan(sqrt((1+e)/(1-e))*tan(E/2))
    sig <- K*(cos(omega+nu)+e*cos(omega))
    dy <- rep(0.4,length(tt))
    list(t=tt,set.id=sid,dy=dy,sig=sig,
         y=sig+ifelse(sid=='set1',10,-5)+rnorm(length(tt),0,dy))
}

kep <- list(P=15.3,K=4,e=0.45,omega=1.1,Mo=2.3)
sim.kep <- make_kepler_data(kep$P,kep$K,kep$e,kep$omega,kep$Mo)
per.kep <- BFP.multiset(t=sim.kep$t,y=sim.kep$y,dy=sim.kep$dy,set.id=sim.kep$set.id,
                        ofac=4,fmin=1/25,fmax=1/10)
fit.kep <- KeplerFit.multiset(per.kep)
expect_true(abs(fit.kep$ParKep$P1-kep$P)<0.1,
            paste('KeplerFit.multiset recovered period',fit.kep$ParKep$P1,'instead of',kep$P))
expect_true(abs(fit.kep$ParKep$K1-kep$K)<0.5,
            paste('KeplerFit.multiset recovered semi-amplitude',fit.kep$ParKep$K1,'instead of',kep$K))
expect_true(abs(fit.kep$ParKep$e1-kep$e)<0.1,
            paste('KeplerFit.multiset recovered eccentricity',fit.kep$ParKep$e1,'instead of',kep$e))
expect_true(abs(fit.kep$ParKep$omega1-kep$omega)<0.3,
            paste('KeplerFit.multiset recovered omega',fit.kep$ParKep$omega1,'instead of',kep$omega))
expect_true(all(c('gamma_set1','gamma_set2')%in%names(fit.kep$ParKep)),
            'KeplerFit.multiset did not return one offset per data set')
expect_true(abs((fit.kep$ParKep$gamma_set1-fit.kep$ParKep$gamma_set2)-15)<1,
            'KeplerFit.multiset did not recover the expected offset separation')
expect_true(cor(fit.kep$ysig,sim.kep$sig)>0.99,
            'KeplerFit.multiset did not reproduce the injected Keplerian signal')
####an eccentric orbit must beat the circular fit that seeded it
expect_true(sd(fit.kep$res)<0.5*sd(per.kep$res),
            'KeplerFit.multiset did not improve on the circular multi-set fit')
####the curve used for the phase plot has to match the fitted signal
expect_true(max(abs(multiset_kepler_curve(fit.kep$ParKep,sim.kep$t-min(sim.kep$t))-fit.kep$ysig))<1e-8,
            'multiset_kepler_curve does not reproduce the fitted Keplerian signal')

####################################################
## Purely stochastic (signal-free) multi-data-set fitting
####################################################
make_red_data <- function(tau,phi=0.8){
    t1 <- sort(runif(70,0,200))
    t2 <- sort(runif(70,10,210))
    tt <- c(t1,t2)
    sid <- c(rep('set1',length(t1)),rep('set2',length(t2)))
    ord <- order(tt)
    tt <- tt[ord]
    sid <- sid[ord]
    dy <- rep(0.5,length(tt))
    red <- rep(0,length(tt))
    for(s in unique(sid)){
        r <- which(sid==s)
        w <- rnorm(length(r),0,2)
        v <- w
        for(i in 2:length(r)){
            v[i] <- w[i]+phi*exp(-abs(tt[r[i]]-tt[r[i-1]])/tau)*w[i-1]
        }
        red[r] <- v
    }
    list(t=tt,set.id=sid,dy=dy,
         y=red+ifelse(sid=='set1',3,-2)+rnorm(length(tt),0,dy))
}

tau.true <- 8
sim.red <- make_red_data(tau.true)
noise <- BFP.multiset(t=sim.red$t,y=sim.red$y,dy=sim.red$dy,set.id=sim.red$set.id,
                      Nma=c(1,1),ofac=4,fmin=1/60,fmax=1/2,noise.only=TRUE)
expect_true(isTRUE(noise$noise_only),
            'BFP.multiset did not return a purely stochastic periodogram')
expect_true(all(noise$ysig==0),
            'the purely stochastic multi-set fit contains a periodic signal')
####the MA time scale is only loosely constrained by ~140 irregularly sampled
####points: over repeated realizations the estimator is unbiased but scatters by
####about 45 per cent, so this bound is roughly three sigma rather than tight
expect_true(abs(log(noise$Popt/tau.true))<log(3),
            paste('the purely stochastic fit recovered a time scale of',noise$Popt,
                  'instead of',tau.true))
expect_true(max(noise$power)>0,
            'the purely stochastic fit found no evidence for correlated noise')
expect_true(all(c('logtau','logtauAR')%in%names(noise$par.opt)),
            'the purely stochastic fit did not report the fitted time scales')
expect_true(sd(noise$res)<sd(noise$base.fit$res),
            'the purely stochastic fit did not reduce the white-noise residual')

####without any AR or MA component there is no stochastic model to fit
flat <- withCallingHandlers(
    BFP.multiset(t=sim.red$t,y=sim.red$y,dy=sim.red$dy,set.id=sim.red$set.id,
                 ofac=2,fmin=1/60,fmax=1/2,noise.only=TRUE),
    warning=function(w) invokeRestart('muffleWarning'))
expect_true(!isTRUE(flat$noise_only),
            'BFP.multiset returned a stochastic fit with no AR or MA component selected')

####################################################
## Keplerian and stochastic multi-set fits through calc.1Dper
####################################################
kep.par <- list(ns=c('RV','Window Function'),ofac=4,frange=c(1/25,1/10),
                per.type='BFP',per.target=c('set1','set2'),sequence=FALSE,
                Nmas=c(0,0),Nars=c(0,0),Inds=list(0,0),Nsig.max=1,
                per.type.seq='BFP',Niter=0,SigType='kepler')
kep.data <- list(set1=data.frame(Time=sim.kep$t[sim.kep$set.id=='set1'],
                                 RV=sim.kep$y[sim.kep$set.id=='set1'],
                                 eRV=sim.kep$dy[sim.kep$set.id=='set1']),
                 set2=data.frame(Time=sim.kep$t[sim.kep$set.id=='set2'],
                                 RV=sim.kep$y[sim.kep$set.id=='set2'],
                                 eRV=sim.kep$dy[sim.kep$set.id=='set2']))
out.kep <- calc.1Dper(Nmax.plots=2,vars='RV',per.par=kep.par,data=kep.data,Ncores=1)
pk <- out.kep$par.list$RV
expect_true(all(c('P1','K1','e1','omega1','Mo1')%in%names(pk)),
            'calc.1Dper did not return Keplerian parameters for a multi-set fit')
expect_true(abs(pk[['e1']]-kep$e)<0.1,
            paste('calc.1Dper recovered eccentricity',pk[['e1']],'instead of',kep$e))
####the simulated curve must line up with the phase-folded signal it is drawn over
ph <- out.kep$phase.list$RV
si <- out.kep$sim.list$RV
curve <- approxfun(si[,'tsim_sig1'],si[,'ysim_sig1'],rule=2)
expect_true(max(abs(curve(ph[,'t_sig1'])-ph[,'ysig_sig1']))<1e-2*sd(ph[,'ysig_sig1']),
            'the phase-folded Keplerian curve does not match the fitted signal')

red.par <- kep.par
red.par$SigType <- 'stochastic'
red.par$Nmas <- c(1,1)
red.par$frange <- c(1/60,1/2)
red.data <- list(set1=data.frame(Time=sim.red$t[sim.red$set.id=='set1'],
                                 RV=sim.red$y[sim.red$set.id=='set1'],
                                 eRV=sim.red$dy[sim.red$set.id=='set1']),
                 set2=data.frame(Time=sim.red$t[sim.red$set.id=='set2'],
                                 RV=sim.red$y[sim.red$set.id=='set2'],
                                 eRV=sim.red$dy[sim.red$set.id=='set2']))
out.red <- calc.1Dper(Nmax.plots=2,vars='RV',per.par=red.par,data=red.data,Ncores=1)
pr <- out.red$par.list$RV
expect_true('tau'%in%names(pr),
            'calc.1Dper did not report the stochastic time scale for a multi-set fit')
expect_true(!any(c('A1','B1','K1')%in%names(pr)),
            'the purely stochastic multi-set fit reported signal parameters')

tmp.two.col <- tempfile(fileext='.dat')
write.table(data.frame(time=c(3,1,2),value=c(9,7,8)),tmp.two.col,
            row.names=FALSE,quote=FALSE)
parsed <- read.agatha.table(tmp.two.col,center.rv=FALSE)
expect_true(identical(colnames(parsed),c('Time','RV','eRV')),
            'read.agatha.table did not normalize two-column upload names')
expect_true(all(parsed$Time==c(1,2,3)),
            'read.agatha.table did not sort uploaded data by time')
expect_true(all(parsed$eRV==1),
            'read.agatha.table did not add unit errors for two-column uploads')

tmp.headerless <- tempfile(fileext='.dat')
writeLines(c('3 9 0.2','1 7 0.2','2 8 0.2'),tmp.headerless)
parsed.headerless <- read.agatha.table(tmp.headerless,center.rv=FALSE)
expect_true(identical(colnames(parsed.headerless),c('Time','RV','eRV')),
            'read.agatha.table did not normalize headerless upload names')

cat('multi-set periodogram tests passed\n')
