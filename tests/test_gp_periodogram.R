source('periodograms.R')
source('periodoframe.R')
source('functions.R')
withProgress <- function(message,value,expr,...) force(expr)
incProgress <- function(...) invisible(NULL)

expect_true <- function(value,message){
    if(!isTRUE(value)){
        stop(message,call.=FALSE)
    }
}

####################################################
## Dense SHO covariance against the independent kernel in mcmc_func.R, and the
## Cholesky likelihood against a direct solve/determinant. (The celerite R port
## previously used by the GP path disagrees with a direct solve by ~15 in lnL
## on this test and is no longer used.)
####################################################
source('mcmc_func.R')
set.seed(11)
t <- sort(runif(60,0,100))
dy <- runif(60,0.3,0.6)
y <- rnorm(60,0,2)
sigmaGP <- 1.7; logProt <- log(12); logtauGP <- log(30)
for(ltau in c(logtauGP,log(1))){   # Q>1/2 (underdamped) and Q<1/2 (overdamped) branches
    Q <- exp(ltau)*pi/exp(logProt)
    K <- gp_sho_cov(t,sigmaGP,logProt,ltau,dy=dy)
    K.ref <- outer(t,t,function(a,b) ker.sho(a,b,sigmaGP,log(Q),logProt))
    diag(K.ref) <- diag(K.ref)+dy^2
    expect_true(max(abs(K-K.ref))<1e-10,paste('gp_sho_cov differs from ker.sho at Q=',Q))
    R <- gp_chol(K)
    lnL.chol <- -0.5*sum(backsolve(R,y,transpose=TRUE)^2)-sum(log(diag(R)))-length(y)/2*log(2*pi)
    lnL.ref <- -0.5*as.numeric(t(y)%*%solve(K.ref,y))-0.5*as.numeric(determinant(K.ref,logarithm=TRUE)$modulus)-length(y)/2*log(2*pi)
    expect_true(abs(lnL.chol-lnL.ref)<1e-8*abs(lnL.ref),paste('Cholesky GP likelihood differs from the direct solve at Q=',Q))
    expect_true(abs(sum(gp_res(t,y,dy,0,sigmaGP,logProt,ltau)^2)-length(t)*off.gp+lnL.ref)<1e-8*abs(lnL.ref),
                'gp_res does not reproduce -logL up to the fixed offset')
}

####################################################
## Multi-set GP periodogram: a circular signal buried in SHO noise shared by
## two data sets. The GP-whitened periodogram must find it; the white one
## is the reference it must beat.
####################################################
set.seed(5)
P <- 23.7; A <- 2.0; B <- -1.2
t1 <- sort(runif(70,0,300)); t2 <- sort(runif(60,40,340))
tt <- c(t1,t2); sid <- c(rep('set1',70),rep('set2',60))
o <- order(tt); tt <- tt[o]; sid <- sid[o]
dy <- rep(0.5,length(tt))
Kgp <- gp_sho_cov(tt,sigmaGP=30,logProt=log(60),logtauGP=log(150),dy=NULL)   # GP sd ~ 5 m/s (S0 w0 Q)
red <- as.numeric(t(chol(Kgp+diag(1e-8,length(tt))))%*%rnorm(length(tt)))
sig <- A*cos(2*pi*tt/P)+B*sin(2*pi*tt/P)
y <- sig+red+ifelse(sid=='set1',10,-5)+rnorm(length(tt),0,dy)
expect_true(sd(red)>1.5*sqrt(A^2+B^2)/sqrt(2),'the simulated red noise is not dominant; the test would not be informative')

gp <- BFP.multiset(t=tt,y=y,dy=dy,set.id=sid,ofac=2,fmin=1/120,fmax=1/8,noise.model='GP')
expect_true(identical(gp$gp$fit,'joint'),'BFP.multiset with GP did not fit the hyperparameters jointly with the signal')
####the hyperparameters must vary along the scan when fitted jointly
expect_true(sd(gp$pars[,'logtauGP'])>0,'joint GP fit did not refit the hyperparameters across trial periods')
gpf <- BFP.multiset(t=tt,y=y,dy=dy,set.id=sid,ofac=2,fmin=1/120,fmax=1/8,noise.model='GP',gp.fit='fixed')
expect_true(abs(gpf$Popt-P)<0.5,paste('fixed-hyperparameter GP periodogram found',gpf$Popt,'instead of',P))
white <- BFP.multiset(t=tt,y=y,dy=dy,set.id=sid,ofac=2,fmin=1/120,fmax=1/8)
expect_true(identical(gp$noise.model,'GP') && all(c('sigmaGP','logProt','logtauGP')%in%names(gp$par.opt)),
            'BFP.multiset did not return a GP fit')
expect_true(abs(gp$Popt-P)<0.5,paste('GP multi-set periodogram found',gp$Popt,'instead of',P))
####lnBF values are not comparable across noise models (the white baseline is
####far worse), so compare each periodogram's peak contrast against its own
####background: the true period must stand out at least as well after whitening
contrast <- function(per){ i <- which.min(abs(per$P-P)); (per$power[i]-median(per$power))/mad(per$power) }
expect_true(contrast(gp)>=0.8*contrast(white),
            paste('the true period stands out less in the GP periodogram (',round(contrast(gp),1),') than in the white one (',round(contrast(white),1),')'))
expect_true(sd(gp$res)<0.6*sd(white$res),'the GP fit did not absorb the red noise')
expect_true(abs(gp$gp$logProt-log(60))<log(2.5),paste('recovered GP rotation period',exp(gp$gp$logProt),'far from 60'))

####Keplerian fit under the same GP covariance, seeded by the harmonics
gp2 <- BFP.multiset(t=tt,y=y,dy=dy,set.id=sid,ofac=2,fmin=1/120,fmax=1/8,noise.model='GP',Nh=2)
expect_true(abs(gp2$Popt-P)<0.5,paste('GP multi-set periodogram with Nh=2 found',gp2$Popt,'instead of',P))
kf <- KeplerFit.multiset(gp2)
expect_true(abs(kf$ParKep$P1-P)<0.5 && kf$ParKep$e1<0.25,
            paste('GP Keplerian multi-set fit gave P=',kf$ParKep$P1,'e=',kf$ParKep$e1,'for a circular signal'))
expect_true(abs(kf$ParKep$K1-sqrt(A^2+B^2))<0.6,paste('GP Keplerian fit gave K=',kf$ParKep$K1))

####################################################
## Stochastic GP periodogram: no signal, scan the rotation period
####################################################
yn <- red+ifelse(sid=='set1',10,-5)+rnorm(length(tt),0,dy)
st <- BFP.multiset(t=tt,y=yn,dy=dy,set.id=sid,ofac=2,fmin=1/200,fmax=1/10,noise.only=TRUE,noise.model='GP')
expect_true(isTRUE(st$noise_only) && all(st$ysig==0),'stochastic GP periodogram contains a signal')
expect_true(max(st$power)>0,'stochastic GP periodogram found no evidence for red noise')
expect_true(abs(log(st$Popt/60))<log(2.5),paste('stochastic GP periodogram recovered Prot=',st$Popt,'instead of 60'))

####################################################
## Through calc.1Dper: multi-set GP with circular and Keplerian signal types,
## and the single-set path with GP
####################################################
renew <- TRUE; Nsamp <- 1
mk <- function(sel) data.frame(Time=tt[sel],RV=y[sel],eRV=dy[sel])
d2 <- list(set1=mk(sid=='set1'),set2=mk(sid=='set2'))
pp <- list(ns=c('RV','Window Function'),ofac=2,frange=c(1/120,1/8),per.type='BFP',
           per.target=c('set1','set2'),sequence=FALSE,Nmas=c(0,0),Nars=c(0,0),
           Inds=list(0,0),Nsig.max=1,per.type.seq='BFP',Niter=0,SigType='circular',Nh=1,noise.model='GP')
out <- calc.1Dper(Nmax.plots=2,vars='RV',per.par=pp,data=d2,Ncores=1)
expect_true(abs(out$par.list$RV[['P1']]-P)<0.5,'calc.1Dper multi-set GP circular fit missed the period')
expect_true(all(c('sigmaGP','logProt','logtauGP')%in%names(out$par.list$RV)),'calc.1Dper did not report the GP hyperparameters')
pp$SigType <- 'kepler'; pp$Nh <- 2
outk <- calc.1Dper(Nmax.plots=2,vars='RV',per.par=pp,data=d2,Ncores=1)
expect_true(all(c('P1','K1','e1')%in%names(outk$par.list$RV)) && abs(outk$par.list$RV[['P1']]-P)<0.5,
            'calc.1Dper multi-set GP Keplerian fit failed')

####single data set, GP noise in BFP (celerite path) with circular and Keplerian signal
d1 <- list(set1=mk(sid=='set1'))
pp1 <- pp; pp1$per.target <- 'set1'; pp1$Nmas <- 0; pp1$Nars <- 0; pp1$Inds <- list(0); pp1$SigType <- 'circular'; pp1$Nh <- 1; pp1$ofac <- 1
out1 <- calc.1Dper(Nmax.plots=2,vars='RV',per.par=pp1,data=d1,Ncores=1)
expect_true(all(c('sigmaGP','logProt','logtauGP')%in%names(out1$par.list$RV)),'single-set BFP with GP did not return the GP hyperparameters')
expect_true(abs(out1$par.list$RV[['P1']]-P)<1,paste('single-set GP BFP found',out1$par.list$RV[['P1']],'instead of',P))
pp1$SigType <- 'kepler'
out1k <- calc.1Dper(Nmax.plots=2,vars='RV',per.par=pp1,data=d1,Ncores=1)
expect_true(abs(out1k$par.list$RV[['P1']]-P)<1 && out1k$par.list$RV[['e1']]<0.3,
            paste('single-set GP Keplerian fit gave P=',out1k$par.list$RV[['P1']],'e=',out1k$par.list$RV[['e1']]))

####################################################
## Fixed time scales (e.g. a photometric rotation period): the fixed value must
## be used as is, the others still fitted, in both the multi-set and single-set paths
####################################################
gpx <- BFP.multiset(t=tt,y=y,dy=dy,set.id=sid,ofac=2,fmin=1/120,fmax=1/8,noise.model='GP',gp.par=c(NA,log(60),NA))
expect_true(abs(gpx$par.opt[['logProt']]-log(60))<1e-12,'multi-set GP did not hold the rotation period fixed')
expect_true(abs(gpx$Popt-P)<0.5,paste('multi-set GP with fixed Prot found',gpx$Popt,'instead of',P))
expect_true(all(abs(gpx$pars[,'logProt']-log(60))<1e-12),'the fixed rotation period drifted along the scan')
gpxx <- BFP.multiset(t=tt,y=y,dy=dy,set.id=sid,ofac=2,fmin=1/120,fmax=1/8,noise.model='GP',gp.par=c(NA,log(60),log(150)))
expect_true(abs(gpxx$par.opt[['logtauGP']]-log(150))<1e-12 && abs(gpxx$Popt-P)<0.5,'multi-set GP with both time scales fixed failed')
ppx <- pp; ppx$SigType <- 'circular'; ppx$Nh <- 1; ppx$gp.Prot <- 60; ppx$gp.tau <- NA
outx <- calc.1Dper(Nmax.plots=2,vars='RV',per.par=ppx,data=d2,Ncores=1)
expect_true(abs(outx$par.list$RV[['logProt']]-log(60))<1e-12,'calc.1Dper did not pass the fixed rotation period to the multi-set GP')
pp1x <- pp1; pp1x$SigType <- 'circular'; pp1x$gp.Prot <- 60; pp1x$gp.tau <- NA
out1x <- calc.1Dper(Nmax.plots=2,vars='RV',per.par=pp1x,data=d1,Ncores=1)
expect_true(!('logProt'%in%names(out1x$par.list$RV)) && 'sigmaGP'%in%names(out1x$par.list$RV),
            'single-set BFP with a fixed rotation period still fitted it (or lost the GP amplitude)')
expect_true(abs(out1x$par.list$RV[['P1']]-P)<1,paste('single-set GP BFP with fixed Prot found',out1x$par.list$RV[['P1']]))

cat('GP periodogram tests passed\n')
