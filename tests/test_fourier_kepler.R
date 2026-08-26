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
## Hansen coefficients against the power series of Delisle et al. 2016, eqs. 10-13
####################################################
e <- 0.1
expect_true(abs(hansen_X(1,e)-(1-e^2))<1e-4,'X_1(e) does not match its series expansion')
expect_true(abs(hansen_X(2,e)-(e-5/4*e^3))<1e-4,'X_2(e) does not match its series expansion')
expect_true(abs(hansen_X(-1,e)-(-e^2/8+e^4/48))<1e-6,'X_-1(e) does not match its series expansion')
expect_true(abs(hansen_X(-2,e)-(-e^3/12+e^5/48))<1e-6,'X_-2(e) does not match its series expansion')

####################################################
## Analytical inversion on exact Fourier coefficients (eqs. 24, 25, 33, 34 + Newton)
####################################################
for(e in c(0.1,0.5,0.8,0.95)) for(w in c(0,45,90,135)*pi/180){
    V <- fourier_V(8,e,w,1.3)
    f <- fourier_to_kepler(V[1],V[2])
    dw <- ((f$omega-w+pi)%%(2*pi)-pi)
    expect_true(abs(f$e-e)<2e-3 && abs(f$K/8-1)<0.05 && abs(dw)<0.02,
                paste('fourier_to_kepler failed on exact coefficients at e=',e,'omega=',w))
}

####################################################
## Single-set BFP and MLP with two harmonics on an eccentric orbit
####################################################
kepler_rv <- function(t,P,K,e,omega,M0){
    M <- (M0+2*pi*t/P)%%(2*pi)
    E <- kep.mt2(M,e)
    nu <- 2*atan(sqrt((1+e)/(1-e))*tan(E/2))
    K*(cos(omega+nu)+e*cos(omega))
}
set.seed(3)
P <- 37.3; K <- 6; e <- 0.8; omega <- 135*pi/180; M0 <- 1.3
t <- sort(runif(120,0,600))
dy <- rep(1,length(t))
y <- kepler_rv(t,P,K,e,omega,M0)+3+rnorm(length(t),0,dy)

b1 <- BFP(t=t,y=y,dy=dy,Nma=0,Nar=0,model.type='man',ofac=2,fmin=1/80,fmax=1/15,quantify=FALSE,renew=TRUE,progress=FALSE,Nh=1)
b2 <- BFP(t=t,y=y,dy=dy,Nma=0,Nar=0,model.type='man',ofac=2,fmin=1/80,fmax=1/15,quantify=FALSE,renew=TRUE,progress=FALSE,Nh=2)
expect_true(all(c('A','B','A2','B2')%in%names(b2$par.opt)),
            'BFP with Nh=2 did not return the first-harmonic amplitudes')
expect_true(abs(b2$Popt[1]-P)<0.5,paste('BFP with Nh=2 found period',b2$Popt[1],'instead of',P))
####the harmonic periodogram must recover the power an eccentric orbit hides from a sinusoid
i1 <- which.min(abs(b1$P-P)); i2 <- which.min(abs(b2$P-P))
expect_true(b2$logLs[i2]-b2$LogLike0>1.2*(b1$logLs[i1]-b1$LogLike0),
            'the two-harmonic BFP did not gain likelihood at the true period of an e=0.8 orbit')

m2 <- MLP(t=t,y=y,dy=dy,Nma=0,Nar=0,ofac=2,fmin=1/80,fmax=1/15,MLP.type='sub',Nh=2)
expect_true(abs(m2$Popt[1]-P)<0.5,paste('MLP with Nh=2 found period',m2$Popt[1],'instead of',P))
expect_true(all(c('A2','B2')%in%names(m2$par.opt)),'MLP with Nh=2 did not return the first-harmonic amplitudes')

####analytical seed from the periodogram fit, then the numerical Keplerian fit
seed <- fourier_kepler_seed(b2$par.opt,b2$Popt[1])
expect_true(!is.null(seed) && abs(seed$e1-e)<0.15,
            paste('the Fourier seed gave e=',if(is.null(seed)) NA else seed$e1,'instead of',e))
fit <- sigfit(per=b2,data=cbind(t,y,dy),SigType='kepler',mcf=FALSE)
pk <- fit$ParSig
expect_true(abs(pk[['e1']]-e)<0.05,paste('KeplerFit seeded by the Fourier solution gave e=',pk[['e1']],'instead of',e))
expect_true(abs(pk[['K1']]/K-1)<0.15,paste('KeplerFit seeded by the Fourier solution gave K=',pk[['K1']],'instead of',K))

####a pure sinusoid must not be reported at twice its period by the harmonic fit
ys <- 3*cos(2*pi*t/P)+rnorm(length(t),0,dy)
bs <- BFP(t=t,y=ys,dy=dy,Nma=0,Nar=0,model.type='man',ofac=2,fmin=1/100,fmax=1/15,quantify=FALSE,renew=TRUE,progress=FALSE,Nh=2)
expect_true(abs(bs$Popt[1]-P)<0.5,paste('BFP with Nh=2 reported a sinusoid at',bs$Popt[1],'instead of',P))

####################################################
## Multi-set: eccentric orbit shared by two data sets, through calc.1Dper
####################################################
t1 <- sort(runif(60,0,600)); t2 <- sort(runif(65,20,620))
mk <- function(tt,off) data.frame(Time=tt,RV=kepler_rv(tt,P,K,e,omega,M0)+off+rnorm(length(tt),0,1),eRV=rep(1,length(tt)))
renew <- TRUE; Nsamp <- 1
pp <- list(ns=c('RV','Window Function'),ofac=2,frange=c(1/80,1/15),per.type='BFP',
           per.target=c('s1','s2'),sequence=FALSE,Nmas=c(0,0),Nars=c(0,0),
           Inds=list(0,0),Nsig.max=1,per.type.seq='BFP',Niter=0,SigType='kepler',Nh=2)
out <- calc.1Dper(Nmax.plots=2,vars='RV',per.par=pp,data=list(s1=mk(t1,10),s2=mk(t2,-5)),Ncores=1)
pm <- out$par.list$RV
expect_true(abs(pm[['e1']]-e)<0.05,paste('multi-set Keplerian fit with Nh=2 gave e=',pm[['e1']],'instead of',e))
expect_true(abs(pm[['P1']]-P)<0.5,paste('multi-set Keplerian fit with Nh=2 gave P=',pm[['P1']],'instead of',P))

cat('Fourier-Keplerian tests passed\n')
