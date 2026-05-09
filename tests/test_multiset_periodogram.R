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

mlp <- MLP.multiset(t=t,y=y,dy=dy,set.id=set.id,ofac=6,fmin=1/20,fmax=1/8)
expect_true(abs(mlp$Popt-sim$P)<0.25,
            paste('MLP.multiset recovered period',mlp$Popt,'instead of',sim$P))

renew <- TRUE
Nsamp <- 1
per.par <- list(ns=colnames(sim$d1),ofac=6,frange=c(1/20,1/8),
                per.type='BFP',per.target=c('set1','set2'),
                SigType='circular',sequence=FALSE,Nmas=c(0,0),Nars=c(0,0),
                Inds=list(0,0),Nsig.max=1,per.type.seq='BFP',Niter=0)
out <- calc.1Dper(Nmax.plots=2,vars='RV',per.par=per.par,
                  data=list(set1=sim$d1,set2=sim$d2),Ncores=1)
expect_true(length(out$per.list)==1,'calc.1Dper did not return one observable')
expect_true(ncol(out$per.list$RV)==2,'calc.1Dper did not return BFP powers')
expect_true(any(grepl('combined',out$tits)),'calc.1Dper did not label the multi-set fit as combined')

cat('multi-set periodogram tests passed\n')
