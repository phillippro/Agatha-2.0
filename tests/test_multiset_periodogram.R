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
