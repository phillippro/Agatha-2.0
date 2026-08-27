source('periodograms.R')
source('periodoframe.R')
source('functions.R')
library(fields)
library(magicaxis)
withProgress <- function(message,value,expr,...) force(expr)
incProgress <- function(...) invisible(NULL)

expect_true <- function(value,message){
    if(!isTRUE(value)){
        stop(message,call.=FALSE)
    }
}

####################################################
## The moving (2D) periodogram across data sets, periodogram types, proxies,
## windows with too little data, and two data sets combined. Each case must
## run through per2D.data() and plotMP() (the app path) without error.
####################################################
data <- list(HIP88962_PFS=read.agatha.table('data/HIP88962_PFS.vels'),
             HD361_PFS=read.agatha.table('data/HD361/HD361_PFS.vels'),
             HD361_HARPSpre=read.agatha.table('data/HD361/HD361_HARPSpre.dat'),
             HD189567=read.agatha.table('data/HD189567_DACE_HARPS03.dat'))

mp.case <- function(targets,per.type,Dtfrac,Nbin,Nma=0,Inds=0,noise.model='ARMA',gp.Prot=NA,gp.tau=NA){
    d <- data[targets]
    tsp <- diff(range(unlist(lapply(d,function(x) x[,1]))))
    pars <- list(ns=c('RV','Window Function'),ofac=1,frange=c(1/tsp,1/2),per.type=per.type,per.target=targets,
                 files=NULL,Niter=0,Nmas=rep(Nma,length(targets)),Nars=rep(0,length(targets)),
                 Inds=rep(list(Inds),length(targets)),Dt=signif(tsp*Dtfrac,3),Nbin=Nbin,
                 alpha=5,scale=TRUE,pmin.zoom=2,pmax.zoom=tsp,show.signal=TRUE,noise.model=noise.model,gp.Prot=gp.Prot,gp.tau=gp.tau)
    v <- per2D.data(vars='RV',per.par=pars,data=d)
    f <- tempfile(fileext='.pdf')
    pdf(f,8,8); plotMP(v,pars); dev.off()
    expect_true(file.info(f)$size>2000,paste('empty 2D figure for',per.type,paste(targets,collapse='+')))
    unlink(f)
    v
}

####all periodogram types on a dense and a sparse data set
for(pt in c('GLS','BGLS','GLST','LS','MLP')){
    v <- mp.case('HD189567',pt,0.3,5)
    expect_true(ncol(v$zz)==5 && nrow(v$zz)==length(v$yy) && all(is.finite(v$zz)),
                paste(pt,'moving periodogram has the wrong shape or NA on dense data'))
    v <- mp.case('HIP88962_PFS',pt,0.4,5)
    expect_true(ncol(v$zz)==5,paste(pt,'moving periodogram has the wrong number of windows on sparse data'))
}

####LS must give real powers (it used to read a field lsp() does not return)
v <- mp.case('HD361_PFS','LS',0.5,6)
expect_true(all(is.finite(v$zz)) && diff(range(v$zz))>0,'LS moving periodogram is empty or constant')

####proxies and a moving-average term in small windows: singular fits must not abort
v <- mp.case('HD361_PFS','MLP',0.5,6,Nma=1,Inds=1)
expect_true(ncol(v$zz)==6,'MLP with proxy and MA(1) did not produce all windows')

####windows with too few points are skipped, not fatal; the others are computed
v <- mp.case('HIP88962_PFS','GLS',0.1,10)
expect_true(ncol(v$zz)==10,'tiny-window case did not return one column per window')
expect_true(any(is.na(v$zz)) && any(is.finite(v$zz)),
            'tiny-window case should have some skipped and some computed windows')

####two data sets combined
for(pt in c('GLS','MLP')){
    v <- mp.case(c('HD361_PFS','HD361_HARPSpre'),pt,0.5,6)
    expect_true(ncol(v$zz)==6 && length(v$idata)==2,paste('two-set moving periodogram failed for',pt))
}

####the moving periodogram of a clean signal peaks at the right period in every window
set.seed(2)
t <- sort(runif(200,0,400)); P <- 21.3
d1 <- data.frame(Time=t,RV=5*sin(2*pi*t/P)+rnorm(200,0,0.5),eRV=rep(0.5,200))
data$sim <- d1
v <- mp.case('sim','GLS',0.5,4)
pk <- apply(v$zz,2,function(z) v$yy[which.max(z)])
expect_true(all(abs(pk-P)<1),paste('window peaks',paste(round(pk,2),collapse=','),'do not track the injected period',P))

####GP red noise inside the windows: free hyperparameters, a fixed oscillation
####period, and two data sets whose GP is removed per set before combining
v <- mp.case('HD361_PFS','MLP',0.6,4,noise.model='GP')
expect_true(ncol(v$zz)==4 && any(is.finite(v$zz)) && grepl('_GP',v$fname),'MLP with a GP in the windows failed')
v <- mp.case('HD361_PFS','MLP',0.6,4,noise.model='GP',gp.Prot=30)
expect_true(ncol(v$zz)==4 && any(is.finite(v$zz)),'MLP with a GP and a fixed oscillation period failed')
v <- mp.case(c('HD361_PFS','HD361_HARPSpre'),'MLP',0.5,4,noise.model='GP')
expect_true(ncol(v$zz)==4 && length(v$idata)==2,'two-set moving periodogram with GP failed')
####proxies are ignored with a GP only if none were selected; with one selected they still enter
v <- mp.case('HD361_PFS','MLP',0.6,4,noise.model='GP',Inds=1)
expect_true(ncol(v$zz)==4,'MLP with GP and a proxy failed')

cat('moving periodogram tests passed\n')
