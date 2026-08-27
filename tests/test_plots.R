source('periodograms.R')
source('periodoframe.R')
source('functions.R')
withProgress <- function(message,value,expr,...) force(expr)
incProgress <- function(...) invisible(NULL)
Nmax.plots <- 50

expect_true <- function(value,message){
    if(!isTRUE(value)){
        stop(message,call.=FALSE)
    }
}

####################################################
## A two-observable result (RV plus one activity proxy) so that the periodogram
## panels of the second observable exercise the level bookkeeping
####################################################
set.seed(9)
P <- 15.3
t <- sort(runif(80,0,300))
rv <- 4*cos(2*pi*t/P)+rnorm(80,0,0.8)
proxy <- 0.3*sin(2*pi*t/40)+rnorm(80,0,0.2)
d <- data.frame(Time=t,RV=rv,eRV=rep(0.8,80),Sindex=proxy)
renew <- TRUE; Nsamp <- 1
pp <- list(ns=c('RV','Sindex','Window Function'),ofac=1,frange=c(1/100,1/5),per.type='BFP',
           per.target='star',sequence=FALSE,Nmas=0,Nars=0,Inds=list(0),Nsig.max=1,
           per.type.seq='BFP',Niter=0,SigType='kepler',Nh=2)
out <- calc.1Dper(Nmax.plots=50,vars=c('RV','Sindex'),per.par=pp,data=list(star=d),Ncores=1)

sp <- list.single.plots(out)
expect_true(is.data.frame(sp) && nrow(sp)>=6,'list.single.plots did not enumerate the panels')
expect_true(all(c('periodogram','phase','fit','residual')%in%sp$kind),'list.single.plots is missing a panel kind')
expect_true(any(sp$ypar=='Sindex' & sp$kind=='periodogram'),'the proxy periodogram is not offered for download')

####every panel in every format must produce a non-empty file
for(k in 1:nrow(sp)) for(fmt in c('pdf','png','jpg')){
    f <- tempfile(fileext=paste0('.',fmt))
    save.single.plot(f,format=fmt,width=6,height=4.5,dpi=120,
                     plot1D.single(out,kind=sp$kind[k],ypar=sp$ypar[k],index=sp$index[k],SigType='kepler'))
    expect_true(file.exists(f) && file.info(f)$size>1000,
                paste('empty or missing',fmt,'for',sp$label[k]))
    unlink(f)
}

####the bundled multi-panel figures still draw with the new panel functions
f <- tempfile(fileext='.pdf')
pdf(f,8,8)
per1D.plot(out$per.list,out$tits,out$pers,out$levels,ylabs=out$ylabs,download=TRUE,SigType='kepler')
phase1D.plot(out$phase.list,out$sim.list,out$tits,download=TRUE,repar=FALSE)
dev.off()
expect_true(file.info(f)$size>5000,'the bundled PDF is empty')

####the periodogram panel reports the peak it annotates
####(with two harmonics the raw maximum may be the 2P alias; the annotation
####must follow the reported period)
pdf(NULL)
pk <- plot1D.single(out,'periodogram','RV',1,SigType='kepler')
dev.off()
expect_true(abs(pk$Popt-P)<0.5,paste('the annotated peak',pk$Popt,'is not the reported period',P))

####titles are human-readable, without the internal ";" coding
expect_true(!grepl(';',pretty.title(out$tits[1])) && grepl('BFP',pretty.title(out$tits[1])),
            paste('pretty.title produced',pretty.title(out$tits[1])))

cat('plot tests passed\n')
