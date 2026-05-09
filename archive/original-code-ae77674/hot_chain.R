tem.up <- 15
ind <- which(!is.na(data[,3] & !is.na(data[,3] & data[,3]!=0)))
if(Niter<1e5){
    tem.low <- 1e-3
}else{
    tem.low <- 1e-6
}
tem <- tem.low
tmp <- run.metropolis.MCMC(startvalue,cov.start,iterations=1e3,tem=tem,bases=bases)
acceptance <-  (1-mean(duplicated(tmp$out)))*100
if(acceptance<50) tem.low <- 1e-3*tem.low

tem.min <- min(tem.low,10^floor(log10(min(data[ind,3]/sd(data[ind,2])))))
verbose <- FALSE
#####Find the optimal tempering parameter: tem
Ntem <- round(log(1/tem.min)/log(2))#tem.min*2^n=tem; tem<=1
i1 <- 0
Ntmp <- Ntmp0 <- 1000
tems <- c()
for(i0 in 0:Ntem){
    tem <- min(tem.min*2^i0,1)
    tems <- c(tems,tem)
    nburn <- floor(Ntmp/2)
    tmp <- run.metropolis.MCMC(startvalue,cov.start,iterations=Ntmp,tem=tem,bases=bases)
    out.mcmc <- tmp$out[-(1:nburn),]
    Npar <- ncol(out.mcmc)-2
    mcmc.out <- out.mcmc[,1:Npar]
    logpost.out <- out.mcmc[,Npar+1]
    loglike.out <- out.mcmc[,Npar+2]
    acceptance <-  (1-mean(duplicated(mcmc.out)))*100
    if(Niter0<1e5 & verbose){
        cat('Adaptive tempering burning: tem=',format(tem,digit=3),'\n')
        cat('maximum likelihood=',max(loglike.out),'\n')
        cat('acceptance:',acceptance,'\n')
    }
    if(acceptance<(tem.up+5)){
        ##        cat('Finding the optimal tempering with longer chain!\n')
        Ntmp <- Ntmp0*10
        if(i1>0 & acceptance<(tem.up-5)) break()
        i1 <- i1+1
    }else{
        Ntmp <- Ntmp0
    }
#    cat('Ntmp=',Ntmp,'\n')
    startvalue <- mcmc.out[which.max(loglike.out),]
}
####Run MCMC with the optimally tempered chain
if(verbose) cat('\nRun hot chain with tem =',tem,'\n')
tmp <- run.metropolis.MCMC(startvalue,cov.start,iterations=floor(Niter/Ncores),tem=tem,bases=bases)
mcmc.out <- tmp$out[-(1:floor(Niter/Ncores/2)),]
startvalue <- mcmc.out[which.max(mcmc.out[,ncol(mcmc.out)]),1:(ncol(mcmc.out)-2)]
####Initial constraint of a signal with cold chain
if(verbose) cat('\nRun short cold chain to constrain the signal\n')
tmp <- run.metropolis.MCMC(startvalue,cov.start,iterations=floor(Niter/Ncores),tem=1,bases=bases)
out.mcmc <- tmp$out[-(1:floor(Niter/Ncores/2)),]
Npar <- ncol(out.mcmc)-2
#mcmc <- out.mcmc[,1:Npar]
logpost <- out.mcmc[,Npar+1]
loglike <- out.mcmc[,Npar+2]
par.hot <- out.mcmc[which.max(loglike),1:Npar]
###use solutions from Np=1 runs? no!
if(length(per.prim)>0 & FALSE){
    ind <- which(abs(par.hot[1]-per.prim)<0.1 & ll.prim>max(loglike))
    if(length(ind)>0){
        ind.opt <- ind[which.max(ll.prim[ind])]
        tmp <- out.mcmc$sig1[[ind.opt]]
        ind.max <- which.max(tmp[,ncol(tmp)])
        par.hot <- out.mcmc$sig1[[ind.opt]][ind.max,1:Npar]
    }
}
Ps <- exp(par.hot[grep('per',names(par.hot))])
if(verbose) cat('Hot and short-cold chains:P=',paste(round(Ps,2),collapse=','),'day; lnL=',round(max(loglike),2),'\n')
#startvalue <- par.hot
