####This file is an example for making PeriodoFrame: the computation part
inds <- 1:2
leg.pos <- 'topright'
############################
####find additional signals
############################
#save(list=ls(all=TRUE),file='test3.Robj')
###a further signal is searched only while the previous one was accepted by
###the model comparison (ln(BF) > lnBF.min); the periodogram of the first
###rejected signal is still reported. With a single data set the accepted
###signals are refined together by a joint MCMC after each addition, and the
###ln(BF) compares the joint models (BIC estimate from the maximum likelihood of
###the MCMC); with several data sets each signal has its own MCMC in the
###residual and the ln(BF) comes from the residual periodogram
joint <- mcf1 && !multi.set
if(!exists('lnBF.min')) lnBF.min <- 5
for(jj in 2:Nsig.max){
    if(!is.null(model.comp) && !all(model.comp$accepted)) break
    cat('\n Find signal ',jj,'!\n')
###the next signal is searched in the data minus the previous signal as fitted
###by sigfit(): with MCMC that is the signal of the maximum-likelihood sample
###(res.s = data - MAP signal), otherwise the periodogram's own optimum
    if(any(names(rv.ls)=='res.s')){
        res <- as.numeric(rv.ls$res.s)
    }else{
        res <- as.numeric(rv.ls$res)
    }
    if(is.matrix(res)){
        rr <- res[1,]
    }else{
        rr <- res
    }
    if(per.type.seq=='BFP'){
        rv.ls <- BFP(t,rr,dy,Nma=Nma,Nar=Nar,Indices=Indices,ofac=ofac,model.type='man',fmin=frange[1],fmax=frange[2],quantify=quantify,renew=renew,progress=progress,noise.only=noise.only,Nh=Nh,GP=GP,gp.par=gp.par)
        ylab <- 'ln(BF)'
        name <- 'logBF'
    }
    if(per.type.seq=='MLP'){
        rv.ls <- MLP(t,rr,dy,Nma=Nma,Nar=Nar,ofac=ofac,mar.type='part',model.type='man',fmin=frange[1],fmax=frange[2],opt.par=NULL,Indices=Indices,MLP.type=MLP.type,noise.only=noise.only,Nh=Nh,GP=GP,gp.par=gp.par)
        ylab <- expression('log(ML/'*ML[max]*')')
        name <- 'logML'
    }
    if(per.type.seq=='GLS'){
        rv.ls <- gls(t,rr,dy,ofac=ofac,fmin=frange[1],fmax=frange[2])
        name <- ylab <- 'Power'
    }
    if(per.type.seq=='BGLS'){
        rv.ls <- bgls(t,rr,dy,ofac=ofac,fmin=frange[1],fmax=frange[2])
        ylab <- expression('log(ML/'*ML[max]*')')
        name <- 'logML'
    }
    if(per.type.seq=='GLST'){
        rv.ls <- glst(t,rr,dy,ofac=ofac,fmin=frange[1],fmax=frange[2])
        name <- ylab <- 'Power'
    }
    if(per.type.seq=='LS'){
        rv.ls <- lsp(times=t,x=rr,ofac=ofac,from=NULL,to=frange[2],alpha=c(0.1,0.01,0.001))
        name <- ylab <- 'Power'
    }
    ylim <- c(min(rv.ls$power),max(rv.ls$power)+0.15*(max(rv.ls$power)-min(rv.ls$power)))
####store data
    if(per.type=='BFP'){
        yy  <- rv.ls$power
    }else if(per.type=='MLP'){
        yy  <- rv.ls$power-max(rv.ls$power)
        rv.ls$sig.level <- NULL#max(yy)-log(c(10,100,1000))
    }else if(per.type=='BGLS'){
        yy  <- rv.ls$power-max(rv.ls$power)
        rv.ls$sig.level <- NULL#max(yy)-log(c(10,100,1000))
    }else{
        yy <- rv.ls$power
    }
    per.data <- cbind(per.data,yy)
    Pmaxs <- c(Pmaxs,format(per.data[which.max(yy),1],digit=1))
#    tit <- paste('Periodogram: BGLS; Target:',instrument,'; Observable',ypar)
    tit <- paste0(per.type.seq,';',instrument,';',ypar,';',jj,' signal')
    f <-  paste0(paste(per.target,collapse='_'),'_',gsub(' ','',ypar),'_',per.type,'_MA',Nma,'proxy',paste(Inds,collapse='.'),'_1sig_',format(rv.ls$P[which.max(rv.ls$power)],digit=1),'d')
    tits <- c(tits,tit)
    fs <- c(fs,f)
    if(length(rv.ls$sig.level)<3){
        sig.levels <- cbind(sig.levels,c(rv.ls$sig.level,rep(NA,3-length(rv.ls$sig.level))))
    }else{
        sig.levels <- cbind(sig.levels,rv.ls$sig.level)
    }
    cnames <- c(cnames,paste0(per.type.seq,jj,'signal:',ypar,':',name))
    ylabs <- c(ylabs,ylab)
    lnBF.per <- if(per.type.seq=='BFP') suppressWarnings(max(rv.ls$power,na.rm=TRUE)) else NA
    if(!is.finite(lnBF.per)) lnBF.per <- NA
###modify the periodogram output for Keplerian fit
    tmp <- sigfit(per=rv.ls,data=cbind(t,rr,dy),SigType=SigType1,basis=basis,mcf=mcf1,Niter=Niter,Ncores=Ncores,Pconv=Pconv,
                  e.prior=if(exists('e.prior')) e.prior else NULL)
    rv.ls <- tmp$per
#    if(!progress) phase.plot(tmp,fold=fold)
    par.new <- addpar(par.data,tmp$ParSig,jj)
    lnBF <- lnBF.per; lnBF.source <- 'periodogram'
    fitj <- NULL
    if(joint){
###joint MCMC of all signals so far, started from the previous joint solution
###plus the new signal; ln(BF) = difference of the maximum log likelihoods
###minus (extra parameters / 2) ln N
        fitj <- mcfit(rv.ls,data=tab[,1:3],tsim=fit$tsim0,Niter=Niter,SigType=SigType,basis=basis,ParSig=par.P2per(par.new),Pconv=TRUE,Ncores=Ncores,
                      e.prior=if(exists('e.prior')) e.prior else NULL)
        if(is.finite(llmax.prev) && is.finite(fitj$llmax)){
            dk <- length(fitj$ParSig)-length(par.prev)
            lnBF <- as.numeric(fitj$llmax-llmax.prev-dk/2*log(length(t)))
            lnBF.source <- 'MCMC'
        }
    }
    accepted <- is.na(lnBF) || lnBF>lnBF.min
    if(!is.null(model.comp)){
        model.comp <- rbind(model.comp,data.frame(n=jj,P=as.numeric(tmp$popt)[1],lnBF=lnBF,source=lnBF.source,accepted=accepted,stringsAsFactors=FALSE))
    }
    cat(' Signal ',jj,': ln(BF) = ',format(lnBF,digits=3),' (',lnBF.source,'); ',if(accepted) 'accepted' else 'rejected','\n',sep='')
    if(!accepted) break

    pp <- cbind(tmp$t,tmp$y,tmp$ysig0)
    colnames(pp) <- paste0(c('t','y','ysig'),'_sig',jj)
    phase.data <- cbind(phase.data,pp)

    qq <- cbind(tmp$tsim,tmp$ysim,tmp$ysim0)
    colnames(qq) <- paste0(c('tsim','ysim','ysim0'),'_sig',jj)
    sim.data <- cbind(sim.data,qq)
    par.data <- par.new
    res.last <- tmp$res
    if(exists('par.stat.data')) par.stat.data <- addpar.stat(par.stat.data,tmp$par.stat,jj,keep=names(tmp$ParSig))
    if(!is.null(fitj)){
        joint.fit <- fitj
        llmax.prev <- fitj$llmax
        par.prev <- fitj$ParSig
        par.data <- fitj$ParSig
###the next signal is searched in the data minus all jointly refined signals
        rv.ls$res.s <- fitj$res.sig
        res.last <- fitj$res
    }
###    phase.plot(tmp)
}

cat('Finished!\n')
