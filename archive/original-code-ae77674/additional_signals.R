####This file is an example for making PeriodoFrame: the computation part
inds <- 1:2
leg.pos <- 'topright'
############################
####find additional signals
############################
#save(list=ls(all=TRUE),file='test3.Robj')
for(jj in 2:Nsig.max){
    cat('\n Find signal ',jj,'!\n')
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
        rv.ls <- BFP(t,rr,dy,Nma=Nma,Nar=Nar,Indices=Indices,ofac=ofac,model.type='man',fmin=frange[1],fmax=frange[2],quantify=quantify,renew=renew,progress=progress,noise.only=noise.only)
        ylab <- 'ln(BF)'
        name <- 'logBF'
    }
    if(per.type.seq=='MLP'){
        rv.ls <- MLP(t,rr,dy,Nma=Nma,Nar=Nar,ofac=ofac,mar.type='part',model.type='man',fmin=frange[1],fmax=frange[2],opt.par=NULL,Indices=Indices,MLP.type=MLP.type,noise.only=noise.only)
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
###modify the periodogram output for Keplerian fit
    tmp <- sigfit(per=rv.ls,data=cbind(t,rr,dy),SigType=SigType,basis=basis,mcf=mcf)
    rv.ls <- tmp$per
#    if(!progress) phase.plot(tmp,fold=fold)

    pp <- cbind(tmp$t,tmp$y,tmp$ysig0)
    colnames(pp) <- paste0(c('t','y','ysig'),'_sig',jj)
    phase.data <- cbind(phase.data,pp)

    qq <- cbind(tmp$tsim,tmp$ysim,tmp$ysim0)
    colnames(qq) <- paste0(c('tsim','ysim','ysim0'),'_sig',jj)
    sim.data <- cbind(sim.data,qq)
    par.data <- addpar(par.data,tmp$ParSig,jj)
###    phase.plot(tmp)
}

cat('Finished!\n')
