library(doMC)
getM0 <- function(e,omega,P,T,T0,type='primary'){
    Tp <- T-getphase(e,omega)*P
    ((T0-Tp)%%P)*2*pi/P
}
kep.mt2 <- function(m,e){
    tol = 1e-8
    E0 <- m
    Ntt <- 1e3
    for(k in 1:Ntt){
        E1 = E0-(E0-e*sin(E0)-m)/(sqrt((1-e*cos(E0))^2-(E0-e*sin(E0)-m)*(e*sin(E0))))
        if(all(abs(E1-E0)<tol)) break()
#        if(k==Ntt) cat('Keplerian solver does not converge:',e,m,E0,E1,'!\n')
        E0 <- E1
    }
    if(k==Ntt){
        cat('Keplerian solver does not converge!\n')
        cat('length(which(abs(E1-E0)>tol))=',length(which(abs(E1-E0)>tol)),'\n')
    }
    return(E1)
}

nrc2 <- function(Nvar){
    nrow <- ceiling(sqrt(Nvar))
    if(Nvar<nrow*(nrow-1)){
        ncol <- nrow-1
    }else{
        ncol <- nrow
    }
    return(c(nrow,ncol))
}

nrc <- function(Nvar){
    nrow <- ceiling(Nvar/2)
    if(Nvar<2){
        ncol <- 1
    }else{
        ncol <- 2
    }
    return(c(nrow,ncol))
}

agatha.dataset.name <- function(filename){
    name <- basename(filename)
    name <- sub('\\.[^.]*$','',name)
    name <- sub('_TERRA.*$','',name)
    name <- gsub('[^[:alnum:]_.-]+','_',name)
    if(is.na(name) || name==''){
        name <- 'dataset'
    }
    return(name)
}

read.agatha.table <- function(path, center.rv=FALSE){
    first.row <- read.table(path,nrows=1,check.names=FALSE)
    has.header <- !is.numeric(first.row[1,1])
    tab <- read.table(path,header=has.header,check.names=FALSE)
    tab <- as.data.frame(tab)
    if(ncol(tab)<2){
        stop('Agatha input tables need at least time and observable columns.',call.=FALSE)
    }
    for(j in 1:ncol(tab)){
        tab[,j] <- as.numeric(tab[,j])
    }
    if(ncol(tab)<3){
        tab$eRV <- rep(1,nrow(tab))
    }
    if(center.rv){
        tab[,2] <- tab[,2]-mean(tab[,2],na.rm=TRUE)
    }
    inds <- sort(tab[,1],index.return=TRUE)$ix
    tab <- tab[inds,,drop=FALSE]
    duplicate.time <- duplicated(tab[,1])
    if(any(duplicate.time)){
        tab <- tab[!duplicate.time,,drop=FALSE]
    }
    if(ncol(tab)>3){
        for(j in 4:ncol(tab)){
            if(is.na(sd(tab[,j])) || sd(tab[,j])==0){
                tab[,j] <- abs(rnorm(nrow(tab),1,0.01))
            }
        }
    }
    if(!has.header || any(grepl('^V[[:digit:]]+$',colnames(tab)[1:min(3,ncol(tab))]))){
        colnames(tab)[1:3] <- c('Time','RV','eRV')
        if(ncol(tab)>3){
            colnames(tab)[4:ncol(tab)] <- paste0('proxy',1:(ncol(tab)-3))
        }
    }else{
        colnames(tab)[1:3] <- c('Time','RV','eRV')
        if(ncol(tab)>3){
            blank <- is.na(colnames(tab)[4:ncol(tab)]) | colnames(tab)[4:ncol(tab)]==''
            colnames(tab)[3+which(blank)] <- paste0('proxy',which(blank))
        }
    }
    rownames(tab) <- NULL
    return(tab)
}

tv.per <- function(targets,ofac,data){
    for(target in targets){
        tab <- data[[target]]
        commandArgs <- function(trailingOnly=TRUE) c('NA',1000,100,ofac,'bgls','res')
        source('time_varying_periodogram.R',local=TRUE)
    }
}

addpar <- function(par.old,par.new,nsig){
    if(nsig==1){
        n0 <- names(par.new)
        if(any(n0=='A')|any(n0=='B')){
            names(par.new)[1:3] <- paste0(n0[1:3],1)
        }
        par <- par.new
    }else{
        n0 <- names(par.new)
        n1 <- n0
        par <- par.old
        if(any(n0=='A')|any(n0=='B')){
            n1[1:3] <- paste0(n0[1:3],nsig)
            names(par.new) <- n1
            Npar.noise <- length(par.new)-3
        }else if(any(n0=='Mo1') & !any(n0=='omega1')){
            n1[1:3] <- gsub('1$',nsig,n0[1:3])
            names(par.new) <- n1
            Npar.noise <- length(par.new)-3
        }else if(any(grepl('omega|Tc',n0))){
            n1[1:5] <- gsub('1$',nsig,n0[1:5])
            names(par.new) <- n1
            Npar.noise <- length(par.new)-5
        }
        par <- c(par.old[-(length(par.old)-(Npar.noise:1)+1)],par.new)
    }
    return(par)
}

calc.1Dper <- function(Nmax.plots, vars,per.par,data,Ncores=4,basis='natural'){
    var <- names(per.par)
    for(k in 1:length(var)){
        assign(var[k],per.par[[var[k]]])
    }
    if(Niter>0){
        mcf <- TRUE
    }else{
        mcf <- FALSE
    }
    if(Ncores>0) {registerDoMC(Ncores)} else {registerDoMC()}
    Nmas <- unlist(Nmas)
    Nars <- unlist(Nars)
    par.list <- sim.list <- phase.list <- per.list <- tits <- mc.list <- list()
    tits <- c()
    fs <- c()
    pars <- list()
    kk <- 1
    if(SigType=='stochastic' & any(per.type=='BFP')){
        noise.only <- TRUE
    }else{
        noise.only <- FALSE
    }
    MLP.type <- 'sub'
    for(j1 in 1:length(vars)){
        for(j2 in 1:length(per.type)){
            if(per.type[j2]=='MLP' | per.type[j2]=='BFP'){
                if(length(per.target)>1){
                    pars[[kk]] <- list(var=vars[j1],per.type=per.type[j2],Inds=0,Nma=Nmas,Nar=Nars)
                }else{
                    pars[[kk]] <- list(var=vars[j1],per.type=per.type[j2],Inds=Inds[[1]],Nma=Nmas[1],Nar=Nars[1])
                }
                kk <- kk+1
            }else{
                pars[[kk]] <- list(var=vars[j1],per.type=per.type[j2],Inds=0,Nma=0,Nar=0)
                kk <- kk+1
            }
        }
    }
    Nvar <- min(length(pars),Nmax.plots)
    sig.levels <- c()
    pers <- c()
    ylabs <- c()
    Pmaxs <- c()
    ypars <- c()
#    lapply(1:Nvar, function(i){
    for(i in 1:Nvar){
        var <- pars[[i]]$var
        Nma <- as.integer(pars[[i]]$Nma)
        Nar <- as.integer(pars[[i]]$Nar)
        Inds <- pars[[i]]$Inds
        if(length(per.target)==1){
            Inds <- as.integer(Inds)
        }
        per.type <- pars[[i]]$per.type
#        instrument <- paste(per.target,collapse='-')
        multi.set <- length(per.target)>1
        if(multi.set){
            instrument <- 'combined'
            subdata <- lapply(1:length(per.target),function(j) data[[per.target[j]]])
            rows <- lapply(1:length(subdata),function(j){
                tabj <- subdata[[j]][,1:3,drop=FALSE]
                cbind(tabj,set.id=per.target[j])
            })
            tab.with.id <- do.call(rbind,rows)
            ord <- sort(as.numeric(tab.with.id[,1]),index.return=TRUE)$ix
            set.id <- tab.with.id[ord,'set.id']
            tab <- tab.with.id[ord,1:3,drop=FALSE]
            tab <- as.data.frame(tab)
            tab[,1] <- as.numeric(tab[,1])
            tab[,2] <- as.numeric(tab[,2])
            tab[,3] <- as.numeric(tab[,3])
            colnames(tab) <- colnames(data[[per.target[1]]])[1:3]
            Inds <- 0
        }else{
            instrument <- per.target
            tab <- data[[per.target]]
        }
###array to store outputs
        per.data <- phase.data <- sim.data <- par.data <- c()
        cnames <- c()
        ypar <- var
#        cat('ypar=',ypar,'\n')
        ypars <- c(ypars,gsub(' ','',ypar))
        Indices <- NULL
        if(ncol(tab)>3){
            Indices <- tab[,4:ncol(tab),drop=FALSE]
        }
        t <- tab[,1]
        if(ypar==ns[1]){
            dy <- tab[,3]
        }else{
            dy <- rep(0.1,nrow(tab))
        }
        if(ypar!='Window Function'){
            y <- tab[,ypar]
            if(ypar!=ns[1]) y <- scale(y)
            if(per.type=='GLST'){
                rv.ls <- glst(t=tab[,1],y=y,err=dy,ofac=ofac,fmin=frange[1],fmax=frange[2])
                ylab <- 'Power'
                name <- 'power'
            }else if(per.type=='GLS'){
                rv.ls <- gls(t=tab[,1]-min(tab[,1]),y=y,err=dy,ofac=ofac,fmin=frange[1],fmax=frange[2])
                ylab <- 'Power'
                name <- 'power'
            }else if(per.type=='BGLS'){
                rv.ls <- bgls(t=tab[,1]-min(tab[,1]),y=y,err=dy,ofac=ofac,fmin=frange[1],fmax=frange[2])
                ylab <- expression('log(ML/'*ML[max]*')')
                name <- 'logML'
            }else if(per.type=='BFP'){
#                if(exists('per.type.seq')){
#                    if(per.type.seq=='BFP'){
#                        quantify <- TRUE
#                    }else{
#                        quantify <- FALSE
#                    }
#                }else{
####quantify is not an important parameter, could either be TRUE or FALSE
                    quantify <- FALSE
#                    quantify <- TRUE
#                }

###preselect Indices according to the value of Inds and Indices
                if(length(Inds)>0){
                    if(all(Inds==0)){
                        Indices <- NULL
                    }else{
                        Inds <- Inds[Inds>0]
                        Indices <- as.matrix(Indices[,Inds,drop=FALSE])
                        for(j in 1:ncol(Indices)){
                            Indices[,j] <- scale(Indices[,j])
                        }
                    }
                }else{
                    Indices <- NULL
                }
#                tmp <- c(Nma=Nma,Nar=Nar,model.type='man',Indices=NULL,
#                                                      ofac=ofac,fmin=frange[1],fmax=frange[2],quantify=quantify)
                if(FALSE){
                    cat('renew=',renew,'\n')
                    cat('t=',head(tab[,1]),'\n')
                    cat('y=',head(y),'\n')
                    cat('dy=',head(dy),'\n')
                    cat('Nma=',Nma,';Nar=',Nar,';model.type=man;Indices=',Indices, ';ofac=',ofac,';fmin=',frange[1],';fmax=',frange[2],';quantify=',quantify, ';renew=',renew,';noise.only=',noise.only,'\n')
                }
                if(multi.set){
                    rv.ls <- BFP.multiset(t=t,y=y,dy=dy,set.id=set.id,Nma=Nma,Nar=Nar,ofac=ofac,fmin=frange[1],fmax=frange[2],progress=FALSE,noise.only=noise.only)
                }else{
                    rv.ls <- BFP(t=t,y=y,dy=dy, Nma=Nma,Nar=Nar,model.type='man',Indices=Indices, ofac=ofac,fmin=frange[1],fmax=frange[2],quantify=quantify, renew=renew,Nsamp=Nsamp,noise.only=noise.only)
                }
###renew: every chi-square minimization start from the initial parameter values
                ylab <- 'ln(BF)'
                name <- 'logBF'
            }else if(per.type=='MLP'){
                if(multi.set){
                    rv.ls <- MLP.multiset(t=t,y=y,dy=dy,set.id=set.id,Nma=Nma,Nar=Nar,ofac=ofac,fmin=frange[1],fmax=frange[2])
                }else{
                    rv.ls <- MLP(t=tab[,1]-min(tab[,1]),y=y,dy=dy,Nma=Nma,Nar=Nar,Indices=Indices,ofac=ofac,fmin=frange[1],fmax=frange[2],MLP.type=MLP.type)
                }
                ylab <- expression('log(ML/'*ML[max]*')')
                name <- 'logML'
            }else if(per.type=='LS'){
                rv.ls <- lsp(times=tab[,1]-min(tab[,1]),x=y,ofac=ofac,from=frange[1],to=frange[2],alpha=c(0.1,0.01,0.001))
                ylab <- 'Power'
                name <- 'power'
            }
#            tit <- paste('Periodogram:',per.type,'; Target:',instrument,'; Observable',ypar)
            tit <- paste0(per.type,'; ',instrument,';', ypar,';1 signal')
            if(!exists('Nma')){
                Nma <- 0
            }
            if(!exists('Nna')){
                Nna <- 0
            }
            if(!exists('Inds')){
                Inds <- 0
            }
#            f <-  paste0(paste(per.target,collapse='_'),'_',gsub(' ','',ypar),'_',per.type,'_MA',paste(Nmas,collapse=''),'proxy',paste(Inds,collapse='.'),'_1sig_',format(rv.ls$P[which.max(rv.ls$power)],digit=2),'d')
            f <-  paste0(paste(per.target,collapse='_'),'_',gsub(' ','',ypar),'_',per.type,'_AR',paste(Nars,collapse=''),'proxy',paste(Inds,collapse='.'),'_1sig_',format(rv.ls$P[which.max(rv.ls$power)],digit=2),'d')
        }else{
            rv.ls <- lsp(times=tab[,1]-min(tab[,1]),x=rep(1,nrow(tab)),ofac=ofac,from=frange[1],to=frange[2],alpha=c(0.1,0.01,0.001))
            tit <- paste0('LS;',instrument,';',ypar)
            pt <- 'LS'
            if(!exists('Nma')){
                Nma <- 0
            }

            if(!exists('Nar')){
                Nar <- 0
            }
            if(!exists('Inds')){
                Inds <- 0
            }
#           f <-  paste0(paste(per.target,collapse='_'),'_',gsub(' ','',ypar),'_',pt,'_MA',paste(Nmas,collapse=''),'proxy',paste(Inds,collapse='.'),'_1sig_',format(rv.ls$P[which.max(rv.ls$power)],digit=2),'d')
            f <-  paste0(paste(per.target,collapse='_'),'_',gsub(' ','',ypar),'_',pt,'_AR',paste(Nmas,collapse=''),'proxy',paste(Inds,collapse='.'),'_1sig_',format(rv.ls$P[which.max(rv.ls$power)],digit=2),'d')
            ylab <- 'Power'
            name <- 'power'
        }
        ylabs <- c(ylabs,ylab)
        tits <- c(tits,tit)
        fs <- c(fs,f)
        pers <- c(pers,per.type)
###plot
#        plotname <- paste("plot", i, sep="")
        if(per.type=='MLP' | per.type=='BGLS'){
            yy  <- rv.ls$power-max(rv.ls$power)
            rv.ls$sig.level <- NULL#max(yy)-log(c(10,100,1000))
        }else{
            yy <- rv.ls$power
        }
        if(!is.null(per.data)){
            if(nrow(per.data)>length(yy)){
                rv.ls$P <- c(rv.ls$P,rv.ls$P[length(rv.ls$P)])
                yy <- c(yy,yy[length(yy)])
            }else if(nrow(per.data)<length(yy)){
                rv.ls$P <- rv.ls$P[-length(rv.ls$P)]
                yy <- yy[-length(yy)]
            }
        }
        if(length(rv.ls$sig.level)<3){
            sig.levels <- cbind(sig.levels,c(rv.ls$sig.level,rep(NA,3-length(rv.ls$sig.level))))
        }else{
            sig.levels <- cbind(sig.levels,rv.ls$sig.level)
        }
#        if(i==1)
        per.data <- cbind(per.data,rv.ls$P)
        per.data <- cbind(per.data,yy)
        Pmaxs <- c(Pmaxs,format(per.data[which.max(yy),1],digit=2))
        inds <- (ncol(per.data)-1):ncol(per.data)
#        if(i==1)
        cnames <- c(cnames,'P')
        cnames <- c(cnames,paste0(pers[i],'1signal:',gsub(' .+','',ypar),':',name))

####calculate the Keplerian fit
	if(Nsig.max>1){
	     Pconv <- FALSE
	}else{
	     Pconv <- TRUE
	}
        fit <- sigfit(per=rv.ls,data=tab,SigType=SigType,basis=basis,Ncores=Ncores,mcf=(mcf & !multi.set),Niter=Niter,Pconv=Pconv)
###update the output from periodogram
        rv.ls <- fit$per

        pp <- cbind(fit$t,fit$y,fit$ysig0)

        colnames(pp) <- paste0(c('t','y','ysig'),'_sig1')
        phase.data <- cbind(phase.data,pp)

        qq <- cbind(fit$tsim,fit$ysim,fit$ysim0)
        colnames(qq) <- paste0(c('tsim','ysim','ysim0'),'_sig1')
        sim.data <- cbind(sim.data,qq)
        par.data <- addpar(c(),fit$ParSig,1)

        if(Nsig.max>1){
            if(per.type==per.type.seq){
                if(length(per.target)>1){
                    Nma <- 0
                    Nar <- 0
                    Inds <- 0
                }
                source('additional_signals.R',local=TRUE)
            }
        }

###use mcmc to update the combined model and output
        mc <- list()
        if(mcf){
            if(Nsig.max>1){
                fit <- mcfit(rv.ls,data=tab[,1:3],tsim=fit$tsim0,Niter=Niter,SigType=SigType,basis=basis,ParSig=par.data,Pconv=TRUE,Ncores=Ncores)
            }
            mc <- fit$mc
            ParSig <- fit$ParSig
            par.data <- fit$par.stat

            pp <- cbind(fit$ysig,fit$ysig0,fit$res)
            colnames(pp) <- paste0(c('y','ysig','res'),'_all')
            phase.data <- cbind(phase.data,pp)
            if(Nsig.max>1){
                qq <- cbind(fit$ysim.sig)
            }else{
                qq <- cbind(fit$ysim0)
            }
            colnames(qq) <- 'ysim_all'
            sim.data <- cbind(sim.data,qq)
        }else{
            res <- fit$res
            if(Nsig.max>1){
                ysig0 <- rowSums(phase.data[,paste0('ysig_sig',1:Nsig.max)])
                ysim <- rowSums(sim.data[,paste0('ysim0_sig',1:Nsig.max)])
            }else{
                ysig0 <- phase.data[,'ysig_sig1']
                ysim <- sim.data[,'ysim0_sig1']
            }
            ysig <- ysig0+res
            tsim0 <- fit$tsim0
            ParSig <- fit$ParSig

            pp <- cbind(ysig,ysig0,res)
            colnames(pp) <- paste0(c('y','ysig','res'),'_all')
            phase.data <- cbind(phase.data,pp)

            qq <- cbind(ysim)
            colnames(qq) <- 'ysim_all'
            sim.data <- cbind(sim.data,qq)
        }

###attach common data
        phase.attach <- cbind(t,y,dy)
        colnames(phase.attach) <- c('t0','y0','ey0')
        sim.attach <- t(t(fit$tsim0))
        colnames(sim.attach) <- 'tsim0'
        phase.data <- cbind(phase.data,phase.attach)
        sim.data <- cbind(sim.data,sim.attach)
        colnames(per.data) <- cnames

###put everything into list
        sim.list[[ypar]] <- sim.data
        phase.list[[ypar]] <- phase.data
        per.list[[ypar]] <- per.data
        par.list[[ypar]] <- par.data
        mc.list[[ypar]] <- mc
    }
    if(!exists('Nsig.max')){
        Nsig.max <- 1
    }
    if(!exists('Nma')){
        Nma <- 0
    }
    if(!exists('Inds')){
        Inds <- 0
    }
    fname <- paste0(paste(per.target,collapse='_'),'_',paste(ypars,collapse='.'),'_',paste(per.type,collapse=''),'_AR',paste(Nar,collapse=''),'MA',paste(Nma,collapse=''),'_proxy',paste(Inds,collapse='.'),'_',Nsig.max,'sig_',paste(Pmaxs,collapse='d'),'d')
#    cat('fname=',fname,'\n')
#    save(list=ls(all=TRUE),file='test1.Robj')
    return(list(per.list=per.list,mc=mc,phase.list=phase.list,sim.list=sim.list,par.list=par.list,tits=tits,pers=pers,levels=sig.levels,ylabs=ylabs,fname=fname,fs=fs,mc.list=mc.list))
}

par.a2m <- function(par,popt,data,SigType='kepler',time.unit=1){
###change parameters from agatha to mcmc
    par <- unlist(par)
    n0 <- names(par)

    startvalue <- c()
#    Nsig <- length(gsub('^A|^Mo',n0))
    if(SigType=='kepler'){
        if(any(grepl('^omega',names(par)))){
            startvalue <- par[1:5]
            if(names(par)[1]=='P1'){
                startvalue[1] <- log(startvalue[1])
            }
        }else{
            phi <- as.numeric(xy2phi(par['A'],par['B']))
            kopt <- as.numeric(sqrt(par['A']^2+par['B']^2))
            startvalue <- c(log(popt),kopt,0,0,phi)
        }
        names(startvalue) <- c('per1','K1','e1','omega1','Mo1')
    }else if(SigType!='stochastic'){
        phi <- as.numeric(xy2phi(par['A'],par['B']))
        kopt <- as.numeric(sqrt(par['A']^2+par['B']^2))
        startvalue <- c(log(popt),kopt,phi)
        names(startvalue) <- c('per1','K1','Mo1')
    }

    par.noise <- c()
    nn <- c()

###fit the trend
    x <- (data[,1]-min(data[,1]))/time.unit
    y <- data[,2]
    fit <- lm(y~x)
    a <- fit$coefficients[2]
    b <- fit$coefficients[1]
    if(any(grepl('beta',n0))){
#        par.noise <- c(par.noise,par['beta']*time.unit)
        par.noise <- c(par.noise,a)
        nn <- c(nn,'a11')
    }
    if(any(grepl('gamma',n0))){
#        par.noise <- c(par.noise,par['gamma'])
        par.noise <- c(par.noise,b)
        nn <- c(nn,'b1')
    }

    if(any(grepl('sj',n0))){
        par.noise <- c(par.noise,par['sj'])
        nn <- c(nn,'s1')
    }else{
        par.noise <- c(par.noise,0)
        nn <- c(nn,'s1')
    }

    if(any(grepl('^l\\d',n0))){
        nar <- length(grep('^l\\d',n0))
        par.noise <- c(par.noise,par[paste0('l',1:nar)])
        nn <- c(nn,paste0('phi1',1:nar))
        par.noise <- c(par.noise,par['logtauAR'])
        nn <- c(nn,'alpha1')
    }

    if(any(grepl('^m\\d',n0))){
        nar <- length(grep('^m\\d',n0))
        par.noise <- c(par.noise,par[paste0('m',1:nar)])
        nn <- c(nn,paste0('w1',1:nar))
        par.noise <- c(par.noise,par['logtau'])
        nn <- c(nn,'beta1')
    }

    if(any(grepl('^d\\d',n0))){
        ii <- grepl('^d\\d',n0)
        par.noise <- c(par.noise,par[ii])
        nn <- c(nn,gsub('d','c',n0[ii]))
    }

    names(par.noise) <- nn

    c(startvalue,par.noise)
}

par.m2a <- function(par.old){
###change parameters from mcmc to agatha
    n0 <- names(par.old)
}

#mcfit <- function(startvalue,Niter,Ncores=1){
mcfit <- function(per,data,tsim,Niter=1e3,SigType='kepler',basis='natural',ParSig=NULL,Pconv=FALSE,Ncores=4,
                  mcmc.method='PT',Ntem=8,tem.min=1e-3,swap.interval=10,mcmc.verbose=FALSE){
###get initial parameters from agatha
#    break()
    time.unit <- 365.25
    par.opt <- unlist(per$par.opt)
    popt <- as.numeric(per$Popt[1])
    if(is.null(ParSig)){
        startvalue <- par.a2m(par.opt,popt,data,SigType=SigType,time.unit=time.unit)
    }else{
        startvalue <- ParSig
    }
####some global parameters for mcmc fit

    tol <- 1e-16
    if(SigType=='kepler'){
        prior.type <- 'mt'
    }else{
        prior.type <- 'e0'
    }
    period.par <- 'logP'
    bases <- rep(basis,10)
    Esd <- 0.1
    phi.min <- wmin <- -1
    phi.max <- wmax <- 1
    ins <- 'none'
    target <- 'TBD'
    offset <- TRUE
    out <- list()
    out$trv.all <- trv.all <- data[,1]
    out$ins <- ins
    out[[ins]] <- list()
    out[[ins]]$RV <- data
    out[[ins]]$index <- 1:nrow(data)
    out$prior.type <- prior.type

    tmin <- min(data[,1])
    tmax <- max(data[,1])
    beta.up <- log(tmax-tmin)#time span of the data
    beta.low <- log(max(1/24,min(1,min(diff(trv.all)))))#1h or minimum separation
    alpha.max <- beta.max <- beta.up#d; limit the range of beta to avoid multimodal or overfitting
    alpha.min <- beta.min <- beta.low#24h
    nqp <- c(length(grep('^c\\d',names(startvalue))),length(grep('^w\\d',names(startvalue))),length(grep('^phi\\d',names(startvalue))))
    out$nqp <- nqp
    out[[ins]]$noise <- list(nqp=nqp)
    Npar <- length(startvalue)
    Sd <- 2.4^2/Npar#hyp
    Dt <- (tmax-tmin)/time.unit
    if(FALSE){
    par.min <- sapply(1:length(startvalue),function(i) startvalue[i]-max(0.1*abs(startvalue[i]),1))
    par.max <- sapply(1:length(startvalue),function(i) startvalue[i]+max(0.1*abs(startvalue[i]),1))
    names(par.min) <- names(par.max) <- names(startvalue)
    }else{
    par.min <- startvalue-1*abs(startvalue)
    par.max <- startvalue+1*abs(startvalue)
    inde <- grep('^e\\d',names(par.min))
    indMo <- grep('^omega|^Mo',names(par.min))
    inda <- grep('^a',names(par.min))
    indb <- grep('^b',names(par.min))
    indK <- grep('^K',names(par.min))
    indc <- grep('^c',names(par.min))
    indphi <- grep('^phi|^w',names(par.min))
    indbeta <- grep('^beta|^alpha',names(par.min))
    inds <- grep('^s\\d',names(par.min))
    indP <- grep('^per',names(par.min))
    if(length(indP)>0){
        par.min[indP] <- startvalue[indP]+log(0.8)
        par.max[indP] <- startvalue[indP]+log(1.2)
    }
    if(length(inds)>0){
        par.min[inds] <- 0
        par.max[inds] <- sd(data[,2])
    }
    if(length(inde)>0){
        par.min[inde] <- 0
        par.max[inde] <- 1
    }
    if(length(indMo)>0){
        par.min[indMo] <- 0
        par.max[indMo] <- 2*pi
    }
    if(length(indb)>0){
        par.min[indb] <- min(par.min[indb],-10*sd(data[,2]))
        par.max[indb] <- max(par.max[indb],10*sd(data[,2]))
    }
    if(length(indK)>0){
        par.min[indK] <- 0.5*startvalue[indK]
        par.max[indK] <- 2*startvalue[indK]
    }
    if(length(inda)>0){
        par.min[inda] <- min(par.min[inda],-10*sd(data[,2])/Dt)
        par.max[inda] <- max(par.max[inda],10*sd(data[,2])/Dt)
    }
    if(length(indc)>0){
        par.min[indc] <- min(par.min[indc],-10*sd(data[,2]))
        par.max[indc] <- max(par.max[indc],10*sd(data[,2]))
    }
    if(length(indphi)>0){
        par.min[indphi] <- min(phi.min,par.min[indphi])
        par.max[indphi] <- max(phi.max,par.max[indphi])
    }
    if(length(indbeta)>0){
        par.min[indbeta] <- min(alpha.min,par.min[indbeta])
        par.max[indbeta] <- max(alpha.max,par.max[indbeta])
    }
    }
    cov.start <- diag(length(startvalue))*1e-6
####mcmc
    source('mcmc_func.R',local=TRUE)
#    mcmc <- foreach(ncore=1:Ncores,.combine='rbind') %dopar% {
    Niter0 <- Niter
    per.prim <- c()
    mcmc <- foreach(ncore=1:Ncores,.errorhandling = 'pass') %dopar% {
        if(mcmc.method=='PT'){
###parallel tempering: the hot replicas explore the aliases of the periodogram
###peak and hand good states to the cold chain through replica exchange, so no
###separate adaptive-tempering burn-in is needed
            tmp <- run.ptmcmc(startvalue,cov.start,iterations=max(as.numeric(Niter),1000),
                              bases=rep(basis,10),Ntem=Ntem,tem.min=tem.min,
                              swap.interval=swap.interval,verbose=mcmc.verbose)
            if(mcmc.verbose){
                cat('PTMCMC tems:',paste(format(tmp$tems,digit=2),collapse=','),'\n')
                cat('PTMCMC per-replica acceptance (%):',paste(round(tmp$acc.all,1),collapse=','),'\n')
                cat('PTMCMC swap acceptance (%):',paste(round(100*tmp$swap.rate,1),collapse=','),'\n')
                cat('PTMCMC out-of-bound proposals (%):',paste(round(tmp$null.rate,1),collapse=','),'\n')
                cat('PTMCMC step scale:',paste(format(tmp$lambda,digit=2),collapse=','),'\n')
            }
        }else{
            cat('use hot_chain.R\n')
            source('hot_chain.R',local=TRUE)
            startvalue <- par.hot
            tmp <- run.metropolis.MCMC(startvalue,cov.start,iterations=max(as.numeric(Niter),1000),tem=1,bases=rep(basis,10))
        }
        tmp$out
    }

#    mcmc  <- list()
#    mcmc[[1]] <- tmp$out
    ind <- which(sapply(1:length(mcmc),function(k) is.null(dim(mcmc[[k]]))))
    if(length(ind)>0) mcmc <- mcmc[-ind]
    mc <- c()
    for(j in 1:length(mcmc)){
        mc <- rbind(mc,mcmc[[j]])
    }

####analyze the MCMC results
    ll <- mc[,'loglike']
    lp <- mc[,'logpost']
    llmax <- max(ll)
    lpmax <- max(lp)
    ind.max <- which.max(mc[,'loglike'])
    par.opt0 <- mc[ind.max,1:Npar]

###derive other parameters
    mc1 <- c()
    if(any(grepl('^omega',colnames(mc))) & Pconv){
        indMo <- grep('^Mo',colnames(mc))
        indP <- grep('^per',colnames(mc))
        t0 <- tmin
        if(tmin<24e5) t0 <- t0+24e5
        Tps <- M02Tp(mc[,indMo],t0,mc[,indP])
        T0 <- t0
        mc1 <- cbind(t0,Tps)
        colnames(mc1) <- c('T0',paste0('Tp',1:length(indP)))
    }

####change per to P
    if(length(indP)>0 & Pconv){
        mc[,indP] <- exp(mc[,indP])
        colnames(mc)[indP] <- gsub('per','P',colnames(mc)[indP])
    }
    ParSig <- mc[ind.max,1:Npar]

    mc.more <- cbind(mc[,1:Npar],mc1)
    par.stat <-  sapply(1:ncol(mc.more),function(i) data.distr(mc.more[,i],ll,plotf=FALSE))

    n <- colnames(mc)[1:Npar]
    if(length(mc1)>0) n <- c(n,colnames(mc1))
    colnames(par.stat) <- n
#    save(list=ls(all=TRUE),file='test0.Robj')

####model prediction
#    rv <- RVsig(ParSig,out=out)
    rv <- RVsig(par.opt0,bases=bases)
#    rv.sig <- RV.kepler(par.opt,bases=bases)[[ins]]
    ysig0 <- rv$ysig
    ytrend <- rv$ytrend
    yproxy <- rv$yproxy
    rv.model <- rv.kep <- rv$y[[ins]]
    yred <- yma <- yar <- 0
    ins <- out$ins
    trv <- out[[ins]]$RV[,1]
    rv.data <- out[[ins]]$RV[,2]
    erv <- out[[ins]]$RV[,3]
    nqp <- out[[ins]]$noise$nqp
    trv <- data[,1]
    if(nqp[2]>0 | nqp[3]>0){
        pp <- arma(t=trv,ymodel=rv.model,ydata=rv.data,pars=par.opt0,ind.set=1,p=nqp[3],q=nqp[2])
        yar <- pp$ar
        yma <- pp$ma
    }
    popt <- NA
    if(any(grepl('^per',names(ParSig)))){
        popt <- exp(ParSig[grepl('^per',names(ParSig))])
    }else if(any(grepl('^P',names(ParSig)))){
        popt <- ParSig[grepl('^P',names(ParSig))]
    }
    y <- ysig0+ytrend+yproxy+yma+yar
    yred <- yma+yar
    res <- data[,2]-y
    res.sig <- data[,2]-ysig0
#    res.sig <- res+ysig0
#
#    break()


#    res <- calc.res(par.opt,bases)[[ins]]

    ysig <- res+ysig0
    ysim.red <- 0
    if(!all(yred==0)){
        redfun <- approxfun(data[,1]-min(data[,1]),yred)
        ysim.red <- redfun(tsim)
    }
    ysim.sig <- RV.kepler(pars.kep=par.opt0,tt=tsim+tmin,kep.only=TRUE,bases=bases)

    if(FALSE){
        ysim.proxy <- 0
        if(!all(yproxy==0)){
            proxyfun <- approxfun(data[,1]-min(data[,1]),yproxy)
            ysim.proxy <- proxyfun(tsim)
        }
        ysim.all <- RV.kepler(pars.kep=par.opt0,tt=tsim+tmin,bases=bases)+ysim.red
    }

    list(mc=mc,llmax=llmax,lpmax=lpmax,ParSig=ParSig,out=out,par.stat=par.stat,yma=yma,yar=yar,yred=yred,ysig=ysig,ysig0=as.numeric(ysig0),ysim.red=ysim.red,ysim.sig=ysim.sig,ytrend=ytrend,yproxy=yproxy,res=res,res.sig=res.sig,popt=popt,tsim0=tsim)#ysim.all=ysim.all
}

sigfit.multiset <- function(per, data, t, tsim, SigType='circular', mcf=FALSE){
###Turn a multi-data-set periodogram into a fitted model, its residual and its
###phase-folded prediction. SigType selects a shared circular signal, a shared
###Keplerian signal or a purely stochastic (signal-free) red-noise model.
    if(mcf){
        warning('MCMC refinement is not implemented for multi-set periodograms; using the weighted linear multi-set fit.')
    }
    if(SigType=='stochastic' & !isTRUE(per$noise_only)){
        warning('The multi-set periodogram was not computed with a purely stochastic model; using the circular signal instead.')
        SigType <- 'circular'
    }
    if(SigType!='stochastic' & isTRUE(per$noise_only)){
        warning('The multi-set periodogram is purely stochastic; reporting the stochastic fit.')
        SigType <- 'stochastic'
    }
    if(SigType=='kepler'){
        fit <- KeplerFit.multiset(per)
        popt <- as.numeric(fit$ParKep$P1)
        ysig0 <- as.numeric(fit$ysig)
        res <- as.numeric(fit$res)
        ysim.sig <- multiset_kepler_curve(fit$ParKep,tsim)
        ParSig <- unlist(fit$ParKep)
        per$par.opt <- fit$ParKep
        per$Popt <- popt
        per$res <- res
        per$ysig <- ysig0
        per$yfull <- fit$yfull
        per$lnBF.kepler <- fit$lnBF
    }else if(SigType=='stochastic'){
        popt <- as.numeric(per$Popt[1])
        ysig0 <- as.numeric(per$yred)
        res <- as.numeric(per$res)
        yred.fun <- approxfun(t,ysig0,rule=2)
        ysim.sig <- as.numeric(yred.fun(tsim))
        ParSig <- c(tau=popt,unlist(per$par.opt))
    }else{
        popt <- as.numeric(per$Popt[1])
        par.opt <- unlist(per$par.opt)
        ysig0 <- as.numeric(per$ysig)
        res <- as.numeric(per$res)
        ysim.sig <- as.numeric(par.opt['A']*cos(2*pi/popt*tsim)+par.opt['B']*sin(2*pi/popt*tsim))
        ParSig <- c(P=popt,par.opt)
    }
    ysig <- ysig0+res
    per$res.s <- res
    if(is.null(popt)){
        popt <- 1e7
    }else if(is.na(popt) | popt<=0){
        popt <- 1e7
    }
    tsim1 <- tsim%%popt
    inds <- sort(tsim1,index.return=TRUE)$ix
    return(list(per=per,t=t%%popt,y=as.numeric(ysig),ey=data[,3],res=as.numeric(res),
                ysig0=as.numeric(ysig0),tsim0=tsim,ysim0=ysim.sig,
                tsim=tsim1[inds],ysim=ysim.sig[inds],ParSig=ParSig,
                par.stat=NULL,popt=popt,mc=c()))
}

sigfit <- function(per,data,SigType='circular',basis='natural',mcf=TRUE,Ncores=4,Niter=1e3,Pconv=FALSE,res.type='sig',
                   mcmc.method='PT',Ntem=8,tem.min=1e-3,swap.interval=10,mcmc.verbose=FALSE){
###This function is to modify the output of various periodograms to give residual, model prediction, and optimal parameters as well as posterior/likelihood samples
    ##x is a list
    ##SigType is either circular or kepler
#    if(any(grepl('gamma',names(per$par.opt)))){
#        per$par.opt['gamma'] <- data[1,2]
#    }
#
    ParSig <- par.opt <- unlist(per$par.opt)
    par.list <- as.list(per$par.opt)
    par.stat <- NULL
    mc <- c()

#    if(any(names(per)=='data')){
#        data <- per$data
#    }else{
        per$data <- data
#    }
    if(!any(names(per$df)=='data')){
        per$df$data <- data
    }
    tmin <- min(data[,1])
    t <- data[,1]-tmin
    tsim <- seq(0,max(t),length.out=1e4)
    if(!any(names(per)=='ysims')) per$ysims <- ysims <- 0
    popt <- per$Popt
    save.data <- FALSE
    if(isTRUE(per$multi_set)){
        return(sigfit.multiset(per=per,data=data,t=t,tsim=tsim,SigType=SigType,mcf=mcf))
    }
#    }else{
        if(SigType=='circular'){
            if(!mcf){
                ysim.sig <- par.opt['A']*cos(2*pi/popt*tsim)+par.opt['B']*sin(2*pi/popt*tsim)
                ysim.all <- ysim.sig
                if(any(names(par.opt)=='gamma')) ysim.all <- ysim.all+par.opt['gamma']
                if(any(names(par.opt)=='beta')) ysim.all <- ysim.all+par.opt['beta']*tsim
                ysig0 <- par.opt['A']*cos(2*pi/popt*t)+par.opt['B']*sin(2*pi/popt*t)
                ysig <- per$res+ysig0
                res <- per$res
                if(any(names(per)=='df') & FALSE){
                    df <- per$df
                    if(any(names(df)=='NI') & any(names(df)=='data') & any(names(df)=='Indices')){
                        if(!any(names(df)=='Nma')) df$Nma <- 0
                        if(!any(names(df)=='NI')) df$NI <- 0
                        if(!any(names(df)=='Nar')) df$Nar <- 0
                        per$df <- df
                        fit <- CircularFit(per,data)
                        ParSig <- unlist(fit$par)
                        per$Popt <- popt <- fit$Popt
                        ysig <- fit$yfull
                        ysig0 <- fit$ysig
                        res <- fit$res
                        cat('0per$res=',sd(per$res),'m/s\n')
                        if(res.type=='sig'){
                            per$res.s <- fit$res.sig
                        }else{
                            per$res.s <- res
                        }
                        cat('per$res=',sd(per$res),'m/s\n')
                        cat('res=',sd(res),'m/s\n')
                        par.list <- as.list(fit$par)
                        sim <- CircularSim(par.list,df,popt,tsim)
                        ysim.sig <- sim$ysig
                        ysim.all <- sim$y
                        cat('popt=',popt,'\n')
                        cat('names(par.opt)=',names(par.opt),'\n')
                        cat('par.opt=',par.opt,'\n')
                        cat('names(ParSig)=',names(ParSig),'\n')
                        cat('ParSig=',ParSig,'\n')
                    }
                }
                ParSig <- c(P=popt,ParSig)
            }
        }else if(SigType=='kepler'){
###Keplerian fitting
            if(!mcf){
                df <- per$df
                if(!any(names(df)=='Nma')) df$Nma <- 0
                if(!any(names(df)=='NI')) df$NI <- 0
                if(!any(names(df)=='Nar')) df$Nar <- 0
                fit <- KeplerFit(per,basis=basis)
                per$par.opt <- fit$ParKep
                per$Popt <- fit$ParKep$P1
                sig <- KeplerSig(fit$ParKep,df,basis=basis)
                sim <- KeplerSim(fit$ParKep,df,tsim,basis=basis)#with trend
                res <- sig$res
                ParSig <- unlist(fit$ParKep)#unlist
                ysig0 <- sig$ysig
                ysig <- sig$res+ysig0
                if(res.type=='sig'){
                    per$res <- per$data[,2]-sig$ysig
                }else{
                    per$res <- sig$res
                }
                popt <- fit$ParKep$P1
                ysim.all <- sim$y
                ysim.sig <- sim$ysig
                per$ysims <- per$ysims+ysim.sig
            }
        }else if(SigType=='stochastic'){
###Stochastic fitting; type=='noise'
            if(!mcf){
                per$Popt <- exp(per$par.opt[grep('logtau',names(per$par.opt))])
                popt <- per$Popt[1]
                df <- per$df
                fit <- CircularSig(par.list,df)
                ysig0 <- fit$yred
                ysig <- fit$res+ysig0
                res <- per$res <- fit$res
                sim <- CircularSim(par.list,df,popt,tsim)
                ysim.sig <- sim$yred
                ysim.all <- sim$y
                per$ysims <- per$ysims+ysim.sig
            }
        }
#    }
#    save(list=ls(all=TRUE),file='test.Robj')
    if(mcf){
        tmp <- mcfit(per=per,data=data,tsim=tsim,Niter=Niter,SigType=SigType,basis=basis,Pconv=Pconv,Ncores=Ncores,
                     mcmc.method=mcmc.method,Ntem=Ntem,tem.min=tem.min,swap.interval=swap.interval,
                     mcmc.verbose=mcmc.verbose)
        startvalue <- ParSig <- tmp$ParSig
        mc <- tmp$mc
	par.stat <- tmp$par.stat
        res <- tmp$res
        if(SigType!='stochastic'){
            ysig0 <- tmp$ysig0
            ysig <- tmp$ysig
            ysim.sig <- tmp$ysim.sig
            popt <- tmp$popt
        }else{
            ysig0 <- tmp$yred
            ysig <- tmp$yred+tmp$res
            ysim.sig <- tmp$ysim.red
            popt <- exp(tmp$ParSig[grep('beta|alpha',names(tmp$ParSig))][1])
        }
        ysim.all <- tmp$ysim.all
        per$ysims <- per$ysims+ysim.sig
        if(res.type=='sig'){
            per$res.s <- tmp$res.sig
        }else{
            per$res.s <- tmp$res
        }
    }

####output
    if(is.na(popt) | is.null(popt)){
        popt <- 1e7
    }
    tsim1 <- tsim%%popt
    inds <- sort(tsim1,index.return=TRUE)$ix
    tsim2 <- tsim1[inds]
    ysim2 <- ysim.sig[inds]
    tsim0 <- tsim
    ysim0 <- ysim.sig
    ts <- t%%popt
    tsims <- tsim2
    ysims <- ysim2
    return(list(per=per,t=ts,y=as.numeric(ysig),ey=data[,3],res=as.numeric(res),ysig0=as.numeric(ysig0),tsim0=tsim0,ysim0=ysim0,tsim=tsims,ysim=ysims,ParSig=ParSig,par.stat=par.stat,popt=popt,mc=mc))
}

phase1D.plot <- function(phase.list,sim.list,tits,download=FALSE,index=NULL,repar=TRUE){
    if(repar){
        if(is.null(index)){
            par(mfrow=c(ceiling(Nmax.plots/2),2),mar=c(5,5,3,1),cex.lab=1.5,cex.axis=1.5,cex=1,cex.main=1.0)
        }
        if(download & is.null(index)){
            par(mfrow=c(2,2),mar=c(5,5,3,1),cex.lab=1.2,cex.main=0.8,cex.axis=1.2,cex=1)
        }
    }
    ypars <- names(phase.list)
    for(ypar in ypars){
        if(!is.null(index)){
            inds <- index
            titles <- rep('',inds)
        }else{
            inds <- 1:length(grep('y_sig',colnames(phase.list[[ypar]])))
            titles <- tits[grepl(paste0(';',ypar,';'),tits)]
        }
        ylab <- unlist(strsplit(titles[1],';'))[3]
        ysig0 <- 0
        ysim0 <- 0
        for(i in inds){
            t <- phase.list[[ypar]][,paste0('t_sig',i)]
            y <- phase.list[[ypar]][,paste0('y_sig',i)]
            ey <- phase.list[[ypar]][,'ey0']

            tsim <- sim.list[[ypar]][,paste0('tsim_sig',i)]
            ysim <- sim.list[[ypar]][,paste0('ysim_sig',i)]

            ysig0 <- ysig0+phase.list[[ypar]][,paste0('ysig_sig',i)]
            ysim0 <- ysim0+sim.list[[ypar]][,paste0('ysim0_sig',i)]

            plot(t,y,xlab='Phase [day]',ylab=ylab,main=titles[i],col='white')
                                        #    lines(df$tsim,df$ysim,col=tcol('red',50))
            lines(tsim,ysim,col='red',lwd=3)
            points(t,y,col='black')
            arrows(t,y-ey,t,y+ey,length=0.03,angle=90,code=3,col=tcol('black',50))
        }
###total fit
        t <- phase.list[[ypar]][,'t0']
        ey <- phase.list[[ypar]][,'ey0']
        res <- phase.list[[ypar]][,'res_all']
                                        #    y <- res+ysig0
        y <- phase.list[[ypar]][,'y_all']
        tsim0 <- sim.list[[ypar]][,'tsim0']
        ysim.all <- sim.list[[ypar]][,'ysim_all']

        plot(t,y,xlab='BJD [day]',ylab=ylab,main=gsub('\\d signal','combined fit',titles[i]),col='white')
        lines(tsim0+min(t),ysim.all,col='red',lwd=3)
        points(t,y,col='black')
        try(arrows(t,y-ey,t,y+ey,length=0.03,angle=90,code=3,col=tcol('black',50)),TRUE)
        legend('topleft',bty='n',legend=paste('RMS =',round(sd(y),1),'\n'),text.col='blue')

###residual
        plot(t,res,xlab='BJD [day]',ylab=ylab,main=gsub('\\d signal','residual',titles[i]))
        try(arrows(t,res-ey,t,res+ey,length=0.03,angle=90,code=3,col=tcol('black',50)),TRUE)
        legend('topleft',bty='n',legend=paste('RMS =',round(sd(res),1),'\n'),text.col='blue')
    }
}

combined.plot <- function(per.list,phase.list,sim.list,tits,pers,levels,ylabs,SigType='circular',download=FALSE,index=NULL){
    per1D.plot(per.list,tits,pers,levels,ylabs,download=download,index=index,SigType=SigType)
    phase1D.plot(phase.list,sim.list,tits=tits,download=download,index=index,repar=FALSE)
}

per1D.plot <- function(per.list,tits,pers,levels,ylabs,download=FALSE,index=NULL,SigType='circular'){
    if(is.null(index)){
        par(mfrow=c(ceiling(Nmax.plots/2),2),mar=c(5,5,3,1),cex.lab=1.5,cex.axis=1.5,cex=1,cex.main=1.0)
    }
    if(download & is.null(index)){
        par(mfrow=c(2,2),mar=c(5,5,3,1),cex.lab=1.2,cex.main=0.8,cex.axis=1.2,cex=1)
    }
    ypars <- names(per.list)
    for(ypar in ypars){
        P <- per.list[[ypar]][,1]
        if(!is.null(index)){
            inds <- index
            titles <- rep('',inds)
        }else{
            inds <- 1:(ncol(per.list[[ypar]])-1)
            titles <- tits[grepl(paste0(';',ypar,';'),tits)]
        }
        for(i in inds){
            power <- per.list[[ypar]][,i+1]
            ylab <- ylabs[i]#gsub('.+:','',colnames(per.list)[i+1])
            per.type <- gsub('[[:digit:]]signal:.+','',colnames(per.list[[ypar]])[i+1])
            f1 <- gsub('signal:.+','',colnames(per.list[[ypar]])[i+1])
            Nsig <- gsub('[A-Z]','',f1)
            if(SigType!='stochastic'){
                ymin <- median(power)
            }else{
                ymin <- max(0,min(power))
            }
                                        #        ymin <- min(power)
            if(grepl('Window',titles[i])){
                ylim <- c(ymin,max(power)+0.15*(max(power)-ymin))
            }else{
                                        #            ylim <- c(ymin,max(max(power)+0.15*(max(power)-ymin),levels[which(!is.na(levels[,i])),i]))
                ylim <- c(ymin,max(power)+0.15*(max(power)-ymin))
            }
            if(SigType=='stochastic'){
                xlab <- 'Time scale [day]'
            }else{
                xlab <- 'Period [day]'
            }
            plot(P,power,xlab=xlab,ylab=ylab,xaxt='n',log='x',type='l',main=titles[i], ylim=ylim)
            magaxis(side=1,tcl=-0.5)
            if(!grepl('Window',titles[i])){
                abline(h=levels[,i],lty=2)
            }
            p <- show.peaks(ps=P,powers=power,levels=levels[,i])
            if(!is.matrix(p)){
                pmaxs <- p[1]
                power.max <- p[2]
            }else{
                pmaxs <- p[,1]
                power.max <- p[,2]
                if(length(pmaxs)>4){
                    pmaxs <- p[1:2,1]
                    power.max <- p[1:2,2]
                }
            }
            if(length(pmaxs)>0){
                par(xpd=TRUE)
                offset <- c(0.08*(max(power)-ymin), 0.02*(max(power)-ymin))
                                        #            pmaxs <- pmaxs[1]
                                        #            power.max <- power.max[1]
                pmaxs <- P[which.max(power)]
                power.max <- power[which.max(power)]
                text(pmaxs,power.max+offset[1],pos=3,labels=format(pmaxs,digit=4),col='red',cex=1.0)
                try(arrows(pmaxs,power.max+offset[1],pmaxs,power.max+offset[2],col='red',length=0.05),TRUE)
                par(xpd=FALSE)
            }
        }
    }
}

per2D.data <- function(vars,per.par,data){
    var <- names(per.par)
    for(k in 1:length(var)){
        assign(var[k],per.par[[var[k]]])
    }
    Nmas <- unlist(Nmas)
    Nars <- unlist(Nars)
    pars <- list()
    kk <- 1
    for(j1 in 1:length(vars)){
        for(j2 in 1:length(per.type)){
            if(per.type[j2]=='MLP' | per.type[j2]=='BFP'){
                pars[[kk]] <- list(var=vars[j1],per.type=per.type[j2],Inds=Inds[[1]],Nma=Nmas[1],Nar=Nars[1])
                kk <- kk+1
            }else{
                pars[[kk]] <- list(var=vars[j1],per.type=per.type[j2],Inds=0,Nma=0,Nar=0)
                kk <- kk+1
            }
        }
    }
    i <- 1
    if(length(per.target)==1){
        Nma <- as.integer(pars[[i]]$Nma)
        Nar <- as.integer(pars[[i]]$Nar)
        Inds <- as.integer(pars[[i]]$Inds)
    }
    Indices <- NULL
    per.type <- pars[[i]]$per.type
    var <- pars[[i]]$var
    if(length(per.target)>1){
        instrument <- 'combined'
        subdata <- lapply(1:length(per.target),function(j) data[[per.target[j]]])
        tmp <- combine.data(data=subdata,Ninds=Inds,Nmas=Nmas)
        tab <- tmp$cdata
        idata <- tmp$idata
        colnames(tab) <- colnames(data[[1]])[1:3]
    }else{
        instrument <- per.target
        tab <- data[[per.target]]
        if(ncol(tab)>3){
            Indices <- as.matrix(tab[,4:ncol(tab),drop=FALSE])
        }
    }
    if(length(Inds)>0){
        if(all(Inds==0)){
            Indices <- NULL
        }else{
            Inds <- Inds[Inds>0]
            Indices <- as.matrix(Indices[,Inds,drop=FALSE])
            for(j in 1:ncol(Indices)){
                Indices[,j] <- scale(Indices[,j])
            }
        }
    }else{
        Indices <- NULL
    }
    ypar <- var
    t <- tab[,1]
    y <- tab[,ypar]
    dy <- tab[,3]
    if(length(per.target)==1){
        mp <- MP(t=t,y=y,dy=dy,Dt=Dt,nbin=Nbin,ofac=ofac,fmin=frange[1],fmax=frange[2],per.type=per.type,sj=0,Nma=Nma,Nar=Nar,Indices=Indices)
    }else{
        mp <- MP(t=t,y=y,dy=dy,Dt=Dt,nbin=Nbin,ofac=ofac,fmin=frange[1],fmax=frange[2],per.type=per.type,sj=0,Nma=0,Nar=Nar,Indices=Indices)
    }
    x2 <- mp$tmid
    y2 <- mp$P
    z2 <- mp$powers
    z2.rel <- mp$rel.powers
    fname <- paste0(paste(per.target,collapse='_'),'_MP_',paste(per.type,collapse=''),'_MA',paste(Nmas,collapse=''),'proxy',paste(Inds,collapse='.'))
    if(length(per.target)==1){
        return(list(t=t,y=y,dy=dy,xx=x2,yy=y2,zz=z2,zz.rel=z2.rel,fname=fname,ypar=ypar))
    }else{
        return(list(t=t,y=y,dy=dy,xx=x2,yy=y2,zz=z2,zz.rel=z2.rel,subdata=subdata,idata=idata,fname=fname,ypar=ypar))
    }
}

plotMP <- function(vals,pars){
    var <- names(pars)
    for(k in 1:length(var)){
        assign(var[k],pars[[var[k]]])
    }
    if(length(per.target)>1){
        subdata <- vals$subdata
        idata <- vals$idata
    }
    ypar <- vals$ypar
    t <- vals$t
    y <- vals$y
    dy <- vals$dy
    xx <- vals$xx
    yy <- vals$yy
    zz <- vals$zz
    zz.rel <- vals$zz.rel
#    save(list=ls(all=TRUE),file='test3.Robj')
    source('MP_plot.R',local=TRUE)
}

calcBF <- function(data,Nbasic,proxy.type,Nma.max,Nar.max,groups=NULL,Nproxy=NULL,Npoly=c(2,0),progress=FALSE){
##add Nar.max
    t <- data[,1]
    y <- data[,2]
    dy <- data[,3]
    NI.max <- ncol(data)-3
    NI.inds <- list(0)
    if(NI.max>0){
        NI.inds <- list()
        Nvary <- NI.max-Nbasic
        if(proxy.type=='cum' & Nvary>0 & Nproxy>0){
            NI.inds[[1]] <- Nbasic
            Indices <- data[,4:ncol(data),drop=FALSE]
            cors <- c()
            ###detrend the data first
            if(Npoly[1]>0){
                x <- t
                p <- lm(y~poly(x,Npoly[1]))
                y1 <- residuals(p)
            }else{
                y1 <- y
            }
            for(j in 1:ncol(Indices)){
                if(sd(Indices[,j])==0){
                    cors <- c(cors,0)
                }else{
                    if(Npoly[2]>0){
                        z <- Indices[,j]
                        p <- lm(z~poly(x,Npoly[2]))
                        z1 <- residuals(p)
                    }else{
                        z1 <- Indices[,j]
                    }
                    cors <- c(cors,abs(cor(z1,y1)))
                }
            }
            inds <- sort(cors,decreasing=TRUE,index.return=TRUE)$ix
            if(Nproxy>Nbasic){
                for(j in 1:(Nproxy-Nbasic)){
                    NI.inds[[j+1]] <- inds[1:j]
                }
            }
        }else if(proxy.type=='group' & Nvary>0){
            NI.inds <- lapply(1:(length(groups)+1),function(i) NI.inds[[i]] <- list())
            if(Nbasic>0){
                NI.inds[[1]] <- 1:Nbasic
            }else{
                NI.inds[[1]] <- 0
            }
            groups <- sort(as.integer(groups))
            for(j in 1:length(groups)){
                if(j==1){
                    NI.inds[[j+1]] <- 1:groups[j]
                }else{
                    if(Nbasic>0){
                        NI.inds[[j+1]] <- c(1:Nbasic,(groups[j-1]+1):groups[j])
                    }else{
                        NI.inds[[j+1]] <- (groups[j-1]+1):groups[j]
                    }
                }
            }
        }else if(proxy.type=='man'){
            NI.inds <- groups
        }else{
            NI.inds <- list(list(Nbasic:NI.max))
        }
    }
    Nmas <- 0:Nma.max
    Nars <- 0:Nar.max
    if(ncol(data)>3){
#        out <- BFP.comp(data, Nmas=0:Nma.max,Nars=0:Nar.max,NI.inds=NI.inds,progress=progress)
        out <- bfp.inf.progress(data,Nmas=Nmas,Nars=Nars,NI.inds=NI.inds)
    }else{
        out <- bfp.inf.progress(data,Nmas=Nmas,Nars=Nars,NI.inds=0)
#        out <- BFP.comp(data, Nmas=0:Nma.max,Nars=0:Nar.max,NI.inds=0,progress=progress)
    }
#    out$logBF
#    if(!is.matrix(out$logBFs)){
    lnBFs <- flatten2d(out$lnBF)
#    }
    return(list(Inds=NI.inds,Inds.opt=out$NI.opt,Nars=0:Nar.max,Nmas=0:Nma.max,Nma.opt=out$Nma.opt,Nar.opt=out$Nar.opt,lnBF=lnBFs,extra=out))
#    return(out)
}

flatten3d <- function(arr){
    dn <- dimnames(arr)
    coln <- gsub('\\d','',c(dn[[1]][1],dn[[2]][1],dn[[3]][1]))
    nn <- outer(outer(gsub('[a-z]|[A-Z]','',dn[[1]]),gsub('[a-z]|[A-Z]','',dn[[2]]),paste),gsub('[a-z]|[A-Z]','',dn[[3]]),paste)
    ns <- t(sapply(1:length(nn),function(i) unlist(strsplit(nn[i],' '))))
    tmp <- data.frame(cbind(ns,flatten(arr)))
    colnames(tmp) <- c(coln,'val')
#    tmp[,1:3] <- gsub('[a-z]|[A-Z]','',tmp[,1:3])
    tmp
}

flatten2d <- function(arr){
    Ncol <- dim(arr)[2]*dim(arr)[3]
    Nrow <- dim(arr)[1]
    dn <- dimnames(arr)
    cn  <- paste0('ARMA(',gsub(' ',',',outer(gsub('[a-z]|[A-Z]','',dn[[3]]),gsub('[a-z]|[A-Z]','',dn[[2]]),paste)),')')
    out <- array(NA,dim=c(Nrow,Ncol))
    colnames(out) <- cn
    rownames(out) <- unlist(dimnames(arr)[1])
    j <- 1
    for(i in 1:Nrow){
        out[i,] <- flatten(arr[i,,])
    }
    out
}

MCMC.panel <- function(){
    id <- 'HD020794_TERRA_1AP1_ervab6ap_ccf'
    Niter <- 1.0e3
    Nbin.per <- 1
    nbin.per <- 1
    tem <- 1
    inicov <- 1e-3
    Pini <- 200#day
    noise.model <- 'ARMA05'#noise.model: white, GP(R), ARMA, TJ(Ntj=1,noise vary with RHK or SA index, the third column of HARPS data), TJ(Ntj=3,vary with FWHM, BIS, RHK), TARMA, TGP, ARMATJ(ARMA+TJ), GPTJ,TJAR(model the RV contributed by index as a AR(p)-like model), ARMATJAR(ARMA+TJAR), PSID (previous-subsequent index dependent, this model is similar to TJAR but without time-varying/index-dependent jitter), ARMAPSID(ARMA+PSID)
    period.par <- 'logp'
    Ncores <- 1
    Np <- 1
    mode <- 'data'#data,sim
    Dtye <- 'D'#Dtype: DE:differential exclusing the target aperture, D: different including all aperture, N: no dependence on differential RV
    Nw <- 1#fit models to multiple wavelength data sets simultaneously
    prior.type <- 'mt'
    calibration <- 0#
    commandArgs <- function(trailingOnly=TRUE){
        cat('args=',c(id,Niter,Nbin.per,nbin.per,tem,inicov,Pini,noise.model,period.par,Ncores,Np,mode,Dtype,Nw,prior.type,calibration)
           ,'\n')
        c(id,Niter,Nbin.per,nbin.per,tem,inicov,Pini,noise.model,period.par,Ncores,Np,mode,Dtype,Nw,prior.type,calibration)
    }
    source('../mcmc_red.R',local=TRUE)
    return(list(folder=folder,pdf=gsub('.+/','',pdf.name)))
}

data.distr <- function(x,xlab,ylab,main='',oneside=FALSE,plotf=TRUE){
    xs <- seq(min(x),max(x),length.out=1e3)
    fitnorm <- fitdistr(x,"normal")
    p <- hist(x,plot=FALSE)
    xfit <- length(x)*mean(diff(p$mids))*dnorm(xs,fitnorm$estimate[1],fitnorm$estimate[2])
    ylim <- range(xfit,p$counts)
    if(plotf){
        plot(p,xlab=xlab,ylab=ylab,main=main,ylim=ylim)
        lines(xs,xfit,col='red')
    }
    x1=Mode(x)
    x2=mean(x)
    x3=sd(x)
    x4=skewness(x)
    x5=kurtosis(x)
    xs = sort(x)
    x1per = max(min(xs),xs[floor(length(xs)*0.01)])
    x99per = min(xs[ceiling(length(xs)*0.99)],max(xs))
#    abline(v=c(x1per,x99per),col='blue')
    if(plotf){
        if(!oneside){
            legend('topleft',legend=c(as.expression(bquote('mode ='~.(format(x1,digit=3)))),as.expression(bquote(mu~'='~.(format(x2,digit=3)))),as.expression(bquote(sigma~'='~.(format(x3,digit=3))))),bty='n')
            legend('topright',legend=c(as.expression(bquote(mu^3~'='~.(format(x4,digit=3)))),as.expression(bquote(mu^4~'='~.(format(x5,digit=3))))),bty='n')
        }else{
            legend('topleft',legend=c(as.expression(bquote('mode ='~.(format(x1,digit=3)))),as.expression(bquote(mu~'='~.(format(x2,digit=3)))),as.expression(bquote(sigma~'='~.(format(x3,digit=3)))),as.expression(bquote(mu^3~'='~.(format(x4,digit=3)))),as.expression(bquote(mu^4~'='~.(format(x5,digit=3))))),bty='n')
        }
    }
    return(c(x1per=x1per,x99per=x99per,mode=x1,mean=x2,sd=x3,skewness=x4,kurtosis=x5))
}
show.peaks <- function(ps,powers,levels=NULL,Nmax=5){
    if(is.null(levels)) levels <- max(max(powers)-log(150),median(powers))
    ind <- which(powers==max(powers) | (powers>(max(powers)-log(100)) & powers>max(levels)))
    if(max(powers)-min(powers)<5) ind <- which.max(powers)
    pmax <- ps[ind]
    ppmax <- powers[ind]
    j0 <- 1
    p0 <- pmax[1]
    pp0 <- ppmax[1]
    pms <- p0
    pos <- pp0
    if(length(pmax)>1){
        for(j in 2:length(pmax)){
            if(abs(pmax[j]-p0) < 0.1*p0){
                if(ppmax[j]>pp0){
                    j0 <- j
                    p0 <- pmax[j]
                    pp0 <- ppmax[j0]
                    pms[length(pms)] <- p0
                    pos[length(pos)] <- pp0
                }
			    }else{
                j0 <- j
                p0 <- pmax[j]
                pp0 <- ppmax[j0]
                pms <- c(pms,p0)
                pos <- c(pos,pp0)
            }
        }
    }else{
        pms <- pmax
        pos <- ppmax
    }
    if(length(pms)>Nmax){
      pms <- pms[1:Nmax]
      pos <- pos[1:Nmax]
    }
    return(cbind(pms,pos))
}
