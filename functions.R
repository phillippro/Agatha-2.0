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
###number of harmonics of the signal in BFP/MLP; 2 or more fits eccentric orbits
    if(!exists('Nh')) Nh <- 1
    Nh <- max(1,as.integer(Nh))
###red-noise model: 'ARMA' (per-set AR/MA orders) or 'GP' (shared SHO Gaussian process)
    if(!exists('noise.model')) noise.model <- 'ARMA'
    GP <- noise.model=='GP'
###gp.par = c(sigmaGP, logProt, logtauGP): NA = free; the rotation period and the
###coherence time scale can be fixed from the panel (e.g. a photometric Prot)
    gp.par <- rep(NA,3)
    if(exists('gp.Prot') && length(gp.Prot)==1 && is.finite(gp.Prot) && gp.Prot>0) gp.par[2] <- log(gp.Prot)
    if(exists('gp.tau') && length(gp.tau)==1 && is.finite(gp.tau) && gp.tau>0) gp.par[3] <- log(gp.tau)
    if(GP){
        Nmas <- rep(0,length(Nmas))
        Nars <- rep(0,length(Nars))
        if(mcf){
            warning('MCMC with a GP noise model is not implemented; using the maximum-likelihood fit.')
            mcf <- FALSE
        }
    }
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
                    rv.ls <- BFP.multiset(t=t,y=y,dy=dy,set.id=set.id,Nma=Nma,Nar=Nar,ofac=ofac,fmin=frange[1],fmax=frange[2],progress=FALSE,noise.only=noise.only,Nh=Nh,noise.model=noise.model,gp.par=gp.par)
                }else{
                    rv.ls <- BFP(t=t,y=y,dy=dy, Nma=Nma,Nar=Nar,model.type='man',Indices=Indices, ofac=ofac,fmin=frange[1],fmax=frange[2],quantify=quantify, renew=renew,Nsamp=Nsamp,noise.only=noise.only,Nh=Nh,GP=GP,gp.par=gp.par)
                }
###renew: every chi-square minimization start from the initial parameter values
                ylab <- 'ln(BF)'
                name <- 'logBF'
            }else if(per.type=='MLP'){
                if(multi.set){
                    rv.ls <- MLP.multiset(t=t,y=y,dy=dy,set.id=set.id,Nma=Nma,Nar=Nar,ofac=ofac,fmin=frange[1],fmax=frange[2],Nh=Nh,noise.model=noise.model,gp.par=gp.par)
                }else{
                    rv.ls <- MLP(t=tab[,1]-min(tab[,1]),y=y,dy=dy,Nma=Nma,Nar=Nar,Indices=Indices,ofac=ofac,fmin=frange[1],fmax=frange[2],MLP.type=MLP.type,Nh=Nh,GP=GP,gp.par=gp.par)
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
            seed <- fourier_kepler_seed(par,popt)
            if(!is.null(seed) && is.finite(seed$K1) && seed$K1>0){
###analytical eccentric start from the fundamental and first harmonic
                startvalue <- c(log(popt),seed$K1,seed$e1,seed$omega1,seed$Mo1)
            }else{
                phi <- as.numeric(xy2phi(par['A'],par['B']))
                kopt <- as.numeric(sqrt(par['A']^2+par['B']^2))
                startvalue <- c(log(popt),kopt,0,0,phi)
            }
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
                  mcmc.method='PT',Ntem=NULL,tem.min=NULL,swap.interval=10,mcmc.verbose=FALSE){
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
###peak and hand good states to the cold chain through replica exchange. With
###Ntem=NULL and tem.min=NULL the ladder is chosen from the data and adapted
###during burn-in, and the chain is extended until the cold chain converges
            tmp <- run.ptmcmc(startvalue,cov.start,iterations=max(as.numeric(Niter),1000),
                              bases=rep(basis,10),Ntem=Ntem,tem.min=tem.min,
                              swap.interval=swap.interval,verbose=mcmc.verbose)
            if(mcmc.verbose){
                cat('PTMCMC tem.min:',format(tmp$tem.min,digit=2),if(tmp$auto.tem) '(automatic)' else '(fixed)',
                    '; rungs:',tmp$Ntem,'; iterations:',tmp$iterations,'; extended:',tmp$extended,'blocks\n')
                cat('PTMCMC initial ladder:',paste(format(tmp$tems.initial,digit=2),collapse=','),'\n')
                cat('PTMCMC adapted ladder:',paste(format(tmp$tems,digit=2),collapse=','),'\n')
                cat('PTMCMC max Rhat:',if(is.null(tmp$Rhat)) NA else round(max(tmp$Rhat),3),'\n')
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
        ysim.sig <- harmonic_signal(par.opt,2*pi/popt,tsim)
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
                   mcmc.method='PT',Ntem=NULL,tem.min=NULL,swap.interval=10,mcmc.verbose=FALSE){
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
                ysim.sig <- harmonic_signal(par.opt,2*pi/popt,tsim)
                ysim.all <- ysim.sig
                if(any(names(par.opt)=='gamma')) ysim.all <- ysim.all+par.opt['gamma']
                if(any(names(par.opt)=='beta')) ysim.all <- ysim.all+par.opt['beta']*tsim
                ysig0 <- harmonic_signal(par.opt,2*pi/popt,t)
                if(isTRUE(per$df$GP)){
                    per$res <- per$res-gp_predict_par(t,data[,3],par.opt,per$df,per$res)
                }
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
                if(isTRUE(df$GP)) df$predict.gp <- TRUE
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
                if(any(names(per$par.opt)=='logProt')){
                    per$Popt <- exp(per$par.opt[['logProt']])
                }else{
                    per$Popt <- exp(per$par.opt[grep('logtau',names(per$par.opt))])
                }
                popt <- per$Popt[1]
                df <- per$df
                fit <- CircularSig(par.list,df)
                ysig0 <- fit$yred
                if(isTRUE(df$GP)){
                    ygp <- gp_predict_par(t,data[,3],par.list,df,fit$res)
                    ysig0 <- ysig0+ygp
                    fit$res <- fit$res-ygp
                }
                ysig <- fit$res+ysig0
                res <- per$res <- fit$res
                sim <- CircularSim(par.list,df,popt,tsim)
                ysim.sig <- sim$yred
                if(isTRUE(df$GP)) ysim.sig <- ysim.sig+approxfun(t,ygp,rule=2)(tsim)
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

###########################################################################
####Publication-quality single-panel plots for the 1D periodogram results.
####Every figure of the app is drawn by one of the panel.* functions below, so
####the on-screen figure, the bundled PDF and an individually downloaded panel
####are the same drawing. Colours follow the Okabe-Ito palette.
###########################################################################
###tcol() lives in mcmc_func.R, which the app sources after this file; keep a
###local fallback so the plotting layer works on its own
if(!exists('tcol')){
    tcol <- function(color, percent=50, name=NULL){
        rgb.val <- col2rgb(color)
        rgb(rgb.val[1],rgb.val[2],rgb.val[3],max=255,alpha=(100-percent)*255/100,names=name)
    }
}
pub.col <- list(model='#D55E00',peak='#0072B2',data='black',err=tcol('black',50),level='grey40',res='grey20')

pub.par <- function(mfrow=c(1,1),cex=1){
    par(mfrow=mfrow,mar=c(4.5,5,2.6,1.2),mgp=c(2.9,0.7,0),las=1,tcl=-0.4,
        cex=cex,cex.lab=1.25,cex.axis=1.05,cex.main=1.1,font.main=1,xaxs='r',yaxs='r')
}

pretty.title <- function(tit){
###'BFP; combined;RV;1 signal' -> 'BFP, RV (combined), signal 1'
    f <- trimws(unlist(strsplit(tit,';')))
    if(length(f)<3) return(tit)
    obs <- f[3]
    if(grepl('Window',obs)) return(paste0('Window function (',f[2],')'))
    sig <- if(length(f)>=4) gsub('(\\d+) signal','signal \\1',f[4]) else ''
    paste0(f[1],', ',obs,' (',f[2],')',if(nzchar(sig)) paste0(', ',sig) else '')
}

ylab.obs <- function(ypar){
    if(ypar=='RV') expression(RV~'[m s'^-1*']') else ypar
}

ylab.power <- function(cname){
###from a per.list column name such as 'BFP1signal:RV:logBF'
    if(grepl('logBF',cname)) 'ln BF' else if(grepl('logML',cname)) expression(ln(ML/ML[max])) else 'Power'
}

pub.axes <- function(xlog=FALSE){
###minor ticks on a log axis via magicaxis when it is installed (the app loads
###it); plain axes otherwise, so the plotting layer has no hard dependency
    if(xlog && requireNamespace('magicaxis',quietly=TRUE)){
        magicaxis::magaxis(side=1,tcl=-0.4,cex.axis=par('cex.axis'))
    }else{
        axis(1)
    }
    axis(2)
    box()
}

panel.periodogram <- function(per.list,ypar,i,title,levels=NULL,SigType='circular',pub=TRUE,Pmark=NULL){
###Pmark: the period reported by the fit; with harmonics the raw maximum can sit
###at the 2P alias, so the reported period is what gets annotated
###per.list[[ypar]]: column 1 = period, columns 2.. = powers; i selects the power column
    P <- per.list[[ypar]][,1]
    power <- per.list[[ypar]][,i+1]
    cname <- colnames(per.list[[ypar]])[i+1]
    ylab <- ylab.power(cname)
    window <- grepl('Window',title)
    ymin <- if(SigType!='stochastic') median(power) else max(0,min(power))
    ylim <- c(ymin,max(power)+0.18*(max(power)-ymin))
    xlab <- if(SigType=='stochastic') 'Time scale [day]' else 'Period [day]'
    plot(P,power,type='n',log='x',xaxt='n',yaxt='n',xlab=xlab,ylab=ylab,ylim=ylim,
         main=if(pub) pretty.title(title) else title)
    pub.axes(xlog=TRUE)
    if(!window && !is.null(levels)){
        lv <- levels[is.finite(levels)]
        if(length(lv)>0) abline(h=lv,lty=2,col=pub.col$level,lwd=1)
    }
    lines(P,power,lwd=if(pub) 1.6 else 1,col=pub.col$data)
    kraw <- which.max(power)
    k <- kraw
    if(!is.null(Pmark) && is.finite(Pmark) && Pmark>0){
        k <- which.min(abs(log(P)-log(Pmark)))
    }
    pmax <- P[k]; wmax <- power[k]
    dyv <- max(power)-ymin
    alias <- abs(log(P[kraw]/pmax))>0.01
    par(xpd=TRUE)
    lab <- if(SigType=='stochastic') paste0('tau = ',format(pmax,digits=4),' d') else paste0('P = ',format(pmax,digits=4),' d')
    if(!alias){
###the reported period is the maximum: arrow and label just above the peak
        text(pmax,wmax+0.08*dyv,pos=3,labels=lab,col=pub.col$peak,cex=1.0)
        try(arrows(pmax,wmax+0.08*dyv,pmax,wmax+0.02*dyv,col=pub.col$peak,length=0.05,lwd=1.5),TRUE)
    }else{
###the reported period is not the maximum (with harmonics the maximum can be
###the 2P alias): a full-height marker so it is visible whatever its power,
###the label in the top margin, and the alias named for what it is
        segments(pmax,ymin,pmax,max(power)+0.05*dyv,col=pub.col$peak,lty=3,lwd=1.5)
        text(pmax,max(power)+0.06*dyv,pos=3,labels=lab,col=pub.col$peak,cex=1.0)
        points(P[kraw],power[kraw],pch=1,col=pub.col$level,cex=1.2)
        ratio <- P[kraw]/pmax
        alab <- if(abs(ratio-2)<0.1) '2P alias' else if(abs(ratio-0.5)<0.03) 'P/2 alias' else 'alias'
        text(P[kraw],power[kraw],labels=paste0(alab,' (',format(P[kraw],digits=3),' d)'),pos=4,col=pub.col$level,cex=0.85)
    }
    par(xpd=FALSE)
    invisible(list(Popt=pmax,power=wmax))
}

panel.phase <- function(phase.list,sim.list,ypar,i,title,pub=TRUE){
    ph <- phase.list[[ypar]]; sm <- sim.list[[ypar]]
    t <- ph[,paste0('t_sig',i)]; y <- ph[,paste0('y_sig',i)]; ey <- ph[,'ey0']
    tsim <- sm[,paste0('tsim_sig',i)]; ysim <- sm[,paste0('ysim_sig',i)]
    ylim <- range(c(y-ey,y+ey,ysim),na.rm=TRUE)
    plot(t,y,type='n',xaxt='n',yaxt='n',xlab='Phase [day]',ylab=ylab.obs(ypar),ylim=ylim,
         main=if(pub) pretty.title(title) else title)
    pub.axes()
    arrows(t,y-ey,t,y+ey,length=0.02,angle=90,code=3,col=pub.col$err)
    lines(tsim,ysim,col=pub.col$model,lwd=2.5)
    points(t,y,pch=21,bg='white',col=pub.col$data,cex=0.9)
}

panel.fit <- function(phase.list,sim.list,ypar,title,pub=TRUE){
    ph <- phase.list[[ypar]]; sm <- sim.list[[ypar]]
    t <- ph[,'t0']; ey <- ph[,'ey0']; y <- ph[,'y_all']
    tsim <- sm[,'tsim0']+min(t); ysim <- sm[,'ysim_all']
    ylim <- range(c(y-ey,y+ey,ysim),na.rm=TRUE)
    plot(t,y,type='n',xaxt='n',yaxt='n',xlab='Time [day]',ylab=ylab.obs(ypar),ylim=ylim,
         main=if(pub) gsub(', signal \\d+',', combined fit',pretty.title(title)) else gsub('\\d signal','combined fit',title))
    pub.axes()
    lines(tsim,ysim,col=pub.col$model,lwd=2)
    try(arrows(t,y-ey,t,y+ey,length=0.02,angle=90,code=3,col=pub.col$err),TRUE)
    points(t,y,pch=21,bg='white',col=pub.col$data,cex=0.9)
    legend('topleft',bg='white',box.col='white',legend=paste0('RMS = ',format(sd(y),digits=3)),text.col=pub.col$peak)
}

panel.residual <- function(phase.list,ypar,title,pub=TRUE){
    ph <- phase.list[[ypar]]
    t <- ph[,'t0']; ey <- ph[,'ey0']; res <- ph[,'res_all']
    ylim <- range(c(res-ey,res+ey),na.rm=TRUE)
    plot(t,res,type='n',xaxt='n',yaxt='n',xlab='Time [day]',ylab=paste('Residual',if(ypar=='RV') '[m/s]' else ''),ylim=ylim,
         main=if(pub) gsub(', signal \\d+',', residual',pretty.title(title)) else gsub('\\d signal','residual',title))
    pub.axes()
    abline(h=0,lty=3,col=pub.col$level)
    try(arrows(t,res-ey,t,res+ey,length=0.02,angle=90,code=3,col=pub.col$err),TRUE)
    points(t,res,pch=21,bg='white',col=pub.col$res,cex=0.9)
    legend('topleft',bg='white',box.col='white',legend=paste0('RMS = ',format(sd(res),digits=3)),text.col=pub.col$peak)
}

list.single.plots <- function(d){
###every panel available from a calc.1Dper() result, as a table used by the
###download selector: label, kind, observable, index
    out <- c()
    for(ypar in names(d$per.list)){
        titles <- d$tits[grepl(paste0(';',ypar,';'),d$tits)]
        np <- ncol(d$per.list[[ypar]])-1
        for(i in 1:np){
            out <- rbind(out,data.frame(label=paste0('Periodogram: ',pretty.title(titles[i])),kind='periodogram',ypar=ypar,index=i,stringsAsFactors=FALSE))
        }
        if(!is.null(d$phase.list[[ypar]])){
            ns <- length(grep('^y_sig',colnames(d$phase.list[[ypar]])))
            for(i in 1:ns){
                out <- rbind(out,data.frame(label=paste0('Phase plot: ',pretty.title(titles[i])),kind='phase',ypar=ypar,index=i,stringsAsFactors=FALSE))
            }
            out <- rbind(out,data.frame(label=paste0('Combined fit: ',ypar),kind='fit',ypar=ypar,index=1,stringsAsFactors=FALSE))
            out <- rbind(out,data.frame(label=paste0('Residual: ',ypar),kind='residual',ypar=ypar,index=1,stringsAsFactors=FALSE))
        }
    }
    out
}

plot1D.single <- function(d,kind,ypar,index=1,SigType='circular',pub=TRUE){
###draw exactly one panel of a calc.1Dper() result on the current device
    pub.par()
    titles <- d$tits[grepl(paste0(';',ypar,';'),d$tits)]
    if(kind=='periodogram'){
        gi <- level.index(d,ypar,index)
        panel.periodogram(d$per.list,ypar,index,titles[index],levels=d$levels[,gi],SigType=SigType,pub=pub,
                          Pmark=reported.period(d$par.list,ypar,index))
    }else if(kind=='phase'){
        panel.phase(d$phase.list,d$sim.list,ypar,index,titles[index],pub=pub)
    }else if(kind=='fit'){
        panel.fit(d$phase.list,d$sim.list,ypar,titles[1],pub=pub)
    }else if(kind=='residual'){
        panel.residual(d$phase.list,ypar,titles[1],pub=pub)
    }
}

level.index <- function(d,ypar,i){
###column of d$levels for the i-th periodogram of observable ypar: levels are
###stored in the order the periodograms were computed, observable by observable
    gi <- 0
    for(yp in names(d$per.list)){
        np <- ncol(d$per.list[[yp]])-1
        if(yp==ypar) return(gi+i)
        gi <- gi+np
    }
    gi+i
}

save.single.plot <- function(file,format='png',width=6,height=4.5,dpi=300,expr){
###open the requested device, evaluate the drawing expression, close it
    format <- tolower(format)
    if(format=='pdf'){
        pdf(file,width=width,height=height,useDingbats=FALSE)
    }else if(format=='png'){
        png(file,width=width,height=height,units='in',res=dpi)
    }else if(format %in% c('jpg','jpeg')){
        jpeg(file,width=width,height=height,units='in',res=dpi,quality=95)
    }else{
        stop('unknown figure format: ',format)
    }
    on.exit(dev.off())
    force(expr)
    invisible(file)
}

phase1D.plot <- function(phase.list,sim.list,tits,download=FALSE,index=NULL,repar=TRUE,pub=TRUE){
    if(repar){
        if(is.null(index)){
            pub.par(mfrow=c(ceiling(Nmax.plots/2),2))
        }
        if(download & is.null(index)){
            pub.par(mfrow=c(2,2),cex=0.9)
        }
    }
    for(ypar in names(phase.list)){
        if(!is.null(index)){
            inds <- index
        }else{
            inds <- 1:length(grep('^y_sig',colnames(phase.list[[ypar]])))
        }
        titles <- tits[grepl(paste0(';',ypar,';'),tits)]
        for(i in inds){
            panel.phase(phase.list,sim.list,ypar,i,titles[i],pub=pub)
        }
        panel.fit(phase.list,sim.list,ypar,titles[max(inds)],pub=pub)
        panel.residual(phase.list,ypar,titles[max(inds)],pub=pub)
    }
}

combined.plot <- function(per.list,phase.list,sim.list,tits,pers,levels,ylabs,SigType='circular',download=FALSE,index=NULL,pub=TRUE,par.list=NULL){
    per1D.plot(per.list,tits,pers,levels,ylabs,download=download,index=index,SigType=SigType,pub=pub,par.list=par.list)
    phase1D.plot(phase.list,sim.list,tits=tits,download=download,index=index,repar=FALSE,pub=pub)
}

per1D.plot <- function(per.list,tits,pers,levels,ylabs,download=FALSE,index=NULL,SigType='circular',pub=TRUE,par.list=NULL){
    if(is.null(index)){
        pub.par(mfrow=c(ceiling(Nmax.plots/2),2))
    }
    if(download & is.null(index)){
        pub.par(mfrow=c(2,2),cex=0.9)
    }
    gi <- 0
    for(ypar in names(per.list)){
        np <- ncol(per.list[[ypar]])-1
        titles <- tits[grepl(paste0(';',ypar,';'),tits)]
        inds <- if(!is.null(index)) index else 1:np
        for(i in inds){
###significance levels are stored per periodogram across all observables
            lv <- if(is.matrix(levels) && ncol(levels)>=gi+i) levels[,gi+i] else NULL
            Pmark <- reported.period(par.list,ypar,i)
            panel.periodogram(per.list,ypar,i,titles[i],levels=lv,SigType=SigType,pub=pub,Pmark=Pmark)
        }
        gi <- gi+np
    }
}

reported.period <- function(par.list,ypar,i){
###period of the i-th signal in the fitted parameters, if available
    if(is.null(par.list) || is.null(par.list[[ypar]])) return(NULL)
    pl <- par.list[[ypar]]
    nm <- paste0('P',i)
    v <- if(is.matrix(pl)) { if(nm%in%colnames(pl)) pl['mode',nm] else NULL } else if(nm%in%names(pl)) pl[[nm]] else NULL
    if(is.null(v) || !is.finite(v)) NULL else as.numeric(v)
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
    Nma <- as.integer(pars[[i]]$Nma)
    Nar <- as.integer(pars[[i]]$Nar)
    Inds.sets <- Inds
    if(length(per.target)==1){
        Inds <- as.integer(pars[[i]]$Inds)
    }else{
###with several data sets the proxies and MA terms are applied per set inside
###combine.data(); the combined residual series is then analysed without them
        Inds <- 0
        Nma <- 0
    }
    Indices <- NULL
    per.type <- pars[[i]]$per.type
    var <- pars[[i]]$var
    if(length(per.target)>1){
        instrument <- 'combined'
        subdata <- lapply(1:length(per.target),function(j) data[[per.target[j]]])
        if(!is.list(Inds.sets)) Inds.sets <- rep(list(Inds.sets),length(per.target))
        tmp <- combine.data(data=subdata,Ninds=Inds.sets,Nmas=Nmas)
        tab <- tmp$cdata
        idata <- tmp$idata
        colnames(tab) <- colnames(data[[per.target[1]]])[1:3]
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
