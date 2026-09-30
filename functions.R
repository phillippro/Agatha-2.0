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

is.signal.par <- function(nm){
###names of the signal parameters (as opposed to the noise, offset and trend
###terms) in either the periodogram or the MCMC naming
    grepl('^(P|per|K|e|omega|Mo|A|B|lnK|sqrecosw|sqresinw|Tc|dP|dTc)[0-9]+$',nm)
}

addpar <- function(par.old,par.new,nsig){
###merge the parameters of the nsig-th signal (par.new: signal terms first,
###then its noise terms) into those of the previous signals: the new signal's
###terms are renamed with its index, the previous signals' terms are kept and
###the noise terms of the latest fit replace the earlier ones
    n0 <- names(par.new)
    if(nsig==1){
        if(any(n0=='A')|any(n0=='B')){
            names(par.new)[1:3] <- paste0(n0[1:3],1)
        }
        return(par.new)
    }
    n1 <- n0
    if(any(n0=='A')|any(n0=='B')){
        n1[1:3] <- paste0(n0[1:3],nsig)
    }else if(any(n0=='Mo1') & !any(n0=='omega1')){
        n1[1:3] <- gsub('1$',nsig,n0[1:3])
    }else if(any(grepl('omega|Tc',n0))){
        n1[1:5] <- gsub('1$',nsig,n0[1:5])
    }
    names(par.new) <- n1
    c(par.old[is.signal.par(names(par.old))],par.new)
}

addpar.stat <- function(stat.old,stat.new,nsig,keep=NULL){
###addpar() for the MCMC statistics of sequential signals: every statistic (row)
###is renamed and merged the same way as the parameter vector, so the table of
###the combined fit keeps the posterior summaries of each signal. keep (the
###names of the fitted parameters) drops derived columns such as reference
###epochs, which the parameter vector does not carry
    if(is.null(stat.new) || !is.matrix(stat.new)) return(NULL)
    if(nsig>1 && (is.null(stat.old) || !is.matrix(stat.old))) return(NULL)
    if(!is.null(keep)){
        keep <- keep[keep%in%colnames(stat.new)]
        if(length(keep)==0) return(NULL)
        stat.new <- stat.new[,keep,drop=FALSE]
    }
    rows <- if(nsig==1) rownames(stat.new) else intersect(rownames(stat.old),rownames(stat.new))
    if(length(rows)==0) return(NULL)
    out <- t(sapply(rows,function(r) addpar(if(nsig==1) c() else stat.old[r,],stat.new[r,],nsig)))
    rownames(out) <- rows
    out
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
    if(!exists('progress')) progress <- FALSE
    Nmas <- unlist(Nmas)
    Nars <- unlist(Nars)
###number of harmonics of the signal in BFP/MLP; 2 or more fits eccentric orbits
    if(!exists('Nh')) Nh <- 1
    Nh <- max(1,as.integer(Nh))
###eccentricity prior of the MCMC (list(type,a,b,sigma); see eprior.check)
    e.prior <- eprior.check(if(exists('e.prior')) e.prior else NULL)
###ln(BF) a further signal must exceed to be accepted by the sequential model comparison
    if(!exists('lnBF.min') || length(lnBF.min)!=1 || !is.finite(lnBF.min)) lnBF.min <- 5
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
        if(mcf && length(per.target)==1){
###the single-set GP path samples with the celerite likelihood, which the MCMC
###wrapper does not support; the multi-set GP MCMC works on the whitened data
            warning('MCMC with a GP noise model is not implemented for a single data set; using the maximum-likelihood fit.')
            mcf <- FALSE
        }
    }
    par.list <- sim.list <- phase.list <- per.list <- tits <- mc.list <- model.list <- list()
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

####calculate the Keplerian fit; periods are reported linearly (P1, P2, ...)
###and converted back to the sampling form by par.P2per() when a fit is chained
        Pconv <- TRUE
###the window function is not a signal to refine: no Keplerian, no MCMC
        mcf1 <- mcf && ypar!='Window Function'
        SigType1 <- if(ypar=='Window Function') 'circular' else SigType
        fit <- sigfit(per=rv.ls,data=tab,SigType=SigType1,basis=basis,Ncores=Ncores,mcf=mcf1,Niter=Niter,Pconv=Pconv,e.prior=e.prior)
###update the output from periodogram
        rv.ls <- fit$per
###model comparison: ln(BF) of the first signal against no signal, from the
###periodogram (BIC estimate; the Keplerian fit's own value with several sets)
        model.comp <- NULL
        if(ypar!='Window Function'){
            lnBF1 <- if(!is.null(rv.ls$lnBF.kepler) && is.finite(rv.ls$lnBF.kepler)) as.numeric(rv.ls$lnBF.kepler) else
                     if(per.type=='BFP') suppressWarnings(max(rv.ls$power,na.rm=TRUE)) else NA
            if(!is.finite(lnBF1)) lnBF1 <- NA
            model.comp <- data.frame(n=1,P=as.numeric(fit$popt)[1],lnBF=lnBF1,source='periodogram',
                                     accepted=is.na(lnBF1) || lnBF1>lnBF.min,stringsAsFactors=FALSE)
        }
        llmax.prev <- if(is.null(fit$llmax)) NA else as.numeric(fit$llmax)
        par.prev <- fit$ParSig
        res.last <- fit$res

        pp <- cbind(fit$t,fit$y,fit$ysig0)

        colnames(pp) <- paste0(c('t','y','ysig'),'_sig1')
        phase.data <- cbind(phase.data,pp)

        qq <- cbind(fit$tsim,fit$ysim,fit$ysim0)
        colnames(qq) <- paste0(c('tsim','ysim','ysim0'),'_sig1')
        sim.data <- cbind(sim.data,qq)
        par.data <- addpar(c(),fit$ParSig,1)
        par.stat.data <- addpar.stat(NULL,fit$par.stat,1,keep=names(fit$ParSig))

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
###with several data sets the sequential signals are each refined by their own
###multi-set MCMC, but there is no joint re-fit of all signals, so the combined
###model is assembled from the per-signal fits in the else branch below
        mcf2 <- mcf1 & !(multi.set & Nsig.max>1)
###number of signals kept in the combined model (those accepted by the model comparison)
        Nacc <- length(grep('^ysig_sig',colnames(phase.data)))
        if(mcf2){
###with a single data set the accepted signals were refined together by the
###joint MCMC of the sequential search (joint.fit); with one signal the fit of
###that signal is already the joint fit
            if(Nsig.max>1 & !multi.set & exists('joint.fit')){
                fit <- joint.fit
            }
            mc <- fit$mc
            ParSig <- fit$ParSig
            par.data <- fit$par.stat

            pp <- cbind(fit$ysig,fit$ysig0,fit$res)
            colnames(pp) <- paste0(c('y','ysig','res'),'_all')
            phase.data <- cbind(phase.data,pp)
            if(!is.null(fit$ysim.sig)){
                qq <- cbind(fit$ysim.sig)
            }else{
                qq <- cbind(fit$ysim0)
            }
            colnames(qq) <- 'ysim_all'
            sim.data <- cbind(sim.data,qq)
        }else{
###the residual after the last accepted (sequentially found and refined) signal
            res <- as.numeric(res.last)
            if(Nacc>1){
                ysig0 <- rowSums(phase.data[,paste0('ysig_sig',1:Nacc)])
                ysim <- rowSums(sim.data[,paste0('ysim0_sig',1:Nacc)])
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
###sequentially refined signals: report the posterior summaries of every signal
            if(!is.null(par.stat.data) && all(names(par.data)%in%colnames(par.stat.data))){
                par.data <- par.stat.data[,names(par.data),drop=FALSE]
            }
        }

###the data-set names travel with the parameters so that tables can label them
        attr(par.data,'sets') <- per.target
        if(!is.null(model.comp)){
            attr(model.comp,'lnBF.min') <- lnBF.min
            attr(model.comp,'Nsig.max') <- Nsig.max
        }
###attach common data; set0 is the index of the data set each point came from
        set0 <- if(multi.set) as.integer(factor(set.id,levels=per.target)) else rep(1L,length(t))
        phase.attach <- cbind(t,y,dy,set0)
        colnames(phase.attach) <- c('t0','y0','ey0','set0')
        sim.attach <- t(t(fit$tsim0))
        colnames(sim.attach) <- 'tsim0'
        phase.data <- cbind(phase.data,phase.attach)
        attr(phase.data,'sets') <- per.target
        sim.data <- cbind(sim.data,sim.attach)
        colnames(per.data) <- cnames

###put everything into list
        sim.list[[ypar]] <- sim.data
        phase.list[[ypar]] <- phase.data
        per.list[[ypar]] <- per.data
        par.list[[ypar]] <- par.data
        mc.list[[ypar]] <- mc
        model.list[[ypar]] <- model.comp
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
    return(list(per.list=per.list,mc=mc,phase.list=phase.list,sim.list=sim.list,par.list=par.list,tits=tits,pers=pers,levels=sig.levels,ylabs=ylabs,fname=fname,fs=fs,mc.list=mc.list,
                model.comp=model.list,Nopt=lapply(model.list,model.Nopt)))
}

model.Nopt <- function(mc){
###the most plausible number of signals: the signals accepted in sequence by
###the model comparison (all of them when no comparison was possible)
    if(is.null(mc) || nrow(mc)==0) return(NA)
    acc <- cumprod(as.logical(mc$accepted))
    sum(acc)
}

par.P2per <- function(par){
###periods reported linearly ('P1') back to the sampling form ('per1' = ln P)
    ip <- grep('^P\\d+$',names(par))
    if(length(ip)>0){
        par[ip] <- log(par[ip])
        names(par)[ip] <- sub('^P','per',names(par)[ip])
    }
    par
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

###########################################################################
####Eccentricity prior of the MCMC fits
eprior.default <- function() list(type='beta',a=0.867,b=3.03,sigma=0.1)

eprior.check <- function(ep=NULL){
###a complete, valid eccentricity prior description from a possibly partial one:
###type 'beta' (Kipping 2013, MNRAS 434, L51: a=0.867, b=3.03 by default),
###'uniform' on [0,1), or 'halfgauss' (half-Gaussian of scale sigma on e>=0)
    d <- eprior.default()
    if(is.null(ep) || !is.list(ep)) return(d)
    ok <- function(v) length(v)==1 && is.finite(v) && v>0
    type <- if(is.null(ep$type)) d$type else as.character(ep$type)[1]
    type <- c(beta='beta',kipping='beta',uniform='uniform',flat='uniform',halfgauss='halfgauss',
              gaussian='halfgauss',gauss='halfgauss')[tolower(type)]
    if(is.na(type)) type <- d$type
    list(type=unname(type),a=if(ok(ep$a)) as.numeric(ep$a) else d$a,b=if(ok(ep$b)) as.numeric(ep$b) else d$b,
         sigma=if(ok(ep$sigma)) as.numeric(ep$sigma) else d$sigma)
}

eprior.log <- function(e,ep=NULL){
###log prior density of the eccentricity e (scalar or vector; the sum is
###returned). Outside [0,1) the density is zero. The beta density with a<1
###diverges at e=0, so e is evaluated no closer than 1e-4 to the edges
    ep <- eprior.check(ep)
    if(length(e)==0) return(0)
    if(any(!is.finite(e)) || any(e<0 | e>=1)) return(-Inf)
    e <- pmin(pmax(e,1e-4),1-1e-4)
    switch(ep$type,
           uniform=0,
           beta=sum(dbeta(e,ep$a,ep$b,log=TRUE)),
           halfgauss=sum(log(2)+dnorm(e,mean=0,sd=ep$sigma,log=TRUE)),
           stop('unknown eccentricity prior: ',ep$type))
}

eprior.describe <- function(ep=NULL){
    ep <- eprior.check(ep)
    switch(ep$type,
           uniform='uniform on [0,1)',
           beta=paste0('beta(a=',format(ep$a,digits=3),', b=',format(ep$b,digits=3),'; Kipping 2013)'),
           halfgauss=paste0('half-Gaussian with sigma=',format(ep$sigma,digits=3)))
}

#mcfit <- function(startvalue,Niter,Ncores=1){
mcfit <- function(per,data,tsim,Niter=1e3,SigType='kepler',basis='natural',ParSig=NULL,Pconv=FALSE,Ncores=4,
                  mcmc.method='PT',Ntem=NULL,tem.min=NULL,swap.interval=10,mcmc.verbose=FALSE,e.prior=NULL){
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
###eccentricity prior seen by prior.func() (mcmc_func.R is sourced locally below)
    e.prior <- eprior.check(e.prior)
    if(mcmc.verbose) cat('eccentricity prior:',eprior.describe(e.prior),'\n')
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
###a forked chain cannot update the interface, so the dialog itself is the
###status indicator; with a single core the cold chain reports its iterations
    prog <- if(Ncores==1) function(k,detail='') incProgress(k,detail=detail) else NULL
    withProgress(message=paste0('Running PT-MCMC (',Ncores,' chain',if(Ncores>1) 's' else '',' x ',
                                max(as.numeric(Niter),1000),' iterations)'), value=0.05, {
    mcmc <- foreach(ncore=1:Ncores,.errorhandling = 'pass') %dopar% {
        if(mcmc.method=='PT'){
###parallel tempering: the hot replicas explore the aliases of the periodogram
###peak and hand good states to the cold chain through replica exchange. With
###Ntem=NULL and tem.min=NULL the ladder is chosen from the data and adapted
###during burn-in, and the chain is extended until the cold chain converges
            tmp <- run.ptmcmc(startvalue,cov.start,iterations=max(as.numeric(Niter),1000),
                              bases=rep(basis,10),Ntem=Ntem,tem.min=tem.min,
                              swap.interval=swap.interval,verbose=mcmc.verbose,progress=prog)
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
    })

#    mcmc  <- list()
#    mcmc[[1]] <- tmp$out
    bad <- sapply(mcmc,function(m) is.null(dim(m)))
    if(all(bad)){
        msg <- tryCatch(conditionMessage(mcmc[[1]]),error=function(e) 'unknown error')
        stop('all MCMC chains failed: ',msg)
    }
    mcmc <- mcmc[!bad]
    mc <- do.call(rbind,mcmc)

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

msmc.stat <- function(x,lp=NULL){
###summary statistics of one MCMC parameter, self-contained so that the result
###does not depend on which data.distr() happens to be visible (functions.R has
###a plotting helper of that name that shadows the mcmc_func.R version when the
###app sources functions.R locally). xopt is the value in the sample with the
###highest lp (the MAP estimate when lp is the log posterior or likelihood).
    q <- as.numeric(quantile(x,c(0.1587,0.8413,0.5,0.01,0.99,0.1,0.9),na.rm=TRUE))
    mode <- tryCatch({ d <- density(x,na.rm=TRUE); d$x[which.max(d$y)] },error=function(e) q[3])
    xopt <- if(is.null(lp)) NA else as.numeric(x[which.max(lp)])
    c(xopt=xopt,xminus.1sig=q[1],xplus.1sig=q[2],mode=as.numeric(mode),med=q[3],
      mean=mean(x,na.rm=TRUE),sd=sd(x),x1per=q[4],x99per=q[5],x10per=q[6],x90per=q[7])
}

mcfit.multiset <- function(per, data, set.id, tsim, Niter=1e3, SigType='kepler', Ncores=4,
                           Ntem=NULL, tem.min=NULL, swap.interval=10, mcmc.verbose=FALSE, e.prior=NULL){
###PT-MCMC refinement of a multi-set fit: a shared signal (Keplerian or
###circular), one offset per data set, a shared linear trend, and one jitter
###per data set. With a GP periodogram the fixed GP covariance from the
###periodogram fit whitens the residuals and replaces the jitters. Per-set
###AR/MA lag terms are not refined here; the jitters absorb what they left.
    t <- as.numeric(data[,1]); t <- t-min(t)
    y <- as.numeric(data[,2])
    dy <- as.numeric(data[,3])
    set.id <- factor(set.id)
    levs <- levels(set.id)
    iset <- as.integer(set.id)
    Nset <- length(levs)
    gp <- !is.null(per$gp.R)
    R <- per$gp.R
    if(!gp && (any(as.integer(per$Nma)>0) || any(as.integer(per$Nar)>0))){
        warning('The per-set AR/MA lag terms are not refined by the multi-set MCMC; the per-set jitters absorb the correlated noise.')
    }
    span <- max(t)-min(t)
    sdy <- sd(y)
    popt <- as.numeric(per$Popt[1])
    gnames <- paste0('gamma_',make.names(levs))
    snames <- paste0('sj_',make.names(levs))

###start from the deterministic fit
    if(SigType=='kepler'){
        kf <- KeplerFit.multiset(per)
        K0 <- max(as.numeric(kf$ParKep$K1),1e-3*sdy)
        start <- c(logP1=log(as.numeric(kf$ParKep$P1)),K1=K0,
                   e1=min(max(as.numeric(kf$ParKep$e1),1e-3),0.9),
                   omega1=as.numeric(kf$ParKep$omega1)%%(2*pi),Mo1=as.numeric(kf$ParKep$Mo1)%%(2*pi))
        low <- c(logP1=log(0.8*popt),K1=0,e1=0,omega1=0,Mo1=0)
        up <- c(logP1=log(1.2*popt),K1=max(5*K0,5*sdy),e1=0.95,omega1=2*pi,Mo1=2*pi)
        gam0 <- unlist(kf$ParKep[gnames])
        beta0 <- if(!is.null(kf$ParKep$beta)) as.numeric(kf$ParKep$beta) else 0
    }else{
        po <- unlist(per$par.opt)
        A0 <- if('A'%in%names(po)) as.numeric(po[['A']]) else 0
        B0 <- if('B'%in%names(po)) as.numeric(po[['B']]) else 0
        amp <- max(sqrt(A0^2+B0^2),1e-3*sdy)
        start <- c(logP1=log(popt),A1=A0,B1=B0)
        low <- c(logP1=log(0.8*popt),A1=-5*max(amp,sdy),B1=-5*max(amp,sdy))
        up <- c(logP1=log(1.2*popt),A1=5*max(amp,sdy),B1=5*max(amp,sdy))
        gam0 <- sapply(levs,function(l) mean(y[set.id==l]))
        if('gamma_1'%in%names(po) || any(gnames%in%names(po))) gam0 <- ifelse(gnames%in%names(po),po[gnames],gam0)
        beta0 <- if('beta'%in%names(po)) as.numeric(po[['beta']]) else 0
    }
    gam0 <- as.numeric(gam0); gam0[!is.finite(gam0)] <- mean(y)
    gwid <- diff(range(y))+10*sdy
    start <- c(start,setNames(gam0,gnames),beta=beta0)
    low <- c(low,setNames(gam0-gwid,gnames),beta=beta0-10*sdy/max(span,1))
    up <- c(up,setNames(gam0+gwid,gnames),beta=beta0+10*sdy/max(span,1))
    if(!gp){
        sj0 <- sapply(levs,function(l) median(dy[set.id==l]))
        start <- c(start,setNames(pmax(sj0,1e-3),snames))
        low <- c(low,setNames(rep(0,Nset),snames))
        up <- c(up,setNames(rep(5*sdy,Nset),snames))
    }
    start <- pmin(pmax(start,low+1e-9),up-1e-9)

###the shared-signal model and its likelihood
    signal.at <- function(par,tt){
        if(SigType=='kepler'){
            multiset_kepler_curve(list(P1=exp(par[['logP1']]),K1=par[['K1']],e1=par[['e1']],
                                       omega1=par[['omega1']],Mo1=par[['Mo1']]),tt)
        }else{
            par[['A1']]*cos(2*pi*tt/exp(par[['logP1']]))+par[['B1']]*sin(2*pi*tt/exp(par[['logP1']]))
        }
    }
    model <- function(par) signal.at(par,t)+par[gnames][iset]+par[['beta']]*t
    loglike <- function(par){
        r <- y-model(par)
        if(gp){
            z <- base::backsolve(R,r,transpose=TRUE)
            -0.5*sum(z^2)-sum(log(diag(R)))-length(y)/2*log(2*pi)
        }else{
            s2 <- dy^2+(par[snames][iset])^2
            -0.5*sum(r^2/s2+log(2*pi*s2))
        }
    }

###a local sampling environment so the generic PT machinery sees this model
    par.min <- low
    par.max <- up
    Npar <- length(start)
    Sd <- 2.4^2/Npar
    tol1 <- 1e-16
###flat priors within the box except for the eccentricity (see eprior.log)
    e.prior <- eprior.check(e.prior)
    if(mcmc.verbose && SigType=='kepler') cat('eccentricity prior:',eprior.describe(e.prior),'\n')
    posterior <- function(param,tem=1,bases='natural'){
        ll <- loglike(param)
        if(!is.finite(ll)) ll <- -1e10
        lp <- if(SigType=='kepler') eprior.log(param[['e1']],e.prior) else 0
        if(!is.finite(lp)) lp <- -1e10
        list(loglike=ll,logprior=lp,post=ll*tem+lp)
    }
    env <- environment()
    rebind <- function(f){ environment(f) <- env; f }
    ptmcmc.proposal <- rebind(ptmcmc.proposal)
    ptmcmc.cov <- rebind(ptmcmc.cov)
    ptmcmc.temmin <- rebind(ptmcmc.temmin)
    covariance.n0 <- rebind(covariance.n0)
    covariance.rep <- rebind(covariance.rep)
    run.pt <- rebind(run.ptmcmc)
    cov.start <- diag((1e-3*(par.max-par.min))^2,nrow=Npar)

###a forked chain cannot update the interface, so the dialog itself is the
###status indicator; with a single core the cold chain reports its iterations
    prog <- if(Ncores==1) function(k,detail='') incProgress(k,detail=detail) else NULL
    withProgress(message=paste0('Running multi-set PT-MCMC (',Ncores,' chain',if(Ncores>1) 's' else '',' x ',
                                max(as.numeric(Niter),1000),' iterations)'), value=0.05, {
    chains <- foreach(ncore=1:Ncores,.errorhandling='pass') %dopar% {
        run.pt(start,cov.start,iterations=max(as.numeric(Niter),1000),bases='natural',
               Ntem=Ntem,tem.min=tem.min,swap.interval=swap.interval,verbose=mcmc.verbose,progress=prog)
    }
    })
    bad <- sapply(chains,function(ch) !(is.list(ch) && !is.null(ch$out)))
    if(all(bad)){
        msg <- tryCatch(conditionMessage(chains[[1]]),error=function(e) 'unknown error')
        stop('all multi-set MCMC chains failed: ',msg)
    }
    chains <- chains[!bad]
    mc <- do.call(rbind,lapply(chains,function(ch) ch$out))
    if(mcmc.verbose){
        cat('multi-set PTMCMC: ',length(chains),'chains;',nrow(mc),'samples; max Rhat:',
            max(unlist(lapply(chains,function(ch) if(is.null(ch$Rhat)) NA else max(ch$Rhat))),na.rm=TRUE),'\n')
    }
###report the period linearly
    mc[,'logP1'] <- exp(mc[,'logP1'])
    colnames(mc)[colnames(mc)=='logP1'] <- 'P1'
    ind.max <- which.max(mc[,'loglike'])
    ParSig <- mc[ind.max,1:Npar]
    pb <- ParSig
    pb[['P1']] <- log(pb[['P1']])
    names(pb)[names(pb)=='P1'] <- 'logP1'
    par.stat <- sapply(1:Npar,function(i) msmc.stat(mc[,i],mc[,'loglike']))
    colnames(par.stat) <- colnames(mc)[1:Npar]
    popt <- as.numeric(ParSig[['P1']])
    ysig0 <- signal.at(pb,t)
    res <- y-model(pb)
    if(gp){
###remove the GP prediction so the residual is white for the plots
        gpm <- gp_predict(t,R,per$gp$sigmaGP,per$gp$logProt,per$gp$logtauGP,res)
        ysig0.gp <- gpm
        res <- res-gpm
    }
    ysim0 <- signal.at(pb,tsim)
    list(mc=mc,ParSig=ParSig,par.stat=par.stat,popt=popt,ysig0=as.numeric(ysig0),
         res=as.numeric(res),ysim0=as.numeric(ysim0))
}

sigfit.multiset <- function(per, data, t, tsim, SigType='circular', mcf=FALSE, Niter=1e3, Ncores=4,
                            Ntem=NULL, tem.min=NULL, swap.interval=10, mcmc.verbose=FALSE, e.prior=NULL){
###Turn a multi-data-set periodogram into a fitted model, its residual and its
###phase-folded prediction. SigType selects a shared circular signal, a shared
###Keplerian signal or a purely stochastic (signal-free) red-noise model.
    if(mcf && SigType=='stochastic'){
        warning('MCMC refinement is not available for the purely stochastic multi-set model; using the maximum-likelihood fit.')
        mcf <- FALSE
    }
    if(SigType=='stochastic' & !isTRUE(per$noise_only)){
        warning('The multi-set periodogram was not computed with a purely stochastic model; using the circular signal instead.')
        SigType <- 'circular'
    }
    if(SigType!='stochastic' & isTRUE(per$noise_only)){
        warning('The multi-set periodogram is purely stochastic; reporting the stochastic fit.')
        SigType <- 'stochastic'
    }
    par.stat <- NULL
    mc <- c()
    if(mcf){
###PT-MCMC over the shared signal, per-set offsets, trend and per-set jitters
        tmp <- mcfit.multiset(per=per,data=data,set.id=per$set.id,tsim=tsim,Niter=Niter,e.prior=e.prior,
                              SigType=SigType,Ncores=Ncores,Ntem=Ntem,tem.min=tem.min,
                              swap.interval=swap.interval,mcmc.verbose=mcmc.verbose)
        popt <- tmp$popt
        ysig0 <- tmp$ysig0
        res <- tmp$res
        ysim.sig <- tmp$ysim0
        ParSig <- tmp$ParSig
        par.stat <- tmp$par.stat
        mc <- tmp$mc
        per$par.opt <- as.list(ParSig)
        per$Popt <- popt
        per$res <- res
        ysig <- ysig0+res
        per$res.s <- res
        if(is.na(popt) | popt<=0) popt <- 1e7
        tsim1 <- tsim%%popt
        inds <- sort(tsim1,index.return=TRUE)$ix
        return(list(per=per,t=t%%popt,y=as.numeric(ysig),ey=data[,3],res=as.numeric(res),
                    ysig0=as.numeric(ysig0),tsim0=tsim,ysim0=ysim.sig,
                    tsim=tsim1[inds],ysim=ysim.sig[inds],ParSig=ParSig,
                    par.stat=par.stat,popt=popt,mc=mc))
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
                   mcmc.method='PT',Ntem=NULL,tem.min=NULL,swap.interval=10,mcmc.verbose=FALSE,e.prior=NULL){
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
        return(sigfit.multiset(per=per,data=data,t=t,tsim=tsim,SigType=SigType,mcf=mcf,Niter=Niter,Ncores=Ncores,e.prior=e.prior,
                               Ntem=Ntem,tem.min=tem.min,swap.interval=swap.interval,mcmc.verbose=mcmc.verbose))
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
        tmp <- mcfit(per=per,data=data,tsim=tsim,Niter=Niter,SigType=SigType,basis=basis,Pconv=Pconv,Ncores=Ncores,e.prior=e.prior,
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
    return(list(per=per,t=ts,y=as.numeric(ysig),ey=data[,3],res=as.numeric(res),ysig0=as.numeric(ysig0),tsim0=tsim0,ysim0=ysim0,tsim=tsims,ysim=ysims,ParSig=ParSig,par.stat=par.stat,popt=popt,mc=mc,
                llmax=if(mcf && !is.null(tmp$llmax)) as.numeric(tmp$llmax) else NA))
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

###A panel grid that, unlike par(mfrow), lets a panel span the full width of
###its row (tables). Panels are placed with par(fig); a full grid starts a new
###page, as mfrow does.
.fig.grid <- new.env()
grid.start <- function(nrow,ncol,cex=1,heights=NULL){
###heights: relative heights of the rows (equal by default), e.g. short rows
###for section headers
    pub.par(cex=cex)
    if(is.null(heights) || length(heights)!=nrow) heights <- rep(1,nrow)
    .fig.grid$nr <- nrow; .fig.grid$nc <- ncol; .fig.grid$k <- 0; .fig.grid$active <- TRUE
    .fig.grid$ytop <- 1-c(0,cumsum(heights))/sum(heights)
    invisible()
}
grid.stop <- function(){ .fig.grid$active <- FALSE; invisible() }
grid.next <- function(span=1){
###move to the next free cell (to the start of the next row for span>1) and
###make it the current figure; no-op when no grid is active
    g <- .fig.grid
    if(!isTRUE(g$active)) return(invisible(FALSE))
    nr <- g$nr; nc <- g$nc; k <- g$k
    span <- min(span,nc)
    if(span>1 && k%%nc!=0) k <- k+(nc-k%%nc)
    if(k>=nr*nc) k <- 0
    row <- k%/%nc+1; col <- k%%nc+1
    yt <- g$ytop
    if(is.null(yt) || length(yt)!=nr+1) yt <- 1-(0:nr)/nr
    fig <- c((col-1)/nc,(col-1+span)/nc,yt[row+1],yt[row])
    if(k==0){ par(fig=c(0,1,0,1),new=FALSE); plot.new() }
    par(fig=fig,new=TRUE)
    g$k <- k+span
    invisible(TRUE)
}

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

###data sets are told apart by a fixed, colour-blind-safe (Okabe-Ito based)
###hue order and, as a second cue, by the point symbol; the vermillion of the
###model curve is deliberately left out
set.palette <- c('#0072B2','#E69F00','#009E73','#CC79A7','#56B4E9','#5D3A9B','#A6761D','#1B9E77')
set.pch <- c(21,22,24,23,25,21,22,24)

set.label <- function(sets){
###display names of the data sets: file names are 'star_instrument', so the
###star is dropped and the instrument kept; the full names stay when that
###would make two sets indistinguishable
    sets <- as.character(sets)
    if(length(sets)==0) return(sets)
    short <- ifelse(grepl('_',sets),sub('^[^_]+_','',sets),sets)
    short[!nzchar(short)] <- sets[!nzchar(short)]
    if(any(duplicated(short))) sets else short
}

set.info <- function(ph){
###which data set each point of a phase.list entry belongs to, with the set
###names, colours and symbols; NULL for a single data set so that single-set
###figures keep their plain look
    if(is.null(ph) || !('set0'%in%colnames(ph))) return(NULL)
    id <- as.integer(ph[,'set0'])
    n <- suppressWarnings(max(id,na.rm=TRUE))
    if(!is.finite(n) || n<2) return(NULL)
    sets <- attr(ph,'sets')
    if(is.null(sets) || length(sets)<n) sets <- paste('Set',1:n)
    k <- (1:n-1)%%length(set.palette)+1
    list(id=id,names=set.label(sets[1:n]),col=set.palette[k],pch=set.pch[k])
}

points.sets <- function(t,y,ey,si,col1=pub.col$data){
###error bars and points; coloured by data set when there are several
    if(is.null(si)){
        try(arrows(t,y-ey,t,y+ey,length=0.02,angle=90,code=3,col=pub.col$err),TRUE)
        points(t,y,pch=21,bg='white',col=col1,cex=0.9)
    }else{
        try(arrows(t,y-ey,t,y+ey,length=0.02,angle=90,code=3,col=adjustcolor(si$col[si$id],alpha.f=0.5)),TRUE)
        points(t,y,pch=si$pch[si$id],bg=si$col[si$id],col='black',cex=0.9,lwd=0.6)
    }
}

legend.sets <- function(si,extra=NULL,pos='topleft'){
###legend naming the data sets, followed by optional text lines (e.g. the RMS)
    n <- length(si$names); m <- length(extra)
    if(n+m==0) return(invisible())
    legend(pos,bg='white',box.col='white',legend=c(si$names,extra),
           pch=c(si$pch,rep(NA,m)),pt.bg=c(si$col,rep(NA,m)),col=c(rep('black',n),rep(NA,m)),
           pt.cex=0.9,text.col=c(rep('black',n),rep(pub.col$peak,m)),cex=if(n>0) 0.9 else 1)
}

ylim.legend <- function(ylim,nleg){
###head room above the data for a legend of nleg lines
    if(nleg<1) return(ylim)
    ylim[2] <- ylim[2]+0.07*(nleg+1)*diff(ylim)
    ylim
}

panel.phase <- function(phase.list,sim.list,ypar,i,title,pub=TRUE){
    ph <- phase.list[[ypar]]; sm <- sim.list[[ypar]]
    t <- ph[,paste0('t_sig',i)]; y <- ph[,paste0('y_sig',i)]; ey <- ph[,'ey0']
    tsim <- sm[,paste0('tsim_sig',i)]; ysim <- sm[,paste0('ysim_sig',i)]
    si <- set.info(ph)
    ylim <- ylim.legend(range(c(y-ey,y+ey,ysim),na.rm=TRUE),length(si$names))
    plot(t,y,type='n',xaxt='n',yaxt='n',xlab='Phase [day]',ylab=ylab.obs(ypar),ylim=ylim,
         main=if(pub) pretty.title(title) else title)
    pub.axes()
    if(is.null(si)) arrows(t,y-ey,t,y+ey,length=0.02,angle=90,code=3,col=pub.col$err)
    lines(tsim,ysim,col=pub.col$model,lwd=2.5)
    if(is.null(si)){
        points(t,y,pch=21,bg='white',col=pub.col$data,cex=0.9)
    }else{
        points.sets(t,y,ey,si)
        legend.sets(si,pos='topright')
    }
}

panel.fit <- function(phase.list,sim.list,ypar,title,pub=TRUE){
    ph <- phase.list[[ypar]]; sm <- sim.list[[ypar]]
    t <- ph[,'t0']; ey <- ph[,'ey0']; y <- ph[,'y_all']
    tsim <- sm[,'tsim0']+min(t); ysim <- sm[,'ysim_all']
    si <- set.info(ph)
    ylim <- ylim.legend(range(c(y-ey,y+ey,ysim),na.rm=TRUE),length(si$names)+1)
    plot(t,y,type='n',xaxt='n',yaxt='n',xlab='Time [day]',ylab=ylab.obs(ypar),ylim=ylim,
         main=if(pub) gsub(', signal \\d+',', combined fit',pretty.title(title)) else gsub('\\d signal','combined fit',title))
    pub.axes()
    lines(tsim,ysim,col=pub.col$model,lwd=2)
    points.sets(t,y,ey,si)
    legend.sets(si,extra=paste0('RMS = ',format(sd(y),digits=3)))
}

panel.residual <- function(phase.list,ypar,title,pub=TRUE){
    ph <- phase.list[[ypar]]
    t <- ph[,'t0']; ey <- ph[,'ey0']; res <- ph[,'res_all']
    si <- set.info(ph)
    ylim <- ylim.legend(range(c(res-ey,res+ey),na.rm=TRUE),length(si$names)+1)
    plot(t,res,type='n',xaxt='n',yaxt='n',xlab='Time [day]',ylab=paste('Residual',if(ypar=='RV') '[m/s]' else ''),ylim=ylim,
         main=if(pub) gsub(', signal \\d+',', residual',pretty.title(title)) else gsub('\\d signal','residual',title))
    pub.axes()
    abline(h=0,lty=3,col=pub.col$level)
    points.sets(t,res,ey,si,col1=pub.col$res)
    legend.sets(si,extra=paste0('RMS = ',format(sd(res),digits=3)))
}

par.table <- function(pl){
###the fitted parameters of one observable as a table. With an MCMC run: the
###MAP value (sample of maximum likelihood), mean, median and the 16 and 84 per
###cent (1 sigma) quantiles; otherwise the maximum-likelihood values alone
    if(is.null(pl) || length(pl)==0) return(NULL)
    if(is.matrix(pl)){
        rn <- rownames(pl)
        get <- function(r) if(r%in%rn) as.numeric(pl[r,]) else rep(NA,ncol(pl))
        map <- if('xopt'%in%rn && any(is.finite(pl['xopt',]))) get('xopt') else get('mode')
        tab <- data.frame(Parameter=colnames(pl),MAP=map,Mean=get('mean'),Median=get('med'),
                          q16=get('xminus.1sig'),q84=get('xplus.1sig'),stringsAsFactors=FALSE)
        if(!('xopt'%in%rn && any(is.finite(pl['xopt',])))) colnames(tab)[2] <- 'Mode'
    }else{
        tab <- data.frame(Parameter=names(pl),MAP=as.numeric(pl),stringsAsFactors=FALSE)
        colnames(tab)[2] <- 'Max. likelihood'
    }
    rownames(tab) <- NULL
###per-set parameters carry the instrument name only (see set.label)
    sets <- attr(pl,'sets')
    if(!is.null(sets) && length(sets)>1){
        full <- make.names(as.character(sets)); short <- make.names(set.label(sets))
        for(k in order(nchar(full),decreasing=TRUE)){
            tab$Parameter <- sub(paste0('_',full[k],'$'),paste0('_',short[k]),tab$Parameter)
        }
    }
    tab
}

draw.table <- function(cells,hdr,main,note=NULL,adj=NULL){
###a text table on the current device: cells is a character matrix, hdr the
###column headers; the first column is left-aligned, the others right-aligned
###unless adj says otherwise. Column widths come from the widest entry; the
###text shrinks when the table is wider than the panel and the columns spread
###when it is narrower
    op <- par(mar=c(1,1,2.6,1),xpd=NA); on.exit(par(op))
    plot.new(); plot.window(xlim=c(0,1),ylim=c(0,1),xaxs='i',yaxs='i')
    title(main=main)
    if(is.null(cells) || nrow(cells)==0){
        text(0.5,0.5,'nothing to report',col=pub.col$level)
        return(invisible(NULL))
    }
    nr <- nrow(cells); nc <- ncol(cells)
    if(is.null(adj)) adj <- c(0,rep(1,nc-1))
    cex <- min(1,20/(nr+3))
    w <- sapply(1:nc,function(j) max(strwidth(c(hdr[j],cells[,j]),cex=cex,font=2)))
    gap <- strwidth('MM',cex=cex)
    tot <- sum(w)+gap*(nc-1)
    if(tot>0.96){
        f <- 0.96/tot
        cex <- cex*f; w <- w*f; gap <- gap*f; tot <- 0.96
    }else if(nc>1){
###spread the columns, but not so far that a narrow table looks scattered
        gap <- gap+min((0.96-tot)/(nc-1),0.22)
    }
    xs <- 0.02+w[1]
    if(nc>1) for(j in 2:nc) xs <- c(xs,xs[j-1]+gap+w[j])
    xs[1] <- 0.02
###columns of adj 0 anchor at their left edge, of adj 1 at their right edge
    xl <- xs-w; xl[1] <- 0.02
    xpos <- ifelse(adj==0,xl,xs)
    lh <- min(0.94/(nr+2.6),2.2*strheight('X',cex=cex))
    y0 <- 0.97-0.5*lh
    for(j in 1:nc){
        text(xpos[j],y0,hdr[j],adj=c(adj[j],0.5),font=2,cex=cex)
        text(rep(xpos[j],nr),y0-(1:nr)*lh,cells[,j],adj=c(adj[j],0.5),cex=cex)
    }
    segments(0,y0-0.55*lh,1,y0-0.55*lh,col=pub.col$level)
    segments(0,y0-(nr+0.55)*lh,1,y0-(nr+0.55)*lh,col=pub.col$level)
    if(!is.null(note)){
        for(k in seq_along(note)){
            cn <- min(0.75*cex,0.96/strwidth(note[k],cex=1))
            text(0.02,y0-(nr+0.7+0.8*k)*lh,note[k],adj=c(0,0.5),cex=cn,col=pub.col$level)
        }
    }
    invisible(NULL)
}

panel.partable <- function(par.list,ypar,title,pub=TRUE){
###a panel listing the fitted parameters (see par.table) so that the table can
###be saved alongside the figures
    tab <- par.table(par.list[[ypar]])
    main <- if(pub) gsub(', signal \\d+',', fitted parameters',pretty.title(title)) else gsub('\\d signal','fitted parameters',title)
    if(is.null(tab)){
        draw.table(NULL,NULL,main)
        return(invisible(NULL))
    }
    hdr <- colnames(tab)
    hdr[hdr=='q16'] <- 'q16%'
    hdr[hdr=='q84'] <- 'q84%'
    fmt <- function(x) ifelse(is.finite(x),formatC(x,digits=4,format='g'),'')
    cells <- cbind(tab$Parameter,sapply(tab[,-1,drop=FALSE],fmt))
    cells <- matrix(as.character(cells),nrow=nrow(tab))
    note <- if(ncol(tab)>2) 'MAP: sample of maximum likelihood; q16%, q84%: 16 and 84 per cent quantiles (1 sigma interval)' else NULL
    draw.table(cells,hdr,main,note=note)
    invisible(tab)
}

panel.modelcomp <- function(model.list,ypar,title,pub=TRUE){
###the sequential model comparison of one observable: for each signal the
###period, the ln(BF) against the model without it, where that number comes
###from, and whether it was accepted; the most plausible number of signals
###follows from the sequence of accepted signals
    mc <- model.list[[ypar]]
    main <- if(pub) gsub(', signal \\d+',', model comparison',pretty.title(title)) else gsub('\\d signal','model comparison',title)
    if(is.null(mc) || nrow(mc)==0){
        draw.table(NULL,NULL,main)
        return(invisible(NULL))
    }
    thr <- attr(mc,'lnBF.min'); if(is.null(thr)) thr <- 5
    nmax <- attr(mc,'Nsig.max')
    nopt <- model.Nopt(mc)
    cells <- cbind(paste('Signal',mc$n),ifelse(is.finite(mc$P),formatC(mc$P,digits=5,format='g'),''),
                   ifelse(is.finite(mc$lnBF),formatC(mc$lnBF,digits=3,format='f'),'n/a'),
                   ifelse(mc$source=='MCMC','joint MCMC (BIC)','periodogram (BIC)'),
                   ifelse(mc$accepted,'yes','no'))
    hdr <- c('Model','Period [d]','ln(BF) vs previous','ln(BF) from','Accepted')
    note <- c(paste0('Most plausible number of signals: ',nopt,
                     if(!is.null(nmax) && nopt>=nmax) paste0(' (the maximum searched, ',nmax,')') else '',
                     '; a signal is accepted when ln(BF) > ',format(thr,digits=3),
                     if(any(is.na(mc$lnBF))) '; n/a: no Bayes factor for this periodogram type (accepted by default)' else ''),
              'ln(BF) = difference of maximum log likelihoods minus (extra parameters/2) ln N, i.e. the BIC estimate.')
    draw.table(cells,hdr,main,note=note,adj=c(0,1,1,0,0))
    invisible(mc)
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
            if(ypar!='Window Function' && !is.null(d$par.list[[ypar]])){
                out <- rbind(out,data.frame(label=paste0('Parameter table: ',ypar),kind='partable',ypar=ypar,index=1,stringsAsFactors=FALSE))
            }
            if(!is.null(d$model.comp[[ypar]])){
                out <- rbind(out,data.frame(label=paste0('Model comparison: ',ypar),kind='modelcomp',ypar=ypar,index=1,stringsAsFactors=FALSE))
            }
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
    }else if(kind=='partable'){
        panel.partable(d$par.list,ypar,titles[1],pub=pub)
    }else if(kind=='modelcomp'){
        panel.modelcomp(d$model.comp,ypar,titles[1],pub=pub)
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

plot1D.grid <- function(download=FALSE){
###the panel grid of the bundled figures: two columns; the on-screen figure is
###one tall page, the download is paginated 2x2
    if(download){
        grid.start(2,2,cex=0.9)
    }else{
        grid.start(ceiling(Nmax.plots/2),2)
    }
}

phase1D.plot <- function(phase.list,sim.list,tits,download=FALSE,index=NULL,repar=TRUE,pub=TRUE,par.list=NULL,model.comp=NULL){
###phase-folded signals, combined fit and residual of every observable; with
###par.list the fitted parameters follow as a table spanning the full width,
###with model.comp the sequential model comparison likewise
    if(repar && is.null(index)) plot1D.grid(download)
    for(ypar in names(phase.list)){
        if(!is.null(index)){
            inds <- index
        }else{
            inds <- 1:length(grep('^y_sig',colnames(phase.list[[ypar]])))
        }
        titles <- tits[grepl(paste0(';',ypar,';'),tits)]
        for(i in inds){
            grid.next()
            panel.phase(phase.list,sim.list,ypar,i,titles[i],pub=pub)
        }
        grid.next()
        panel.fit(phase.list,sim.list,ypar,titles[max(inds)],pub=pub)
        grid.next()
        panel.residual(phase.list,ypar,titles[max(inds)],pub=pub)
        if(ypar!='Window Function' && !is.null(par.list[[ypar]])){
            grid.next(2)
            panel.partable(par.list,ypar,titles[max(inds)],pub=pub)
        }
        if(!is.null(model.comp[[ypar]])){
            grid.next(2)
            panel.modelcomp(model.comp,ypar,titles[max(inds)],pub=pub)
        }
    }
}

combined.plot <- function(per.list,phase.list,sim.list,tits,pers,levels,ylabs,SigType='circular',download=FALSE,index=NULL,pub=TRUE,par.list=NULL,model.comp=NULL){
    per1D.plot(per.list,tits,pers,levels,ylabs,download=download,index=index,SigType=SigType,pub=pub,par.list=par.list)
    phase1D.plot(phase.list,sim.list,tits=tits,download=download,index=index,repar=FALSE,pub=pub,par.list=par.list,model.comp=model.comp)
}

per1D.plot <- function(per.list,tits,pers,levels,ylabs,download=FALSE,index=NULL,SigType='circular',pub=TRUE,par.list=NULL,repar=TRUE){
    if(repar && is.null(index)) plot1D.grid(download)
    gi <- 0
    for(ypar in names(per.list)){
        np <- ncol(per.list[[ypar]])-1
        titles <- tits[grepl(paste0(';',ypar,';'),tits)]
        inds <- if(!is.null(index)) index else 1:np
        for(i in inds){
###significance levels are stored per periodogram across all observables
            lv <- if(is.matrix(levels) && ncol(levels)>=gi+i) levels[,gi+i] else NULL
            Pmark <- reported.period(par.list,ypar,i)
            grid.next()
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

per2D.data <- function(vars,per.par,data,fit1D=NULL){
###fit1D: a calc.1Dper() result; with per.par$use.fit TRUE the moving
###periodogram is computed on its signal-only series (the data minus the
###per-set offsets, the trend and the red-noise model of that fit, i.e. the
###MCMC solution when MCMC was run), which combines the data sets on a common
###zero point and so extends the baseline for testing the time consistency of
###long-period signals
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
###red-noise model of the windows: ARMA (default) or a shared SHO Gaussian process
    if(!exists('noise.model')) noise.model <- 'ARMA'
    GP <- noise.model=='GP'
    gp.par <- rep(NA,3)
    if(exists('gp.Prot') && length(gp.Prot)==1 && is.finite(gp.Prot) && gp.Prot>0) gp.par[2] <- log(gp.Prot)
    if(exists('gp.tau') && length(gp.tau)==1 && is.finite(gp.tau) && gp.tau>0) gp.par[3] <- log(gp.tau)
    if(GP){
        Nma <- 0
        Nar <- 0
        Nmas <- rep(0,length(Nmas))
        Nars <- rep(0,length(Nars))
    }
###the signal-only series of the 1D fit replaces the data and the noise model
    use.fit <- exists('use.fit') && isTRUE(use.fit) && !is.null(fit1D)
    if(use.fit){
        var <- pars[[i]]$var
        per.type <- pars[[i]]$per.type
        ph <- fit1D$phase.list[[var]]
        if(is.null(ph) || !('y_all'%in%colnames(ph))){
            stop('the 1D fit has no combined model for ',var,'; compute the 1D periodogram for this observable first')
        }
        sets <- attr(ph,'sets')
        if(is.null(sets)) sets <- per.target
        if(!setequal(sets,per.target)){
            stop('the 1D fit was made for ',paste(sets,collapse=', '),' but the moving periodogram is asked for ',
                 paste(per.target,collapse=', '),'; select the same data sets or recompute the 1D periodogram')
        }
        ord <- order(as.numeric(ph[,'t0']))
        t <- as.numeric(ph[ord,'t0']); y <- as.numeric(ph[ord,'y_all']); dy <- as.numeric(ph[ord,'ey0'])
        set0 <- if('set0'%in%colnames(ph)) as.integer(ph[ord,'set0']) else rep(1L,length(t))
        adaptive <- exists('adaptive') && isTRUE(adaptive)
        mp <- MP(t=t,y=y,dy=dy,Dt=Dt,nbin=Nbin,ofac=ofac,fmin=frange[1],fmax=frange[2],per.type=per.type,sj=0,Nma=0,Nar=0,Indices=NULL,adaptive=adaptive)
        fname <- paste0(paste(per.target,collapse='_'),'_MP_',paste(per.type,collapse=''),'_signalonly',if(adaptive) '_adaptive' else '')
        out <- list(t=t,y=y,dy=dy,xx=mp$tmid,yy=mp$P,zz=mp$powers,zz.rel=mp$rel.powers,fname=fname,ypar=var,
                    set=per.target[match(sets[set0],per.target)],signal.only=TRUE,tstart=mp$tstart,tend=mp$tend,adaptive=adaptive)
        if(length(per.target)>1){
            out$subdata <- lapply(1:length(per.target),function(j) data[[per.target[j]]])
            out$idata <- lapply(per.target,function(s){ k <- set0==match(s,sets); cbind(t[k],y[k],dy[k]) })
        }
        return(out)
    }
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
        tmp <- combine.data(data=subdata,Ninds=Inds.sets,Nmas=Nmas,GP=GP,gp.par=gp.par)
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
    adaptive <- exists('adaptive') && isTRUE(adaptive)
    if(length(per.target)==1){
        mp <- MP(t=t,y=y,dy=dy,Dt=Dt,nbin=Nbin,ofac=ofac,fmin=frange[1],fmax=frange[2],per.type=per.type,sj=0,Nma=Nma,Nar=Nar,Indices=Indices,GP=GP,gp.par=gp.par,adaptive=adaptive)
    }else{
###the per-set noise (ARMA or GP) was removed in combine.data(); the combined residual is analysed white
        mp <- MP(t=t,y=y,dy=dy,Dt=Dt,nbin=Nbin,ofac=ofac,fmin=frange[1],fmax=frange[2],per.type=per.type,sj=0,Nma=0,Nar=0,Indices=Indices,adaptive=adaptive)
    }
    x2 <- mp$tmid
    y2 <- mp$P
    z2 <- mp$powers
    z2.rel <- mp$rel.powers
    fname <- paste0(paste(per.target,collapse='_'),'_MP_',paste(per.type,collapse=''),if(GP) '_GP' else paste0('_MA',paste(Nmas,collapse='')),'proxy',paste(Inds,collapse='.'),if(adaptive) '_adaptive' else '')
    if(length(per.target)==1){
        return(list(t=t,y=y,dy=dy,xx=x2,yy=y2,zz=z2,zz.rel=z2.rel,fname=fname,ypar=ypar,tstart=mp$tstart,tend=mp$tend,adaptive=adaptive))
    }else{
        return(list(t=t,y=y,dy=dy,xx=x2,yy=y2,zz=z2,zz.rel=z2.rel,subdata=subdata,idata=idata,fname=fname,ypar=ypar,tstart=mp$tstart,tend=mp$tend,adaptive=adaptive))
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

###########################################################################
####Overall diagnosis of the signals in several RV data sets, after Feng et
####al. (2020, ApJS 250, 29, Figs. 5-16): BFPs of the combined data with the
####signals subtracted in sequence for several noise models (black), BFPs of
####each data set (grey), BFPs of the activity indices and the window function
####(blue), and the moving periodogram of the signal-only combined series. A
####Keplerian signal should be significant, robust to the noise model, absent
####from the activity indices, and consistent over time.
diag.noise <- list(W=c(Nma=0,Nar=0),MA=c(Nma=1,Nar=0),AR=c(Nma=0,Nar=1),GP=c(Nma=0,Nar=0))
diag.col <- list(combined='black',individual='grey45',proxy='#0072B2',window='#56B4E9',signal='#D55E00',rotation='#009E73')

diagnose.signals <- function(data,sets,noise.models=c('W','MA','AR'),Nsig.max=3,ofac=1,frange=NULL,lnBF.min=5,
                             SigType='circular',Nh=1,Ncores=1,Nwin=5,individual=TRUE,proxies=TRUE,moving=TRUE,progress=NULL,
                             gp.Prot=NA,gp.tau=NA,adaptive=TRUE){
###data: the named list of data tables; sets: the RV data sets to diagnose;
###noise.models: any of 'W' (white), 'MA' (MA(1)), 'AR' (AR(1)) and 'GP' (the
###SHO Gaussian process, with gp.Prot and gp.tau optionally fixed); the first
###one also serves the individual sets, the activity indices are analysed
###white. Returns the calc.1Dper() results of every part plus the accepted periods
    sets <- sets[sets%in%names(data)]
    if(length(sets)==0) stop('no data set selected')
    noise.models <- noise.models[noise.models%in%names(diag.noise)]
    if(length(noise.models)==0) noise.models <- 'W'
    rv <- colnames(data[[sets[1]]])[2]
    tspan <- diff(range(unlist(lapply(data[sets],function(x) x[,1]))))
    if(is.null(frange)) frange <- c(1/tspan,1/2)
    step <- function(msg) if(is.function(progress)) progress(msg)
    base <- function(target,model,sequence,Nsig){
        n <- length(target)
        nm <- diag.noise[[model]]
        list(ns=c(rv,'Window Function'),ofac=ofac,frange=frange,per.type='BFP',per.target=target,sequence=sequence,
             Nmas=rep(as.integer(nm['Nma']),n),Nars=rep(as.integer(nm['Nar']),n),Inds=rep(list(0),n),Nsig.max=as.integer(Nsig),
             per.type.seq='BFP',Niter=0,SigType=SigType,Nh=Nh,lnBF.min=lnBF.min,
             noise.model=if(model=='GP') 'GP' else 'ARMA',gp.Prot=gp.Prot,gp.tau=gp.tau)
    }
    out <- list(sets=sets,rv=rv,noise.models=noise.models,lnBF.min=lnBF.min,Nsig.max=Nsig.max,tspan=tspan,frange=frange,
                combined=list(),individual=list(),proxies=list(),moving=NULL)
    for(m in noise.models){
        step(paste0('combined data, ',m,' noise'))
        pp <- base(sets,m,sequence=Nsig.max>1,Nsig=Nsig.max)
        out$combined[[m]] <- calc.1Dper(Nmax.plots=50,vars=c(rv,'Window Function'),per.par=pp,data=data,Ncores=Ncores)
    }
    m1 <- noise.models[1]
    if(individual && length(sets)>1){
        for(s in sets){
            step(paste0('data set ',s))
            pp <- base(s,m1,sequence=FALSE,Nsig=1)
            out$individual[[s]] <- calc.1Dper(Nmax.plots=50,vars=rv,per.par=pp,data=data,Ncores=Ncores)
        }
    }
    if(proxies){
        for(s in sets){
            idx <- colnames(data[[s]])[-(1:3)]
            if(length(idx)==0) next
            step(paste0('activity indices of ',s))
            pp <- base(s,'W',sequence=FALSE,Nsig=1)
            pp$ns <- c(rv,idx,'Window Function')
            out$proxies[[s]] <- calc.1Dper(Nmax.plots=50,vars=idx,per.par=pp,data=data,Ncores=Ncores)
        }
    }
    mc <- out$combined[[m1]]$model.comp[[rv]]
    out$model.comp <- mc
    out$periods <- if(is.null(mc)) numeric(0) else mc$P[cumprod(as.logical(mc$accepted))==1]
    if(moving && length(out$periods)>0){
###one moving periodogram per accepted signal, on the signal-only combined
###series with the other accepted signals subtracted (Feng et al. 2020, Fig. 5)
        pp <- base(sets,'W',sequence=FALSE,Nsig=1)
        pp$gp.Prot <- NA; pp$gp.tau <- NA
        pp <- c(pp,list(files=NULL,Dt=signif(tspan/2,3),Nbin=as.integer(Nwin),alpha=5,scale=TRUE,
                        pmin.zoom=1/frange[2],pmax.zoom=1/frange[1],show.signal=TRUE,use.fit=TRUE,adaptive=adaptive))
        out$moving.par <- pp
        out$moving <- list()
        fit <- out$combined[[m1]]
        ph <- fit$phase.list[[rv]]
        for(k in seq_along(out$periods)){
            step(paste0('moving periodogram of signal ',k))
            others <- setdiff(seq_along(out$periods),k)
            fitk <- fit
            if(length(others)>0){
                cols <- paste0('ysig_sig',others)
                cols <- cols[cols%in%colnames(ph)]
                if(length(cols)>0) fitk$phase.list[[rv]][,'y_all'] <- ph[,'y_all']-rowSums(ph[,cols,drop=FALSE])
            }
            out$moving[[k]] <- tryCatch(per2D.data(vars=rv,per.par=pp,data=data,fit1D=fitk),
                                        error=function(e){ warning('moving periodogram of signal ',k,' failed: ',conditionMessage(e)); NULL })
        }
    }
    out
}

diag.peak <- function(per,P0,tol=0.1){
###the highest power of a periodogram matrix (P, power) within tol*P0 of P0
    k <- which(abs(per[,1]-P0)<tol*P0)
    if(length(k)==0) return(NA)
    max(per[k,2],na.rm=TRUE)
}

diagnosis.table <- function(diag,tol=0.1){
###one row per accepted signal of the first noise model: ln(BF) under every
###noise model, robustness, overlap with activity indices and the window
###function, time consistency from the moving periodogram, and a verdict
    P <- diag$periods
    if(length(P)==0) return(NULL)
    thr <- diag$lnBF.min; rv <- diag$rv
    rows <- list()
    for(k in seq_along(P)){
        row <- list(Signal=k,`Period [d]`=signif(P[k],5))
###ln(BF) of the signal under each noise model (from the sequential search)
        robust <- TRUE
        for(m in diag$noise.models){
            mc <- diag$combined[[m]]$model.comp[[rv]]
            if(m==diag$noise.models[1]){
                j <- k
            }else{
###the accepted signal of this noise model closest to the period, within tol
                acc <- if(is.null(mc)) integer(0) else which(cumprod(as.logical(mc$accepted))==1 & abs(mc$P-P[k])<tol*P[k])
                j <- if(length(acc)>0) acc[which.min(abs(mc$P[acc]-P[k]))] else integer(0)
            }
            v <- if(length(j)>0 && j<=nrow(mc)) mc$lnBF[j] else NA
            row[[paste0('ln(BF) ',m)]] <- if(is.finite(v)) round(v,1) else NA
            if(!(length(j)>0 && (is.na(v) || v>thr))) robust <- FALSE
        }
        row$`Robust to noise model` <- if(robust) 'yes' else 'no'
###a period within tol of an earlier accepted signal is likely its alias or residual
        near <- which(seq_along(P)<k & abs(P-P[k])<tol*P[k])
###activity indices whose periodogram peaks within tol of the period
        act <- c()
        for(s in names(diag$proxies)){
            pl <- diag$proxies[[s]]$per.list
            for(idx in names(pl)){
                if(idx=='Window Function') next
                v <- diag.peak(pl[[idx]][,1:2],P[k],tol)
                if(is.finite(v) && v>thr) act <- c(act,paste0(set.label(s),': ',idx))
            }
        }
        row$`Activity overlap` <- if(length(act)>0) paste(act,collapse='; ') else 'none'
###window function: is the sampling pattern's strongest peak at this period?
        wf <- diag$combined[[1]]$per.list[['Window Function']]
        wfp <- if(is.null(wf)) NA else wf[which.max(wf[,2]),1]
        row$`Window function peak` <- if(is.finite(wfp) && abs(wfp-P[k])<tol*P[k]) 'yes' else 'no'
###time consistency: windows of the moving periodogram in which the signal exceeds the threshold
        tc <- NA
        mv <- if(is.list(diag$moving) && length(diag$moving)>=k) diag$moving[[k]] else NULL
        if(!is.null(mv) && is.matrix(mv$zz)){
            zz <- mv$zz
            if(nrow(zz)!=length(mv$yy) && ncol(zz)==length(mv$yy)) zz <- t(zz)
            sel <- abs(mv$yy-P[k])<tol*P[k]
###a malformed matrix (no window computed) leaves the test as n/a
            if(nrow(zz)==length(mv$yy) && any(sel)){
                pk <- apply(zz[sel,,drop=FALSE],2,function(z) suppressWarnings(max(z,na.rm=TRUE)))
                ok <- is.finite(pk)
                tc <- paste0(sum(pk[ok]>thr),'/',sum(ok),' windows')
                row$tc.frac <- if(sum(ok)>0) sum(pk[ok]>thr)/sum(ok) else NA
            }
        }
        row$`Consistent over time` <- if(is.na(tc)) 'n/a' else tc
###verdict
        v1 <- is.finite(row[[paste0('ln(BF) ',diag$noise.models[1])]]) && row[[paste0('ln(BF) ',diag$noise.models[1])]]>thr
        verdict <- if(!v1) 'not significant' else if(!robust) 'depends on the noise model' else if(length(act)>0) 'possible activity' else 'candidate'
        if(row$`Window function peak`=='yes') verdict <- paste0(verdict,'; at a window-function peak')
        if(length(near)>0) verdict <- paste0(verdict,'; within ',round(100*tol),'% of signal ',near[1],' (alias or residual?)')
        if(!is.null(row$tc.frac) && is.finite(row$tc.frac) && row$tc.frac<0.5 && P[k]<diag$tspan/2) verdict <- paste0(verdict,'; not consistent over time')
        row$Verdict <- verdict
        row$tc.frac <- NULL
        rows[[k]] <- as.data.frame(row,check.names=FALSE,stringsAsFactors=FALSE)
    }
    do.call(rbind,rows)
}

diag.sections <- c(combined='Combined RV data: signals subtracted in sequence, one row per noise model (black)',
                   individual='Individual RV data sets (grey)',
                   proxy='Activity indices of each data set (blue)',
                   window='Window function of the combined data: sampling aliases (light blue)')

diagnosis.panels <- function(diag){
###the list of periodogram panels of the diagnosis figure, grouped as in Feng
###et al. (2020): the combined data per noise model and signal, the individual
###data sets, the activity indices, the window function
    rv <- diag$rv
    panels <- list()
    add <- function(per,title,kind,ylab='ln BF'){ panels[[length(panels)+1]] <<- list(per=per,title=title,kind=kind,ylab=ylab) }
    for(m in diag$noise.models){
        pl <- diag$combined[[m]]$per.list[[rv]]
        if(is.null(pl)) next
        for(j in 2:ncol(pl)){
            k <- j-1
            add(pl[,c(1,j)],paste0(m,' noise, ',if(k==1) 'all data' else paste0(k-1,' signal',if(k>2) 's' else '',' subtracted')),'combined')
        }
    }
    for(s in names(diag$individual)){
        pl <- diag$individual[[s]]$per.list[[rv]]
        if(!is.null(pl)) add(pl[,1:2],paste0(set.label(s),' (',diag$noise.models[1],' noise)'),'individual')
    }
    for(s in names(diag$proxies)){
        pl <- diag$proxies[[s]]$per.list
        for(idx in names(pl)){
            if(idx=='Window Function') next
            add(pl[[idx]][,1:2],paste0(set.label(s),': ',idx),'proxy')
        }
    }
    wf <- diag$combined[[1]]$per.list[['Window Function']]
    if(!is.null(wf)) add(wf[,1:2],'window function (combined data)','window','Power')
    panels
}

diagnosis.layout <- function(diag,ncol=3){
###rows of the diagnosis figure: a short header row before each section and
###the panel rows of that section; returns the relative row heights and the
###sections present
    panels <- diagnosis.panels(diag)
    kinds <- sapply(panels,function(p) p$kind)
    secs <- names(diag.sections)[names(diag.sections)%in%kinds]
    heights <- c()
    for(sname in secs){
        ng <- sum(kinds==sname)
        heights <- c(heights,0.22,rep(1,ceiling(ng/ncol)))
    }
    list(heights=heights,sections=secs,n=length(panels))
}

panel.header <- function(label,note=NULL){
###a section header spanning the figure width
    op <- par(mar=c(0.1,0.5,0.1,0.5),xpd=NA); on.exit(par(op))
    plot.new(); plot.window(xlim=c(0,1),ylim=c(0,1))
    rect(0,0.05,1,0.95,col='grey93',border=NA)
    text(0.01,0.55,label,adj=c(0,0.5),font=2,cex=1.15)
    if(!is.null(note)) text(0.99,0.55,note,adj=c(1,0.5),cex=0.85,col=pub.col$level)
    invisible()
}

panel.diag <- function(p,periods=NULL,lnBF.min=5,Prot=NA,number=NULL,xlim=NULL){
###one panel of the diagnosis figure
    P <- p$per[,1]; pw <- p$per[,2]
    ok <- is.finite(P) & is.finite(pw)
    P <- P[ok]; pw <- pw[ok]
    col <- diag.col[[p$kind]]
    ylim <- range(c(pw,if(p$kind!='window') lnBF.min),na.rm=TRUE)
    ylim[2] <- ylim[2]+0.12*diff(ylim)
    if(is.null(xlim)) xlim <- range(P)
    plot(P,pw,type='n',log='x',xaxt='n',yaxt='n',xlab='Period [day]',ylab=p$ylab,xlim=xlim,ylim=ylim,
         main=paste0(if(!is.null(number)) paste0('(',number,') ') else '',p$title),cex.main=0.95)
    pub.axes(xlog=TRUE)
    if(length(periods)>0) abline(v=periods,col=diag.col$signal,lwd=1.2)
    if(length(Prot)==1 && is.finite(Prot) && Prot>0) abline(v=Prot,col=diag.col$rotation,lty=3,lwd=1.5)
    if(p$kind!='window') abline(h=lnBF.min,lty=2,col=pub.col$level)
    lines(P,pw,col=col,lwd=if(p$kind=='combined') 1.4 else 1.1)
    invisible()
}

diagnosis.plot <- function(diag,ncol=3,Prot=NA,pub=TRUE){
###the diagnosis figure: the panels of each section under a header row, on a
###grid of ncol columns; returns the number of periodogram panels drawn
    panels <- diagnosis.panels(diag)
    n <- length(panels)
    if(n==0) return(invisible(0))
    lay <- diagnosis.layout(diag,ncol)
    xlim <- 1/rev(diag$frange)
    grid.start(length(lay$heights),ncol,cex=if(ncol>2) 0.85 else 1,heights=lay$heights)
    marks <- paste0('red: accepted periods',if(length(Prot)==1 && is.finite(Prot) && Prot>0) '; green dotted: rotation period' else '',
                    '; dashed: ln(BF) = ',format(diag$lnBF.min,digits=3))
    i <- 0
    for(sname in lay$sections){
        grid.next(ncol)
        panel.header(diag.sections[[sname]],note=if(sname==lay$sections[1]) marks else NULL)
        for(p in panels[sapply(panels,function(p) p$kind)==sname]){
            i <- i+1
            grid.next()
            panel.diag(p,periods=diag$periods,lnBF.min=diag$lnBF.min,Prot=Prot,number=i,xlim=xlim)
        }
###a partly filled panel row: skip to the next row for the next header
        k <- .fig.grid$k
        if(k%%ncol!=0) .fig.grid$k <- k+(ncol-k%%ncol)
    }
    invisible(n)
}
