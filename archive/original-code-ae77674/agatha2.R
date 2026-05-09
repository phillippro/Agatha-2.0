library(magicaxis)
library(foreach)
library(doMC)
library(parallel)
source('periodograms.R')
source('periodoframe.R')
source('mcmc_func.R')
source('sofa.R')
source('orbit.R')
options(warn=2)
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    PerType <- args[1]
    SigType <- args[2]
    Nmax <- as.integer(args[3])
    ofac <- as.numeric(args[4])
    noise <- args[5]
    dir.in <- args[6]
    fs <- args[7:length(args)]
}else{
    PerType <- 'BFP'
    SigType <- 'kepler'
    Nmax <- 1
    ofac <- 1
    noise <- 'MA'
#    dir.in <- '../data/combined/HD39060/'
    dir.in <- '~/Desktop/'
#    dir.in <- '../dwarfs/output/betaPic/'
    fs <- 'HD217987_rv_red.res'
#    fs <- 'HD39060_LagBin.dat'
#    fs <- 'betaPic_Nw1_N1e+05_logLmax-229.9_Np2_dec_res.dat'
}
#example: Rscript agatha2.R BFP kepler 2 0.1 MA data HD210193_PFS.vels HD103949_PFS.vels
if(!grepl('\\/$',dir.in)) dir.in <- paste0(dir.in,'/')
dir <- 'results/'
save.data <- TRUE
targets <- gsub('_.+','',fs)
stars <- gsub('\\..+','',fs)
fs <- paste0(dir.in,fs)
Ncores <- 1
Nsamp <- 1
Pmin <- 1.1
Pmax <- 1e5
fmin <- 1/Pmax
fmax <- 1/Pmin
lnBF <- pars <- Popt <- list()
if(Ncores>0) {registerDoMC(Ncores)} else {registerDoMC()}
labels <-c(paste0(PerType,' for raw RV'),rep(paste0(PerType,' for raw RV subtracted by signal'),Nmax-1),paste0(PerType,' for Sindex'),paste0(PerType,' for Halpha'),paste0(PerType,' for Photon Count'),paste0(PerType,' for Observation Time'),'Lomb-Scargle window function',paste0(PerType,' for raw RV subtracted by signal'))
for(jj in 1:length(targets)){
#for(jj in 1){
#out <- foreach(jj = 1:length(targets)) %dopar% {
    f <- fs[jj]
    if(FALSE){
        target <- gsub('HD','',targets[jj])
        target <- tolower(target)
        target <- gsub('a$','A',target)
        target <- gsub('b$','B',target)
        target <- gsub('c$','C',target)
    }else{
        target <- targets[jj]
    }
    cat('\ntarget:',target,'\n')
    syntax <- paste0(target,'_',PerType,'_',noise)
    tab <- read.table(f,header=TRUE)
    tmin <- min(tab[,1])
    if(tmin<2400000) tmin <- tmin+2400000
    t <- tab[,1]-tmin
    lnBF1 <- Popt1 <- par1 <- out1 <- list()
    Nper <- Nmax+ncol(tab)-2
    ysims <- 0
    gamma <- 0
    gammadot <- 0
    Popts <- c()
    tsim <- seq(min(t),max(t),length.out=1e3)
    for(k in 1:Nper){
        tit <- labels[k]
        if(k==1){
            y <- tab[,2]
            ey <- tab[,3]
            n <- paste0('sig',k)
        }else if(k<=Nmax){
            y <- per$res
            ey <- tab[,3]
            n <- paste0('sig',k)
            tit <- gsub('signal',paste0(paste(round(Popts,2),collapse='and'),'-day signal'),tit)
        }else if(k<Nper){
            y <- scale(tab[,k-Nmax+3])
            w <- sqrt(tab[,'PhotonCount']/max(tab[,'PhotonCount']))
            ey <- 0.01/w
            n <- colnames(tab)[k-Nmax+3]
        }else{
            y <- rep(1,length(y))
#            ey <- rnorm(length(y),0.1,0.001)
            n <- 'window'
        }
        if(length(y)>10){
            SigType <- 'kepler'
        }else{
            SigType <- 'circular'
        }

####save data
        Indices <- NULL
        if(n!='window'){
            if(PerType=='BFP'){
                p <- q <- 0
                if(noise=='MA' & length(y)>10) q <- 1
                if(noise=='AR' & length(y)>10) p <- 1
                per <- BFP(t,y,ey,Nma=q,Nar=p,Indices=Indices,ofac=ofac,model.type='man',fmin=fmin,fmax=fmax,quantify=TRUE,progress=FALSE,GP=FALSE,gp.par=rep(NA,3),noise.only=FALSE,Nsamp=Nsamp,sampling='combined',par.opt=NULL,renew=TRUE)
                if(SigType=='kepler'){
                    kep <- Circ2kep(per,basis='natural')
                }
            }else if(PerType=='GLST'){
                per <- glst(t,y,ey,fmax=fmax,ofac=ofac,fmin=fmin)
            }
        }else{
            per <- lsp(x=y,times=t,from=fmin,to=fmax,ofac=ofac)
        }
        if(any(names(per)=='lnBFs')){
            pp <- per$lnBFs
            cat('annual signal logL:',per$logLs[which.min(abs(per$ps-365.25))],'\n')
#            pp <- per$llmax
        }else{
            pp <- per$power
        }
        lnBF1[[n]] <- lnbf <- pp
        if(PerType=='GLST'){
            par1[[n]] <- par.opt <- unlist(per$par.opt)
        }else{
            par1[[n]] <- par.opt <- per$par.opt[1,]
        }
        Popt1[[n]] <- popt <- per$P[which.max(pp)]
        if(k<=Nmax){
            Popts <- c(Popts,popt)
        }

####plot periodograms
        fout <- paste0(dir,syntax,'_periodogram_',n,'.pdf')
        cat('\n',fout,'\n')
        pdf(fout,6,6)
        par(mar=c(4,4,4,1),mgp=c(2,1,0))
        if(n!='window'){
            ylab <- 'ln(BF)'
        }else{
            ylab <- 'Power'
        }
        ps <- per$P
        isort <- sort(ps,index.return=TRUE)$ix
#        Pact <- per$P[which.max(pp)]
        plot(ps[isort],pp[isort],type='l',log='x',ylim=range(0,median(pp),1.1*max(pp)),xaxt='n',yaxt='n',xlab='Period [day]',ylab=ylab,main=tit)
        abline(v= Popts,col=tcol('red',50),lwd=2)
        if(k>Nmax){
            abline(v= popt,col=tcol('green',50),lwd=2)
        }
#        lines(ps,pp)
        magaxis(side=c(1,2),tcl=-0.5)
        if(n!='window'){
            abline(h=5,lty=2,lwd=3,col='grey')
        }

        if(k<=Nmax){
            legend('topright',legend=paste0('Psig=',round(popt,2),'d'),bty='n',col='red')
        }else if(k<Nper){
            legend('topright',legend=paste0('Pact=',round(popt,2),'d'),bty='n',col='red')
        }else{
            legend('topright',legend=paste0('Pwindow=',round(popt,2),'d'),bty='n',col='red')
        }
        if(save.data){
            f1 <- gsub('pdf','txt',fout)
            cat(f1,'\n')
            write.table(cbind(ps[isort],pp[isort]),file=f1,row.names=FALSE,quote=FALSE,col.names=c('Period.day','lnBF'))
        }
        dev.off()

####plot phase curve
        if(k<=Nmax){
            fout <- paste0(dir,syntax,'_phase_',n,'.pdf')
            cat('\n',fout,'\n')
            pdf(fout,6,6)
            ##                    par(mfrow=c(2,1),mar=c(0,4,4,1))
            layout(matrix(data=c(1,2), nrow=2, ncol=1), widths=c(1,1), heights=c(2,1))
            par(mar=c(0,4,4,1))
            if(SigType!='kepler'){
                ysim0 <- par.opt[1]*cos(2*pi/popt*tsim)+par.opt[2]*sin(2*pi/popt*tsim)
                ysim <- ysim0+par.opt[3]+par.opt[4]*tsim
                ysims <- ysims+ysim
                y1 <- y-(par.opt[3]+par.opt[4]*t)
                if(save.data){
                    f0 <- gsub('\\.pdf','_OptPar.txt',fout)
                    cat(f0,'\n')
                    pars <- t(c(Popt=popt,par.opt))
                    write.table(pars,file=f0,quote=FALSE,row.names=FALSE)
                }
            }else{
                y1 <- y-kep$Pred$rv.red-kep$Pred$rv.trend
                df <- kep$df
                per$res <- y-kep$Pred$rv.kep-kep$Pred$rv.trend#preserve red noise component
                popt <- kep$ParKep$P1
                Popts[length(Popts)] <- popt
                df$tsim <- tsim
                par.sim0 <- par.sim <- kep$ParKep
                gamma <- gamma+par.sim$gamma
                gammadot <- gammadot+par.sim$beta
                par.sim0$beta <- par.sim0$gamma  <- 0
                ysim0 <- KeplerRv(par.sim0,df)$rv#without trend
                ysim <- KeplerRv(par.sim,df)$rv#with trend
                ysims <- ysims+ysim
                if(save.data){
                    f0 <- gsub('\\.pdf','_OptPar.txt',fout)
                    cat(f0,'\n')
                    write.table(t(unlist(kep$ParKep)),file=f0,quote=FALSE,row.names=FALSE)
                }
            }
            tt <- t%%popt
            ylim <- 1.1*range(y1,ysim0)
            plot(tt,y1,xlab=paste0('Orbital Phase [days]'),ylab='RV (data-model) [m/s]',main=paste(k,'signal'),xaxt='n',ylim=ylim)
            axis(side=1,label=FALSE)
            inds <- sort(tsim%%popt,index.return=TRUE)$ix
            lines(tsim[inds]%%popt,ysim0[inds],col='red')
            if(SigType!='kepler'){
                legend('topright',legend=paste0('P=',round(popt,2),'d'),bty='n',col='red')
            }else{
                legend('top',legend=paste0('P=',round(popt,2),'d; K=',round(kep$ParKep$K1,2),'m/s; e=',round(kep$ParKep$e1,2)),bty='n',col='red',horiz=TRUE,inset=c(0,-0.1),xpd=NA)
            }
            try(arrows(tt,y1+ey,tt,y1-ey,length=0.05,angle=90,code=3,col='grey'),TRUE)
        if(save.data){
            f1 <- gsub('\\.pdf','_DataPhase.txt',fout)
            cat(f1,'\n')
            write.table(cbind(tt,y1,ey),file=f1,row.names=FALSE,quote=FALSE,col.names=c('phase.day','RV-model','eRV'))
            f2 <- gsub('\\.pdf','_SimPhase.txt',fout)
            cat(f2,'\n')
            write.table(cbind(tsim[inds]%%popt,ysim0[inds]),file=f2,row.names=FALSE,quote=FALSE,col.names=c('phase.day','RVsim'))
}
####
            par(mar=c(4,4,0,1))
            plot(tt,per$res,xlab='Orbital Phase [day]',ylab='O-C [m/s]',col='grey50',ylim=range(min(per$res),max(per$res)+0.4*(max(per$res)-min(per$res))))
            abline(h=0,lty=2,lwd=3,col='darkgrey')
            legend('topright',legend=paste0('RMSraw=',format(round(sd(per$res),1),1),' m/s'),bty='n',text.col='black',cex=1.2,inset=c(0,-0.0))
            try(arrows(tt,per$res+tab[,3],tt,per$res-tab[,3],length=0.05,angle=90,code=3,col='grey'),TRUE)
### Add bined data
            data.bin <- wtb(t=t,x=per$res,ex=tab[,3],dt=0.5)
            tbin <- data.bin[,1]%%popt
            points(tbin,data.bin[,2],col='blue',pch=20,cex=0.8)
            legend('topleft',legend=paste0('RMSbin=',format(round(sd(data.bin[,2]),1),1),' m/s'),bty='n',text.col='blue',cex=1.2,inset=c(0,-0.0))
            try(arrows(tbin,data.bin[,2]-data.bin[,3],tbin,data.bin[,2]+data.bin[,3],length=0.03,angle=90,code=3,col='blue'),TRUE)
        if(save.data){
            f1 <- gsub('\\.pdf','_RawRes.txt',fout)
            cat(f1,'\n')
            write.table(cbind(tt,per$res,tab[,3]),file=f1,row.names=FALSE,quote=FALSE,col.names=c('phase.day','RVres','eRV'))
            f2 <-gsub('\\.pdf','_bin.txt',fout)
            cat(f2,'\n')
            write.table(cbind(tbin,data.bin[,2],data.bin[,3]),file=f2,row.names=FALSE,quote=FALSE,col.names=c('BinPhase.day','BinRV','eBinRV'))
}
            dev.off()
        }
        if(k==Nmax){
            fout <- paste0(dir,target,'_',PerType,'_fit_allsig.pdf')
            cat('\n',fout,'\n')
            pdf(fout,6,6)
            layout(matrix(data=c(1,2), nrow=2, ncol=1), widths=c(1,1), heights=c(2,1))
            par(mar=c(0,4,4,1))
            ylim0 <- range(tab[,2],ysims)
            ylim <- c(ylim0[1],ylim0[2]+0.2*(ylim0[2]-ylim0[1]))
            plot(t,tab[,2],xlab=paste0('orbital phase [days]'),ylab='RV (data-model) [m/s]',main='Combined fit',xaxt='n',ylim=ylim)
            tyr <- jd2yr(cbind(tmin,t))
            for(Ndig in 3:0){
                labs <- seq(min(tyr),max(tyr),by=0.001)
                labs <- unique(round(labs,Ndig))
                if(length(labs)<=20) break()
            }
            tjd <- rowSums(yr2jd(labs))-tmin
            par(mgp=c(2,0.5,0))
            axis(side=3,at=tjd,label=labs)
            par(mgp=c(3,1,0))
            axis(side=1,label=FALSE)
            lines(tsim,ysims,col='red')
            points(t,tab[,2])
#            legend('topright',legend=c(paste0('gamma dot=',round(gammadot*365.25,1),'m/s/year',paste0('gamma=',round(gamma),'m/s'))),,bty='n',col='red')
            if(SigType=='kepler' & (max(tab[,1])-min(tab[,1]))>10){
                legend('topleft',legend=bquote(gamma==.(round(gamma,1))~'m/s'),bty='n',col='red')
                legend('topright',legend=bquote(dot(gamma)==.(round(gammadot*365.25,1))~'m/s/year'),bty='n',col='red')
            }
            try(arrows(t,tab[,2]+tab[,3],t,tab[,2]-tab[,3],length=0.05,angle=90,code=3,col='grey'),TRUE)
        if(save.data){
            f1 <- gsub('\\.pdf','_data.txt',fout)
            cat(f1,'\n')
            write.table(cbind(t,tab[,2],tab[,3]),file=f1,row.names=FALSE,quote=FALSE,col.names=c('t-tmin','RV','eRV'))
            f2 <-gsub('\\.pdf','_fit.txt',fout)
            cat(f2,'\n')
            write.table(cbind(tsim,ysim),file=f2,row.names=FALSE,quote=FALSE,col.names=c('tsim','RVsim'))
}
####
            par(mar=c(4,4,0,1))
            plot(t,per$res,xlab=paste0('BJD - ',min(tab[,1])),ylab='O-C [m/s]',col='grey50',ylim=range(min(per$res),max(per$res)+0.2*(max(per$res)-min(per$res))))
            abline(h=0,lty=2,lwd=3,col='darkgrey')
            legend('topright',legend=paste0('RMSraw=',format(round(sd(per$res),1),1),' m/s'),bty='n',text.col='black',cex=1.2,inset=c(0,-0.0))
            try(arrows(t,per$res+tab[,3],t,per$res-tab[,3],length=0.05,angle=90,code=3,col='grey'),TRUE)
### Add bined data
            data.bin <- wtb(t=t,x=per$res,ex=tab[,3],dt=0.5)
            tbin <- data.bin[,1]
            points(tbin,data.bin[,2],col='blue',pch=20,cex=0.8)
            legend('topleft',legend=paste0('RMSbin=',format(round(sd(data.bin[,2]),1),1),' m/s'),bty='n',text.col='blue',cex=1.2,inset=c(0,-0.0))
            try(arrows(tbin,data.bin[,2]-data.bin[,3],tbin,data.bin[,2]+data.bin[,3],length=0.03,angle=90,code=3,col='blue'),TRUE)
        if(save.data){
            f1 <- gsub('\\.pdf','_RawRes.txt',fout)
            cat(f1,'\n')
            write.table(cbind(t,per$res,tab[,3]),file=f1,row.names=FALSE,quote=FALSE,col.names=c('t-tmin','RVres','eRVres'))
            f2 <- gsub('\\.pdf','_BinRes.txt',fout)
            cat(f2,'\n')
            write.table(cbind(tbin,data.bin[,2],data.bin[,3]),file=f2,row.names=FALSE,quote=FALSE,col.names=c('tbin-tmin','RVbinRes','eRVbinRes'))
}
            dev.off()

####residual periodogram
            if(PerType=='BFP'){
                p <- q <- 0
                if(noise=='MA') q <- 1
                if(noise=='AR') p <- 1
                per <- BFP(t,per$res,ey,Nma=q,Nar=p,Indices=Indices,ofac=ofac,model.type='man',fmin=fmin,fmax=fmax,quantify=TRUE,progress=FALSE,GP=FALSE,gp.par=rep(NA,3),noise.only=FALSE,Nsamp=Nsamp,sampling='combined',par.opt=NULL,renew=TRUE)
            }else if(PerType=='GLST'){
                per <- glst(t,per$res,ey,fmax=fmax,ofac=ofac,fmin=fmin)
            }
            fout <- paste0(dir,syntax,'_periodogram_res.pdf')
            cat('\n',fout,'\n')
            pdf(fout,6,6)
            par(mar=c(4,4,4,1),mgp=c(2,1,0))
            if(n!='window'){
                ylab <- 'ln(BF)'
            }else{
                ylab <- 'Power'
            }
            ps <- per$P
            if(!any(names(per)=='lnBFs')){
                pp <- per$power
            }else{
                pp <- per$lnBFs
            }
            lnBF1[['BFPres']] <- lnbf <- pp
            par1[['BFPres']] <- par.opt <- unlist(per$par.opt)[1,]
            Popt1[['BFPres']] <- popt <- per$P[which.max(pp)]
            tit <- paste0(PerType,' for raw RV subtracted by ',paste(round(Popts,2),collapse='and'),'-day signal')
            isort <- sort(ps,index.return=TRUE)$ix
            plot(ps[isort],pp[isort],type='l',log='x',ylim=range(0,median(pp),1.1*max(pp)),xaxt='n',yaxt='n',xlab='Period [day]',ylab=ylab,main=tit)
            abline(v= Popts,col=tcol('red',50))
            lines(ps,pp)
            magaxis(side=c(1,2),tcl=-0.5)
            if(n!='window'){
                abline(h=5,lty=2,lwd=3,col='grey')
            }
            abline(v= popt,col=tcol('orange',50),lwd=2)
            legend('topright',legend=paste0('Pres=',round(popt,2),'d'),bty='n',col='red')
        if(save.data){
            f1 <- gsub('\\.pdf','.txt',fout)
            cat(f1,'\n')
            write.table(cbind(ps[isort],pp[isort]),file=f1,row.names=FALSE,quote=FALSE,col.names=c('Period.day','lnBF'))
}
            dev.off()
        }
    }
#    out[[jj]] <- list(lnBF=lnBF1,)
#    lnBF[[target]] <- lnBF1
#    pars[[target]] <- par1
#    Popt[[target]] <- Popt1
}
