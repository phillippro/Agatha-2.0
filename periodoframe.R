library(minpack.lm)
#library(lomb)
tol11 <- 1e-15
tol22 <- 1e-20
tol33 <- 1e-30
off <- 0
#trend <- FALSE
solve.try <- function(lin.mat,vec.rh){
    white.par <- try(solve(lin.mat,vec.rh),TRUE)#gamma,beta,dj
    if(class(white.par)=='try-error')     white.par <- try(solve(lin.mat,vec.rh,tol=tol33),TRUE)
#    if(class(white.par)=='try-error')     white.par <- try(solve(lin.mat,vec.rh,tol=tol3),TRUE)
    if(class(white.par)=='try-error')     white.par <- try(solve(nearPD(lin.mat)$mat,vec.rh,tol=tol33),TRUE)
    if(class(white.par)=='try-error'){
        cat('lin.mat reversion error!\n')
        cat('lin.mat=',lin.mat,'\n')
        cat('vec.rh=',vec.rh,'\n')
    }
    return(white.par)
}

xy2phi <- function(x,y){
    phi <- atan(y/x)
    inds <- which(x<0)
    phi[inds] <- phi[inds]+pi
    inds <- which(x>=0 & y<0)
    phi[inds] <- phi[inds]+2*pi
    return(phi)
}

KeplerSig <- function(par,df,basis='natural'){
    data <- df$data
    t <- data[,1]-min(data[,1])
    y <- data[,2]
    dy <- data[,3]

    Np <- length(grep('^K',names(par)))
    yred <- yma <- yar <- ytrend <- yproxy <- ysig <- ypred <- rep(0,length(t))

###signal
    if(Np>0){
        for(h in 1:Np){
            P <- par[[paste0('P',h)]]
            K <- par[[paste0('K',h)]]
            if(basis=='natural'){
                e <- par[[paste0('e',h)]]
                Mo <- par[[paste0('Mo',h)]]
                omega <- par[[paste0('omega',h)]]
            }else if(basis=='linear'){
                Tc <- par[[paste0('Tc',h)]]
                if(any(names(par)==paste0('sqresinw',h))){
                    yy <- sqresin <- par[[paste0('sqresinw',h)]]
                    xx <- sqrecos <- par[[paste0('sqrecosw',h)]]
                    e <- sqresin^2+sqrecos^2
                }else{
                    yy <- esin <- par[[paste0('esinw',h)]]
                    xx <- ecos <- par[[paste0('ecosw',h)]]
                    e <- sqrt(esin^2+ecos^2)
                }
                if(xx!=0){
                    omega <- atan(yy/xx)
                }else{
                    omega <- atan(yy/1e-6)
                }
                Mo <- getM0(e=e,omega=omega,P=P,T=Tc,T0=tmin,type='primary')
            }
#            cat('e=',e,'\n')
            m <- (Mo+2*pi*t/P)%%(2*pi)#mean anomaly
#            cat('m=',head(m),'\n')
#            cat('head(m)=',head(m),'\n')
#            cat('e=',e,'\n')
            E <- kep.mt2(m,e)
            T <- 2*atan(sqrt((1+e)/(1-e))*tan(E/2))#true anomaly
            ysig <- ysig+K*(cos(omega+T)+e*cos(omega))
        }
    }

###trend
    if(any(names(par)=='gamma')){
        ytrend <- ytrend+par[['gamma']]
    }
    if(any(names(par)=='beta')){
        ytrend <- ytrend+par[['beta']]*t
    }

###get noise model parameters
    m <- unlist(par[grep('^m\\d$',names(par))])
    l <- unlist(par[grep('^l\\d$',names(par))])
    d <- unlist(par[grep('^d\\d$',names(par))])
    Nma <- length(m)
    Nar <- length(l)
    NI <- length(d)

###proxy
    if(NI>0){
        yproxy <- unlist(d%*%t(df$Indices))
    }

    ypred <- ysig+yproxy+ytrend
###MA
    if(Nma>0){
        es <- c()
        rs <- c()
        for(i in 1:Nma){
            ei <- c(rep(0,i),exp(-abs(t[-(length(t)+1-(1:i))]-t[-(1:i)])/exp(par[['logtau']])))
            es <- rbind(es,ei)
            ri <- c(rep(0,i),ypred[-(length(t)+1-(1:i))])
            rs <- rbind(rs,ri)
        }
        yma <- as.numeric(m%*%(es*(df$ys-rs)))
    }

###AR
    if(Nar>0){
        ea <- c()
        for(i in 1:Nar){
            ei <- c(rep(0,i),exp(-abs(t[-(length(t)+1-(1:i))]-t[-(1:i)])/exp(par[['logtauAR']])))
            ea <- rbind(ea,ei)
        }
        yar <- l%*%(ea*df$ya)
    }

###sum
    yred <- yma+yar
    ypred <- ypred+yma+yar
    res <- y-ypred

    return(list(res=res,y=ypred,ysig=ysig,ytrend=ytrend,yproxy=yproxy,yred=yred,yma=yma,yar=yar))
}

KeplerRes <- function(par,df){
    y <- df$data[,2]
    dy <- df$data[,3]
    v <- KeplerSig(par,df,df$basis)$y
    neglnLs <- (y-v)^2/(2*(dy^2+par[['sj']]^2)) + 0.5*log(dy^2+par[['sj']]^2)+log(sqrt(2*pi))+off
    if(any(neglnLs<0)) neglnLs[neglnLs<0] <- 0
    sqrneglnLs <- sqrt(neglnLs)
    return(sqrneglnLs)
}

partI <- function(d,Is){
    if(length(d)==1){
        d*Is[,1]
    }else{
        d%*%t(Is)
    }
}
wI <- function(w,Is,NI){
    unlist(lapply(1:NI,function(i) sum(w*Is[,i])))
}
me <- function(m,t,logtau,x){
    if(length(m)==1){
        m*exp(-abs(t[-1]-t[-length(t)])/exp(logtau))*x[-length(x)]
    }else{
#        unlist(lapply(1:length(m), function(i) m[i]*exp(-abs(t[-c(1:i)]-t[-(length(t)+1-c(1:i))])/exp(logtau))*x[-(length(x)+1-c(1:i))]))
    }
}

CircularSig <- function(par,df){
#####setting up
    if(!is.null(df$par.fix)){
        phi <- df$par.fix$phi
        omega <- df$par.fix$omega
#        if(any(names(par.fix)=='logtau') | any(names(par.fix)=='logtauAR')){
#            logtau <- logtauAR <- df$par.fix$logtau
#        }
    }
    Indices <- df$Indices
    type <- df$type
    NI <- df$NI
    Nma <- df$Nma
    Nar <- df$Nar
    t <- df$data[,1]
    y <- df$data[,2]
    dy <- df$data[,3]
    if(!is.null(df$omega)){
        if(!is.na(df$omega)){
            if(Nma>0) logtau <- log(2*pi/df$omega)
            if(Nar>0) logtauAR <- log(2*pi/df$omega)
        }else{
            logtau <- par$logtau
            logtauAR <- par$logtauAR
        }
    }else{
        logtau <- par$logtau
        logtauAR <- par$logtauAR
    }
    m <- unlist(par[grepl('m\\d',names(par))])#
    l <- unlist(par[grepl('l\\d',names(par))])#
    sj <- par$sj

#####notations
    W <- sum(1/(sj^2+dy^2))
    w <- 1/(sj^2+dy^2)/W

###different components of the y value
    yma <- yar <- yred <- ytrend <- yproxy <- ysig <- 0

###MA
    es <- c()
    if(Nma>0){
        for(i in 1:Nma){
            ei <- c(rep(0,i),exp(-abs(t[-(length(t)+1-(1:i))]-t[-(1:i)])/exp(logtau)))
            es <- rbind(es,ei)
        }
        yp <- y-m%*%(es*df$ys)
        wp <- 1-m%*%es
        tp <- t-m%*%(es*df$ts)
        if(type=='period'){
            cp <- cos(omega*t)-m%*%(es*df$cs)
            sp <- sin(omega*t)-m%*%(es*df$ss)
        }
    }else{
        yp <- y
        wp <- 1
        tp <- t
        if(type=='period'){
            cp <- cos(omega*t)
            sp <- sin(omega*t)
        }
    }
#    cat('head(m%*%(es*df$ys))=',head(m%*%(es*df$ys)),'\n')
#    cat('head(m)=',head(m),'\n')

###AR
    ea <- c()
    if(Nar>0){
        for(i in 1:Nar){
            ei <- c(rep(0,i),exp(-abs(t[-(length(t)+1-(1:i))]-t[-(1:i)])/exp(logtauAR)))
            ea <- rbind(ea,ei)
        }
        yar <- l%*%(ea*df$ya)
        yp <- yp-yar
    }

    WIp <- ip <- TIp <- YIp <- CIp <- SIp <- c()
    if(NI>0){
        IIp <- array(data=NA,dim=c(NI,NI))
        for(j in 1:NI){
            if(Nma>0){
                ip <- rbind(ip,t(Indices[,j,drop=FALSE])-m%*%(es*df$Is[,j,]))
            }else{
                ip <- rbind(ip,t(Indices[,j,drop=FALSE]))
            }
            YIp <- c(YIp,sum(w*yp*ip[j,]))
            WIp <- c(WIp,sum(w*wp*ip[j,]))
            TIp <- c(TIp,sum(w*tp*ip[j,]))
            if(type=='period'){
                CIp <- c(CIp,sum(w*cp*ip[j,]))
                SIp <- c(SIp,sum(w*sp*ip[j,]))
            }
        }
        for(j in 1:NI){
            for(i in j:NI){
                IIp[j,i] <- IIp[i,j] <- sum(w*ip[i,]*ip[j,])
            }
        }
    }

    WWp <- sum(w*wp*wp)
    YWp <- sum(w*wp*yp)
    YTp <- sum(w*yp*tp)
    Wp <- sum(w*wp)
    WTp <- sum(w*wp*tp)
    TTp <- sum(w*tp^2)
    if(type=='period'){
        Cp <- sum(w*cp)
        Sp <- sum(w*sp)
        YCp <- sum(w*yp*cp)
        YSp <- sum(w*yp*sp)
        CCp <- sum(w*cp^2)
        SSp <- sum(w*sp^2)
        CSp <- sum(w*cp*sp)
        CTp <- sum(w*cp*tp)
        CWp <- sum(w*cp*wp)
        STp <- sum(w*sp*tp)
        SWp <- sum(w*sp*wp)
    }
    if(type=='noise'){
        if(NI>0){
            lin.mat <- matrix(c(WWp,WTp,WIp,WTp,TTp,TIp),byrow=TRUE,ncol=2+NI)
            for(j in 1:NI){
                lin.mat <- rbind(lin.mat,c(WIp[j],TIp[j],IIp[j,]))
            }
            vec.rh <- c(YWp,YTp,YIp)
        }else{
            dI <- 0
            d <- 0
            lin.mat <- matrix(c(WWp,WTp,WTp,TTp),byrow=TRUE,ncol=2)
            vec.rh <- c(YWp,YTp)
        }
###optimized parameterse for the trend model
    }else if(type=='period'){
#####calculate the optimal white noise model parameters:gamma,beta,dj
        lin.mat <- c()
        if(NI>0){
            for(j in 1:(4+NI)){
                if(j==1){
                    lin.mat <- rbind(lin.mat,c(CCp,CSp,CWp,CTp,CIp))
                }
                if(j==2){
                    lin.mat <- rbind(lin.mat,c(CSp,SSp,SWp,STp,SIp))
                }
                if(j==3){
                    lin.mat <- rbind(lin.mat,c(CWp,SWp,WWp,WTp,WIp))
                }
                if(j==4){
                    lin.mat <- rbind(lin.mat,c(CTp,STp,WTp,TTp,TIp))
                }
                if(j>4){
                    lin.mat <- rbind(lin.mat,c(CIp[j-4],SIp[j-4],WIp[j-4],TIp[j-4],IIp[j-4,]))
                }
            }
            vec.rh <- c(YCp,YSp,YWp,YTp,YIp)
        }else{
            lin.mat <- rbind(lin.mat,c(CCp,CSp,CWp,CTp))
            lin.mat <- rbind(lin.mat,c(CSp,SSp,SWp,STp))
            lin.mat <- rbind(lin.mat,c(CWp,SWp,WWp,WTp))
            lin.mat <- rbind(lin.mat,c(CTp,STp,WTp,TTp))
            vec.rh <- c(YCp,YSp,YWp,YTp)
        }
    }
    white.par <- solve.try(lin.mat,vec.rh)
    ind0 <- length(white.par)-NI-1
    r <- ytrend <- white.par[ind0]+white.par[ind0+1]*t
    if(NI>0){
        yproxy <- white.par[(ind0+2):length(white.par)]%*%t(Indices)
        r <- r+yproxy
    }
    if(type=='period'){
        ysig <- white.par[1]*cos(omega*t)+white.par[2]*sin(omega*t)
        r <- r+ysig
    }

####the residual in the MA model is not subtracted by the AR component so that the MA and AR components could be independent
    v <- as.numeric(r)
    if(Nma>0){
        r <- as.numeric(r)
        rs <- c()
        for(i in 1:Nma){
            ri <- c(rep(0,i),r[-(length(t)+1-(1:i))])
            rs <- rbind(rs,ri)
        }
        yma <- as.numeric(m%*%(es*(df$ys-rs)))
        v <- v+yma
    }

###the AR component is added to the total mdoel prediction after the MA component is added
    if(Nar>0) v <- v+yar
    yred <- yma+yar
    opt.par <- white.par
    return(list(v=v,ysig=ysig,ytrend=ytrend,yproxy=yproxy,yred=yred,yma=yma,yar=yar,par=opt.par,white.par=white.par,NI=NI,Indices=Indices,res=y-v))
}

CircularSim <- function(par,df,popt,tsim){
###red component
    tmp <- CircularSig(par,df)
    data <- df$data
    yma <- yar <- yproxy <- ytrend <- ysig <- 0
    t <- data[,1]-min(data[,1])

###MA component
    if(df$Nma>0){
        yma.fun <- approxfun(t,tmp$yma)
        yma <- yma.fun(tsim)
    }

###AR component
    if(df$Nar>0){
        yar.fun <- approxfun(t,tmp$yar)
        yar <- yar.fun(tsim)
    }

###proxy component
    if(df$NI>0){
        yproxy.fun <- approxfun(t,tmp$yproxy)
        yproxy <- yproxy.fun(tsim)
    }

##signal
    if(any(names(par)=='A')) ysig <- par$A*cos(2*pi/popt*tsim)+par$B*sin(2*pi/popt*tsim)

##trend
    ytrend <- 0
    if(any(names(par)=='gamma')){
        ytrend <- par$gamma
    }
    if(any(names(par)=='beta')){
        ytrend <- ytrend+par$beta*tsim
    }
    y <- ysig+ytrend+yproxy+yma+yar

    return(list(y=y,ysig=ysig,ytrend=ytrend,yproxy=yproxy,yred=yma+yar))
}

KeplerSim <- function(par,df,tsim,basis='natural'){
###simulate Keplerian signal and various noise components
    yma <- yar <- yproxy <- ytrend <- ysig <- 0
    Np <- length(grep('^K',names(par)))
    if(Np>0){
        for(h in 1:Np){
            P <- par[[paste0('P',h)]]
            K <- par[[paste0('K',h)]]
            if(basis=='natural'){
                e <- par[[paste0('e',h)]]
                Mo <- par[[paste0('Mo',h)]]
                omega <- par[[paste0('omega',h)]]
            }else if(basis=='linear'){
                Tc <- par[[paste0('Tc',h)]]
                if(any(names(par)==paste0('sqresinw',h))){
                    yy <- sqresin <- par[[paste0('sqresinw',h)]]
                    xx <- sqrecos <- par[[paste0('sqrecosw',h)]]
                    e <- sqresin^2+sqrecos^2
                }else{
                    yy <- esin <- par[[paste0('esinw',h)]]
                    xx <- ecos <- par[[paste0('ecosw',h)]]
                    e <- sqrt(esin^2+ecos^2)
                }
                if(xx!=0){
                    omega <- atan(yy/xx)
                }else{
                    omega <- atan(yy/1e-6)
                }
                Mo <- getM0(e=e,omega=omega,P=P,T=Tc,T0=tmin,type='primary')
            }
            m <- (Mo+2*pi*tsim/P)%%(2*pi)#mean anomaly
            E <- kep.mt2(m,e)
            T <- 2*atan(sqrt((1+e)/(1-e))*tan(E/2))#true anomaly
            ysig <- ysig+K*(cos(omega+T)+e*cos(omega))
        }
    }

###red component
    tmp <- KeplerSig(par,df,basis=basis)
    data <- df$data
    t <- data[,1]-min(data[,1])

###MA component
    if(df$Nma>0){
        yma.fun <- approxfun(t,tmp$yma)
        yma <- yma.fun(tsim)
    }

###AR component
    if(df$Nar>0){
        yar.fun <- approxfun(t,tmp$yar)
        yar <- yar.fun(tsim)
    }

###proxy component
    if(df$NI>0){
        yproxy.fun <- approxfun(t,tmp$yproxy)
        yproxy <- yproxy.fun(tsim)
    }

###trend
    if(any(names(par)=='gamma')){
        ytrend <- par$gamma
    }
    if(any(names(par)=='beta')){
        ytrend <- ytrend+par$beta*tsim
    }

###sum
    yred <- yma + yar
    y <- ysig+ytrend+yproxy+yma+yar

    return(list(y=y,ysig=ysig,ytrend=ytrend,yproxy=yproxy,yred=yma+yar))
}

rv.white <- function(par,df){
#####setting up
    var <- names(df)
    for(k in 1:length(var)){
        assign(var[k],df[[var[k]]])
    }
    if(!is.null(par.fix)){
        phi <- par.fix$phi
        omega <- par.fix$omega
    }
    t <- data[,1]
    y <- data[,2]
    dy <- data[,3]
    sj <- par$sj
#####notations
    W <- sum(1/(sj^2+dy^2))
    w <- 1/(sj^2+dy^2)/W
    if(type=='period'){
        c <- cos(omega*t)
        s <- sin(omega*t)
    }
    Y <- sum(w*y)
    YT <- sum(w*y*t)
    T <- sum(w*t)
    TT <- sum(w*t^2)
    if(NI>0){
        I <- TI <- YI <- c()
        II <- array(data=NA,dim=c(NI,NI))
        for(j in 1:NI){
            I <- c(I,sum(w*Indices[,j]))
            TI <- c(TI,sum(w*t*Indices[,j]))
            YI <- c(YI,sum(w*y*Indices[,j]))
            for(i in j:NI){
                II[j,i] <- II[i,j] <- sum(w*Indices[,i]*Indices[,j])
            }
        }
        if(type=='period'){
            CI <- SI <- c()
            for(j in 1:NI){
                CI <- c(CI,sum(w*cos(omega*t)*Indices[,j]))
                SI <- c(SI,sum(w*sin(omega*t)*Indices[,j]))
            }
        }
    }
####
#    if w is normalized, 1 is used; otherwise, W is used.
    if(type=='noise'){
        if(NI>0){
            lin.mat <- matrix(c(1,T,I,T,TT,TI),byrow=TRUE,ncol=2+NI)
            for(j in 1:NI){
                lin.mat <- rbind(lin.mat,c(I[j],TI[j],II[j,]))
            }
            vec.rh <- c(Y,YT,YI)
            white.par <- solve.try(lin.mat,vec.rh)
            indd <- (length(white.par)-NI+1):length(white.par)
        }else{
            dI <- 0
            d <- 0
            lin.mat <- matrix(c(1,T,T,TT),byrow=TRUE,ncol=2)
            vec.rh <- c(Y,YT)
            white.par <- solve.try(lin.mat,vec.rh)

        }
###optimized parameterse for the trend model
    }else if(type=='period'){
        C <- sum(w*cos(omega*t))
        S <- sum(w*sin(omega*t))
        CC <- sum(w*cos(omega*t)^2)
        SS <- sum(w*sin(omega*t)^2)
        CS <- sum(w*sin(omega*t)*cos(omega*t))
        CT <- sum(w*t*cos(omega*t))
        ST <- sum(w*t*sin(omega*t))
        YC <- sum(w*y*cos(omega*t))
        YS <- sum(w*y*sin(omega*t))
#####calculate the optimal white noise model parameters:gamma,beta,dj
        lin.mat <- c()#matrix(c(CC,CS,C,CT,CIW,T,I,T,TT,TI),byrow=TRUE,ncol=2+NI)
        if(NI>0){
            for(j in 1:(4+NI)){
                if(j==1){
                    lin.mat <- rbind(lin.mat,c(CC,CS,C,CT,CI))
                }
                if(j==2){
                    lin.mat <- rbind(lin.mat,c(CS,SS,S,ST,SI))
                }
                if(j==3){
                    lin.mat <- rbind(lin.mat,c(C,S,W,T,I))
                }
                if(j==4){
                    lin.mat <- rbind(lin.mat,c(CT,ST,T,TT,TI))
                }
                if(j>4){
                    lin.mat <- rbind(lin.mat,c(CI[j-4],SI[j-4],I[j-4],TI[j-4],II[j-4,]))
                }
            }
            vec.rh <- c(YC,YS,Y,YT,YI)
            white.par <- solve.try(lin.mat,vec.rh)
            indd <- (length(white.par)-NI+1):length(white.par)
        }else{
            lin.mat <- rbind(lin.mat,c(CC,CS,C,CT))
            lin.mat <- rbind(lin.mat,c(CS,SS,S,ST))
            lin.mat <- rbind(lin.mat,c(C,S,W,T))
            lin.mat <- rbind(lin.mat,c(CT,ST,T,TT))
            vec.rh <- c(YC,YS,Y,YT)
            white.par <- solve.try(lin.mat,vec.rh)
        }
    }
    ind0 <- length(white.par)-NI-1
    r <- white.par[ind0]+white.par[ind0+1]*t
    if(NI>0){
        r <- r+white.par[(ind0+2):length(white.par)]%*%Indices
    }
    if(type=='period'){
        r <- r+white.par[1]*cos(omega*t)+white.par[2]*sin(omega*t)
    }
    if(NI>0){
        indd <- (length(white.par)-NI+1):length(white.par)
        opt.par <- c(white.par[-indd],d=white.par[indd])
    }else{
        opt.par <- white.par
    }
    return(list(v=r,par=opt.par))
}
celerite <- function(t,y,dy,s,term){
    a <- term[,1]
    b <- term[,2]
    c <- term[,3]
    d <- term[,4]
###a, b, c, d are for complex terms
###ar, cr are for real terms
#    Jc <- length(a)
#    Jr <- length(ar)
#    J <- Jc+Jr
#    R <- 2*J-Jr
    R <- 2*length(a)
    J <- length(a)
    N <- length(t)
    phi <- Wp <- Up <- Vp <- array(NA,dim=c(N,R))
    A <- D <- array(0,dim=c(N,N))
    S <- array(NA,dim=c(N,R,R))
#    diag(D) <- s^2+sum(a)+sum(ar)
    js <- 1:J
    ns <- 1:N
    dt <- outer(d,t,'*')
    ct <- outer(c,t,'*')
    cdt <- cos(dt)
    sdt <- sin(dt)
    nct <- exp(-ct)
    pct <- exp(ct)
    diag(A) <- dy^2+s^2+sum(a)#add white jitter
####pre-conditioned variables
    Up[,2*js-1] <- a%*%cdt+b%*%sdt
    Up[,2*js] <- a%*%sdt-b%*%cdt
    Wp[,2*js-1] <- cdt
    Wp[,2*js] <- sdt
    ect <- outer(c,t[2:N]-t[1:(N-1)],'*')
    phi[2:N,2*js-1] <- phi[2:N,2*js] <- exp(-ect)
    phi[1,] <- 0
####calculate S, D and W
    S[1,,] <- 0
    D[1,1] <- A[1,1]
    Wp[1,] <- 1/D[1,1]*Wp[1,]
#    cat('A[N,N]=',A[N,N],'\n')

    for(n in 2:N){
        S[n,,] <- outer(phi[n-1,],phi[n-1,],'*')*(S[n-1,,]+D[n-1,n-1]*outer(Wp[n-1,],Wp[n-1,],'*'))
        D[n,n] <- A[n,n]-Up[n,]%*%(S[n,,]%*%Up[n,])
        Wp[n,] <- 1/D[n,n]*(Wp[n,]-Up[n,]%*%S[n,,])
    }
    diagD <- diag(D)
    diagD[diagD<0] <- min(abs(diagD))
    lndetKs <- log(diagD)#sometimes NAs generated
    lndetK <- sum(lndetKs)
    f <- array(NA,dim=c(N,R))
    z <- rep(NA,N)
    z[1] <- y[1]
    f[1,] <- 0
    ykys <- y[1]^2/D[1,1]
    for(n in 2:N){
        f[n,] <- phi[n,]*(f[n-1,]+Wp[n-1,]*z[n-1])
        z[n] <- y[n]-sum(Up[n,]*f[n,])
        ykys <- c(ykys,z[n]^2/D[n,n])
    }
    yky <- sum(ykys)
    return(list(lndetK=lndetK,yky=yky,ykys=ykys,lndetKs=lndetKs))
}
sho.term <- function(S0,Q,w0){
    a <- b <- c <- d<- 0
    if(Q<0.5){
        f <- sqrt(1-4*Q^2)
        a <- 0.5*S0*w0*Q*c(1+1/f,1-1/f)
        c <- 0.5*w0/Q*c(1-f,1+f)
    }else{
        f <- sqrt(4.0*Q^2-1)
        a <- S0*w0*Q
        b <- S0*w0*Q/f
        c <- 0.5*w0/Q
        d <- 0.5*w0*f/Q
    }
    return(cbind(a,b,c,d))
}
CircularRes <- function(par,df){
###par is free variables
###data and fixed parameters are in the df list
    y <- df$data[,2]
    dy <- df$data[,3]
#    cat('df$omega=',df$omega,'\n')
    v <- CircularSig(par,df)$v
    sj <- par$sj
    if(!df$GP){
        neglnLs <- (y-v)^2/(2*(dy^2+sj^2)) + 0.5*log(dy^2+sj^2)+log(sqrt(2*pi))+off
        if(any(neglnLs<0)) neglnLs[neglnLs<0] <- 0
        sqrneglnLs <- sqrt(neglnLs)
    }else{
        t <- df$data[,1]
        logProt <- c(df$logProt,par$logProt)
        logtauGP <- c(df$logtauGP,par$logtauGP)
        sigmaGP <- c(df$sigmaGP,par$sigmaGP)
        Q <- exp(logtauGP)*pi/exp(logProt)
#cat('Q=',Q,'\n')
#cat('w0=',2*pi/exp(logProt),'\n')
        term <- sho.term(S0=sigmaGP,Q=Q,w0=2*pi/exp(logProt))
        out <- celerite(t=t,y=y-v,dy=dy,s=sj,term=term)
#        cat('out$lndetKs=',out$lndetKs,'\n')
#        cat('out$ykys=',out$ykys,'\n')
        neglogl <- 1/2*out$lndetKs+1/2*out$ykys+length(t)/2*log(2*pi)
#        cat('min(neglogl)=',min(neglogl),'\n')
        neglogl <- neglogl+off
#        if(any(neglogl<0)) cat('neglogl=',neglogl,'\n')
        if(any(neglogl<0)) neglogl[neglogl<0] <- 0
        sqrneglnLs <- sqrt(neglogl)#for likelihood maximization
    }
    return(sqrneglnLs)
}

######express other parameters as functions of correlated noise parameters, optimize correlated noise parameters using the LM algorithm
sopt <- function(omega,phi,Nma,Nar,Indices,data,type='noise',par.low,par.up,start,noise.only,GP,gp.par,Nrep=1){
    par.fix <- list()
    if(type=='period'){
        par.fix <- list(omega=omega,phi=phi)
    }
#else{
#        par.fix <- NULL
#    }
    if(noise.only & !all(is.na(omega))){
        if(GP){
            gp.par[2] <- log(2*pi/omega)
        }else{
###assume that the MA and AR time scale is the same; this could be generalized
#            par.fix$logtauAR <- par.fix$logtau <- log(2*pi/omega)
        }
    }
    NI <- 0
    if(!is.null(Indices)){
        NI <- ncol(Indices)
    }
    ts <- cs <- ys <- ya <- ss <- c()
    t <- data[,1]
    y <- data[,2]
    dy <- data[,3]
    if(NI>0){
        Is <- array(data=NA,dim=c(Nma,NI,length(t)))
    }else{
        Is <- c()
    }
    if(Nma>0){
        for(i in 1:Nma){
            yi <- c(rep(0,i),y[-(length(t)+1-(1:i))])
            ys <- rbind(ys,yi)
            ti <- c(rep(0,i),t[-(length(t)+1-(1:i))])
            ts <- rbind(ts,ti)
            if(type=='period'){
                ci <- c(rep(0,i),cos(omega*t[-(length(t)+1-(1:i))]))
                cs <- rbind(cs,ci)
                si <- c(rep(0,i),sin(omega*t[-(length(t)+1-(1:i))]))
                ss <- rbind(ss,si)
            }
            if(NI==1){
                Is[i,,] <- c(rep(0,i*NI),Indices[-(length(t)+1-(1:i)),1])
            }else if(NI>1){
                Is[i,,] <- c(matrix(rep(0,i*NI),nrow=NI),t(Indices[-(length(t)+1-(1:i)),]))
            }
        }
    }
    if(Nar>0){
        for(i in 1:Nar){
            yi <- c(rep(0,i),y[-(length(t)+1-(1:i))])
            ya <- rbind(ya,yi)
        }
    }
    df <- list(data=data,Indices=Indices,par.fix=par.fix,Nma=Nma,Nar=Nar,NI=NI,type=type,ts=ts,cs=cs,ss=ss,Is=Is,ys=ys,GP=GP,ya=ya)
    if(GP){
        if(!is.na(gp.par[1])){
            df$sigmaGP <- gp.par[1]
        }
        if(!is.na(gp.par[2])){
            df$logProt <- gp.par[2]
        }
        if(!is.na(gp.par[3])){
            df$logtauGP <- gp.par[3]
        }
    }

    if(noise.only){
        df$omega <- omega
    }
##########numerical fit
    Ntry <- 10
#    cat('start=',unlist(start),'\n')
    start0 <- start
    Ls <- par.ini <- c()
    for(k in 1:Nrep){
        if(k>1){
            start <- detIni(start0,par.low,par.up)
        }
        out0 <- nls.lm(par = start,lower=par.low,upper=par.up,fn = CircularRes,df=df,control=nls.lm.control(maxiter=1024))
        Ls <- c(Ls,-sum(out0$fvec^2-off))
        par.ini <- rbind(par.ini,unlist(start))
    }
    start <- as.list(par.ini[which.max(Ls),])
    names(start) <- names(par.low)
    out <- nls.lm(par = start,lower=par.low,upper=par.up,fn = CircularRes,df=df,control=nls.lm.control(maxiter=1024))

#####retrieve parameters
    logL <- -sum(out$fvec^2-off)
    lnls <- -(out$fvec^2-off)
    pars <- opt.par0 <- opt.par <- as.list(coef(out))
    tmp <- CircularSig(par=opt.par,df = df)
    yp.full <- tmp$v
    nams <- c('gamma','beta')
    if(NI>0) nams <- c(nams,paste0('d',1:NI))
#if(type=='period' | (GP & noise.only)) nams <- c('A','B',nams)
    if(type=='period') nams <- c('A','B',nams)
    names(tmp$par) <- nams
    opt.par <- c(unlist(tmp$par),opt.par)
####model prediction and chi2
    res <- as.numeric(y-yp.full)
####
    pars <- opt.par
    yp.noise <- tmp$ytrend+tmp$yred
    ysig <- tmp$ysig
    res.sig <- as.numeric(y-ysig)
    return(list(df=df,res=res,res.sig=res.sig,noise=yp.noise,logL=logL,lnls=lnls,par=opt.par,par0=opt.par0,par.low=par.low,par.up=par.up,ysig=ysig,yfull=yp.full))
}

par.integral <- function(data,Indices=NULL,sj,m,d,type='noise',logtau=NULL,omega=NULL,Nma=0, Nar=0){
    if(is.null(Indices)){
        dI <- 0
        NI <- 0
    }else{
        NI <- ncol(Indices)
        dI <- partI(d,Indices)
    }
    t <- data[,1]
    y <- data[,2]
    dy <- data[,3]
    ind <- which(dy==0)
    dy[ind] <- 1e-6
    W <- sum(1/(dy^2+sj^2))#new weight sum is 1.
    w <- 1/(dy^2+sj^2)/W#normalized weighting
    if(Nma>0){
        Is <- c()
        vs <- c()
        ts <- c()
        trep <- c()
        for(k in 1:Nma){
            if(NI>1){
                Is <- rbind(Is,c(rep(0,k),partI(d,Indices[-(length(t)+1-(1:k)),]) ))
            }else{
                Is <- 0
            }
            ts <- rbind(ts,c(rep(0,k),t[-(length(t)+1-(1:k))]))#shift forward by k points
            trep <- rbind(trep,c(rep(0,k),t[-(1:k)]))#unshifted but with the previous points chopped
            vs <- rbind(vs,c(rep(0,k),y[-(length(t)+1-(1:k))]))
        }
        c <- m*exp(-abs(trep-ts)/exp(logtau))
    }else{
        c <- 0
        ts <- 0
        vs <- 0
        Is <- 0
    }
    if(Nma>0){
        if(is.matrix(Is) | is.data.frame(Is)){
            yp <- y-colSums(c*vs)-dI-colSums(c*Is)
        }else if(length(Is)>1){
            yp <- y-colSums(c*vs)-dI-c*Is
        }else{
            yp <- y-colSums(c*vs)-dI
        }
        tp <- t+colSums(c*ts)
        wp <- 1+colSums(c)
    }else{
        yp <- y-dI
        tp <- t
        wp <- 1
    }
    YWp <- sum(w*wp*yp)
    Wp <- sum(w*wp)
    WWp <- sum(w*wp^2)
    WTp <- sum(w*wp*tp)
    TTp <- sum(w*tp^2)
    WTp <- sum(w*wp*tp)
    YTp <- sum(w*yp*tp)
    YWp <- sum(w*yp*wp)
    YYp <- sum(w*yp^2)
    logL0 <- log(2*pi/sqrt(WWp*TTp-WTp^2))-log(W)+W*(((YTp*WWp-WTp*YWp)^2/(WWp*TTp-WTp^2)+YWp^2-YYp*WWp)/(2*WWp))
    logLp <- logL <- logL0
    if(type=='period'){
###determine phase to eliminate CS
        if(Nma>0){
            cc <- colSums(c*cos(omega*ts))
            ss <- colSums(c*sin(omega*ts))
        }else{
            ss <- cc <- 0
        }
        CCpp <- sum(w*(cos(omega*t)-cc)^2)
        SSpp <- sum(w*(sin(omega*t)-ss)^2)
        CSpp <- sum(w*(sin(omega*t)-ss)*(cos(omega*t)-cc))
        phi <- 0.5*atan(2*CSpp/(CCpp-SSpp))
        phip <- 0.5*atan(sum(w*sin(2*omega*t))/sum(w*cos(2*omega*t)))
###define notations using the determined phase
        if(Nma>0){
            cc <- colSums(c*cos(omega*ts-phi))
            ss <- colSums(c*sin(omega*ts-phi))
        }else{
            ss <- cc <- 0
        }
        cp <- cos(omega*t-phi)-cc
        sp <- sin(omega*t-phi)-ss
        Cp <- sum(w*cp)
        Sp <- sum(w*sp)
        CCp <- sum(w*cp^2)
        SSp <- sum(w*sp^2)
        CSp <- sum(w*sp*cp)
        CTp <- sum(w*cp*tp)
        STp <- sum(w*sp*tp)
        CWp <- sum(w*cp*wp)
        SWp <- sum(w*sp*wp)
        YCp <- sum(w*yp*cp)
        YSp <- sum(w*yp*sp)
        U <- WWp-CWp^2/CCp-SWp^2/SSp
        R <- CWp*CTp/CCp+SWp*STp/SSp-WTp
        Q <- YWp-CWp*YCp/CCp-SWp*YSp/SSp
        V <- CCp*SSp*TTp*U-SSp*CTp^2*U-CCp*STp^2*U-CCp*SSp*R^2
        if(V<=0){
            cat('V=',V,'is not positive for f=',omega/(2*pi),'!\n')
#            V <- V+1e-30
            cat('Q=',Q,'\n')
            cat('U=',U,'\n')
            cat('R=',R,'\n')
        }
        X <- SSp*YCp^2*U+CCp*YSp^2*U-YYp*CCp*SSp*U+CCp*SSp*Q^2
        logL <- log((2*pi)^2/sqrt(abs(V)))-2*log(W)+W*((X+(CCp*SSp*YTp*U-YCp*CTp*SSp*U-YSp*STp*CCp*U+CCp*SSp*Q*R)^2/V)/(2*CCp*SSp*U))
    }
    return(list(logL0=logL0,logL=logL,logLp=logLp))
}

#notations depending on frequency omega
local.notation <- function(t,y,dy,Indices,NI,omega,phi){
    W <- sum(1/dy^2)
    w <- 1/dy^2/W
    C <- sum(w*cos(omega*t-phi))
    S <- sum(w*sin(omega*t-phi))
    YC <- sum(w*y*cos(omega*t-phi))
    YS <- sum(w*y*sin(omega*t-phi))
    CC <- sum(w*cos(omega*t-phi)^2)
    SS <- sum(w*sin(omega*t-phi)^2)
    CS <- sum(w*sin(omega*t-phi)*cos(omega*t-phi))
    ST <- sum(w*sin(omega*t-phi)*t)
    CT <- sum(w*cos(omega*t-phi)*t)
    if(NI>0){
        CI <- (w*cos(omega*t-phi))%*%Indices[,1:NI]
        SI <- (w*sin(omega*t-phi))%*%Indices[,1:NI]
    }else{
        CI <- SI <- rep(0,NI)
    }
#    vars <- unique(c('omega','phi','C','S','CC','SS','CS','CT','ST','SI','CI','YC','YS'))
    vars <- unique(c('phi','C','S','CC','SS','CS','CT','ST','SI','CI','YC','YS'))
    pars <- list()
    for(k in 1:length(vars)){
        if(exists(vars[k])){
            pars[[vars[k]]] <- eval(parse(text = vars[k]))
        }
    }
    return(pars)
}

#####definition of notations
global.notation <- function(t,y,dy,Nma,Nar,Indices,GP,gp.par){
#######notations
    data <- cbind(t,y,dy)
    W <- sum(1/dy^2)
    w <- 1/dy^2/W
    T <- sum(w*t)
    Y <- sum(w*y)
    Inds <- NI <- 0
    TI <- YI <- I <- 0
    if(!is.null(Indices)){
        if(!any(is.na(Indices))){
            NI <- ncol(Indices)
            Inds <- 1:NI
            I <- w%*%Indices
            YI <- (w*y)%*%Indices
            TI <- (w*t)%*%Indices
        }
    }
    YY <- sum(w*y^2)
    YT <- sum(w*y*t)
    TT <- sum(w*t^2)
    II <- array(data=NA,dim=c(NI,NI))
    if(NI>0){
        for(i in 1:NI){
                II[i,] <- (w*Indices[,i])%*%Indices
        }
    }
#####parameter boundaries
    gamma.min <- min(y)
    gamma.max <- max(y)
    gamma.ini <- (gamma.min+gamma.max)/2
    if(NI>0){
        dmin <- -2*(max(y)-min(y))/(max(Indices)-min(Indices))
        dmax <- 2*(max(y)-min(y))/(max(Indices)-min(Indices))
        dini <- rep((dmin+dmax)/2,1)
    }else{
        dini <- dmin <- dmax <- 0#rep(0,NI)
    }
    beta.min <- -(max(y)-min(y))/(max(t)-min(t))
    beta.max <- (max(y)-min(y))/(max(t)-min(t))
    beta.ini <- (beta.min+beta.max)/2
    trend <- TRUE
    if(!trend){
        beta.min <- -1e-6
        beta.max <- 1e-6
        beta.ini <- 0
    }
    sigmaGP.min <- sj.min <- 0
    sj.max <- 10*sd(y)
    sigmaGP.max <- max(1e4,100*sd(y))
    sigmaGP.ini <- sj.ini <- max(sj.min,sj.max/100)
    lmin <- mmin <- -1
    lmax <- mmax <- 1
    rr <- 1
    lini <- mini <- rep((mmin+rr*mmax)/(1+rr),1)
    logProt.min <- logtauGP.min <- logtauAR.min <- logtau.min <- min(-10,log(max(min(diff(t)),1e-3)))
    logProt.max <- logtauGP.max <- logtauAR.max <- logtau.max <- max(20,log(1e4*(max(t)-min(t))))
    logProt.ini <- logtauGP.ini <- logtauAR.ini <- logtau.ini <- (logtau.min+rr*logtau.max)/(1+rr)
    gp.min <- c(sigmaGP=sigmaGP.min,logProt=logProt.min,logtauGP=logtauGP.min)
    gp.max <- c(sigmaGP=sigmaGP.max,logProt=logProt.max,logtauGP=logtauGP.max)
    gp.ini <- c(sigmaGP=sigmaGP.ini,logProt=logProt.ini,logtauGP=logtauGP.ini)
    Amin <- Bmin <- -2*(max(y)-min(y))
    Amax <- Bmax <- 2*(max(y)-min(y))
    Aini <- Bini <- (Amin+Amax)/2
###start
    start <- list()
    par.low <- par.up <- c()
    start$sj <- sj.ini
    par.low <- c(par.low,sj=sj.min)
    par.up <- c(par.up,sj=sj.max)
    if(GP){
        if(is.na(gp.par[1])){
            start$sigmaGP <- sigmaGP.ini
            par.low <- c(par.low,sigmaGP=sigmaGP.min)
            par.up <- c(par.up,sigmaGP=sigmaGP.max)
        }
        if(is.na(gp.par[2])){
            start$logProt <- logProt.ini
            par.low <- c(par.low,logProt=logProt.min)
            par.up <- c(par.up,logProt=logProt.max)
        }
        if(is.na(gp.par[3])){
            start$logtauGP <- logtauGP.ini
            par.low <- c(par.low,logtauGP=logtauGP.min)
            par.up <- c(par.up,logtauGP=logtauGP.max)
        }
    }
    if(Nma>0 ){
        for(k in 1:Nma){
            start[[paste0('m',k)]] <- mini
            par.low <- c(par.low,mmin)
            par.up <- c(par.up,mmax)
        }
        start <- c(start,logtau=logtau.ini)
        par.low <- c(par.low,logtau=logtau.min)
        par.up <- c(par.up,logtau=logtau.max)
    }
    if(Nar>0){
        for(k in 1:Nar){
            start[[paste0('l',k)]] <- lini
            par.low <- c(par.low,lmin)
            par.up <- c(par.up,lmax)
        }
        start <- c(start,logtauAR=logtauAR.ini)
        par.low <- c(par.low,logtauAR=logtauAR.min)
        par.up <- c(par.up,logtauAR=logtauAR.max)
    }
    names(par.up) <- names(par.low) <- names(start)
    start0 <- start
    vars <- unique(c('err2','II','T','TT','TI','II','w','W','I','Y','YT','YI','data','Indices','logtau','phi','m','xs','Amin','Bmin','Amax','Bmax','logtau.min','logtau.max','lmin','lmax','lini','mmin','mmax','beta.min','beta.max','dmin','dmax','gamma.max','gamma.min','Aini','Bini','gamma.ini','beta.ini','dini','mini','logtau.ini','logtauAR.ini','logtauAR.min','logtauAR.max','sj.ini','sj.max','sj.min','logProt.ini','logProt.min','logProt.max','logtauGP.ini','logtauGP.min','logtauGP.max','sigmaGP.ini','sigmaGP.min','sigmaGP.max','GP','gp.par','gp.min','gp.max','gp.ini','start','par.low','par.up','NI','Nma','Nar','start0'))
    pars <- list()
    for(k in 1:length(vars)){
        if(exists(vars[k])){
            pars[[vars[k]]] <- eval(parse(text = vars[k]))
        }
    }
    return(pars)
}

multiset_lag_terms <- function(y, set.id, Nar=0, Nma=0, residuals=NULL, t=NULL, tau=NULL, tauAR=NULL){
###tau/tauAR are the moving-average/auto-regressive time scales of the exponential
###kernel exp(-|dt|/tau) applied to each lag term. They default to NULL, which
###reproduces the plain lagged terms used by the multi-set signal periodogram.
    set.id <- factor(set.id)
    levels.id <- levels(set.id)
    normalize.order <- function(order){
        order <- as.integer(order)
        if(length(order)==0){
            order <- 0
        }
        if(length(order)<length(levels.id)){
            order <- rep(order,length.out=length(levels.id))
        }
        order[is.na(order)] <- 0
        return(order)
    }
    Nar <- normalize.order(Nar)
    Nma <- normalize.order(Nma)
    terms <- c()
    add.lags <- function(source, prefix, orders, terms, timescale){
        decay <- !is.null(t) && !is.null(timescale) && is.finite(timescale) && timescale>0
        for(j in 1:length(levels.id)){
            rows <- which(set.id==levels.id[j])
            if(orders[j]>0 && length(rows)>1){
                for(lag in 1:orders[j]){
                    col <- rep(0,length(y))
                    if(length(rows)>lag){
                        target <- rows[(lag+1):length(rows)]
                        origin <- rows[1:(length(rows)-lag)]
                        vals <- source[origin]
                        if(decay){
                            vals <- vals*exp(-abs(t[target]-t[origin])/timescale)
                        }
                        col[target] <- vals
                    }
                    terms <- cbind(terms,col)
                    colnames(terms)[ncol(terms)] <- paste0(prefix,'_',make.names(levels.id[j]),'_lag',lag)
                }
            }
        }
        return(terms)
    }
    terms <- add.lags(y,'ar',Nar,terms,tauAR)
    if(is.null(residuals)){
        residuals <- rep(0,length(y))
    }
    terms <- add.lags(residuals,'ma',Nma,terms,tau)
    return(terms)
}

multiset_design <- function(t, set.id, omega=NULL, trend=TRUE, lag.terms=NULL){
    set.id <- factor(set.id)
    levs <- levels(set.id)
###built directly rather than with model.matrix so that a single data set works too
    offsets <- matrix(unlist(lapply(levs,function(l) as.numeric(set.id==l))),
                      nrow=length(set.id))
    colnames(offsets) <- paste0('gamma_', make.names(levs))
    design <- offsets
    if(trend){
        design <- cbind(design, beta=t)
    }
    if(!is.null(lag.terms)){
        design <- cbind(design, lag.terms)
    }
    if(!is.null(omega)){
        design <- cbind(A=cos(omega*t), B=sin(omega*t), design)
    }
    return(design)
}

multiset_wls <- function(t, y, dy, set.id, omega=NULL, trend=TRUE, Nar=0, Nma=0, residuals=NULL,
                         tau=NULL, tauAR=NULL, signal=NULL){
    lag.terms <- multiset_lag_terms(y=y,set.id=set.id,Nar=Nar,Nma=Nma,residuals=residuals,
                                    t=t,tau=tau,tauAR=tauAR)
    X <- multiset_design(t=t,set.id=set.id,omega=omega,trend=trend,lag.terms=lag.terms)
    if(!is.null(signal)){
        X <- cbind(signal,X)
    }
    w <- 1/dy^2
    fit <- lm.wfit(x=X,y=y,w=w)
    coef <- fit$coefficients
    coef[is.na(coef)] <- 0
    yfit <- as.numeric(X%*%coef)
    res <- y-yfit
    logL <- sum(-0.5*(res^2/dy^2+log(2*pi*dy^2)))
    return(list(coef=coef,yfit=yfit,res=res,logL=logL,design=X))
}

multiset_periodogram <- function(t, y, dy, set.id, Nma=0, Nar=0, ofac=1, fmax=NULL, fmin=NA,
                                 tspan=NULL, sampling='combined', section=1,
                                 trend=TRUE, max.iter=3){
    unit <- 1
    t <- (t-min(t))/unit
    set.id <- factor(set.id)
    data <- cbind(t,y,dy)
    if(is.null(tspan)){
        tspan <- max(t)-min(t)
    }
    if(is.na(fmin)){
        fmin <- 1/(tspan*ofac)
    }
    fnyq <- 0.5*length(y)/tspan
    if(is.null(fmax)){
        fmax <- fnyq
    }
    f <- fsample(fmin,fmax,sampling,section,ofac,unit)
    omegas <- 2*pi*f
    base <- multiset_wls(t=t,y=y,dy=dy,set.id=set.id,trend=trend,Nar=Nar,Nma=0)
    if(any(as.integer(Nma)>0)){
        for(iter in 1:max.iter){
            base <- multiset_wls(t=t,y=y,dy=dy,set.id=set.id,trend=trend,Nar=Nar,Nma=Nma,residuals=base$res)
        }
    }
    logLs <- rep(NA,length(omegas))
    opt.pars <- c()
    yfits <- matrix(NA,nrow=length(t),ncol=length(omegas))
    residuals <- matrix(NA,nrow=length(t),ncol=length(omegas))
    for(kk in 1:length(omegas)){
        fit <- multiset_wls(t=t,y=y,dy=dy,set.id=set.id,omega=omegas[kk],trend=trend,Nar=Nar,Nma=Nma,residuals=base$res)
        if(any(as.integer(Nma)>0)){
            for(iter in 1:max.iter){
                fit <- multiset_wls(t=t,y=y,dy=dy,set.id=set.id,omega=omegas[kk],trend=trend,Nar=Nar,Nma=Nma,residuals=fit$res)
            }
        }
        logLs[kk] <- fit$logL
        opt.pars <- rbind(opt.pars,fit$coef)
        yfits[,kk] <- fit$yfit
        residuals[,kk] <- fit$res
    }
    Ndata <- length(t)
    Nextra <- 2
    logBF <- logLs-base$logL-Nextra/2*log(Ndata)
    inds <- sort(logBF,decreasing=TRUE,index.return=TRUE)$ix
    P <- unit/f
    Popt <- P[inds[1]]
    opt.par <- opt.pars[inds[1],]
    opt.par.top <- opt.pars[inds[1:min(10,length(inds))],,drop=FALSE]
    omega.opt <- 2*pi/Popt
    ysig <- as.numeric(opt.par['A']*cos(omega.opt*t)+opt.par['B']*sin(omega.opt*t))
    yfull <- yfits[,inds[1]]
    res <- residuals[,inds[1]]
    ps <- P[inds[1:min(5,length(inds))]]
    power.opt <- logBF[inds[1:min(5,length(inds))]]
    df <- list(data=data,set.id=set.id,multi_set=TRUE,Nma=Nma,Nar=Nar,NI=0,type='period')
    return(list(data=data,logBF=logBF,logLs=logLs,llmax=max(logLs),lnbfs=matrix(NA,nrow=length(P),ncol=Ndata),
                P=P,Popt=Popt,Popts=P[inds[1:min(10,length(inds))]],logBF.opt=logBF[inds[1:min(10,length(inds))]],
                par.opt=opt.par,opt.par=opt.par.top,res=res,res.nst=res,res.n=res,res.s=y-ysig,
                res.st=res,res.nt=res,ysig=ysig,yfull=yfull,base.fit=base,pars=opt.pars,df=df,
                power=logBF,ps=ps,power.opt=power.opt,sig.level=c(-Nextra/2*log(Ndata),0,log(150)),
                ParLow=rep(NA,length(opt.par)),ParUp=rep(NA,length(opt.par)),LogLike0=base$logL,
                multi_set=TRUE,set.id=set.id,trend=trend,Nma=Nma,Nar=Nar))
}

BFP.multiset <- function(t, y, dy, set.id, Nma=0, Nar=0, ofac=1, fmax=NULL, fmin=NA,
                         tspan=NULL, sampling='combined', section=1, progress=FALSE,
                         noise.only=FALSE){
    if(noise.only & !any(as.integer(Nma)>0) & !any(as.integer(Nar)>0)){
        warning('A purely stochastic multi-set fit needs at least one AR or MA component; fitting the periodic model instead.')
        noise.only <- FALSE
    }
    if(noise.only){
        return(multiset_noise_periodogram(t=t,y=y,dy=dy,set.id=set.id,Nma=Nma,Nar=Nar,ofac=ofac,
                                          fmax=fmax,fmin=fmin,tspan=tspan,sampling=sampling,
                                          section=section,trend=TRUE))
    }
    multiset_periodogram(t=t,y=y,dy=dy,set.id=set.id,Nma=Nma,Nar=Nar,ofac=ofac,fmax=fmax,fmin=fmin,
                         tspan=tspan,sampling=sampling,section=section,trend=TRUE)
}

MLP.multiset <- function(t, y, dy, set.id, Nma=0, Nar=0, ofac=1, fmax=NULL, fmin=NULL,
                         tspan=NULL, sampling='combined', section=1){
    if(is.null(fmin)){
        fmin <- NA
    }
    out <- multiset_periodogram(t=t,y=y,dy=dy,set.id=set.id,Nma=Nma,Nar=Nar,ofac=ofac,fmax=fmax,fmin=fmin,
                                tspan=tspan,sampling=sampling,section=section,trend=TRUE)
    out$power <- out$power-max(out$power)
    out$logBF <- out$power
    return(out)
}

###########################################################################
####Keplerian and purely stochastic fitting of multi-data-set time series
###########################################################################
multiset_kepler_basis <- function(t, P, e, Mo){
###Basis of a Keplerian signal in the form for which the model is linear:
###K*(cos(omega+nu)+e*cos(omega)) = C1*(cos(nu)+e)+C2*sin(nu)
###with C1=K*cos(omega) and C2=-K*sin(omega). The nonlinear parameters are
###therefore only P, e and Mo, which keeps the multi-set optimization robust.
    e <- min(max(e,0),0.99)
    M <- (Mo+2*pi*t/P)%%(2*pi)
    E <- kep.mt2(M,e)
    nu <- 2*atan(sqrt((1+e)/(1-e))*tan(E/2))
    return(cbind(kc=cos(nu)+e,ks=sin(nu)))
}

multiset_kepler_amp <- function(C1, C2){
###Convert the linear Keplerian coefficients back into K and omega
    K <- sqrt(C1^2+C2^2)
    if(K==0){
        omega <- 0
    }else{
        omega <- as.numeric(xy2phi(C1,-C2))
    }
    return(c(K=as.numeric(K),omega=omega))
}

multiset_kepler_curve <- function(par, t){
###Keplerian radial velocity of a multi-set fit at arbitrary times
    P <- par[['P1']]
    K <- par[['K1']]
    e <- par[['e1']]
    omega <- par[['omega1']]
    Mo <- par[['Mo1']]
    M <- (Mo+2*pi*t/P)%%(2*pi)
    E <- kep.mt2(M,e)
    nu <- 2*atan(sqrt((1+e)/(1-e))*tan(E/2))
    return(as.numeric(K*(cos(omega+nu)+e*cos(omega))))
}

multiset_kepler_wls <- function(t, y, dy, set.id, P, e, Mo, trend=TRUE, Nar=0, Nma=0,
                                residuals=NULL, max.iter=3){
###Weighted least squares for a shared Keplerian signal, one offset per data set,
###a shared trend and per-data-set AR/MA lag terms
    kep <- multiset_kepler_basis(t=t,P=P,e=e,Mo=Mo)
    fit <- multiset_wls(t=t,y=y,dy=dy,set.id=set.id,omega=NULL,trend=trend,Nar=Nar,Nma=Nma,
                        residuals=residuals,signal=kep)
    if(any(as.integer(Nma)>0)){
        for(iter in 1:max.iter){
            fit <- multiset_wls(t=t,y=y,dy=dy,set.id=set.id,omega=NULL,trend=trend,Nar=Nar,Nma=Nma,
                                residuals=fit$res,signal=kep)
        }
    }
    fit$ysig <- as.numeric(kep%*%fit$coef[c('kc','ks')])
    fit$kep <- kep
    return(fit)
}

KeplerFit.multiset <- function(per, frac=0.2, emax=0.95, max.iter=3){
####################################################
## Keplerian fitting of a multi-data-set periodogram
## Input:
##   per - output of BFP.multiset or MLP.multiset
##
## Output:
##   ParKep - shared Keplerian parameters plus the per-data-set nuisance parameters
####################################################
###sigfit() replaces per$data with the raw input table, so the fit has to start
###from the data the periodogram itself used
    if(!is.null(per$df) && !is.null(per$df$data)){
        data <- per$df$data
    }else{
        data <- per$data
    }
    t <- as.numeric(data[,1])
    t <- t-min(t)
    y <- as.numeric(data[,2])
    dy <- as.numeric(data[,3])
    set.id <- factor(per$set.id)
    trend <- if(is.null(per$trend)) TRUE else isTRUE(per$trend)
    Nar <- if(is.null(per$Nar)) 0 else per$Nar
    Nma <- if(is.null(per$Nma)) 0 else per$Nma
    Popt <- as.numeric(per$Popt[1])
    par.opt <- unlist(per$par.opt)
    A <- if(any(names(par.opt)=='A')) as.numeric(par.opt['A']) else 0
    B <- if(any(names(par.opt)=='B')) as.numeric(par.opt['B']) else 0
    if(A==0 & B==0){
        phi <- 0
    }else{
        phi <- as.numeric(xy2phi(A,B))
    }
    Plow <- max((1-frac)*Popt,1e-6)
    Pup <- (1+frac)*Popt
    fit.at <- function(v){
        multiset_kepler_wls(t=t,y=y,dy=dy,set.id=set.id,P=v[1],e=v[2],Mo=v[3],
                            trend=trend,Nar=Nar,Nma=Nma,max.iter=max.iter)
    }
    obj <- function(v){
        out <- try(fit.at(v),TRUE)
        if(class(out)[1]=='try-error') return(1e10)
        if(!is.finite(out$logL)) return(1e10)
        return(-out$logL)
    }
###the circular solution fixes the phase, so only e and Mo need to be scanned
    e.start <- c(0,0.1,0.3,0.5,0.7)
    Mo.start <- ((-phi)+seq(0,2*pi,length.out=5)[-5])%%(2*pi)
    starts <- as.matrix(expand.grid(P=Popt,e=e.start,Mo=Mo.start))
    lower <- c(Plow,0,0)
    upper <- c(Pup,emax,2*pi)
    lls <- c()
    pars <- c()
    best <- NULL
    best.ll <- -Inf
    for(k in 1:nrow(starts)){
        v0 <- pmin(pmax(as.numeric(starts[k,]),lower+1e-8),upper-1e-8)
        out <- try(optim(par=v0,fn=obj,method='L-BFGS-B',lower=lower,upper=upper,
                         control=list(maxit=500)),TRUE)
        if(class(out)[1]=='try-error') next
        ll <- -out$value
        if(!is.finite(ll) | ll<=-1e9) next
        lls <- c(lls,ll)
        pars <- rbind(pars,out$par)
        if(ll>best.ll){
            best.ll <- ll
            best <- out$par
        }
    }
    if(is.null(best)){
        warning('Keplerian multi-set optimization failed; falling back to the circular solution.')
        best <- c(Popt,0,(-phi)%%(2*pi))
    }
    fit <- fit.at(best)
    coef <- fit$coef
    amp <- multiset_kepler_amp(as.numeric(coef['kc']),as.numeric(coef['ks']))
    ParKep <- list(P1=as.numeric(best[1]),K1=as.numeric(amp['K']),e1=as.numeric(best[2]),
                   omega1=as.numeric(amp['omega']),Mo1=as.numeric(best[3]))
    extra <- coef[!(names(coef)%in%c('kc','ks'))]
    ParKep <- c(ParKep,as.list(extra))
    logL0 <- if(is.null(per$LogLike0)) NA else as.numeric(per$LogLike0)
###five extra free parameters relative to the no-signal multi-set model
    Nextra <- 5
    lnBF <- fit$logL-logL0-Nextra/2*log(length(y))
    return(list(ParKep=ParKep,ysig=fit$ysig,yfull=fit$yfit,res=fit$res,logL=fit$logL,
                lnBF=lnBF,ll0=logL0,lls=lls,pars=pars,coef=coef,Popt=ParKep$P1,
                data=data,set.id=set.id,trend=trend,Nar=Nar,Nma=Nma))
}

multiset_noise_periodogram <- function(t, y, dy, set.id, Nma=0, Nar=0, ofac=1, fmax=NULL, fmin=NA,
                                       tspan=NULL, sampling='combined', section=1,
                                       trend=TRUE, max.iter=3){
###Purely stochastic multi-set fit: no periodic signal is included and the
###periodogram scans the time scale of the exponential AR/MA kernel instead of
###the period of a signal. The baseline is the white-noise multi-set model.
    unit <- 1
    t <- (t-min(t))/unit
    set.id <- factor(set.id)
    data <- cbind(t,y,dy)
    if(is.null(tspan)){
        tspan <- max(t)-min(t)
    }
    if(is.na(fmin)){
        fmin <- 1/(tspan*ofac)
    }
    fnyq <- 0.5*length(y)/tspan
    if(is.null(fmax)){
        fmax <- fnyq
    }
    f <- fsample(fmin,fmax,sampling,section,ofac,unit)
    taus <- unit/f
    base <- multiset_wls(t=t,y=y,dy=dy,set.id=set.id,trend=trend,Nar=0,Nma=0)
    logLs <- rep(NA,length(taus))
    opt.pars <- c()
    yfits <- matrix(NA,nrow=length(t),ncol=length(taus))
    residuals <- matrix(NA,nrow=length(t),ncol=length(taus))
    for(kk in 1:length(taus)){
        fit <- multiset_wls(t=t,y=y,dy=dy,set.id=set.id,trend=trend,Nar=Nar,Nma=Nma,
                            residuals=base$res,tau=taus[kk],tauAR=taus[kk])
        if(any(as.integer(Nma)>0)){
            for(iter in 1:max.iter){
                fit <- multiset_wls(t=t,y=y,dy=dy,set.id=set.id,trend=trend,Nar=Nar,Nma=Nma,
                                    residuals=fit$res,tau=taus[kk],tauAR=taus[kk])
            }
        }
        logLs[kk] <- fit$logL
        opt.pars <- rbind(opt.pars,c(fit$coef,logtau=log(taus[kk]),logtauAR=log(taus[kk])))
        yfits[,kk] <- fit$yfit
        residuals[,kk] <- fit$res
    }
    Ndata <- length(t)
    lag.terms <- multiset_lag_terms(y=y,set.id=set.id,Nar=Nar,Nma=Nma,
                                    residuals=rep(0,length(y)),t=t,tau=1,tauAR=1)
    Nlag <- if(is.null(lag.terms)) 0 else ncol(lag.terms)
###one lag coefficient per term plus the shared MA and AR time scales
    Nextra <- Nlag+sum(any(as.integer(Nma)>0),any(as.integer(Nar)>0))
    logBF <- logLs-base$logL-Nextra/2*log(Ndata)
    inds <- sort(logBF,decreasing=TRUE,index.return=TRUE)$ix
    Popt <- taus[inds[1]]
    opt.par <- opt.pars[inds[1],]
    opt.par.top <- opt.pars[inds[1:min(10,length(inds))],,drop=FALSE]
    yfull <- yfits[,inds[1]]
    res <- residuals[,inds[1]]
###the red-noise component is the full prediction minus the offsets and trend
    white <- multiset_design(t=t,set.id=set.id,omega=NULL,trend=trend,lag.terms=NULL)
    yred <- as.numeric(yfull-white%*%opt.par[colnames(white)])
    ysig <- rep(0,length(t))
    ps <- taus[inds[1:min(5,length(inds))]]
    power.opt <- logBF[inds[1:min(5,length(inds))]]
    df <- list(data=data,set.id=set.id,multi_set=TRUE,Nma=Nma,Nar=Nar,NI=0,type='noise')
    return(list(data=data,logBF=logBF,logLs=logLs,llmax=max(logLs),lnbfs=matrix(NA,nrow=length(taus),ncol=Ndata),
                P=taus,Popt=Popt,Popts=taus[inds[1:min(10,length(inds))]],logBF.opt=logBF[inds[1:min(10,length(inds))]],
                par.opt=opt.par,opt.par=opt.par.top,res=res,res.nst=res,res.n=res,res.s=res,
                res.st=res,res.nt=res,ysig=ysig,yred=yred,yfull=yfull,base.fit=base,pars=opt.pars,df=df,
                power=logBF,ps=ps,power.opt=power.opt,sig.level=c(-Nextra/2*log(Ndata),0,log(150)),
                ParLow=rep(NA,length(opt.par)),ParUp=rep(NA,length(opt.par)),LogLike0=base$logL,
                multi_set=TRUE,noise_only=TRUE,set.id=set.id,trend=trend,Nma=Nma,Nar=Nar))
}

####Marginalized likelihood periodogram
MLP <- function(t, y, dy, Nma=0, Nar=0,mar.type='part',sj=0,logtau=NULL,ofac=1,fmax=NULL,fmin=NULL,tspan=NULL,model.type='MA',opt.par=NULL,Indices=NULL,MLP.type='sub',sampling='combined',section=1, GP=FALSE,gp.par=NULL,noise.only=FALSE,Nsamp=5){
#    save(list=ls(all=TRUE),file='test.Robj')
    if(Nma==0 & Nar==0 & !GP) noise.only <- FALSE
    if(noise.only) quantify <- FALSE
    unit <- 1
    t <- (t-min(t))/unit#rescale time
    if(is.null(Indices)){
        NI <- 0
    }else{
        NI <- ncol(Indices)
        Indices <- as.matrix(Indices)
    }
    data <- cbind(t,y,dy)
    if(is.null(tspan)){
        tspan <- max(t)-min(t)
    }
    if(is.null(fmin)){
        fmin <- 1/(tspan*ofac)
    }
    fnyq <- 0.5*length(y)/tspan
    if(is.null(fmax)){
        fmax <- fnyq
    }
    f <- fsample(fmin,fmax,sampling,section,ofac,unit)
    nout <- length(f)
    t <- (t-min(t))/unit#
    Ndata <- length(t)
    omegas <- 2*pi*f
    phi <- 0
#######define notations and variables
    vars <- global.notation(t,y,dy,Nma,Nar,Indices,GP,gp.par)

############################################################
#####optimization; select the initial condition
############################################################
####fix the non-marginzed parameters at their optimal values
    if(is.null(opt.par) | MLP.type=='sub'){
#        tmp <- sopt(omega=NA,phi=NA,Nma=Nma,Nar=Nar,Indices=Indices,data=data,type='noise',par.low=vars$par.low,par.up=vars$par.up,start=vars$start,noise.only=noise.only,GP=GP,gp.par=gp.par,Nrep=Nsamp)
        tmp <- sopt(omega=NA,phi=NA,Nma=Nma,Nar=Nar,Indices=NULL,data=data,type='noise',par.low=vars$par.low,par.up=vars$par.up,start=vars$start,noise.only=noise.only,GP=GP,gp.par=gp.par,Nrep=Nsamp)
        opt.par <- tmp$par
        if(MLP.type=='sub'){
            y <- tmp$res
            data[,2] <- y
            vars$y <- y
            vars$data <- data
        }
    }
    if(Nma>0 & MLP.type=='assign'){
        ind <- grep('^m',names(opt.par))
        m <- unlist(lapply(ind,function(i) opt.par[[i]]))
        logtau <- opt.par$logtau
    }else{
        m <- 0
        logtau <- 1
    }
    if(NI>0 & MLP.type=='assign'){
        ind <- grep('^d',names(opt.par))
        d <- unlist(lapply(ind,function(i) opt.par[[i]]))
    }else{
        d <- rep(0,NI)
    }
    if(sj==0 & any(sj==names(opt.par))){
        sj <- opt.par$sj
    }
########################################
########marginalized posterior
########################################
    if(!exists('m')) m <- 0
    tmp <- par.integral(data=data,Indices=Indices,sj=sj,m=m,d=d,logtau=logtau,Nma=Nma,Nar=Nar,type='noise')
    logL0 <- tmp$logL0
    logBFmax <- 0
    p <- pn <- rep(NA,length(omegas))
    logBF.noise <- rep(NA,length(omegas))
    logBF <- rep(NA,length(omegas))
    logLp <- rep(NA,length(omegas))
    for(kk in 1:length(f)){
        omega <- omegas[kk]
        tmp <- par.integral(data,Indices,sj=sj,m=m,d=d,type='period',logtau=logtau,omega=omega,Nma=Nma,Nar=Nar)
        logL1 <- tmp$logL0#signal dependent noise model log likelihood
        logL <- tmp$logL#likelihood for the full model for f=f[k]
        logLp[kk] <- tmp$logLp
        logBF.noise[kk] <- logL1-logL0
        logBF[kk] <- logL-logL0
    }
    P <- unit/f
    ind <- which(logBF>(max(logBF)+log(0.01)))
    Popt <- (unit/f)[ind]
    omega.opt <- 2*pi*(f[which.max(logBF)]/unit)
####optimized parameter
    vars <- global.notation(t,y,dy,Nma=Nma,Nar=Nar,Indices=Indices,GP=GP,gp.par=gp.par)
    tmp <- sopt(omega=omega.opt,phi=0,Nma=Nma,Nar=Nar,Indices=Indices,data=data,type='period',par.low=vars$par.low,par.up=vars$par.up,start=vars$start,noise.only=FALSE,GP=FALSE,gp.par=rep(NA,3),Nrep=Nsamp)#
    opt.par <- tmp$par
    inds <- sort(logBF,decreasing=TRUE,index.return=TRUE)$ix
    ps <- P[inds[1:5]]
    power.opt <- logBF[inds[1:5]]
    return(list(data=data,P=unit/f, Popt=P[inds[1]],logp=power, ps=ps,power.opt=power.opt,logBF.opt=logBF[ind],logBF=logBF,power=logBF,logBF.noise=logBF.noise,Nma=Nma,Nar=Nar,NI=NI,par.opt=opt.par,res=tmp$res,sig.level=log(c(10,100,1000)),ParLow=vars$par.low,ParUp=vars$par.up,df=tmp$df))
}

#####BFP-based model inference/selection
model.infer.combined <- function(data,proxy=NULL,Nma.max=6,Nar.max=0,Nrep=5,GP=FALSE){
    Ndata <- nrow(data)
    data[,1] <- data[,1]-min(data[,1])
    t <- data[,1]
    y <- data[,2]
    dy <- data[,3]
    vars <- global.notation(t,y,dy,Indices=NULL,Nma=0,Nar=0,GP=FALSE,gp.par=rep(NA,3))
    tmp <- sopt(omega=NA,phi=NA,Nma=0,Nar=0,Indices=NULL,data=data,type='noise',par.low=vars$par.low,par.up=vars$par.up,start=vars$start,noise.only=FALSE,GP=FALSE,gp.par=rep(NA,3),Nrep=Nrep)#
    L0 <- tmp$logL
    lnLmaxs <- list(base=L0)
    lnBFs <- list()
    Inds.opt <- 0
    for(k in 1:Nma.max){
        lnBFs <- rep(NA,Nma.max)
        if(!is.null(proxy)){
            lnBFs <- rep(NA,Nma.max)
            lnLmax.proxy <- lnBF.proxy <- c()
            for(j in 1:ncol(proxy)){
                vars <- global.notation(t,y,dy,Indices=proxy[,1:j,drop=FALSE],Nma=0,Nar=0,GP=FALSE,gp.par=rep(NA,3))
                tmp <- sopt(omega=NA,phi=NA,Nma=0,Nar=0,Indices=proxy[,1:j,drop=FALSE],data=data,type='noise',par.low=vars$par.low,par.up=vars$par.up,start=vars$start,noise.only=FALSE,GP=FALSE,gp.par=rep(NA,3),Nrep=Nrep)#
                lnBF <- tmp$logL-L0-0.5*log(Ndata)
                lnBF.proxy <- c(lnBF.proxy,lnBF)
                lnLmax.proxy <- c(lnLmax.proxy,tmp$logL)
                if(lnBF<5){
                    break()
                }else{
                    L0 <- tmp$logL
                }
            }
            if(j>1){
                Inds.opt <- 1:(j-1)
            }
            indice <- proxy[,Inds.opt,drop=FALSE]
            lnBFs[['proxy']][['MA']] <- lnBF.proxy
            lnLmaxs[['proxy']][['']] <- lnLmax.proxy
        }else{
            indice <- NULL
        }
    }
    for(j in 1:Nma.max){
        lnLmax.MA <- lnBF.MA <- c()
        vars <- global.notation(t,y,dy,Indices=indice,Nma=j,Nar=0,GP=FALSE,gp.par=rep(NA,3))
          tmp <- sopt(omega=NA,phi=NA,Nma=j,Nar=0,Indices=indice,data=data,type='noise',par.low=vars$par.low,par.up=vars$par.up,start=vars$start,noise.only=FALSE,GP=FALSE,gp.par=rep(NA,3),Nrep=Nrep)
      if(j==1){
            lnBF <- tmp$logL-L0-1*log(Ndata)
        }else{
            lnBF <- tmp$logL-L0-0.5*log(Ndata)
        }
        lnBF.MA <- c(lnBF.MA,lnBF)
        lnLmax.MA <- c(lnLmax.MA,tmp$logL)
        if(lnBF<5){
            break()
        }else{
            L0 <- tmp$logL
        }
    }
    lnBFs[['MA']] <- lnBF.MA
    lnLmaxs[['MA']] <- lnLmax.MA
    Nma.opt <- j-1
    if(Nar.max>0){
        lnLmax.AR <- lnBF.AR <- c()
        for(j in 1:Nar.max){
            vars <- global.notation(t,y,dy,Indices=indice,Nma=Nma.opt,Nar=j,GP=FALSE,gp.par=rep(NA,3))
            tmp <- sopt(omega=NA,phi=NA,Nma=Nma.opt,Nar=0,Indices=indice,data=data,type='noise',par.low=vars$par.low,par.up=vars$par.up,start=vars$start,noise.only=FALSE,GP=FALSE,gp.par=rep(NA,3),Nrep=Nrep)
            if(j==1){
                lnBF <- tmp$logL-L0-1*log(Ndata)
            }else{
                lnBF <- tmp$logL-L0-0.5*log(Ndata)
            }
            lnBF.AR <- c(lnBF.AR,lnBF)
            lnLmax.AR <- c(lnLmax.AR,tmp$logL)
            if(lnBF<5){
                break()
            }else{
                L0 <- tmp$logL
            }
        }
        Nar.opt <- j-1
        lnBFs[['AR']] <- lnBF.AR
    }else{
        Nar.opt <- 0
    }
    GPf <- FALSE
    if(GP){
        vars <- global.notation(t,y,dy,Indices=indices,Nma=Nmar.opt,Nar=Nar.opt,GP=TRUE,gp.par=rep(NA,3))
        tmp <- sopt(omega=NA,phi=NA,Nma=0,Nar=0,Indices=NULL,data=data,type='noise',par.low=vars$par.low,par.up=vars$par.up,start=vars$start,noise.only=FALSE,GP=TRUE,gp.par=vars$gp.par,Nrep=Nrep)
        lnBF <- tmp$logL-L0-1.5*log(Ndata)
        if(lnBF>=5) GPf <- TRUE
        lnBFs[['GP']] <- lnBF
        lnLmaxs[['GP']] <- tmp$logL
    }
    cat('Best model: ARMA(p=',Nar.opt,',q=',Nma.opt,') ')
    if(any(Inds.opt!=0)){
        cat('+',Inds.opt,'proxy ')
    }
    if(GP){
        cat('+GP')
    }
    cat('\n')
###
    cat('Optimal noise model hyper parameter: Nar=',Nar.opt,'; Nma=',Nma.opt,'; Inds=',Inds.opt,';GP',GPf,'\n')
    return(list(Nar=Nar.opt,Nma=Nma.opt,Inds=Inds.opt,lnLmaxs=lnLmaxs,lnBFs=lnBFs,GP=GPf))
}

model.infer <- function(data,proxy=NULL,Nma.max=6,Nar.max=0,Nrep=5,GP=FALSE){
    Ndata <- nrow(data)
    data[,1] <- data[,1]-min(data[,1])
    t <- data[,1]
    y <- data[,2]
    dy <- data[,3]
    vars <- global.notation(t,y,dy,Indices=NULL,Nma=0,Nar=0,GP=FALSE,gp.par=rep(NA,3))
    tmp <- sopt(omega=NA,phi=NA,Nma=0,Nar=0,Indices=NULL,data=data,type='noise',par.low=vars$par.low,par.up=vars$par.up,start=vars$start,noise.only=FALSE,GP=FALSE,gp.par=rep(NA,3),Nrep=Nrep)#
    L0 <- tmp$logL
    lnLmaxs <- list(base=L0)
    lnBFs <- list()
    Inds.opt <- 0
    if(!is.null(proxy)){
        lnLmax.proxy <- lnBF.proxy <- c()
        for(j in 1:ncol(proxy)){
            vars <- global.notation(t,y,dy,Indices=proxy[,1:j,drop=FALSE],Nma=0,Nar=0,GP=FALSE,gp.par=rep(NA,3))
            tmp <- sopt(omega=NA,phi=NA,Nma=0,Nar=0,Indices=proxy[,1:j,drop=FALSE],data=data,type='noise',par.low=vars$par.low,par.up=vars$par.up,start=vars$start,noise.only=FALSE,GP=FALSE,gp.par=rep(NA,3),Nrep=Nrep)#
            lnBF <- tmp$logL-L0-0.5*log(Ndata)
            lnBF.proxy <- c(lnBF.proxy,lnBF)
            lnLmax.proxy <- c(lnLmax.proxy,tmp$logL)
            if(lnBF<5){
                break()
            }else{
                L0 <- tmp$logL
            }
        }
        if(j>1){
            Inds.opt <- 1:(j-1)
        }
        indice <- proxy[,Inds.opt,drop=FALSE]
        lnBFs[['proxy']] <- lnBF.proxy
        lnLmaxs[['proxy']] <- lnLmax.proxy
    }else{
        indice <- NULL
    }
    for(j in 1:Nma.max){
        lnLmax.MA <- lnBF.MA <- c()
        vars <- global.notation(t,y,dy,Indices=indice,Nma=j,Nar=0,GP=FALSE,gp.par=rep(NA,3))
          tmp <- sopt(omega=NA,phi=NA,Nma=j,Nar=0,Indices=indice,data=data,type='noise',par.low=vars$par.low,par.up=vars$par.up,start=vars$start,noise.only=FALSE,GP=FALSE,gp.par=rep(NA,3),Nrep=Nrep)
      if(j==1){
            lnBF <- tmp$logL-L0-1*log(Ndata)
        }else{
            lnBF <- tmp$logL-L0-0.5*log(Ndata)
        }
        lnBF.MA <- c(lnBF.MA,lnBF)
        lnLmax.MA <- c(lnLmax.MA,tmp$logL)
        if(lnBF<5){
            break()
        }else{
            L0 <- tmp$logL
        }
    }
    lnBFs[['MA']] <- lnBF.MA
    lnLmaxs[['MA']] <- lnLmax.MA
    Nma.opt <- j-1
    if(Nar.max>0){
        lnLmax.AR <- lnBF.AR <- c()
        for(j in 1:Nar.max){
            vars <- global.notation(t,y,dy,Indices=indice,Nma=Nma.opt,Nar=j,GP=FALSE,gp.par=rep(NA,3))
            tmp <- sopt(omega=NA,phi=NA,Nma=Nma.opt,Nar=0,Indices=indice,data=data,type='noise',par.low=vars$par.low,par.up=vars$par.up,start=vars$start,noise.only=FALSE,GP=FALSE,gp.par=rep(NA,3),Nrep=Nrep)
            if(j==1){
                lnBF <- tmp$logL-L0-1*log(Ndata)
            }else{
                lnBF <- tmp$logL-L0-0.5*log(Ndata)
            }
            lnBF.AR <- c(lnBF.AR,lnBF)
            lnLmax.AR <- c(lnLmax.AR,tmp$logL)
            if(lnBF<5){
                break()
            }else{
                L0 <- tmp$logL
            }
        }
        Nar.opt <- j-1
        lnBFs[['AR']] <- lnBF.AR
    }else{
        Nar.opt <- 0
    }
    GPf <- FALSE
    if(GP){
        vars <- global.notation(t,y,dy,Indices=indices,Nma=Nmar.opt,Nar=Nar.opt,GP=TRUE,gp.par=rep(NA,3))
        tmp <- sopt(omega=NA,phi=NA,Nma=0,Nar=0,Indices=NULL,data=data,type='noise',par.low=vars$par.low,par.up=vars$par.up,start=vars$start,noise.only=FALSE,GP=TRUE,gp.par=vars$gp.par,Nrep=Nrep)
        lnBF <- tmp$logL-L0-1.5*log(Ndata)
        if(lnBF>=5) GPf <- TRUE
        lnBFs[['GP']] <- lnBF
        lnLmaxs[['GP']] <- tmp$logL
    }
    cat('Best model: ARMA(p=',Nar.opt,',q=',Nma.opt,') ')
    if(any(Inds.opt!=0)){
        cat('+',Inds.opt,'proxy ')
    }
    if(GP){
        cat('+GP')
    }
    cat('\n')
###
    cat('Optimal noise model hyper parameter: Nar=',Nar.opt,'; Nma=',Nma.opt,'; Inds=',Inds.opt,';GP',GPf,'\n')
    return(list(Nar=Nar.opt,Nma=Nma.opt,Inds=Inds.opt,lnLmaxs=lnLmaxs,lnBFs=lnBFs,GP=GPf))
}

bfp.inf.combined <- function(data,Nmas=NULL,Nars=NULL,NI.inds=list(0),Nrep=5,GP=FALSE){
    Ndata <- nrow(data)
    data[,1] <- data[,1]-min(data[,1])
    t <- data[,1]
    y <- data[,2]
    dy <- data[,3]
    Indices <- NULL
    if(ncol(data)>3){
        Indices <- as.matrix(data[,4:ncol(data),drop=FALSE])
    }
    if(is.null(Nmas)){
        Nmas <- 0:2
    }
    if(is.null(Nars)){
        Nars <- 0:2
    }
    nis <- c()
    for(j in 1:length(NI.inds)){
            ni.ind <- NI.inds[[j]]
            nis <- c(nis,length(ni.ind[ni.ind!=0]))
    }
    Ndata <- nrow(data)
    logLmaxs <- Npars <- logBFs <- array(data=NA,dim=c(length(NI.inds),length(Nmas),length(Nars)),dimnames=list(paste0('NI',nis),paste0('Nma',Nmas),paste0('Nar',Nars)))
    ind.opt <- 1
    for(i in 1:length(Nars)){
        for(k in 1:length(Nmas)){
            for(j in 1:length(NI.inds)){
                if(!all(NI.inds[[j]]==0)){
                    ni <- length(NI.inds[[j]][NI.inds[[j]]!=0])
                    Inds <- NI.inds[[j]]
                    indices <- Indices[,Inds,drop=FALSE]
                }else{
                    ni <- 0
                    indices <- NULL
                    Inds <- 0
                }
                vars <- global.notation(t,y,dy,Indices=indices,Nma=Nmas[k],Nar=0,GP=FALSE,gp.par=rep(NA,3))
                tmp <- sopt(omega=NA,phi=NA,Nma=Nmas[k],Nar=0,Indices=indices,data=data,type='noise',par.low=vars$par.low,par.up=vars$par.up,start=vars$start,noise.only=FALSE,GP=FALSE,gp.par=vars$gp.par,Nrep=Nrep)
                Npar.ma <- Nmas[k]+1
                if(Nmas[k]==0){
                    Npar.ma <- 0
                }
                Npar.ar <- Nars[i]+1
                if(Nars[i]==0){
                    Npar.ar <- 0
                }
                Npar <- ni+Npar.ma+Npar.ar
                logLmaxs[j,k,i] <- tmp$logL
                logBFs[j,k,i] <- tmp$logL-logLmaxs[1,1,1]-0.5*Npar*log(Ndata)
                Npars[j,k,i] <- Npar
            }
        }
    }
###
    ind.opt <- c(1,1,1)
    for(i in 1:length(Nars)){
        for(k in 1:length(Nmas)){
            for(j in 1:length(NI.inds)){
#                lnBF <- logLmaxs-logLmaxs[j,k,i]-0.5*(Npars-Npars[j,k,i])*log(Ndata)
#                cat('i=',i,';k=',k,';j=',j,'lnBF',lnBF,'\n')
                if(!any((Npars>Npars[j,k,i]) & (logBFs-logBFs[j,k,i])>=5) & (logBFs[j,k,i]-logBFs[ind.opt[1],ind.opt[2],ind.opt[3]])>=5){
                    ind.opt <- c(j,k,i)
                }
            }
        }
    }
    Nma.opt=Nmas[ind.opt[2]]
    Nar.opt=Nars[ind.opt[3]]
    NI.opt=NI.inds[[ind.opt[1]]]
    NI.opt <- length(which(NI.opt!=0))
    if(length(length(NI.opt)>0)){
        proxy.opt <- Indices[,NI.opt,drop=FALSE]
    }else{
        proxy.opt <- NULL
    }
#    ind <- which(logBFs>5,arr.ind=TRUE)
    return(list(lnL=logLmaxs,lnBF=logBFs,Npars=Npars,Npar.opt=Npars[ind.opt[1],ind.opt[2],ind.opt[3]],ind.opt=ind.opt,nqp=c(NI.opt=NI.opt,Nma.opt=Nma.opt,Nar.opt=Nar.opt),proxy.opt=proxy.opt))
}
bfp.inf.progress <- function(data,Nmas=NULL,Nars=NULL,NI.inds=list(0),Nrep=5,GP=FALSE){
    Ndata <- nrow(data)
    data[,1] <- data[,1]-min(data[,1])
    t <- data[,1]
    y <- data[,2]
    dy <- data[,3]
    Indices <- NULL
    if(ncol(data)>3){
        Indices <- as.matrix(data[,4:ncol(data),drop=FALSE])
    }
    if(is.null(Nmas)){
        Nmas <- 0:2
    }
    if(is.null(Nars)){
        Nars <- 0:2
    }
    nis <- c()
    for(j in 1:length(NI.inds)){
            ni.ind <- NI.inds[[j]]
            nis <- c(nis,length(ni.ind[ni.ind!=0]))
    }
    Ndata <- nrow(data)
    logLmaxs <- Npars <- logBFs <- array(data=NA,dim=c(length(NI.inds),length(Nmas),length(Nars)),dimnames=list(paste0('NI',nis),paste0('Nma',Nmas),paste0('Nar',Nars)))
    ind.opt <- 1

    withProgress(message = 'Calculating log(BF) table\n', value = 0, {
    for(i in 1:length(Nars)){
        for(k in 1:length(Nmas)){
            for(j in 1:length(NI.inds)){
                incProgress(1/(length(Nmas)*length(NI.inds)*length(Nars)), detail = paste("ARMA(", Nars[i],',',Nmas[k],')+Proxy',paste(NI.inds[[j]],collapse=',')))
                if(!all(NI.inds[[j]]==0)){
                    ni <- length(NI.inds[[j]][NI.inds[[j]]!=0])
                    Inds <- NI.inds[[j]]
                    indices <- Indices[,Inds,drop=FALSE]
                }else{
                    ni <- 0
                    indices <- NULL
                    Inds <- 0
                }
                vars <- global.notation(t,y,dy,Indices=indices,Nma=Nmas[k],Nar=0,GP=FALSE,gp.par=rep(NA,3))
                tmp <- sopt(omega=NA,phi=NA,Nma=Nmas[k],Nar=0,Indices=indices,data=data,type='noise',par.low=vars$par.low,par.up=vars$par.up,start=vars$start,noise.only=FALSE,GP=FALSE,gp.par=vars$gp.par,Nrep=Nrep)
                Npar.ma <- Nmas[k]+1
                if(Nmas[k]==0){
                    Npar.ma <- 0
                }
                Npar.ar <- Nars[i]+1
                if(Nars[i]==0){
                    Npar.ar <- 0
                }
                Npar <- ni+Npar.ma+Npar.ar
                logLmaxs[j,k,i] <- tmp$logL
                logBFs[j,k,i] <- tmp$logL-logLmaxs[1,1,1]-0.5*Npar*log(Ndata)
                Npars[j,k,i] <- Npar
            }
        }
    }
    })
###
    ind.opt <- c(1,1,1)
    for(i in 1:length(Nars)){
        for(k in 1:length(Nmas)){
            for(j in 1:length(NI.inds)){
#                lnBF <- logLmaxs-logLmaxs[j,k,i]-0.5*(Npars-Npars[j,k,i])*log(Ndata)
#                cat('i=',i,';k=',k,';j=',j,'lnBF',lnBF,'\n')
                if(!any((Npars>Npars[j,k,i]) & (logBFs-logBFs[j,k,i])>=5) & (logBFs[j,k,i]-logBFs[ind.opt[1],ind.opt[2],ind.opt[3]])>=5){
                    ind.opt <- c(j,k,i)
                }
            }
        }
    }
    Nma.opt=Nmas[ind.opt[2]]
    Nar.opt=Nars[ind.opt[3]]
    NI.opt=NI.inds[[ind.opt[1]]]
    NI.opt <- length(which(NI.opt!=0))
    if(length(length(NI.opt)>0)){
        proxy.opt <- Indices[,NI.opt,drop=FALSE]
    }else{
        proxy.opt <- NULL
    }
#    ind <- which(logBFs>5,arr.ind=TRUE)
    return(list(lnL=logLmaxs,lnBF=logBFs,Npars=Npars,Npar.opt=Npars[ind.opt[1],ind.opt[2],ind.opt[3]],ind.opt=ind.opt,Nma.opt=Nma.opt,NI.opt=NI.opt,Nar.opt=Nar.opt,nqp=c(NI.opt=NI.opt,Nma.opt=Nma.opt,Nar.opt=Nar.opt),proxy.opt=proxy.opt))
}

bfp.inf.norm <- function(data,Nmas=NULL,Nars=NULL,NI.inds=NULL,Nrep=5,GP=FALSE){
    Ndata <- nrow(data)
    data[,1] <- data[,1]-min(data[,1])
    t <- data[,1]
    y <- data[,2]
    dy <- data[,3]
    Indices <- NULL
    if(ncol(data)>3){
        Indices <- as.matrix(data[,4:ncol(data),drop=FALSE])
    }
    if(is.null(Nmas)){
        Nmas <- 0:2
    }
    if(is.null(Nars)){
        Nars <- 0:2
    }
    if(is.null(NI.inds)){
        NI.inds <- list(0,1:3,1:5,c(1:3,6:10),c(1:3,11:18))
    }
    nis <- c()
    for(j in 1:length(NI.inds)){
            ni.ind <- NI.inds[[j]]
            nis <- c(nis,length(ni.ind[ni.ind!=0]))
    }
    ni.min <- min(nis)
    Nma.opt <- Nmas[1]
    Nar.opt <- Nars[1]
    Inds.opt <- NI.inds[[1]]
    NI.opt <- length(Inds.opt)
    Ndata <- nrow(data)
    logLmaxs <- array(data=NA,dim=c(length(NI.inds),length(Nmas),length(Nars)))
    logBFs <- array(data=NA,dim=c(length(NI.inds),length(Nmas),length(Nars)),dimnames=list(paste0('NI',nis),paste0('Nma',Nmas),paste0('Nar',Nars)))
    ind.opt <- 1
    for(j in 1:length(NI.inds)){
        if(!all(NI.inds[[j]]==0)){
            ni <- nis[j]
            Inds <- NI.inds[[j]]
            indices <- Indices[,Inds,drop=FALSE]
        }else{
            ni <- 0
            indices <- NULL
            Inds <- 0
        }
        if(j==1){
            Inds.opt <- Inds
            indices.opt <- indices
        }
#cat('indices=',indices,'\n')
        vars <- global.notation(t,y,dy,Indices=indices,Nma=0,Nar=0,GP=FALSE,gp.par=rep(NA,3))
        tmp <- sopt(omega=NA,phi=NA,Nma=0,Nar=0,Indices=indices,data=data,type='noise',par.low=vars$par.low,par.up=vars$par.up,start=vars$start,noise.only=FALSE,GP=FALSE,gp.par=rep(NA,3),Nrep=Nrep)#
        logLmaxs[j,1,1] <- tmp$logL
        logBFs[j,1,1] <- tmp$logL-logLmaxs[1,1,1]-0.5*(ni-ni.min)*log(Ndata)
        if(j>1){
            if(logBFs[j,1,1]>(logBFs[ind.opt,1,1]+5)){
                NI.opt <- ni
                Inds.opt <- Inds
                indices.opt <- indices
                ind.opt <- j
            }
        }
        if(NI.opt<(ni-1)) break()
    }
    ind.ma <- ind.ar <- 1
    for(k in 1:length(Nars)){
        nar <- Nars[k]
        for(i in 1:length(Nmas)){
            nma <- Nmas[i]
            if(k>1 | i>1){
                vars <- global.notation(t,y,dy,Indices=indices.opt,Nma=nma,Nar=nar,GP=FALSE,gp.par=rep(NA,3))
                tmp <- sopt(omega=NA,phi=NA,Nma=nma,Nar=nar,Indices=indices.opt,data=data,type='noise',par.low=vars$par.low,par.up=vars$par.up,start=vars$start,noise.only=FALSE,GP=FALSE,gp.par=rep(NA,3),Nrep=Nrep)#
                dN <- nis[ind.opt]-ni.min#number of additional free parameters
                if(nma>0){
                    dN <- dN+nma+1
                }
                if(nar>0){
                    dN <- dN+nar+1
                }
                logLmaxs[ind.opt,i,k] <-  tmp$logL
                logBFs[ind.opt,i,k] <-  tmp$logL-logLmaxs[1,1,1]-0.5*dN*log(Ndata)
                if(i>1 | k>1){
                    if(logBFs[ind.opt,i,k]>logBFs[ind.opt,ind.ma,ind.ar]+5){
                        Nar.opt <- nar
                        Nma.opt <- nma
                        ind.ma <- i
                        ind.ar <- k
                    }
                }
            }
            if(Nma.opt<nma) break()
            if(Nar.opt<nar) break()
        }
    }
    best.model <- paste0('ARMA(',Nar.opt,',',Nma.opt,')')
###compare the optimal noise model with GP
#cat('GP loglike\n')
    if(GP){
        vars <- global.notation(t,y,dy,Indices=indices,Nma=0,Nar=0,GP=TRUE,gp.par=rep(NA,3))
        tmp <- sopt(omega=NA,phi=NA,Nma=0,Nar=0,Indices=NULL,data=data,type='noise',par.low=vars$par.low,par.up=vars$par.up,start=vars$start,noise.only=FALSE,GP=TRUE,gp.par=vars$gp.par,Nrep=Nrep)
        logLmax.gp <- tmp$logL
        logBF.gp <- logLmax.gp-logLmaxs[1,1,1]-1.5*log(Ndata)
#cat('logBF.gp=',logBF.gp,'\n')
        if(logBF.gp>logBFs[ind.opt,ind.ma,ind.ar]) best.model <- 'GP'
    }else{
        logLmax.gp <- logBF.gp <- NULL
    }
    cat('best logBF=',logBFs[ind.opt,ind.ma,ind.ar],'\n')
#    cat('ind.opt=',ind.opt,';ind.ma=',ind.ma,';ind.ar=',ind.ar,'\n')
###
    cat('The optimal Nar=',Nar.opt,'; Nma=',Nma.opt,'; Inds=',Inds.opt,'\n')
    return(list(Nar=Nar.opt,Nma=Nma.opt,Inds=Inds.opt,logBFs=logBFs,logLmaxs=logLmaxs,logLmax.gp=logLmax.gp,logBF.gp=logBF.gp,best.model=best.model))
}

combine.data <- function(data,Ninds,Nmas,GP=FALSE,gp.par=NULL){
###data is a list of matrices
    out <- c()
    idata <- data
    for(j in 1:length(data)){
        Inds <- Ninds[[j]]
        tab <- data[[j]]
        t0 <- tab[,1]
        t <- tab[,1]-min(tab[,1])
        y <- tab[,2]
        dy <- tab[,3]
        if(ncol(tab)>4){
            Indices <- as.matrix(tab[,4:ncol(tab)])
        }else if(ncol(tab)==4){
            Indices <- matrix(tab[,4],ncol=1)
        }else{
            Indices <- NA
        }
        if(all(Inds==0)){
            NI <- 0
            Indices <- NA
        }else{
            NI <- length(Inds)
            if(NI>0){
                if(NI<2){
                    Indices <- matrix(Indices[,Inds],ncol=1)
                }else{
                    Indices <- Indices[,Inds]
                }
            }
        }
        vars <- global.notation(t,y,dy,Indices,Nmas[j],NI,GP,gp.par)
        tmp <- sopt(data=tab,Indices=Indices,NI=NI,Nma=Nmas[j],type='noise',pars=vars)
        val <- cbind(t0,tmp$res,dy)
        idata[[j]] <- val
        out <- rbind(out,val)
    }
    inds <- sort(out[,1],index.return=TRUE)$ix
    return(list(cdata=out[inds,],idata=idata))
}
fsample <- function(fmin,fmax,sampling,section=1,ofac=1,unit=1){
    logP.max <- log(1/fmin)
    logP.min <- log(1/fmax)
    dlogP <- (logP.max-logP.min)/section
    if(sampling=='logP'){
        f <- c()
        for(j in 1:section){
            f <- c(f,1/exp(seq(logP.min+dlogP*(j-1),logP.min+dlogP*j,length.out=1000*ofac))*unit)
        }
    }else if(sampling=='freq'){
        f <- c()
        for(j in 1:section){
            f1 <- 1/exp(logP.max-dlogP*(j-1))
            f2 <- 1/exp(logP.max-dlogP*j)
            f <- c(f,seq(f1,f2,length.out=1000/section*min(5^(j-1),20)*ofac)*unit)
        }
    }else if(sampling=='combined'){
        fl <- ff <- c()
        for(j in 1:section){
            fl <- c(fl,1/exp(seq(logP.min+dlogP*(j-1),logP.min+dlogP*j,length.out=1000*ofac))*unit)
        }
        for(j in 1:section){
            f1 <- 1/exp(logP.max-dlogP*(j-1))
            f2 <- 1/exp(logP.max-dlogP*j)
            ff <- c(ff,seq(f1,f2,length.out=1000/section*min(5^(j-1),20)*ofac)*unit)
        }
        f <- sort(c(fl,ff))
    }else if(sampling=='random'){
        f <- sort(1/runif(1000*ofac,exp(logP.min),exp(logP.max)))*unit
    }
    f <- unique(f)
    return(f)
}
####Bayes factor periodogram
#####moving periodogram
MP <- function(t, y, dy,Dt,nbin,fmax=1,ofac=1,fmin=1/1000,tspan=NULL,Indices=NA,per.type='MLP',...){
    n <- nbin-1
    dt <- (max(t)-min(t)-Dt)/n
    tstart <- min(t)+(0:n)*dt
    tend <- min(t)+(0:n)*dt+Dt
    tmid <- (tstart+tend)/2
    df <- 1/(Dt*ofac)
    rel.powers <- powers <- c()
    ndata <- rep(NA,nbin)
    withProgress(message = 'Calculating moving periodogram', value = 0, {
        for(j in 0:n){
            incProgress(1/nbin, detail = paste0('window ',j+1,'/',nbin))
            inds <- which(t>=min(t)+j*dt & t<min(t)+j*dt+Dt)
            index <- NULL
            if(!is.null(Indices) & length(inds)>0){
                index <- as.matrix(Indices[inds,,drop=FALSE])
            }
            if(per.type=='BGLS'){
                tmp <- bgls(t=t[inds],y=y[inds],err=dy[inds],fmax=fmax,ofac=ofac,fmin=fmin,tspan=Dt)
            }else if(per.type=='GLS'){
                tmp <- gls(t=t[inds],y=y[inds],err=dy[inds],fmax=fmax,ofac=ofac,fmin=fmin,tspan=Dt)
            }else if(per.type=='GLST'){
                tmp <- glst(t=t[inds],y=y[inds],err=dy[inds],fmax=fmax,ofac=ofac,fmin=fmin,tspan=Dt)
            }else if(per.type=='MLP'){
                tmp <- MLP(t=t[inds],y=y[inds],dy=dy[inds],fmax=fmax,ofac=ofac,fmin=fmin,tspan=Dt,Indices=index,...)
            }else if(per.type=='BFP'){
                tmp <- BFP(t=t[inds],y=y[inds],dy=dy[inds],fmax=fmax,ofac=ofac,fmin=fmin,tspan=Dt,Indices=index,...)
            }else if(per.type=='LS'){
                tmp <- lsp(times=t[inds],x=y[inds],ofac=ofac,from=fmin,to=fmax,tspan=Dt,alpha=c(0.1,0.01,0.001))
                ps <- 1/tmp$scanned
                pps <- tmp$power
                inds <- sort(ps,index.return=TRUE,decreasing=TRUE)$ix
                tmp[['power']] <- pps[inds]
                tmp[['P']] <- ps[inds]
            }
            index <- sort(tmp$P,index.return=TRUE)$ix
            relpower <- (tmp$power[index]-mean(tmp$power[index]))/(max(tmp$power[index])-mean(tmp$power[index]))
            rel.powers <- cbind(rel.powers,relpower)
            if(per.type=='MLP' | per.type=='BGLS'){
                pp <- tmp$power[index]-max(tmp$power[index])
            }else{
                pp <- tmp$power[index]
            }
            powers <- cbind(powers,pp)
            ndata[j+1] <- length(inds)
        }
    })
    return(list(tmid=tmid,P=tmp$P[index],powers=powers,rel.powers=rel.powers,ndata=ndata))
}

####MP without progress
MP.norm <- function(t, y, dy,Dt,nbin,fmax=1,ofac=1,fmin=1/1000,per.type='MLP',Indices=NA,...){
    n <- nbin-1
    dt <- (max(t)-min(t)-Dt)/n
    tstart <- min(t)+(0:n)*dt
    tend <- min(t)+(0:n)*dt+Dt
    tmid <- (tstart+tend)/2
    df <- 1/(Dt*ofac)
    rel.powers <- powers <- c()
    ndata <- rep(NA,nbin)
        for(j in 0:n){
#            cat(paste0(Dt,'d window'),j+1,'/',nbin,'\n')
            inds <- which(t>=min(t)+j*dt & t<min(t)+j*dt+Dt)
            index <- NULL
            if(!is.null(Indices) & length(inds)>0){
                index <- as.matrix(Indices[inds,,drop=FALSE])
            }
            if(per.type=='BGLS'){
                tmp <- bgls(t=t[inds],y=y[inds],err=dy[inds],fmax=fmax,ofac=ofac,fmin=fmin,tspan=Dt)
            }else if(per.type=='GLS'){
                tmp <- gls(t=t[inds],y=y[inds],err=dy[inds],fmax=fmax,ofac=ofac,fmin=fmin,tspan=Dt)
            }else if(per.type=='GLST'){
                tmp <- glst(t=t[inds],y=y[inds],err=dy[inds],fmax=fmax,ofac=ofac,fmin=fmin,tspan=Dt)
            }else if(per.type=='MLP'){
                tmp <- MLP(t=t[inds],y=y[inds],dy=dy[inds],fmax=fmax,ofac=ofac,fmin=fmin,tspan=Dt,Indices=index,...)
            }else if(per.type=='BFP'){
                tmp <- BFP(t=t[inds],y=y[inds],dy=dy[inds],fmax=fmax,ofac=ofac,fmin=fmin,tspan=Dt,progress=FALSE,...)
            }else if(per.type=='LS'){
                tmp <- lsp(times=t[inds],x=y[inds],ofac=ofac,from=fmin,to=fmax,tspan=Dt,alpha=c(0.1,0.01,0.001))
                ps <- 1/tmp$scanned
                pps <- tmp$power
                inds <- sort(ps,index.return=TRUE,decreasing=TRUE)$ix
                tmp[['power']] <- pps[inds]
                tmp[['P']] <- ps[inds]
            }
            index <- sort(tmp$P,index.return=TRUE)$ix
#            relpower <- (tmp$power[index]-mean(tmp$power[index]))/(max(tmp$power[index])-mean(tmp$power[index]))
            relpower <- (tmp$power[index]-min(tmp$power[index]))/length(inds)
            rel.powers <- cbind(rel.powers,relpower)
            if(per.type=='MLP' | per.type=='BGLS'){
                pp <- tmp$power[index]-max(tmp$power[index])
            }else{
                pp <- tmp$power[index]
            }
            powers <- cbind(powers,pp)
            ndata[j+1] <- length(inds)
        }
    return(list(tmid=tmid,P=tmp$P[index],powers=powers,rel.powers=rel.powers,ndata=ndata))
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

BFP <- function(t, y, dy, Nma=0, Nar=0,Indices=NULL,ofac=1, fmax=NULL,fmin=NA,tspan=NULL,sampling='combined',model.type='man',progress=TRUE,quantify=FALSE,dP=0.1,section=1,GP=FALSE,gp.par=rep(NA,3),noise.only=FALSE,Nsamp=1,par.opt=NULL,renew=TRUE,sj=0){
###gp.par is the free parameters of SHO-GP
    if(FALSE){
        cat('head(t)=',head(t),'\n')
        cat('head(y)=',head(y),'\n')
        cat('head(dy)=',head(dy),'\n')
        cat('Nma=',Nma,';Nar=',Nar,';model.type=man;Indices=',head(Indices), ';ofac=',ofac,';fmin=',fmin,';fmax=',fmax,';quantify=',quantify, ';renew=',renew,';noise.only=',noise.only,'\n')
    }
#    if(Nma==0 & Nar==0 & !GP) noise.only <- FALSE
    if(noise.only) quantify <- FALSE
    unit <- 1
    t <- (t-min(t))/unit#rescale time
    if(is.null(Indices)){
        NI <- 0
    }else{
        NI <- ncol(Indices)
        Indices <- as.matrix(Indices)
    }
    data <- cbind(t,y,dy)
    if(is.null(tspan)){
        tspan <- max(t)-min(t)
    }
    step <- 1/(tspan*ofac)
    if(is.na(fmin)){
        fmin <- 1/(tspan*ofac)
    }
    NN <- 1
    if(quantify) NN <- 2
    f <- fsample(fmin,fmax,sampling,section,ofac,unit)
    nout <- length(f)
    Ndata <- length(t)
    omegas <- 2*pi*f
    phi <- 0
    #######define notations and variables
    vars1 <- global.notation(t,y,dy,Indices=Indices,Nma=Nma,Nar=Nar,GP=GP,gp.par=gp.par)
#    if(noise.only & GP){
    if(noise.only){
        vars0 <- global.notation(t,y,dy,Indices=Indices,Nma=0,Nar=0,GP=FALSE,gp.par=gp.par)
    }else{
        vars0 <- vars1
    }
#    var <- names(vars)
####shor
#    for(k in 1:length(var)){
#        assign(var[k],vars[[var[k]]])
#    }
##
##########################################
#####optimizing the noise parameters
#########################################
    if(model.type=='auto'){
        t1 <- proc.time()
        out <- bfp.inf(vars1,Indices)
        t2 <- proc.time()
        dur <- format((t2-t1)[3],digit=3)
        cat('model comparison computation time:',dur,'s\n\n')
        NI <- out$NI
        Nma <- out$Nma
        Nar <- out$Nar
        vars1 <- out$vars1
        Indices <- out$Indices
    }

########################################################
##############baseline model
########################################################
    if(!is.null(par.opt)){
        vars0$start <- as.list(par.opt)
    }
    tmp <- sopt(omega=NA,phi=NA,Nma=vars0$Nma,Nar=vars0$Nar,Indices=Indices,data=data,type='noise',par.low=vars0$par.low,par.up=vars0$par.up,start=vars0$start,noise.only=noise.only,GP=vars0$GP,gp.par=vars0$gp.par,Nrep=5)#Nrep=1
#    pars <- unlist(c(omega=NA,phi=NA,Nma=vars0$Nma,Nar=vars0$Nar,Indices=Indices,data=data[1,],type='noise',par.low=vars0$par.low,par.up=vars0$par.up,start=vars0$start,noise.only=noise.only,GP=vars0$GP,gp.par=vars0$gp.par,Nrep=1))
    logL0 <- tmp$logL
    lnl0 <- tmp$lnls
    if(FALSE){
        ind <- match(names(vars1$start),names(tmp$par))
        vars1$start <- tmp$par[ind]
    }
    ##
    logLs <- rep(NA,length(omegas))
    lnls <- array(NA,dim=c(length(omegas),Ndata))
    Npar <- length(tmp$par)
    if(noise.only){
        Nextra <- 0
        if(GP){
            ind <- which(is.na(gp.par[-2]))
            Nextra <- Nextra+length(ind)
        }
        if(Nma>0){
            Nextra <- Nextra+2*Nma
        }
        if(Nar>0){
            Nextra <- Nextra+2*Nar
        }
    }else{
        Nextra <- 2
    }
    omegas.all <- c()
    opt.pars <- array(NA,dim=c(length(f),Npar+Nextra))
    P <- c()
#    if(!GP){
#        Nrep <- Nsamp
#    }else{
        Nrep <- Nsamp
#    }
#    cat('Nrep=',Nrep,'\n')
    for(nn in 1:NN){
        P <- c(unit/f,P)
        if(progress){
            withProgress(message = 'Calculating BFP', value = 0, {
                for(kk in 1:length(f)){
#####for periodogram, set sj=0
                    incProgress(1/length(f), detail = paste0(round(kk/length(f)*100),'%'))
                    omega <- omegas[kk]
################################################
####I optimization
####################################################
                    t0 <- proc.time()
                    if(noise.only){
                        if(GP){
                            gp.par[2] <- log(2*pi/omega)
                            vars1$gp.par <- gp.par
                        }
                        opt <- sopt(omega=omega,phi=phi,Nma=Nma,Nar=Nar,Indices=Indices,data=data,type='noise',par.low=vars1$par.low,par.up=vars1$par.up,start=vars1$start,noise.only=noise.only,GP=GP,gp.par=gp.par,Nrep=Nrep)
                    }else{
                        tmp <- local.notation(t,y,dy,Indices,NI,omega,phi)
                        vars1 <- c(vars1,tmp)
                        opt <- sopt(omega=omega,phi=phi,Nma=Nma,Nar=Nar,Indices=Indices,data=data,type='period',par.low=vars1$par.low,par.up=vars1$par.up,start=vars1$start,noise.only=noise.only,GP=GP,gp.par=gp.par,Nrep=Nrep)
                    }
                    if(renew){
                        vars1$start <- opt$par0
                    }
                    logLs[kk] <- opt$logL
                    lnls[kk,] <- opt$lnls
                    opt.pars[kk,] <- unlist(opt$par)
                }
            })
        }else{
            if(GP & noise.only){
#                if(renew){
                    vars1$start <- vars1$start[-grep('logProt',names(vars1$start))]
#                }
                    vars1$par.low <- vars1$par.low[-grep('logProt',names(vars1$par.low))]
                    vars1$par.up <- vars1$par.up[-grep('logProt',names(vars1$par.up))]
            }
            for(kk in 1:length(f)){
#            for(kk in length(f):1){
#####for periodogram, set sj=0
#                cat(kk,'/',length(f),'\n')

                omega <- omegas[kk]
################################################
####I optimization
####################################################
                t0 <- proc.time()
                if(noise.only){
                    if(GP){
                        gp.par[2] <- log(2*pi/omega)
                        vars1$gp.par <- gp.par
                    }
#                    cat('omega=',omega,'\n')
#                    opt <- sopt(omega=omega,phi=phi,Nma=Nma,Nar=Nar,Indices=Indices,data=data,type='noise',par.low=vars1$par.low,par.up=vars1$par.up,start=vars1$start,noise.only=noise.only,GP=GP,gp.par=gp.par,Nrep=Nrep)
                    opt <- sopt(omega=omega,phi=phi,Nma=Nma,Nar=Nar,Indices=Indices,data=data,type='noise',par.low=vars1$par.low,par.up=vars1$par.up,start=vars1$start,noise.only=noise.only,GP=GP,gp.par=gp.par,Nrep=10)
                }else{
                    tmp <- local.notation(t,y,dy,Indices,NI,omega,phi)
                    vars1 <- c(vars1,tmp)
                    opt <- sopt(omega=omega,phi=phi,Nma=Nma,Nar=Nar,Indices=Indices,data=data,type='period',par.low=vars1$par.low,par.up=vars1$par.up,start=vars1$start,noise.only=noise.only,GP=GP,gp.par=gp.par,Nrep=Nrep)
                }
                if(renew){
                    vars1$start <- opt$par0
                }
                logLs[kk] <- opt$logL
                lnls[kk,] <- opt$lnls
                opt.pars[kk,] <- unlist(opt$par)
            }
        }

####signals
        Pmax <- P[which.max(logLs)]
        if(quantify & nn==1){
            P1 <- P
            logLs1 <- logLs
            lnls1 <- lnls
            opt.pars1 <- opt.pars
            ind.max <- which.max(logLs)
            fmin <- (1-dP)*f[ind.max]
            fmax <- (1+dP)*f[ind.max]
###oversampling
            f <- seq(fmin,fmax,by=step/20)*unit
                                        #        f <- seq(fmin,fmax,length.out=1000)*unit
            Nf <- length(f)
            omegas <- f*2*pi
            logLs <- c(rep(NA,length(f)),logLs)
            lnls <- rbind(array(NA,dim=c(length(f),Ndata)),lnls)
            opt.pars <- rbind(array(data=NA,dim=c(length(f),ncol(opt.pars))),opt.pars)
        }else if(quantify & nn==2){
            ind.max <- which.max(logLs)
            ind1 <- which.max(logLs1)
            if(ind.max<Nf){
                P1[ind1] <- P[ind.max]
                logLs1[ind1] <- logLs[ind.max]
                lnls1[ind1,] <- lnls[ind.max,]
                opt.pars1[ind1,] <- opt.pars[ind.max,]
            }
            P <- P1
            logLs <- logLs1
            lnls <- lnls1
            opt.pars <- opt.pars1
        }
        omegas.all <- c(omegas.all,omegas)
    }
    colnames(opt.pars) <- names(opt$par)
    if(noise.only){
        name0 <- names(opt$par)
#        cat('names(opt$par)=',names(opt$par),'\n')
        if(Nma>0){
            opt.pars[,'logtau'] <- log(2*pi/omegas.all)
#            opt.pars <- cbind(opt.pars,)
#            name0 <- colnames(opt.pars) <- c(name0,'logtau')
        }
        if(Nar>0){
            opt.pars[,'logtauAR'] <- log(2*pi/omegas.all)
#            opt.pars <- cbind(opt.pars,log(2*pi/omegas.all))
#            colnames(opt.pars) <- c(name0,'logtauAR')
        }
    }
###calculate BF
    logBF <- logLs-logL0-Nextra/2*log(Ndata)#BIC-estimated BF; the extra free parameter n=2, which are {A, B}
    lnbfs <- lnls-lnl0-(Nextra/2*log(Ndata))/Ndata#BIC-estimated BF; the extra free parameter n=2, which are {A, B}
####signals
    ind.max <- which.max(logBF)
    llmax <- max(logLs)
    inds <- sort(logBF,decreasing=TRUE,index.return=TRUE)$ix[1:10]
    Popt <- P[inds[1]]
    Popts <- P[inds]
    opt.par <- opt.pars[inds,]
    par.opt <- opt.pars[inds[1],]
    logBF.opt <- logBF[inds]

####calculate the residual
    par.fix <- list(omega=2*pi/Popt[1],phi=0)
    df <- opt$df
    df$par.fix <- par.fix
    if(noise.only){
        df$type <- 'noise'
    }else{
        df$type <- 'period'
    }
    if(is.matrix(opt.par) | is.data.frame(opt.par)){
        pp <- as.list(opt.par[1,])
    }else{
        pp <- as.list(opt.par)
    }
    names(pp) <- colnames(opt.pars)
    yall <- as.numeric(CircularSig(pp,df)$v)
    if(!noise.only){
        ysig <- pp$A*cos(2*pi/Popt[1]*t)+pp$B*sin(2*pi/Popt[1]*t)
        ysigt <- ysig+pp$gamma+pp$beta*t
    }else{
        ysig <- ysigt <- rep(0,length(t))
    }
    yred <- yall-ysigt
    yredt <- yall-ysig
    res.nst <- y-yall
    res.n <- y-yred
    res.nt <- y-yredt
    res.s <- y-ysig
    res.st <- y-ysigt
    inds <- sort(logBF,decreasing=TRUE,index.return=TRUE)$ix
    ps <- P[inds[1:5]]
    power.opt <- logBF[inds[1:5]]
    return(list(data=data,logBF=logBF,logLs=logLs,llmax=llmax,lnbfs=lnbfs,P=P,Popt=Popt,Popts=Popts,logBF.opt=logBF.opt,par.opt=par.opt,opt.par=opt.par,res.nst=res.nst,res=res.nst,res.n=res.n,res.s=res.s,res.st=res.st,res.nt=res.nt,sig.level=c(-Nextra/2*log(length(y)),0,log(150)), power=logBF,ps=ps,power.opt=power.opt,ysig=ysig,pars=opt.pars,df=df,ParLow=vars1$par.low,ParUp=vars1$par.up,LogLike0=logL0))
}
Circ2kep <- function(per,basis='linear'){
####################################################
## Find the Keplerian signal corresponding to the circular signal
## Input:
##   per - Output of a peridogram
##
## Output:
##   ParKep - Parameters of a Keplerian model
####################################################
    frac <- 0.1
    t <- per$df$data[,1]
    Popt <- as.numeric(per$Popt[1])#day
#    if(is.null(dim(per$par.opt))) per$par.opt <- t(per$par.opt)
    PhiOpt <- as.numeric(xy2phi(per$par.opt['A'],per$par.opt['B']))
    Kopt <- as.numeric(sqrt(per$par.opt['A']^2+per$par.opt['B']^2))
    ind <- which(names(per$par.opt)!='A' & names(per$par.opt)!='B')
    ParLow <- per$ParLow
    ParUp <- per$ParUp
    ParOpt <- per$par.opt[ind]
    if(any(names(per$par.opt)=='beta')){
        betaLow <- as.numeric(per$par.opt['beta']-5*sd(per$par.opt['beta']))
        betaUp <- as.numeric(per$par.opt['beta']+5*sd(per$par.opt['beta']))
        ParLow <- c(beta=betaLow,ParLow)
        ParUp <- c(beta=betaUp,ParUp)
    }
    if(any(colnames(per$par.opt)=='gamma')){
        gammaLow <- as.numeric(per$par.opt['gamma']-5*sd(per$pars[,'gamma']))
        gammaUp <- as.numeric(per$par.opt['gamma']+5*sd(per$pars[,'gamma']))
        ParLow <- c(gamma=gammaLow,ParLow)
        ParUp <- c(gamma=gammaUp,ParUp)
    }
    ParLow0 <- ParLow
    ParUp0 <- ParUp
    if(basis=='natural'){
        ParLow <- c(P1=(1-frac)*Popt,K1=max((1-frac)*Kopt,0),e1=0,omega1=0,Mo1=0,ParOpt-frac*(ParUp0-ParLow0))
        ParUp <- c(P1=(1+frac)*Popt,K1=(1+frac)*Kopt,e1=1,omega1=2*pi,Mo1=2*pi,ParOpt+frac*(ParUp0-ParLow0))
    }else{
        ParLow <- c(P1=(1-frac)*Popt,K1=(1-frac)*Kopt,esinw1=-sqrt(2)/2,ecosw1=-sqrt(2)/2,Tc1=min(t),ParOpt-frac*(ParUp0-ParLow0))
        ParUp <- c(P1=(1+frac)*Popt,K1=(1+frac)*Kopt,esinw1=sqrt(2)/2,ecosw1=sqrt(2)/2,Tc1=min(t)+(1+frac)*Popt,ParOpt+frac*(ParUp0-ParLow0))
    }
    Ntry <- 100
    pars <- c()
    ll <- c()
    n <- 1
    for(j in 1:Ntry){
        ParIni <- runif(length(ParLow),ParLow,ParUp)
        names(ParIni) <- names(ParLow)
        ParIni <- as.list(ParIni)
        df <- c(per$df,basis=basis)
        eps <- 1e-16
        out <- try(nls.lm(par = ParIni,lower=ParLow,upper=ParUp,fn = KepRes,data=df,control=nls.lm.control(maxiter=1024,ptol=eps,gtol=eps,ftol=eps)),TRUE)
#        out <- nls.lm(par = ParIni,lower=ParLow,upper=ParUp,fn = KepRes,data=df,control=nls.lm.control(maxiter=1024,ptol=eps,gtol=eps,ftol=eps))
        if(class(out)!='try-error'){
            pars <- rbind(pars,coef(out))
            ll <- c(ll,-sum(out$fvec^2))#logLike
            n <- n+1
        }
#cat('n=',n,'\n\n')
    }
    ind <- which.max(ll)
    llmax <- ll[ind]
    ParKep <- as.list(pars[ind,])
    if(basis=='linear' & FALSE){
        e  <- sqrt(ParKep$ecosw^2+ParKep$esinw^2)
    }
    DlogLike <-llmax-per$LogLike0
    lnBF3 <- DlogLike-1.5*log(length(per$df$data))
    lnBF5 <- DlogLike-2.5*log(length(per$df$data))
    Pred <- KeplerRv(ParKep,df)
    return(list(ParKep=ParKep,ParLow=ParLow,ParUp=ParUp,LogLike=ll,lnBF3=lnBF3,lnBF5=lnBF5,ll0=per$LogLike0,lls=ll,pars=pars,Pred=Pred,df=df))
}

LMkepler <- function(ParIni,ParLow,ParUp,df){
####################################################
## Using LM algorithm to constrain the Keplerian parameters for a sinusoidal/circular signal found by BFP or GLST or other periodograms
## Input:
##   ParIni - Parameters of a circular signal
##
## Output:
##   ParKep - Parameters of a Keplerian signal
####################################################
#Set the lower and upper boundary
    return(ParKep)
}
detIni <- function(par,par.low,par.up){
    Ntry <- 100
    for(k3 in 1:Ntry){
#        start <- rnorm(length(par),unlist(par),as.numeric((par.up-par.low)/2))
        start <- runif(length(par),par.low,par.up)
        if(all(start>par.low & start<par.up)) break()
    }
    names(start) <- names(par)
    as.list(start)
}

scale.proxy <- function(proxy){
    if(!is.null(proxy)){
        if(!is.null(dim(proxy))){
            for(j in 1:ncol(proxy)){
                proxy[,j] <- scale(proxy[,j])
            }
        }else{
            proxy <- scale(proxy)
        }
    }
    if(is.null(dim(proxy))) proxy <- t(t(proxy))
    return(proxy)
}

CircularFit <- function(per,data,Nsamp=10){
    out <- list()
    llmax <- c()
    ps <- c()
    for(j in 1:Nsamp){
        if(j==1){
            popt <- per$Popt
        }else{
            popt <- rnorm(1,per$Popt,1e-3*per$Popt)
        }
        ps <- c(ps,popt)
        vars <- global.notation(data[,1],data[,2],data[,3],Nma=per$df$Nma,Nar=per$df$Nar,Indices=per$df$Indices,GP=per$df$GP,gp.par=rep(NA,3))
        tmp <- sopt(omega=2*pi/popt,phi=0,Nma=per$df$Nma,Nar=per$df$Nar,Indices=per$df$Indices,data=data,type='period',par.low=vars$par.low,par.up=vars$par.up,start=vars$start,noise.only=FALSE,GP=FALSE,gp.par=rep(NA,3),Nrep=1)#
        cat('tmp$logL=',tmp$logL,'\n')
        llmax <- c(llmax,tmp$logL)
        out[[j]] <- tmp
    }
    ind <- which.max(llmax)
    c(out[[ind]],list(Popt=ps[ind]))
}

KeplerFit <- function(per,data,basis='natural'){
####################################################
## Keplerian fitting for a given circular fitting
## Input:
##   per - Output of a peridogram
##
## Output:
##   ParKep - Parameters of a Keplerian model
####################################################
    frac <- 0.1
    t <- per$data[,1]
    Popt <- as.numeric(per$Popt[1])#day
    if(is.list(per$par.opt)) per$par.opt <- unlist(per$par.opt)
    if(is.null(dim(per$par.opt))){
        per$par.opt <- t(replicate(2,per$par.opt))
    }
###original prediction
    PhiOpt <- as.numeric(xy2phi(per$par.opt[1,'A'],per$par.opt[1,'B']))
    Kopt <- as.numeric(sqrt(per$par.opt[1,'A']^2+per$par.opt[1,'B']^2))
    ind <- which(colnames(per$par.opt)!='A' & colnames(per$par.opt)!='B')
    if(is.null(per$ParLow)){
#        vars <- global.notation(per$data[,1],per$data[,2],per$data[,3],Nma=0,Nar=0,Indices=NA,GP=FALSE,gp.par=FALSE)
#        ParLow <- vars$par.low
#        ParUp <- vars$par.up
        ParLow <- c(sj=-1e-3)
        ParUp <-  c(sj=1e-3)
    }else{
        ParLow <- per$ParLow
        ParUp <- per$ParUp
    }
    ParOpt <- per$par.opt[1,ind]
    if(!any(names(ParOpt)=='sj')){
        ParOpt <- c(ParOpt,sj=0)
    }
    ii <- grep('^d\\d',colnames(per$par.opt))
    if(length(ii)>0){
        dlow <- dup <- c()
        for(i in ii){
            dlow <- c(dlow,per$par.opt[1,i]-max(5*sd(per$par.opt[,ii]),10))
            dup <- c(dup,per$par.opt[1,i]+max(5*sd(per$par.opt[,ii]),10))
        }
        names(dlow) <- names(dup) <- colnames(per$par.opt)[ii]
        ParLow <- c(dlow,ParLow)
        ParUp <- c(dup,ParUp)
    }
    if(any(colnames(per$par.opt)=='beta')){
        betaLow <- as.numeric(per$par.opt[1,'beta']-5*sd(per$par.opt[,'beta']))
        betaUp <- as.numeric(per$par.opt[1,'beta']+5*sd(per$par.opt[,'beta']))
        ParLow <- c(beta=betaLow,ParLow)
        ParUp <- c(beta=betaUp,ParUp)
    }
    if(any(colnames(per$par.opt)=='gamma')){
        gammaLow <- as.numeric(per$par.opt[1,'gamma']-max(5*sd(per$par.opt[,'gamma']),100))
        gammaUp <- as.numeric(per$par.opt[1,'gamma']+max(5*sd(per$par.opt[,'gamma']),100))
        ParLow <- c(gamma=gammaLow,ParLow)
        ParUp <- c(gamma=gammaUp,ParUp)
    }
    jj <- match(names(ParOpt),names(ParLow))
    ParLow0 <- ParLow <- ParLow[jj]
    ParUp0 <- ParUp <- ParUp[jj]
    if(basis=='natural'){
        ParLow <- c(P1=(1-frac)*Popt,K1=max((1-frac)*Kopt,0),e1=0,omega1=0,Mo1=0,ParOpt-frac*(ParUp0-ParLow0))
        ParUp <- c(P1=(1+frac)*Popt,K1=(1+frac)*Kopt,e1=1-1e-3,omega1=2*pi,Mo1=2*pi,ParOpt+frac*(ParUp0-ParLow0))
    }else{
        ParLow <- c(P1=(1-frac)*Popt,K1=(1-frac)*Kopt,esinw1=-sqrt(2)/2,ecosw1=-sqrt(2)/2,Tc1=min(t),ParOpt-frac*(ParUp0-ParLow0))
        ParUp <- c(P1=(1+frac)*Popt,K1=(1+frac)*Kopt,esinw1=sqrt(2)/2,ecosw1=sqrt(2)/2,Tc1=min(t)+(1+frac)*Popt,ParOpt+frac*(ParUp0-ParLow0))
    }
    Ntry <- 100
    pars <- c()
    ll <- c()
    n <- 1
    for(j in 1:Ntry){
        ParIni <- runif(length(ParLow),ParLow,ParUp)
        names(ParIni) <- names(ParLow)
        ParIni <- as.list(ParIni)
#    df <- list(data=tmp,Indices=Indices,par.fix=par.fix,Nma=Nma,NI=NI,GP=GP)
        if(!is.null(per$df)){
            df <- per$df
            df$basis <- basis
        }else{
            df <- list(data=per$data,Indices=NULL, basis=basis)
        }
#        df <- list(data=cbind(t,y,dy),Indices=Indices,par.fix=par.fix,Nma=Nma,Nar=Nar,NI=NI,GP=GP)
        eps <- 1e-16
        out <- try(nls.lm(par = ParIni,lower=ParLow,upper=ParUp,fn = KeplerRes,df=df,control=nls.lm.control(maxiter=1024,ptol=eps,gtol=eps,ftol=eps)),TRUE)
#        cat('ParIni=',unlist(ParIni),'\n')
#        out <- nls.lm(par = ParIni,lower=ParLow,upper=ParUp,fn = KeplerRes,df=df,control=nls.lm.control(maxiter=1024,ptol=eps,gtol=eps,ftol=eps))
        if(class(out)!='try-error'){
            pars <- rbind(pars,coef(out))
            ll <- c(ll,-sum(out$fvec^2))#logLike
            n <- n+1
        }
#cat('j=',j,'\n')
    }
    ind <- which.max(ll)
    llmax <- ll[ind]
    ParKep <- as.list(pars[ind,])
    DlogLike <-llmax-per$LogLike0
    lnBF3 <- DlogLike-1.5*log(nrow(per$data))
    lnBF5 <- DlogLike-2.5*log(nrow(per$data))
    return(list(ParKep=ParKep,ParLow=ParLow,ParUp=ParUp,LogLike=ll,lnBF3=lnBF3,lnBF5=lnBF5,ll0=per$LogLike0,lls=ll,pars=pars,df=df))
}
