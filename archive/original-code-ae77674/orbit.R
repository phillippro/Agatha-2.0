######################################################################
#####This file contains functions for the orbits of target system.
######################################################################
#####Solve Kepler's equation
solveKepler <- function(m,e){
    tol = 1e-8
    E0 <- m
    Ntt <- 1e2
    for(k in 1:Ntt){
        E1 = E0-(E0-e*sin(E0)-m)/(sqrt((1-e*cos(E0))^2-(E0-e*sin(E0)-m)*(e*sin(E0))))
        if(all(abs(E1-E0)<tol)) break()
        if(k==Ntt){
            cat('Warning: Keplerian solver does not converge!\n')
            cat('Warning: length(which(abs(E1-E0)>tol))=',length(which(abs(E1-E0)>tol)),'\n')
        }
        E0 <- E1
    }
    return(E1)
}

kepler.classic <- function(te,tpos,pars){
##input:m1(target mass, Msun); m2(companion/planet mass, Msun); P(orbital period, yr); e(eccentricity); Omega(longitude of ascending node,rad); omega(argument of periastron, rad); Mo(mean anomaly and typically a vector, rad)
##output: motion in the plane coordinates, {X,Y,Z,VX,VY,VZ} [au; au/yr]
##output: geometry delays
###auxiliary variables
#    m1 <- (pars['mtot']-pars['m2'])*Tsun
    m1 <- (pars['mtot']-pars['m2'])
    m <- pars['mtot']
#    m2 <- pars['m2']*Tsun
    m2 <- pars['m2']
    mu <-m1*m2/m
    Py <- pars['pb']#yr
    Pd <- Py*DJY#yr
    e <- pars['ecc']
    I <- asin(pars['sini'])
    Omega <- pars['Omega']
    omega <- pars['omz']
    dt <- (te[,1]-pars['Tp'])+te[,2]
    Mt <- 2*pi*(dt%%Pd)/Pd#mean anomaly
    n <- 2*pi/Py#1/yr
    m <- m1+m2#Msun
    arr <- (m*Py^2)^{1/3}#au
    ar <- arr*m2/m
    E <- solveKepler(Mt,e)%%(2*pi)
####Keplerian motion in the orbital plane
    x <- ar*(cos(E)-e)#au
    y <- ar*(sqrt(1-e^2)*sin(E))
    vx <- -ar*n*sin(E)/(1-e*cos(E))#au/yr
    vy <- ar*n*sqrt(1-e^2)*cos(E)/(1-e*cos(E))
###Thiele Innes constants; ref. Wright et al. 2009
###Definition:(ex:North; ey:East; ez:target to observer)
    A <- cos(Omega)*cos(omega)-sin(Omega)*sin(omega)*cos(I)
    B <- sin(Omega)*cos(omega)+cos(Omega)*sin(omega)*cos(I)
    F <- -cos(Omega)*sin(omega)-sin(Omega)*cos(omega)*cos(I)
    G <- -sin(Omega)*sin(omega)+cos(Omega)*cos(omega)*cos(I)
    C <- sin(omega)*sin(I)#this is negative in Catanzarite 2010
    H <- cos(omega)*sin(I)#this is negative in Catanzarite 2010
####sky plane coordinates [pb, qb, ub]
    Y <- A*x+F*y
    X <- B*x+G*y
    Z <- C*x+H*y#Z is from observer to target; opposite in Catanzarite 2010
    VY <- A*vx+F*vy
    VX <- B*vx+G*vy
    VZ <- C*vx+H*vy
####return state vector;[{X,Y,Z}]=au;[{VX,VY,VZ}]=au/yr
    list(state=cbind(X,Y,Z,VX,VY,VZ),u=E)
}

###ref: DD86; Damour & Taylor 1991; Taylor & Weisberg 1989
kepler.PN <- function(tB,pars){
##input parameter set: bbat(binary barycentric time), m1 and m2 are fittable only for astrometry and RV data; PK (Keplerian parameters); PP (separately measurable post-Keplerian parameters); PN (not separately measurable post-Keplerian parameters)
##output parameters: position and location of target star at a given epoch
##Note 1: all orbital parameters are derive in proper time T
###load parameters
#Keplerian parameters
    dt <-rowSums(tB-tpos)
    pb <- pars['pb']
    e0 <- pars['e0']
    omega0 <- pars['omega']
    Omega0 <- pars['Omega']
    x0 <- pars['x0']#combined function of ar and I

#measurable Post-Keplerian parameters
    Pdot <- pars['Pdot']#10^-18 year year^-1; or 10^-18 s s^-1; variation of period per year or per s
    wdot <- pars['wdot']#degree yr^-1
    gamma <- pars['gamma']#ms or 10^-3 s
    r <- pars['r']#mus; 10^-6 s
    s <- pars['s']#sin(I)
    xdot <- pars['xdot']#10^-13
    edot <- pars['edot']#10^-14 s^-1
    dt <- pars['dt']##fittable in DT91, but not in TW89
    dr <- pars['dr']##fittable in high precision astrometry and RV models
####parameters A and B are only fittable in pulsar timing and is negnigible
    A <- pars['A']
    B <- pars['B']

###rotation matrix
    R1 <- rbind(c(sin(Omega),-cos(Omega),0),c(cos(Omega),sin(Omega),0),c(0,0,1))
    R2 <- rbind(c(1,0,0),c(0,-cos(I),-sin(I)),c(0,sin(I),-cos(I)))
    R <- R1%*%R2
##new variables
    P <- pb^2/(pb+0.5*pi*Pdot*(tB-tpos))
#    n <- 2*pi/P0+pi*Pdot*(tB-tpos)/P0
    m <- m1+m2#Msun
    arr <- (m*P^2)^{1/3}#au
    ar <- a*m2/m#au
    e <- e0+edot*dt
    n <- 2*pi/P0+pi*Pdot*dt/P0^2
    Mo <- (n*dt)%%(2*pi)
    u <- solveKepler(Mo,e)%%(2*pi)
    er <- e*(1+dr)
    M <- u-e*sin(u)
    et <- e*(1+dt)
    Ae <- 2*atan(((1+e)/(1-e))^2*tan(u/2))
    Aet <- 2*atan(((1+et)/(1-et))^2*tan(u/2))
    k <- wdot/n
    omega <- omega0+k*Ae
    theta <- omega+Ae
    x0 <- ar*sin(I)*AULT#Keplerian timing parameter
    x <- x0+xdot*dt
###according to Damour & Taylor 1992, PRD, 45, 6
    C <- cos(omega+Ae)+e*cos(omega)
    S <- sin(omega+Ae)+e*sin(omega)
    beta <- n*x*(1-e^2)^{-0.5}
###position
    b <- ar*(1-er*cos(u))
    r0 <- t(cbind(b*cos(theta),b*sin(theta),0))
    r1 <- R%*%r0
###velocity
    v0 <- t(beta/sin(I)*cbind(-S,C,0)*CKMPS)#km/s
    v1 <- R%*%v0
###Einstein delay
    dET <- einstein.target(gamma,u)
###Shapiro delay
    dST <- shapiro.target(omega,u,r,e,w,s)
    list(state=cbind(t(r1),t(v1)),dET=dET,dST=dST)
}

calMag <- function(dv,vref){
###calculate magnitude or length of a vector
###dv: vector increment (au) vref: reference vector (pc)
###Taylor expansion to second order of dv, residual O(dv^3)
    Vref <- as.numeric(sqrt(vref%*%vref))
    vref <- t(replicate(nrow(dv),vref))
    term1 <- rowSums(vref*dv)/Vref#au
    term2 <- 0.5*sapply(1:nrow(dv), function(i) dv[i,]%*%dv[i,]/Vref)/pc2au#au
    term3 <- -0.5*sapply(1:nrow(dv), function(i) (vref%*%dv[i,])^2/Vref^3)/pc2au#au
    return(term1+term2+term3)
}

####SB states
calSB <- function(par.astro,tB,tpos){
##input: par.astro
##output: SB (r[pc],v[au/yr])
###coordinate frame formed by the [p q u] triad
    d0 <- 1/par.astro['plx']#kpc
    pqu <- cbind(par.astro[paste0('p',1:3)],par.astro[paste0('q',1:3)],par.astro[paste0('u',1:3)])

    v1 <- c(d0*c(par.astro['pmra'],par.astro['pmdec']),par.astro['rv']/auyr2kms)#au/yr in p,q,u system, assuming zero galactic acceleration, vSB(t)=vSB(t0)
    vSB <- as.numeric(pqu%*%v1)#au/yr
    dt <- as.numeric((tB[,1,drop=FALSE]-tpos[1])+(tB[,2,drop=FALSE]-tpos[2]))/DJY#year
#    if(UNITS=='TDB'){
#        dt <- dt/IFTE.K
#    }
    drSB <- cbind(dt*vSB[1],dt*vSB[2],dt*vSB[3])/pc2au#pc
                                        #initial rSB
    rSB0 <- d0*1e3*c(cos(dec0)*cos(ra0),cos(dec0)*sin(ra0),sin(dec0))#pc
cat('dec0=',dec0,'\n')
cat('drSB=',drSB,'\n')
    rSB <- cbind(drSB[,1]+rSB0[1],drSB[,2]+rSB0[2],drSB[,3]+rSB0[3])
#    if(UNITS=='TDB'){
#        rSB <- rSB/IFTE.K
#        rSB0 <- rSB0/IFTE.K
#    }
    RSB <- sqrt(rowSums(rSB^2))
    RSB0 <- sqrt(sum(rSB0^2))
#    SB <- cbind(asNumeric(rSB),t(replicate(nrow(t),vSB)))
    SB <- cbind(rSB,t(replicate(nrow(tB),vSB)))

###compared with the TEMPO2 approach
    if(TRUE){
        p <- pqu[,1]
        q <- pqu[,2]
        u <- pqu[,3]
        mu.perp <- par.astro['pmra']*p+par.astro['pmdec']*q/IFTE.K#mas/yr; from ephemeris to coordinate
        mu.para <- par.astro['rv']/auyr2kms/d0/IFTE.K#mas/yr
        u1 <- outer(dt,mu.perp*DMAS2R,FUN='*')#rad
        u2 <- -outer(dt^2*DMAS2R^2,0.5*sum(mu.perp^2)*u,'*')#rad
        u3 <- -outer(dt^2*DMAS2R^2,mu.para*mu.perp,'*')#rad
        u0 <- t(replicate(nrow(tB),u))
        du1 <- t(replicate(length(dt),mu.perp*DMAS2R))#rad/yr
        du2 <- -outer(2*dt*DMAS2R^2,0.5*sum(mu.perp^2)*u,'*')#rad/yr
        du3 <- -outer(2*dt*DMAS2R^2,mu.para*mu.perp,'*')#rad/yr
        uSBt <- u0+u1+u2+u3
        tempo <- list(uSBt=uSBt,u0=u0,u1=u1,u2=u2,u3=u3,du1=du1,du2=du2,du3=du3)
    }else{
        tempo <- NULL
    }
    ##uSBt <- uSBt/sqrt(rowSums(uSB1^2))
    list(SB=SB,dRSB=RSB-RSB0,tempo=tempo)
}

calBT <- function(tB,tpos,par.astro,pars,model){
    if(model=='kepler'){
##classical Keplerian orbit
        out <- kepler.classic(tB,tpos,pars)
        BT <- out$state
    }else if(model=='DD'){
        out <- DDmodel(tB=tB,tpos=tpos,pars)
        BT <- cbind(out$rvec,out$vvec)
    }else if(model=='DDGR'){
##Einstein theory-based Keplerian orbit
        out <- DDGRmodel(tB=tB,tpos,pars)#bbat!=te; change
        BT <- cbind(out$rvec,out$vvec)
    }else if(model=='SIM'){
##simulated orbit
        BT <- array(0,c(nrow(tB),6))
        out <- c()
    }
##transform BT from sky-plane coordinates, [p,q,u] to ICRS coordinates
    ra0 <- par.astro['ra']*pi/180
    dec0 <- par.astro['dec']*pi/180

#perspective acceleration
    p <- c(-sin(ra0),cos(ra0),0)
    q <- c(-sin(dec0)*cos(ra0),-sin(dec0)*sin(ra0),cos(dec0))
    u <- c(cos(dec0)*cos(ra0),cos(dec0)*sin(ra0),sin(dec0))
    x <- BT[,1]%o%p
    y <- BT[,2]%o%q
    z <- BT[,3]%o%u
    vx <- BT[,4]%o%p
    vy <- BT[,5]%o%q
    vz <- BT[,6]%o%u
    rBT <- x+y+z
    vBT <- vx+vy+vz
#    if(UNITS=='TDB') rBT <- rBT/IFTE.K
    c(list(BT=cbind(rBT,vBT)),out)
}


##calculate classical RV
RV.classic <- function(t,pars){
    tt <- t-min(t)
    rv <- c()
    for(j in 1:Np){
        P <- exp(pars[paste0(k,'logP',j)])#day
        K <- pars[paste0(k,'K',j)]#m/s
        e <- pars[paste0(k,'e',j)]#
        omega <- pars[paste0(k,'omega',j)]#
        M0 <- pars[paste0(k,'Mo',j)]#
        M <- (Mo+2*pi*tt/P)%%(2*pi)#mean anomaly
        E <- solveKepler(M,e)
        T <- 2*atan(sqrt((1+e)/(1-e))*tan(E/2))
        rv <- rv+K*(cos(omega+T)+e*cos(omega))
    }
###noise model parameters for different sets
###ind.set contain all the indices for different data sets
    for(k in 1:Nset){
        ind <- ind.set[[k]]
##trend
        b <- pars[grep(paste0('^',k,'b'),names(pars))]
        if(Npoly>0){
            rv[ind] <- rv[ind]+a*tt[ind]+b
        }
    }
    return(rv)
}

RV.relativistic <- function(pars,uOT,vST,vSO){
    VSO <- rowSums(vSO^2)
    VST <- rowSums(vST^2)
}

one2two <- function(var){
    var <- gsub('\\.',' 0\\.',var)
    var <- unlist(strsplit(var,' '))
    as.numeric(var[var!=''])
}

numDeriv <- function(x,y,N=2){
    if(N==2){
        out <- (y[-1]-y[-length(y)])/(x[-1]-x[-length(x)])
    }else if(N==3){
        if(length(y)>2){
            out <- (y[-(1:2)]-y[-(length(y)-(0:1))])/(x[-(1:2)]-x[-(length(y)-(0:1))])
        }else{
            out <- c()
        }
        out <- c((y[2]-y[1])/(x[2]-x[1]),out)#add the edge
    }else{
        out <- NULL
    }
    return(out)
}

numDeriv.2part <- function(x,y,N=2){
    out <- c()
    if(N==2){
        out <- rowSums(y[-1,]-y[-nrow(y),])/rowSums(x[-1,]-x[-nrow(x),])
    }else if(N==3){
        if(nrow(y)>2){
            if(nrow(y)==3){
                out <- sum(y[-(1:2),]-y[-(nrow(y)-(0:1)),])/sum(x[-(1:2),]-x[-(nrow(y)-(0:1)),])
            }else{
                out <- rowSums(y[-(1:2),]-y[-(nrow(y)-(0:1)),])/rowSums(x[-(1:2),]-x[-(nrow(y)-(0:1)),])
            }
        }
#        out <- c(sum(y[2,]-y[1,])/sum(x[2,]-x[1,]),out)#add the edge
    }
    return(out)
}

numZ.2part <- function(x,y,N=2){
    out <- c()
    if(N==2){
        if(nrow(y)>2){
            out <- (rowSums(y[-1,]-x[-1,])-rowSums(y[-nrow(y),]-x[-nrow(x),]))/rowSums(x[-1,]-x[-nrow(x),])
        }else{
            out <- sum((y[-1,]-x[-1,])-(y[-nrow(y),]-x[-nrow(x),]))/sum(x[-1,]-x[-nrow(x),])
        }
    }else if(N==3){
        if(nrow(y)>2){
            if(nrow(y)==3){
                out <- (sum(y[-(1:2),]-x[-(1:2),])-sum(y[-(nrow(y)-(0:1)),]-x[-(nrow(y)-(0:1)),]))/sum(x[-(1:2),]-x[-(nrow(y)-(0:1)),])
            }else{
                out <- (rowSums(y[-(1:2),]-x[-(1:2),])-rowSums(y[-(nrow(y)-(0:1)),]-x[-(nrow(y)-(0:1)),]))/rowSums(x[-(1:2),]-x[-(nrow(y)-(0:1)),])
            }
        }
#        out <- c(sum(y[2,]-y[1,])/sum(x[2,]-x[1,]),out)#add the edge
    }
    return(out)
}

##not accurate: https://www.mathworks.com/matlabcentral/fileexchange/15285-geodetic-toolbox
jd2yr2 <- function(jd){
    if(!is.null(dim(jd))){
        mjd <- jd[,1]-DJM0
        yr <- 1970+((mjd-40587)+0.5+jd[,2])/DJY
    }else{
        mjd <- jd-DJM0
        yr <- 1970+((mjd-40587)+0.5)/DJY
    }
    return(yr)
}

###https://www.mathworks.com/matlabcentral/fileexchange/15285-geodetic-toolbox
doy2jd <- function(yr,doy){
    jd <- Cal2jd(cbind(yr,1,0))
    cbind(jd[,1]+doy,jd[,2])
}

###accurate: https://www.mathworks.com/matlabcentral/fileexchange/15285-geodetic-toolbox
jd2yr <- function(jd){
    cal <- Jd2cal(jd)
    jd0  <-  Cal2jd(cbind(cal[,1],1,1))
    jd1  <-  Cal2jd(cbind(cal[,1]+1,1,1))
    if(!is.null(jd0)){
        return(cal[,1] + rowSums(jd-jd0)/(rowSums(jd1-jd0)))
    }else{
        return(cal[,1] + (jd-jd0)/(jd1-jd0))
    }
}
###https://www.mathworks.com/matlabcentral/fileexchange/15285-geodetic-toolbox
yr2jd <- function(yr){
##    jd1 <- Cal2jd(cbind(yr%/%1,0,0))
    ##    change.base(cbind(jd1[,1],jd1[,2]+(yr%%1)*DJY),denom=1)
    iyr <-  floor(yr)
    jd0  <-  Cal2jd(cbind(iyr,1,1))
    days  <-  rowSums(Cal2jd(cbind(iyr+1,1,1)) - jd0)
    doy  <-  (yr-iyr)*days + 1
    doy2jd(iyr,doy)
}

change.base <- function(t,denom=1){
    t1 <- (rowSums(t)%/%denom)*denom
    t2 <- (t[,1,drop=FALSE]-t1)+t[,2,drop=FALSE]
    return(cbind(t1,t2))
}

decompose.vector <- function(s){
#amplitude
    R <- rowSums(s[,1:3]^2)
    V <- rowSums(s[,4:6]^2)
#unit direction vector
    ur <- s[,1:3]/R
    uv <- s[,4:6]/V
    return(list(S=cbind(R,V),uS=cbind(ur,uv)))
}
