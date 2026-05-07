dnint <- function(x){
    x[x>0] <- floor(x[x>0]+0.5)
    x[x<=0] <- floor(x[x<=0]-0.5)
    x
}

Tt2tcg <- function(tt){
    t77t <- DJM77 + TTMTAI/DAYSEC

    ## TT to TCG rate
    elgg <- ELG/(1.0-ELG)

    tt1 <- tt[,1]
    tt2 <- tt[,2]
    ##Result, safeguarding precision.
    tcg1 <- tt1
    tcg2 <- tt2 + ( ( tt1 - DJM0 ) + ( tt2 - t77t ) ) * elgg
    cbind(tcg1,tcg2)
}

Tcbtdb <- function(tcb){
## 1977 Jan 1 00:00:32.184 TT, as two-part JD
   t77td <- DJM0 + DJM77
   t77tf <- TTMTAI/DAYSEC

## TDB (days) at TAI 1977 Jan 1.0
   tdb0 <- TDB0/DAYSEC
   tcb1 <- tcb[,1]
   tcb2 <- tcb[,2]
##Result, safeguarding precision.
   d <- tcb1 - t77td
   tdb1 <- tcb1
   tdb2 <- tcb2 + tdb0 - ( d + ( tcb2 - t77tf ) ) * ELB
   cbind(tdb1,tdb2)
}

Tttdb <- function(tt){
##Result, safeguarding precision.
    tt1 <- tt[1]
    tt2 <- tt[2]
###only an approximation; absolute value should be calculated according to JPL ephemeris
    g <- 6.24 + 0.017202*(sum(tt)-2451545)
    dtr <- 0.001657*sin(g)#s
    dtrd = dtr / DAYSEC
    if ( tt1 > tt2 ) {
        tdb1 = tt1
        tdb2 = tt2 + dtrd
    } else {
        tdb1 <- tt1 + dtrd
        tdb2 <- tt2
    }
    c(tdb1,tdb2)
}

Tdbtcb <- function(tdb){
    ##1977 Jan 1 00:00:32.184 TT, as two-part JD
    t77td = DJM0 + DJM77
    t77tf = TTMTAI/DAYSEC

    ##TDB (days) at TAI 1977 Jan 1.0
    tdb0 = TDB0/DAYSEC

    ## TDB to TCB rate
    elbb = ELB/(1.0-ELB)

    tdb1 <- tdb[,1]
    tdb2 <- tdb[,2]
    ##Result, preserving date format but safeguarding precision.
    d <- t77td - tdb1
    f <- tdb2 - tdb0
    tcb1 <- tdb1
    tcb2 <- f - ( d - ( f - t77tf ) ) * elbb
    cbind(tcb1,tcb2)
}

#ref: https://uk.mathworks.com/matlabcentral/fileexchange/15285-geodetic-toolbox?focused=6513469&tab=function
#ref: Algorithm is from Fliegel and van Flandern (1968); also see julian.py
####Julian Date to Gregorian year, month, day, fraction
## dj1             dj2
##  2450123.7           0.0       (JD method)
##  2451545.0       -1421.3       (J2000 method)
##  2400000.5       50123.2       (MJD method)
##  2450123.5           0.2       (date & time method)
##  Returned (arguments):
##     iy        int      year
##     im        int      month
##     id        int      day
##     fd        double   fraction of day
Jd2cal <- function(jd){
## Minimum and maximum allowed JD */
   DJMIN = -68569.5
   DJMAX = 1e9

## Verify date is acceptable. */
   if(any(rowSums(jd) < DJMIN) | any(rowSums(jd) > DJMAX)) return(NULL)

## Copy the date, big then small, and re-align to midnight. */
   jd1 <- jd[,1]+jd[,2]-0.5+1
#   f <- jd1-jd
   f <- (jd[,1]-as.integer(jd1))+jd[,2]-0.5+1
   jd <- as.integer(jd1)

## Express day in Gregorian calendar.
   l <- jd + 68569
   n <- (4 * l) %/% 146097
   l <- l-((146097 * n + 3) %/% 4)
   i <- (4000 * (l + 1)) %/% 1461001
   l <- l-(1461 * i) %/% 4 + 31
   k <- (80 * l) %/% 2447
   id <- l - (2447 * k) %/% 80
   l <- k %/% 11
   im <- k + 2 - 12* l
   iy <- 100 * (n - 49) + i + l
   fd <- f
   return(cbind(iy,im,id,fd))
}

Jd2cal.hms <- function(jd){
    tmp <- Jd2cal(jd)
    fd <- tmp[,4]*DAYSEC
    h <- fd%/%3600
    m <- (fd-h*3600)%/%60
    s <- fd-h*3600-m*60
    return(cbind(tmp[,1:3,drop=FALSE],cbind(h,m,s)))
}

Cal2jd.hms <- function(cal){
    cal2 <- cal[,1:3,drop=FALSE]
    jd2 <- Cal2jd(cal2)
    fd <- cal[,4]/24+cal[,5]/(24*60)+cal[,6]/(24*3600)
    cbind(jd2[,1],jd2[,2]+fd)
}

#Gregorian Calendar to Julian Date.
#ref: https://uk.mathworks.com/matlabcentral/fileexchange/15285-geodetic-toolbox?focused=6513448&tab=function
Cal2jd <- function(cal){
                                        # CAL2JD  Converts calendar date to Julian date using algorithm
                                        #   from "Practical Ephemeris Calculations" by Oliver Montenbruck
                                        #   (Springer-Verlag, 1989). Uses astronomical year for B.C. dates
                                        #   (2 BC = -1 yr). Non-vectorized version. See also DOY2JD, GPS2JD,
                                        #   JD2CAL, JD2DOW, JD2DOY, JD2GPS, JD2YR, YR2JD.
                                        # Version: 2011-11-13
                                        # Usage:   jd=Cal2jd(yr,mn,dy)
                                        # Input:   yr - calendar year (4-digit including century)
                                        #          mn - calendar month
                                        #          dy - calendar day (including factional day)
                                        # Output:  jd - jJulian date
                                        # Copyright (c) 2011, Michael R. Craymer
                                        # All rights reserved.
                                        # EMAIL: mike@craymer.com
    y <- yr <- cal[,1]
    m <- mn <- cal[,2]
    dy <- cal[,3]
    ind <- which(mn<=2)
    y[ind] <- y[ind]-1
    m[ind] <- m[ind]+12
    date1 <- 4.5+31*(10+12*1582)# Last day of Julian calendar (1582.10.04 Noon)
    date2 <- 15.5+31*(10+12*1582)# First day of Gregorian calendar (1582.10.15 Noon)
    date <- dy+31*(mn+12*yr)
    ind1 <- which(date<=date1)
    ind2 <- which(date>=date2)
    b <- y
    b[ind1] <- -2
    b[ind2] <- trunc(y/400) - trunc(y/100)
    if(length(ind1)==0 & length(ind2)==0){
        cat('Dates between October 5 & 15, 1582 do not exist!\n')
    }
    ind1 <- which(y>0)
    ind2 <- which(y<0)
    jd <- y
    jd[ind1] <- trunc(365.25*y[ind1]) + trunc(30.6001*(m[ind1]+1)) + b[ind1] + 1720996.5 + dy[ind1]
    jd[ind2] <- trunc(365.25*y[ind2]-0.75) + trunc(30.6001*(m[ind2]+1)) + b[ind2] + 1720996.5 + dy[ind2]
#    return(cbind(DJM0,jd-DJM0))
    return(cbind(jd%/%1,jd%%1))
}

##Delta(AT) (=TAI-UTC) for a given UTC date
Dat <- function(g){
#####
# Given:
#      iy           UTC:  year (Notes 1 and 2)
#      im           month (Note 2)
#      id           day (Notes 2 and 3)
#      fd           fraction of day (Note 4)
#   Returned:
#      Dt       TAI minus UTC, seconds
##http://www.iausofa.org/2018_0130_C/sofa/dat.c
##Release year for this version of iauDat
    IYV <- 2017
    iy <- g[,1]
    im <- g[,2]
    id <- g[,3]
    fd <- g[,4]
##Reference dates (MJD) and drift rates (s/day), pre leap seconds
   drift <-rbind(
      c(37300.0, 0.0012960 ),
      c( 37300.0, 0.0012960 ),
      c( 37300.0, 0.0012960 ),
      c( 37665.0, 0.0011232 ),
      c( 37665.0, 0.0011232 ),
      c( 38761.0, 0.0012960 ),
      c( 38761.0, 0.0012960 ),
      c( 38761.0, 0.0012960 ),
      c( 38761.0, 0.0012960 ),
      c( 38761.0, 0.0012960 ),
      c( 38761.0, 0.0012960 ),
      c( 38761.0, 0.0012960 ),
      c( 39126.0, 0.0025920 ),
      c( 39126.0, 0.0025920 )
   )

##Number of Delta(AT) expressions before leap seconds were introduced */
   NERA1 <- nrow(drift)

##Number of Delta(AT) changes
   NDAT <- nrow(delat)

## Miscellaneous local variables
#   int j, i, m;
#   double da, djm0, djm;

##Initialize the result to zero.
   Dt <- rep(0,nrow(g))

##Convert the date into an MJD.
   jd <- Cal2jd(cbind(iy, im, id))
   djm0 <- jd[,1]
   djm <- jd[,2]

#If pre-UTC year, set warning status and give up.
   if(any(iy < delat[1,1])){
       cat('pre-UTC year ',delat[1,1],'!\n')
   }

# If suspiciously late year, set warning status but proceed.
   if(any(iy > IYV + 5)){
       cat('Warning: suspiciously late years after ',IYV + 5,'!\n')
   }

#Combine year and month to form a date-ordered integer...
#the original SOFA routine only use integer months, I add fraction months to account for the recent dates in the delat table.
   m <- 12*iy + im

# Get the Delta(AT)
    ms <- 12*delat[,1] + delat[,2]
    mmax <- max(ms)
    mmin <- min(ms)
    Dt[m>mmax] <- delat[nrow(delat),3]
    Dt[m<mmin] <- 0#or NA

#...and use it to find the preceding table entry.
    ind0 <- which(m>=mmin & m<=mmax)
    inds <- unlist(sapply(m[ind0], function(x) which(x<ms)[1]-1))
    inds[is.na(inds)] <- length(ms)
    Dt[ind0] <- delat[inds,3]

# If pre-1972, adjust for drift
# m: ind0; ms: inds
    ind1 <- which(m[ind0]<ms[NERA1+1] & m[ind0]>=mmin)
    if(length(ind1)>0){
        ind.m <- ind0[ind1]
        ind.ms <- inds[ind1]
        Dt[ind.m] <- Dt[ind.m] + (djm[ind.m] + fd[ind.m] - drift[ind.ms,1])*drift[ind.ms,2]
    }
   return(Dt)
}

##  Given:
##     tt1,tt2    double    TT as a 2-part Julian Date
##
##  Returned:
##     tai1,tai2  double    TAI as a 2-part Julian Date
Tttai <- function(utc,tt){
##tt
    if(TT.standard!='TAI'){
        dt <- bipm.corr(utc)
    }else{
        dt <- 0
    }
##TT minus TAI (days).
    dtat <- TTMTAI/DAYSEC

##Result, safeguarding precision.
    tai1 <- tt1
    tai2 <- tt2 - dtat-dt

    return(cbind(tai1,tai2))
}

###Taitt is a modified version of iauTaitt in SOFA to account for the BIPMXX realization of TT.
Taitt <- function(utc,tai){
    if(TT.standard!='TAI'){
        dt <- bipm.corr((utc[,1]-DJM0)+utc[,2])/DAYSEC
    }else{
        dt <- 0
    }
    dtat = TTMTAI/DAYSEC
#Result, safeguarding precision.
    cbind(tai[,1],tai[,2] + dtat+dt)
}

Utctai <- function(utc,tai.type='scale'){
#utc1/u1 should always be larger than utc2/u2
    u1 <- utc[,1]
    u2 <- utc[,2]
    j <- Jd2cal(utc)
    iy <- j[,1]
    im <- j[,2]
    id <- j[,3]
    fd <- j[,4]
    if(tai.type=='instant'){
        dat0 <- Dat(j[,1:4,drop=FALSE])
        dat.past <- Dat(cbind(j[,1:3,drop=FALSE],j[,4]-1e-6))
        dat.future <- Dat(cbind(j[,1:3,drop=FALSE],j[,4]+1e-6))
        dlod <- 2.0 * (dat0 - dat.past)
        dleap <- dat.future - (dat.past + dlod)
        fd <- fd+dleap/DAYSEC
        fd <- fd+dlod/DAYSEC
#        tai <- cbind(u1,fd)
    }else{
        ## Get TAI-UTC at 0h at given epochs.
        dat0 <- Dat(cbind(j[,1:3,drop=FALSE],t(t(rep(0,nrow(j))))))

        ##Get TAI-UTC at 12h (to detect drift).
        dat12 <- Dat(cbind(j[,1:3,drop=FALSE],t(t(rep(0.5,nrow(j))))))

###Get TAI-UTC at 0h tomorrow (to detect jumps).
###Note that the original sofa code at http://www.iausofa.org/2018_0130_C/sofa/utctai.c has a typo the following commented function would yield 12h tomorrow
        ##    j <- Jd2cal(cbind(u1+1.5, u2-fd))
        j <- Jd2cal(cbind(u1, u2-fd+1))
        iyt <- j[,1]
        imt <- j[,2]
        idt <- j[,3]
        w <- j[,4]
        dat24 <- Dat(cbind(j[,1:3,drop=FALSE],t(t(rep(0,nrow(j))))))

        ## Separate TAI-UTC change into per-day (DLOD) and any jump (DLEAP). */
        dlod <- 2.0 * (dat12 - dat0)
        dleap <- dat24 - (dat0 + dlod)

        ##Remove any scaling applied to spread leap into preceding day.
        fd <- fd*(DAYSEC+dleap)/DAYSEC

        ## Scale from (pre-1972) UTC seconds to SI seconds.
        fd <- fd*(DAYSEC+dlod)/DAYSEC
    }
    ##given epoch calendar date to 2-part JD.
    tmp <- Cal2jd(cbind(iy, im, id))
    z1 <- tmp[,1]
    z2 <- tmp[,2]

    ## Assemble the TAI result, preserving the UTC split and order.
    a2 <- z1 - u1
    a2 <- a2+z2
    a2 <- a2+fd + dat0/DAYSEC
    tai <- cbind(u1,a2)
    return(list(tai=tai,leap=dleap))
}

#Taiut1
##  Given:
##    tai1,tai2  double    TAI as a 2-part Julian Date
##     dta        double    UT1-TAI in seconds
##  Returned:
##     ut11,ut12  double    UT1 as a 2-part Julian Date
##
##  Returned (function value):
##                int       status:  0 = OK
Taiut1 <- function(tai1,tai2,dta){
   dtad <- dta/DAYSEC
   ut11 <- tai1
   ut12 <- tai2 + dtad
   return(cbind(ut11,ut12))
}

##Given:
##     tai1,tai2  double   TAI as a 2-part Julian Date (Note 1)
##  Returned:
##     utc1,utc2  double   UTC as a 2-part quasi Julian Date (Notes 1-3)
Taiutc <- function(tai){
## Put the two parts of the TAI into big-first order.
    a1 <- tai[,1]
    a2 <- tai[,2]
##Initial guess for UTC.
   u1 <- a1
   u2 <- a2

##Iterate (though in most cases just once is enough).
   for(i in 1:3){
##Guessed UTC to TAI.
      g <- Utctai(cbind(u1, u2))
##Adjust guessed UTC.
      u2 <- u2+a1 - g[,1]
      u2 <- u2+a2 - g[,2]
   }

##Return the UTC result, preserving the TAI order.
      utc1 <- u1
      utc2 <- u2
##return
   return(c(u1,u2))
}

##http://maia.usno.navy.mil/ser7/ser7.dat
##
pred.dut1 <- function(MJD){
    A <- 2*pi*(MJD-58367)/365.25
    C <- 2*pi*(MJD-58367)/435
##ref. http://www.iausofa.org/2018_0130_C/sofa/epb.c
    T <- 1900.0 + ((DJM0-DJ00)+(MJD+36524.68648))/DTY
##ref. https://datacenter.iers.org/eop/-/somos/5Rgv/getTX/6/bulletina-xxxi-036.txt
    dUT2.UT1 <- 0.022*sin(2*pi*T) - 0.012*cos(2*pi*T) -0.006*sin(4*pi*T) + 0.007*cos(4*pi*T)
    x <- 0.1156 + 0.1029*cos(A) - 0.0276*sin(A) - 0.0048*cos(C) + 0.0217*sin(C)
    y <- 0.3554 - 0.0209*cos(A) - 0.0992*sin(A) + 0.0217*cos(C) + 0.0048*sin(C)
    dut1 <- 0.0239 - 0.00067*(MJD - 58375) - dUT2.UT1
    return(dut1)
}
###UTC to UT1
### Given:
##    utc1,utc2  double   UTC as a 2-part quasi Julian Date (Notes 1-4)
##     dut1       double   Delta UT1 = UT1-UTC in seconds (Note 5)
##
##  Returned:
##     ut11,ut12  double   UT1 as a 2-part Julian Date (Note 6)
Utcut1 <- function(utc,tempo=TRUE){
    utc1 <- utc[,1]
    utc2 <- utc[,2]
    ##Look up TAI-UTC.
    tmp <- Jd2cal(utc)
    iy <- tmp[,1]
    im <- tmp[,2]
    id <- tmp[,3]
    w <- tmp[,4]
    dat <- Dat(cbind(iy, im, id, 0.0))

    ## Form UT1-TAI.
    mjd <- (utc1-DJM0)+utc2
    if(tempo){
        dut1 <- calDut1dot(mjd)$dut1
    }else{
        dut1 <- find.dut(mjd)
    }
    dta <- dut1 - dat

    ## UTC to TAI to UT1. */
    tai <- Utctai(utc)
    Taiut1(tai[,1], tai[,2], dta)
}

Eform <- function(n){
## Given:
##     n    int         ellipsoid identifier (Note 1)
##
##  Returned:
##     a    double      equatorial radius (meters, Note 2)
##     f    double      flattening (Note 2)
##
##  n    ellipsoid
##        1     WGS84
##        2     GRS80
##        3     WGS72
    if(n==1){
        a <- 6378137.0
        f <- 1.0 / 298.257223563
    }

    if(n==2){
        a = 6378137.0
        f = 1.0 / 298.257222101
    }

    if(n==3){
        a = 6378135.0
        f = 1.0 / 298.26
    }
    return(c(a,f))
   }

Gd2gce <- function(a,f,llh){
##  Given:
##     a       double     equatorial radius (Notes 1,4)
##     f       double     flattening (Notes 2,4)
##     elong   double     longitude (radians, east +ve)
##     phi     double     latitude (geodetic, radians, Note 4)
##     height  double     height above ellipsoid (geodetic, Notes 3,4)
##
##  Returned:
##     xyz     double[3]  geocentric vector (Note 3)
    if(is.null(dim(llh))) llh <- t(llh)
    elong <- llh[,1]
    phi <- llh[,2]
    height <- llh[,3]
    sp = sin(phi)
   cp = cos(phi)
   w = 1.0 - f
   w = w * w
   d = cp*cp + w*sp*sp
   if( d <= 0.0 ) return(NULL)
   ac = a / sqrt(d)
   as = w * ac

## Geocentric vector.
    r = (ac + height)*cp
    x = r * cos(elong)
    y = r * sin(elong)
    z = (as + height)*sp
    return(cbind(x,y,z))
}

Gd2gc <- function(n, llh){
#  Given:
#     n            ellipsoid identifier (Note 1)
#     elong        longitude (radians, east +ve)
#     phi          latitude (geodetic, radians, Note 3)
#     height       height above ellipsoid (geodetic, Notes 2,3)
#
#  Returned:
#     xyz          geocentric vector (Note 2)
#
# Obtain reference ellipsoid parameters.
    af <- Eform(n)

# transform longitude, geodetic latitude, height to x,y,z.
    Gd2gce(af[1], af[2], llh)
}

Tdbtcb <- function(tdb){
    #1977 Jan 1 00:00:32.184 TT, as two-part JD */
    tdb1 <- tdb[,1]
    tdb2 <- tdb[,2]
    t77td = DJM0 + DJM77
    t77tf = TTMTAI/DAYSEC

    ##TDB (days) at TAI 1977 Jan 1.0 */
    tdb0 = TDB0/DAYSEC

    ##TDB to TCB rate
    elbb = ELB/(1.0-ELB)

    ##Result, preserving date format but safeguarding precision.
    d <- t77td - tdb1
    f <- tdb2 - tdb0
    tcb1 <- tdb1
    tcb2 <- f - ( d - ( f - t77tf ) ) * elbb
    return(cbind(tcb1,tcb2))
}
Fal03 <- function(t){
          ((485868.249036  +
             t * ( 1717915923.2178 +
             t * (         31.8792 +
             t * (          0.051635 +
             t * (        -0.00024470 ) ) ) ))%%TURNAS ) * DAS2R
}

Falp03 <- function(t){
        ( ( 1287104.793048 +
             t * ( 129596581.0481 +
             t * (       -0.5532 +
             t * (         0.000136 +
             t * (       -0.00001149 ) ) ) ))%%TURNAS)*DAS2R
}

Fad03 <- function(t){
          ((1072260.703692 +
             t * ( 1602961601.2090 +
             t * (        - 6.3706 +
             t * (          0.006593 +
             t * (        - 0.00003169 ) ) ) ))%%TURNAS)*DAS2R
}

Faf03 <- function(t){
    ((335779.526232 +
             t * ( 1739527262.8478 +
             t * (       - 12.7512 +
             t * (        - 0.001037 +
             t * (          0.00000417 ) ) ) ))%%TURNAS)*DAS2R
}

Fave03 <- function(t){
    (3.176146697 + 1021.3285546211 * t)%%D2PI
}

Faom03 <- function(t){
    (( 450160.398036 +
             t * ( - 6962890.5431 +
             t * (         7.4722 +
             t * (         0.007702 +
             t * (       - 0.00005939 ) ) ) ))%%TURNAS ) * DAS2R;
}

Fae03 <- function(t){
    (1.753470314 + 628.3075849991 * t)%%D2PI
}

Fapa03 <- function(t){
    (0.024381750 + 0.00000538691 * t) * t
}
##
S06 <- function(tt,x,y){
    ## Given:
    ##     tt   double    TT as a 2-part Julian Date (Note 1)
    ##     x,y           double    CIP coordinates (Note 3)
    ##  Returned (function value):
    ##                   double    the CIO locator s in radians (Note 2)
    sp <-   c(94.00e-6,  3808.65e-6,  -122.68e-6,-72574.11e-6, 27.98e-6,15.62e-6)
    s0 <-  rbind(
        c( 0,  0,  0,  0,  1,  0,  0,  0, -2640.73e-6,   0.39e-6 ),
        c( 0,  0,  0,  0,  2,  0,  0,  0,   -63.53e-6,   0.02e-6 ),
        c( 0,  0,  2, -2,  3,  0,  0,  0,   -11.75e-6,  -0.01e-6 ),
        c( 0,  0,  2, -2,  1,  0,  0,  0,   -11.21e-6,  -0.01e-6 ),
        c( 0,  0,  2, -2,  2,  0,  0,  0,     4.57e-6,   0.00e-6 ),
        c( 0,  0,  2,  0,  3,  0,  0,  0,    -2.02e-6,   0.00e-6 ),
        c( 0,  0,  2,  0,  1,  0,  0,  0,    -1.98e-6,   0.00e-6 ),
        c( 0,  0,  0,  0,  3,  0,  0,  0,     1.72e-6,   0.00e-6 ),
        c( 0,  1,  0,  0,  1,  0,  0,  0,     1.41e-6,   0.01e-6 ),
        c( 0,  1,  0,  0, -1,  0,  0,  0,     1.26e-6,   0.01e-6 ),
        c( 1,  0,  0,  0, -1,  0,  0,  0,     0.63e-6,   0.00e-6 ),
        c( 1,  0,  0,  0,  1,  0,  0,  0,     0.63e-6,   0.00e-6 ),
        c( 0,  1,  2, -2,  3,  0,  0,  0,    -0.46e-6,   0.00e-6 ),
        c( 0,  1,  2, -2,  1,  0,  0,  0,    -0.45e-6,   0.00e-6 ),
        c( 0,  0,  4, -4,  4,  0,  0,  0,    -0.36e-6,   0.00e-6 ),
        c( 0,  0,  1, -1,  1, -8, 12,  0,     0.24e-6,   0.12e-6 ),
        c( 0,  0,  2,  0,  0,  0,  0,  0,    -0.32e-6,   0.00e-6 ),
        c( 0,  0,  2,  0,  2,  0,  0,  0,    -0.28e-6,   0.00e-6 ),
        c( 1,  0,  2,  0,  3,  0,  0,  0,    -0.27e-6,   0.00e-6 ),
        c( 1,  0,  2,  0,  1,  0,  0,  0,    -0.26e-6,   0.00e-6 ),
        c( 0,  0,  2, -2,  0,  0,  0,  0,     0.21e-6,   0.00e-6 ),
        c( 0,  1, -2,  2, -3,  0,  0,  0,    -0.19e-6,   0.00e-6 ),
        c( 0,  1, -2,  2, -1,  0,  0,  0,    -0.18e-6,   0.00e-6 ),
        c( 0,  0,  0,  0,  0,  8,-13, -1,     0.10e-6,  -0.05e-6 ),
        c( 0,  0,  0,  2,  0,  0,  0,  0,    -0.15e-6,   0.00e-6 ),
        c( 2,  0, -2,  0, -1,  0,  0,  0,     0.14e-6,   0.00e-6 ),
        c( 0,  1,  2, -2,  2,  0,  0,  0,     0.14e-6,   0.00e-6 ),
        c( 1,  0,  0, -2,  1,  0,  0,  0,    -0.14e-6,   0.00e-6 ),
        c( 1,  0,  0, -2, -1,  0,  0,  0,    -0.14e-6,   0.00e-6 ),
        c( 0,  0,  4, -2,  4,  0,  0,  0,    -0.13e-6,   0.00e-6 ),
        c( 0,  0,  2, -2,  4,  0,  0,  0,     0.11e-6,   0.00e-6 ),
        c( 1,  0, -2,  0, -3,  0,  0,  0,    -0.11e-6,   0.00e-6 ),
        c( 1,  0, -2,  0, -1,  0,  0,  0,    -0.11e-6,   0.00e-6 ))

    s1 <- rbind(c( 0,  0,  0,  0,  2,  0,  0,  0,    -0.07e-6,   3.57e-6 ),
                c( 0,  0,  0,  0,  1,  0,  0,  0, 1.73e-6,  -0.03e-6 ),
                c( 0,  0,  2, -2,  3,  0,  0,  0,     0.00e-6,   0.48e-6 ))

    s2 <- rbind(
        c( 0,  0,  0,  0,  1,  0,  0,  0,   743.52e-6,  -0.17e-6 ),
        c( 0,  0,  2, -2,  2,  0,  0,  0,    56.91e-6,   0.06e-6 ),
        c( 0,  0,  2,  0,  2,  0,  0,  0,     9.84e-6,  -0.01e-6 ),
        c( 0,  0,  0,  0,  2,  0,  0,  0,    -8.85e-6,   0.01e-6 ),
        c( 0,  1,  0,  0,  0,  0,  0,  0,    -6.38e-6,  -0.05e-6 ),
        c( 1,  0,  0,  0,  0,  0,  0,  0,    -3.07e-6,   0.00e-6 ),
        c( 0,  1,  2, -2,  2,  0,  0,  0,     2.23e-6,   0.00e-6 ),
        c( 0,  0,  2,  0,  1,  0,  0,  0,     1.67e-6,   0.00e-6 ),
        c( 1,  0,  2,  0,  2,  0,  0,  0,     1.30e-6,   0.00e-6 ),
        c( 0,  1, -2,  2, -2,  0,  0,  0,     0.93e-6,   0.00e-6 ),
        c( 1,  0,  0, -2,  0,  0,  0,  0,     0.68e-6,   0.00e-6 ),
        c( 0,  0,  2, -2,  1,  0,  0,  0,    -0.55e-6,   0.00e-6 ),
        c( 1,  0, -2,  0, -2,  0,  0,  0,     0.53e-6,   0.00e-6 ),
        c( 0,  0,  0,  2,  0,  0,  0,  0,    -0.27e-6,   0.00e-6 ),
        c( 1,  0,  0,  0,  1,  0,  0,  0,    -0.27e-6,   0.00e-6 ),
        c( 1,  0, -2, -2, -2,  0,  0,  0,    -0.26e-6,   0.00e-6 ),
        c( 1,  0,  0,  0, -1,  0,  0,  0,    -0.25e-6,   0.00e-6 ),
        c( 1,  0,  2,  0,  1,  0,  0,  0,     0.22e-6,   0.00e-6 ),
        c( 2,  0,  0, -2,  0,  0,  0,  0,    -0.21e-6,   0.00e-6 ),
        c( 2,  0, -2,  0, -1,  0,  0,  0,     0.20e-6,   0.00e-6 ),
        c( 0,  0,  2,  2,  2,  0,  0,  0,     0.17e-6,   0.00e-6 ),
        c( 2,  0,  2,  0,  2,  0,  0,  0,     0.13e-6,   0.00e-6 ),
        c( 2,  0,  0,  0,  0,  0,  0,  0,    -0.13e-6,   0.00e-6 ),
        c( 1,  0,  2, -2,  2,  0,  0,  0,    -0.12e-6,   0.00e-6 ),
        c( 0,  0,  2,  0,  0,  0,  0,  0,    -0.11e-6,   0.00e-6 ))

    s3 <- rbind(c( 0,  0,  0,  0,  1,  0,  0,  0,     0.30e-6, -23.42e-6 ),
                c( 0,  0,  2, -2,  2,  0,  0,  0,    -0.03e-6,  -1.46e-6 ),
                c( 0,  0,  2,  0,  2,  0,  0,  0,    -0.01e-6,  -0.25e-6 ),
                c( 0,  0,  0,  0,  2,  0,  0,  0,     0.00e-6,   0.23e-6 ))

    s4 <- rbind(c( 0,  0,  0,  0,  1,  0,  0,  0,    -0.26e-6,  -0.01e-6 ))

    NS0 <- nrow(s0)
    NS1 <- nrow(s1)
    NS2 <- nrow(s2)
    NS3 <- nrow(s3)
    NS4 <- nrow(s4)

    t <- ((tt[,1] - DJ00) + tt[,2]) / DJC
    ## Fundamental Arguments (from IERS Conventions 2003) */

    fa <- array(NA,dim=c(length(t),8))
    ## Mean anomaly of the Moon. */
    fa[,1] = Fal03(t)

    ## Mean anomaly of the Sun. */
    fa[,2] = Falp03(t)

    ## Mean longitude of the Moon minus that of the ascending node. */
    fa[,3] = Faf03(t)

    ## Mean elongation of the Moon from the Sun. */
    fa[,4] = Fad03(t)

    ## Mean longitude of the ascending node of the Moon. */
    fa[,5] = Faom03(t)

    ## Mean longitude of Venus. */
    fa[,6] = Fave03(t)

    ## Mean longitude of Earth. */
    fa[,7] = Fae03(t)

    ##General precession in longitude. */
    fa[,8] = Fapa03(t)

    ## Evaluate s. */
    w0 <- sp[1]
    w1 <- sp[2]
    w2 <- sp[3]
    w3 <- sp[4]
    w4 <- sp[5]
    w5 <- sp[6]

    for(i in NS0:1){
        a <- 0
        for(j in 1:8){
            a <- a+ s0[i,j]*fa[,j]
        }
        w0 <- w0+s0[i,9]*sin(a) + s0[i,10] * cos(a)
    }

    for (i in NS1:1){
        a <- 0.0
        for (j in 1:8) {
            a <- a + s1[i,j]*fa[,j]
        }
        w1 <- w1+s1[i,9] * sin(a) + s1[i,10]*cos(a)
    }

    for (i in NS2:1){
        a <- 0.0
        for(j in 1:8){
            a <- a+s2[i,j] * fa[,j]
        }
        w2 <- w2+s2[i,9]*sin(a) + s2[i,10]*cos(a)
    }

    for (i in NS3:1){
        a <- 0.0
        for(j in 1:8){
            a <- a+s3[i,j]*fa[,j]
        }
        w3 <- w3+s3[i,9]*sin(a)+s3[i,10] * cos(a);
    }

    for (i in NS4:1){
        a <- 0.0
        for(j in 1:8){
            a <- a+s4[i,j]*fa[,j]
        }
        w4 <- w4+s4[i,9]*sin(a) + s4[i,10]*cos(a)
    }
    s <- (w0 + (w1 + (w2 + (w3 + (w4 + w5 * t) * t) * t) * t) * t)*DAS2R - x*y/2.0
    return(s)
}

Pom00 <- function(xp, yp, sp){
##  Given:
##     xp,yp    double    coordinates of the pole (radians, Note 1)
##     sp       double    the TIO locator s' (radians, Note 2)
##
##  Returned:
##     rpom     double[3][3]   polar-motion matrix (Note 3)
#Construct the matrix W in https://www.iers.org/SharedDocs/Publikationen/EN/IERS/Publications/tn/TechnNote32/tn32.pdf?__blob=publicationFile&v=1
    rpom <- array(dim=c(3,3,length(xp)))
    for(j in 1:length(xp)){
        rpom[,,j] <- R1(-yp[j])%*%R2(-xp[j])%*%R3(sp[j])
    }
    return(rpom)
}

Era00 <- function(ut){
## Days since fundamental epoch.
   t = (ut[,1]-DJ00) + ut[,2]
## Fractional part of T (days).
   f = rowSums(ut)%%1

##Earth rotation angle at this UT1.
   theta <- 2*pi* ((f + 0.7790572732640 + 0.00273781191135448 * t)%%1)
   return(theta)
}
C2t00b <- function(tt, ut, xp, yp,rpom){
    ##  Given:
    ##     tta,ttb  double         TT as a 2-part Julian Date (Note 1)
    ##     uta,utb  double         UT1 as a 2-part Julian Date (Note 1)
    ##     xp,yp    double         coordinates of the pole (radians, Note 2)
    ##
    ##  Returned:
    ##     rc2t     double[3][3]   celestial-to-terrestrial matrix (Note 3)
    ## Form the celestial-to-intermediate matrix for this TT (IAU 2000B).
    rc2i <- C2i00b(tt)
    ## Predict the Earth rotation angle for this UT1.
    era <- Era00(ut)
    ## Form the polar motion matrix (neglecting s').
#    rpom <- Pom00(xp, yp, rep(0,length(xp)))
    ## Combine to form the celestial-to-terrestrial matrix.
    rc2t <- array(NA,dim=c(3,3,nrow(tt)))
    for(j in 1:length(xp)){
        rc2t[,,j] <- rpom[,,j]%*%R3(era[j])%*%rc2i[,,j]
    }
    rc2t
#    return(list(Mcio=rc2i,Mc2t=rc2t,era=era,rpom=rpom))
}

##
C2i00a <- function(tt){
##  Given:
##     date1,date2 double       TT as a 2-part Julian Date (Note 1)
##
##  Returned:
##     rc2i        double[3][3] celestial-to-intermediate matrix (Note 2)
    rbpn <- Pnm00a(tt)
    rc2i <- C2ibpn(tt, rbpn)
}

C2ibpn <- function(tt,rbpn){
##Extract the X,Y coordinates.
    x <- rbpn[3,1,]
    y <- rbpn[3,2,]
##Form the celestial-to-intermediate matrix (n.b. IAU 2000 specific).
    rc2i <- C2ixy(tt, x, y)
    return(rc2i)
}

Pn00b <- function(tt){
##  Given:
##     date1,date2  double          TT as a 2-part Julian Date (Note 1)
##
##  Returned:
##     dpsi,deps    double          nutation (Note 2)
##     epsa         double          mean obliquity (Note 3)
##     rb           double[3][3]    frame bias matrix (Note 4)
##     rp           double[3][3]    precession matrix (Note 5)
##     rbp          double[3][3]    bias-precession matrix (Note 6)
##     rn           double[3][3]    nutation matrix (Note 7)
##     rbpn         double[3][3]    GCRS-to-true matrix (Notes 8,9)

##Nutation.
    tmp1 <- Nut00b(tt)

##Remaining results.
    tmp2 <- Pn00(tt,tmp1$dpsi,tmp1$deps)
    return(c(tmp1,tmp2))
}

C2i00b <- function(tt){
##Obtain the celestial-to-true matrix (IAU 2000B).
#   rbpn <- Pnm00b(tt)
   rbpn <- Pn00b(tt)$rbpn

##Form the celestial-to-intermediate matrix.
   rc2i <- C2ibpn(tt, rbpn)

   return(rc2i)
}
Bi00 <- function(){
    DPBIAS = -0.041775  * DAS2R
    DEBIAS = -0.0068192 * DAS2R
    DRA0 = -0.0146 * DAS2R
    dpsibi = DPBIAS
    depsbi = DEBIAS
    dra = DRA0
    return(c(dpsibi, depsbi, dra))
}

Bp00 <- function(tt){
##J2000.0 obliquity (Lieske et al. 1977)
   EPS0 = 84381.448 * DAS2R
   Nt <- nrow(tt)

##Interval between fundamental epoch J2000.0 and current date (JC).
   t = ((tt[,1] - DJ00) + tt[,2]) / DJC

##Frame bias.
   tmp <- Bi00()
   dpsibi <- tmp[1]
   depsbi <- tmp[2]
   dra0 <- tmp[3]

##Precession angles (Lieske et al. 1977)
   psia77 = (5038.7784 + (-1.07259 + (-0.001147) * t) * t) * t * DAS2R
   oma77  =       EPS0 + ((0.05127 + (-0.007726) * t) * t) * t * DAS2R
   chia   = (  10.5526 + (-2.38064 + (-0.001125) * t) * t) * t * DAS2R

##Apply IAU 2000 precession corrections.
   out <- Pr00(tt)
   psia = psia77 + out$dpsipr
   oma  = oma77  + out$depspr

   rb <- rp <- rbp <- array(NA,dim=c(3,3,Nt))
   for(j in 1:Nt){
       ##Frame bias matrix: GCRS to J2000.0.
       rbw <- R1(-depsbi)%*%R2(dpsibi*sin(EPS0))%*%R3(dra0)
       rb[,,j] <- rbw

       ##Precession matrix: J2000.0 to mean of date.
       rp[,,j] <- R3(chia[j])%*%R1(-oma[j])%*%R3(-psia[j])%*%R1(EPS0)

       ##Bias-precession matrix: GCRS to mean of date.
       rbp[,,j] <- rp[,,j]%*%rbw
   }
##return
   return(list(rb=rb,rp=rp,rbp=rbp))
}

Pr00 <- function(tt){
    ##Precession and obliquity corrections (radians per century)
    PRECOR = -0.29965 * DAS2R
    OBLCOR = -0.02524 * DAS2R

    ##Interval between fundamental epoch J2000.0 and given date (JC).
    t = ((tt[,1] - DJ00) + tt[,2]) / DJC

    ##Precession rate contributions with respect to IAU 1976/80.
    dpsipr <- PRECOR * t
    depspr <- OBLCOR * t
    return(list(dpsipr=dpsipr,depspr=depspr))
}

Pn00 <- function(tt,dpsi,deps){
    Nt <- nrow(tt)
    ## IAU 2000 precession-rate adjustments. , &dpsipr, &depspr)
    tmp <- Pr00(tt)
    dpsipr <- tmp$dpsipr
    depspr <- tmp$depspr

    ##Mean obliquity, consistent with IAU 2000 precession-nutation.
    epsa <- Obl80(tt) + depspr

    ##Frame bias and precession matrices and their product.
    rbpw <- Bp00(tt)$rbp
    rb <- Bp00(tt)$rb
    rp <- Bp00(tt)$rp
    rbp <- rbpw

    ##Nutation matrix.
    rn <- rbpn <- array(NA,dim=c(3,3,Nt))
    for(j in 1:Nt){
        rnw <- R1(-(epsa[j]+deps[j]))%*%R3(-dpsi[j])%*%R1(epsa[j])
        rn[,,j] <- rnw
        ##Bias-precession-nutation matrix (classical).
        rbpn[,,j] <- rnw%*%rbpw[,,j]
    }
    return(list(rb=rb,rp=rp,rbp=rbp,rn=rn,rbpn=rbpn))
}

Obl80 <- function(tt){
##Interval between fundamental epoch J2000.0 and given date (JC).
   t = ((tt[,1] - DJ00) + tt[,2]) / DJC

##Mean obliquity of date.
   eps0 = DAS2R * (84381.448  +
                  (-46.8150   +
                  (-0.00059   +
                  ( 0.001813) * t) * t) * t)
return(eps0)
}

Pnm06a <- function(tt){
##Fukushima-Williams angles for frame bias and precession.
   gppe <- Pfw06(tt)

##Nutation components.
   dpe <- Nut06a(tt)

##Equinox based nutation x precession x bias matrix. return rnpb
   R1(-gppe[4])%*%R3(-(gppe[3] + dpe[1]))%*%R1(gppe[2])%*%R3(gppe[1])
}

Pfw06 <- function(tt){
    t <- ((tt[1] - DJ00) + tt[2]) / DJC

    ##P03 bias+precession angles.
    gamb <- (    -0.052928     +
                 (    10.556378     +
                      (     0.4932044    +
                           (    -0.00031238   +
                                (    -0.000002788  +
                                     (     0.0000000260 )
                                 * t) * t) * t) * t) * t) * DAS2R
    phib <- ( 84381.412819     +
                 (   -46.811016     +
                      (     0.0511268    +
                           (     0.00053289   +
                                (    -0.000000440  +
                                     (    -0.0000000176 )
                                 * t) * t) * t) * t) * t) * DAS2R
    psib <- (    -0.041775     +
                 (  5038.481484     +
                      (     1.5584175    +
                           (    -0.00018522   +
                                (    -0.000026452  +
                                     (    -0.0000000148 )
                                 * t) * t) * t) * t) * t) * DAS2R
    epsa <- iauObl06(date1, date2)
    return(c(gamb,phib,psib,epsa))
}

##calculate classical NPB matrix
Pnm00a <- function(tt){
    tmp <- Nut00a(tt)/DAS2R
    Pn00a(tt)
}

#Nut06a
Nut06a <- function(tt){
    date1 <- tt[,1]
    date2 <- tt[,2]
    ##Interval between fundamental date J2000.0 and given date (JC).
    t = ((date1 - DJ00) + date2) / DJC;

    ##Factor correcting for secular variation of J2.
    fj2 = -2.7774e-6 * t;

    ##Obtain IAU 2000A nutation.
    dpe <- iauNut00a(tt)

    ##Apply P03 adjustments (Wallace & Capitaine, 2006, Eqs.5).
    dpsi <- dp + dp * (0.4697e-6 + fj2)
    deps <- de + de * fj2
    return(c(epsi,deps))
}

Gc2gd <- function(n,xyz){
## Given:
##     n       int        ellipsoid identifier (Note 1)
##     xyz     double[3]  geocentric vector (Note 2)
##
##  Returned:
##     elong   double     longitude (radians, east +ve, Note 3)
##     phi     double     latitude (geodetic, radians, Note 3)
##     height  double     height above ellipsoid (geodetic, Notes 2,3)
    af <- Eform( n )
    Gc2gde ( af[1], af[2], xyz)
}


Gc2gde <- function(a,f,xyz){
##Functions of ellipsoid parameters (with further validation of f).
   aeps2 = a*a * 1e-32
   e2 = (2.0 - f) * f
   e4t = e2*e2 * 1.5
   ec2 = 1.0 - e2
   if ( ec2 <= 0.0 ) return(-1)
   ec = sqrt(ec2)
   b = a * ec

##Cartesian components.
   x = xyz[1]
   y = xyz[2]
   z = xyz[3]

##Distance from polar axis squared.
   p2 <- x*x + y*y

##Longitude.
   elong <- atan2(y, x)

##Unsigned z-coordinate.
   absz <- abs(z)

##Proceed unless polar case.
   if ( p2 > aeps2 ){
##Distance from polar axis.
      p <- sqrt(p2)
##Normalization.
      s0 = absz/a
      pn = p/a
      zc = ec*s0

##Prepare Newton correction factors.
      c0 = ec * pn
      c02 = c0 * c0
      c03 = c02 * c0
      s02 = s0 * s0
      s03 = s02 * s0
      a02 = c02 + s02
      a0 = sqrt(a02)
      a03 = a02 * a0
      d0 = zc*a03 + e2*s03
      f0 = pn*a03 - e2*c03

##Prepare Halley correction factor.
      b0 <- e4t * s02 * c02 * pn * (a0 - ec)
      s1 <- d0*f0 - b0*s0
      cc <- ec * (f0*f0 - b0*c0)

##Evaluate latitude and height.
      phi <- atan(s1/cc)
      s12 <- s1 * s1
      cc2 <- cc * cc
      height <- (p*cc + absz*s1 - a * sqrt(ec2*s12 + cc2)) /sqrt(s12 + cc2)
   }else{
##Exception: pole
      phi <- DPI / 2.0
      height <- absz - b
   }

##Restore sign of latitude.
   if ( z < 0 ) phi <- -phi
   return(c(long=elong,lat=phi,h=height))
}

C2ixy <- function(tt,x,y){
    Nt <- nrow(tt)
    rc2i <- array(NA,c(3,3,Nt))
    s <- S06(tt, x, y)
    ## Obtain the spherical angles E and d. */
    r2 <- x*x + y*y
    e <- atan2(y, x)
    d <- atan(sqrt(r2/(1.0 - r2)))

    ## Form the matrix.
    for(j in 1:Nt){
        rc2i[,,j] <- R3(-(e[j]+s[j]))%*%R2(d[j])%*%R3(e[j])
    }
    return(rc2i)
}

