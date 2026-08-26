source('mcmc_func.R')

expect_true <- function(value,message){
    if(!isTRUE(value)){
        stop(message,call.=FALSE)
    }
}

set.seed(101)

####################################################
## A two-mode target that a single-temperature chain cannot cross: the modes are
## far apart compared with their width, and the higher mode is the one the chain
## does not start in. This is the situation a periodogram peak leaves the MCMC in
## when an alias or a higher eccentricity solution fits better.
####################################################
mu1 <- -4
mu2 <- 4
sig <- 0.25
dlogL <- 6

par.min <- c(x=-10,y=-10)
par.max <- c(x=10,y=10)
Npar <- length(par.min)
Sd <- 2.4^2/Npar

loglikelihood.test <- function(param){
    x <- param[1]
    y <- param[2]
    l1 <- -0.5*((x-mu1)/sig)^2
    l2 <- -0.5*((x-mu2)/sig)^2+dlogL
    log(exp(l1)+exp(l2))-0.5*(y/2)^2
}

posterior <- function(param,tem=1,bases='natural'){
    llike <- loglikelihood.test(param)
    pr <- -log(prod(par.max-par.min))
    list(loglike=llike,logprior=pr,post=llike*tem+pr)
}

startvalue <- c(x=mu1,y=0)
cov.start <- diag(c(sig,1)^2)

Niter <- 20000
pt <- run.ptmcmc(startvalue,cov.start,iterations=Niter,Ntem=8,tem.min=1e-3,swap.interval=10)
single <- run.ptmcmc(startvalue,cov.start,iterations=Niter,tems=1)

frac.pt <- mean(pt$out[,'x']>0)
frac.single <- mean(single$out[,'x']>0)

expect_true(frac.pt>0.5,
            paste('parallel tempering spent only',round(100*frac.pt,1),
                  '% of the cold chain in the dominant mode'))
expect_true(frac.single<0.05,
            paste('the single-temperature control chain unexpectedly crossed the barrier (',
                  round(100*frac.single,1),'%); the test target is no longer separating the samplers'))

####the dominant mode carries exp(dlogL) times the mass of the other one
expected <- exp(dlogL)/(1+exp(dlogL))
expect_true(abs(frac.pt-expected)<0.1,
            paste('parallel tempering recovered a mode weight of',round(frac.pt,3),
                  'instead of',round(expected,3)))
expect_true(abs(mean(pt$out[pt$out[,'x']>0,'x'])-mu2)<0.1,
            'parallel tempering did not recover the location of the dominant mode')

####every replica should be near the target acceptance rate
expect_true(all(pt$acc.all>5 & pt$acc.all<60),
            paste('replica acceptance rates are out of range:',
                  paste(round(pt$acc.all,1),collapse=',')))
####neighbouring replicas have to actually exchange for tempering to do anything
expect_true(all(pt$swap.rate>0.01),
            paste('some neighbouring replicas never exchanged:',
                  paste(round(100*pt$swap.rate,1),collapse=',')))

####the ladder runs from the cold chain down to tem.min
expect_true(abs(pt$tems[1]-1)<1e-12 && abs(pt$tems[length(pt$tems)]-1e-3)<1e-12,
            'the tempering ladder does not span 1 to tem.min')
expect_true(all(diff(pt$tems)<0),'the tempering ladder is not monotonic')

####the returned chain must carry the columns the rest of the code expects
expect_true(all(c('x','y','logpost','loglike')%in%colnames(pt$out)),
            'run.ptmcmc did not return the expected chain columns')
expect_true(nrow(pt$out)==Niter-floor(Niter/2),
            'run.ptmcmc did not discard the burn-in half of the cold chain')

cat('PT-MCMC tests passed\n')
