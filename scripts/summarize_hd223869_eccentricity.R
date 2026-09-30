args <- commandArgs(trailingOnly=TRUE)
outdir <- if(length(args)) args[1] else 'results/HD223869_20260919/eccentricity_mcmc'
fits <- readRDS(file.path(outdir,'fits.rds'))
chains <- lapply(fits,`[[`,'mc')
split_chains <- function(vs){
  n <- min(vapply(vs,length,integer(1))) %/% 2
  do.call(cbind,lapply(vs,function(v)cbind(head(v,n),tail(v,n))))
}
rhat <- function(z){
  n <- nrow(z); W <- mean(apply(z,2,var)); B <- n*var(colMeans(z))
  sqrt(((n-1)/n*W+B/n)/W)
}
rank_rhat <- function(z){
  rankz <- function(v) matrix(qnorm((rank(v,ties.method='average')-3/8)/(length(v)+1/4)),nrow=nrow(z))
  max(rhat(rankz(z)),rhat(rankz(abs(z-median(z)))))
}
ess <- function(v){
  ac <- as.numeric(acf(v,lag.max=min(2000,length(v)-1),plot=FALSE)$acf)
  pairs <- ac[seq(1,length(ac)-1,2)]+ac[seq(2,length(ac),2)]
  k <- which(pairs<0)[1]
  if(!is.na(k)) pairs <- head(pairs,k-1)
  length(v)/max(1,-1+2*sum(cummin(pairs)))
}
pars <- setdiff(colnames(chains[[1]]),c('omega1','Mo1','logpost','loglike'))
checks <- do.call(rbind,lapply(pars,function(p){
  vs <- lapply(chains,function(z)z[,p]); z<-split_chains(vs)
  data.frame(parameter=p,split_Rhat=rhat(z),rank_normalized_Rhat=rank_rhat(z),ESS_approx=sum(vapply(vs,ess,numeric(1))))
}))
write.csv(checks,file.path(outdir,'convergence.csv'),row.names=FALSE)
print(checks)
for(i in seq_along(chains)) {
  cat('Chain',i,'retained samples',nrow(chains[[i]]),'\n')
  print(quantile(chains[[i]][,'e1'],c(.025,.16,.5,.84,.95,.975)))
}
mc <- do.call(rbind,chains)
print(quantile(mc[,'e1'],c(.025,.16,.5,.84,.95,.975)))
cat('P interval:',quantile(exp(mc[,'per1']),c(.16,.5,.84)),'\n')
cat('e<0.1:',mean(mc[,'e1']<.1),'\n')
