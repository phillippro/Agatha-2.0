args <- commandArgs(trailingOnly=TRUE)
outdir <- args[1]; title <- args[2]
fits <- readRDS(file.path(outdir,'fits.rds'))
vs <- lapply(fits,function(f)f$mc[,'e1'])
e <- unlist(vs); qs <- quantile(e,c(.16,.5,.84))
png(file.path(outdir,'eccentricity_posterior.png'),width=1500,height=850,res=160)
par(mfrow=c(1,2),mar=c(4,4,3,1))
hist(e,breaks=50,freq=FALSE,xlab='Eccentricity',main=paste(title,'posterior'),col='lightblue',border='white')
abline(v=qs,lty=c(2,1,2),col='firebrick')
plot(vs[[1]],type='l',col='#0072B280',xlab='Post-burn-in draw',ylab='Eccentricity',
  main='Independent cold chains',ylim=c(0,max(e)*1.05),xlim=c(1,max(lengths(vs))))
lines(vs[[2]],col='#D55E0080')
legend('topright',c('Run 1','Run 2'),col=c('#0072B2','#D55E00'),lty=1,bty='n')
dev.off()
