library(foreach)
library(doParallel)
registerDoParallel(cores=8)

alphab = function(p, e1, e2)
	plogis(p[1] + p[2]*e1 + p[3]*e2 + p[4]*e1^2 + p[5]*e2^2)

alphat = function(p, e1, e2)
	plogis(p[6] + p[7]*e1 + p[8]*e2 + p[9]*e1^2 + p[10]*e2^2)

rtob = function(pars, e1, e2, T, B, M)
	alphab(pars, e1, e2) * (M + B) * (1 - alphat(pars, e1, e2) * (T + M))
	
rtot = function(pars, e1, e2, T, B, M)
	alphat(pars, e1, e2) * (T + M) * (1 - alphab(pars, e1, e2) * (B + M))
	
rtom = function(pars, e1, e2, T, B, M)
	alphab(pars, e1, e2) * (M + B) * alphat(pars, e1, e2) * (T + M)

rtor = function(pars, e1, e2, T, B, M)
	1 - rtot(pars, e1, e2, T, B, M) - rtob(pars, e1, e2, T, B, M) - rtom(pars, e1, e2, T, B, M)
	
inits1 = readRDS("dat/inits_r1.rds")
oldSamples = readRDS("res/old results/full/samples_full.rds")
oldSamples = oldSamples[seq(1, nrow(oldSamples), length.out=1000),]
## params=as.numeric(inits1$full[2,])


temp = seq(-3, 3, 0.01)
precip = 0
TT = c(0.90, 0.50, 0.05, 0.00)
MM = c(0.10, 0.45, 0.45, 0.10)
BB = c(0.00, 0.05, 0.50, 0.90)
 
lwd=2
par(mfrow=c(2,2))

cols=c(T='#1b9e77', M = '#d95f02', B='#7570b3', R='#e7298a')
for(i in 1:4)
{
	T = TT[i]
	B = BB[i]
	M = MM[i]


	yb = foreach(params = iter(oldSamples, by='row'), .combine=rbind) %dopar%
		rtob(params, temp, precip, T, B, M)
	yt = foreach(params = iter(oldSamples, by='row'), .combine=rbind) %dopar%
		rtot(params, temp, precip, T, B, M)
	ym = foreach(params = iter(oldSamples, by='row'), .combine=rbind) %dopar%
		rtom(params, temp, precip, T, B, M)
	yr = foreach(params = iter(oldSamples, by='row'), .combine=rbind) %dopar%
		rtor(params, temp, precip, T, B, M)

	plot(0,0,xlim=range(temp), ylim=c(0,1), xlab="Scaled Temperature", ylab="Transition Prob", type='n', main=paste0("T = ", T, "; M = ", M, "; B = ", B))
	polygon(c(temp, rev(temp)), c(apply(yb, 2, quantile, 0.05), rev(apply(yb, 2, quantile, 0.95))), border=NA, col=paste0(cols['B'], '44'))
	lines(temp, colMeans(yb), col=cols['B'], lwd=lwd)

	polygon(c(temp, rev(temp)), c(apply(yt, 2, quantile, 0.05), rev(apply(yt, 2, quantile, 0.95))), border=NA, col=paste0(cols['T'], '44'))
	lines(temp, colMeans(yt), col=cols['T'], lwd=lwd)

	polygon(c(temp, rev(temp)), c(apply(ym, 2, quantile, 0.05), rev(apply(ym, 2, quantile, 0.95))), border=NA, col=paste0(cols['M'], '44'))
	lines(temp, colMeans(ym), col=cols['M'], lwd=lwd)

	polygon(c(temp, rev(temp)), c(apply(yr, 2, quantile, 0.05), rev(apply(yr, 2, quantile, 0.95))), border=NA, col=paste0(cols['R'], '44'))
	lines(temp, colMeans(yr), col=cols['R'], lwd=lwd)

	legend(-3, 1, cex=0.8, col=cols, legend=c('R->T', 'R->M', 'R->B', 'R->R'), lwd=1.5, bg='white')
}
