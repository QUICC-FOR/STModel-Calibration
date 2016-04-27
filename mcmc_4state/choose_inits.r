#!/usr/bin/env Rscript

# prep_mcmc.r
# March 1 2016
# Matt Talluto
#
# Selects initial values for 4 chains from a single initial run
#
# following this script, re-run the prep_mcmc script, using -i=file.rds (where file.rds
# is the output file from this script) and using the -nocalib option

# arguments
# -r5: uses a basal area cutoff of 5 for state R

library(parallel)

arg = commandArgs(trailingOnly = TRUE)

r = 'r1'
ch.names = c('full', 'int', 'ab')
if('-r5' %in% arg) {
	r = 'r5'
} else if('-g1' %in% arg) {
	r = 'g1' 
	ch.names = c('full', 'ab')
} else if(-'g5' %in% arg) {
	r = 'g5'
	ch.names = c('full', 'ab')
}



shan = function(init, orig,length.out=3)
{
	while(nrow(init) < length.out)
	{
		dm = as.matrix(dist(rbind(init, orig)))[,1:nrow(init), drop=FALSE]
		sh = apply(dm, 1, function(x) -sum(x * log(x)))
		sh[is.nan(sh)] = 0
		select = ((which(sh == min(sh))) - nrow(init))[1]
		init = rbind(init, orig[select,])
	}
	init
}

get.inits = function(ch.name, r.val)
{
	samples = read.csv(paste0('res/', ch.name, '_', r.val, '_1/posterior.csv'))
	# take the final value of the chain as one initial value
	inits = samples[nrow(samples),,drop=FALSE]

	# get initial values
	shan(inits, samples, length.out=4)
}

inits = mclapply(ch.names, get.inits, r.val=r, mc.cores=detectCores())
names(inits) = ch.names

outfile=paste0('dat/inits_', r, '.rds')
saveRDS(inits, outfile)