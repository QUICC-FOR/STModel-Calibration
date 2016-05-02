#!/usr/bin/env Rscript

# process_samples.r
# May 2 2016
# Matt Talluto
#
# set up thinning and process raw .csv samples into coda objects; also grab DIC from text file

# arguments:
# the only argument is the name of the model, something like:
# ab_r1_1
# where ab is the model type, r1 is the r threshold, _1 indicates the chain number to use

library(coda)

arg = commandArgs(trailingOnly = TRUE)
if(length(arg) != 1) stop("Enter one and only one argument, the model name")

burnin = 0.5
length.out = 10000
tol=1e-8


samples = read.csv(file.path('res', arg, 'posterior.csv'))
dic = read.table(file.path('res', arg, 'dic.csv'), header=FALSE, sep=':', row.names=1)
colnames(dic) = "value"

sampsize = nrow(samples)

# drop columns where all values are zero
# apply thin and burnin
strt = ceiling(nrow(samples) * burnin) + 1
keep = seq(strt, nrow(samples), length.out=length.out)
samples = samples[keep,!apply(samples, 2, function(x) all(abs(x) < tol))]
samples = mcmc(samples, thin = as.integer(keep[2] - keep[1]), start=strt)

modname = substr(arg, 1, nchar(arg)-2)
saveRDS(list(dic=dic, samples=samples), paste0('res/', modname, '_samples.rds'))