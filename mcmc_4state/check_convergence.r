#!/usr/bin/env Rscript

# on mammouth:
# module load bioinformatics/R

# on froggy
# source /applis/ciment/v2/env.bash
# module load R/3.2.3_gcc-4.6.2

# check_convergence.r
# March 1 2016
# Matt Talluto
#
# examine posterior distributions and report convergence. prints to the console

# arguments:
# the only argument is the name of the model, something like:
# ab_r1
# where ab is the model type, r1 is the r threshold
# chains will be opened automagically

library(coda)

arg = commandArgs(trailingOnly = TRUE)
if(length(arg) != 1) stop("Enter one and only one argument, the model name")

models = grep(paste0('^', arg, '_[0-9]'), list.files('res'), value=TRUE)
if(length(models) < 2) stop("No model or only one model found; check model name")

samples = lapply(models, function(mod) read.csv(file.path('res', mod, 'posterior.csv')))
sampsize = min(sapply(samples, nrow))
samples = lapply(samples, function(sam) sam[1:sampsize,])
cat("Computing convergence diagnostic for", sampsize, "samples\n")

# drop columns where all values are zero
samples = lapply(samples, function(sam, tol=1e-8)
	sam[,!apply(sam, 2, function(x) all(abs(x) < tol))])
	
samples.mcmc = as.mcmc.list(lapply(samples, mcmc))
print(gelman.diag(samples.mcmc))

