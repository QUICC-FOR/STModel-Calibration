#!/usr/bin/env Rscript

# prep_mcmc.r
# March 1 2016
# Matt Talluto
#
# Prepare 4-state data for the MCMC runs.
#
# Inputs:
#	calibration_data.rds

# outputs:
#	dat/mcmc_calib.txt  (transition dataset for the mcmc)


# arguments
# -r5: uses a basal area cutoff of 5 for state R
# -ifile=filename.rds: rds file containing initial values to use; if missing set to 0

arg = commandArgs(trailingOnly = TRUE)

r.level = 1
if('-r5' %in% arg) r.level = 5

# handle initial values
ifile = "dat/inits_zero.rds"
if(length(iarg <- grep('-i=', arg, value=T)) > 0)
	ifile = substr(iarg, 4, nchar(iarg))

initialVals = readRDS(ifile)

calib.name = paste0("dat/calibration_data_r", r.level, ".rds")
calib.output.name = paste0("dat/mcmc_calib_r", r.level, ".txt")
calib = readRDS(calib.name)


# set up the transition file:
mcmcData = data.frame(
	env1 = calib[,"ENV1"],
	env2 = calib[,"ENV2"],
	initial = calib[,"st0"],
	final = calib[,"st1"],
	interval = calib[,"itime"],
	plot_id = calib[,"plot_id"],
	prevalenceB = calib[,"B"],
	prevalenceT = calib[,"T"],
	prevalenceM = calib[,"M"],
	prevalenceR = calib[,"R"]
)

write.csv(mcmcData, calib.output.name, row.names=FALSE)


# prep the inits - full model
inits = initialVals$full
parNames = colnames(inits)
isConstant = rep(0,length(parNames))
isConstant[grep("[a-z][a-z]?[5-6]", parNames)] = 1
for(i in 1:nrow(inits))
{
	mcmcIFile = paste0("dat/mcmcInits_full_r", r.level, "_", i, ".txt")
	mcmcInits = data.frame(
		name = parNames,
		initialValue = as.numeric(inits[i,]),
		samplerVariance = 0.5,
		priorMean = 0,
		priorSD = 2.5,
		priorDist = "Cauchy",
		isConstant = isConstant)
	mcmcInits$priorSD[grep("[a-z][a-z]?[0]", parNames)] = 10
	write.csv(mcmcInits, mcmcIFile, row.names=FALSE)
}


# prep the ab constant inits
isConstant = rep(0,length(parNames))
isConstant[grep("[a-z][a-z]?[5-6]", parNames)] = 1
isConstant[grep("ab[1-4]", parNames)] = 1
inits = initialVals$ab
for(i in 1:nrow(inits))
{
	mcmcIFile = paste0("dat/mcmcInits_ab_r", r.level, "_", i, ".txt")
	mcmcInits = data.frame(
		name = parNames,
		initialValue = as.numeric(inits[i,]),
		samplerVariance = 0.5,
		priorMean = 0,
		priorSD = 2.5,
		priorDist = "Cauchy",
		isConstant = isConstant)
	mcmcInits$priorSD[grep("[a-z][a-z]?[0]", parNames)] = 10
	write.csv(mcmcInits, mcmcIFile, row.names=FALSE)
}


# prep the intercept only inits
isConstant = rep(1,length(parNames))
isConstant[grep("[a-z][a-z]?[0]", parNames)] = 0
inits = initialVals$int
for(i in 1:nrow(inits)) {
	mcmcIFile = paste0("dat/mcmcInits_int_r", r.level, "_", i, ".txt")
	mcmcInits = data.frame(
		name = parNames,
		initialValue = 0,
		samplerVariance = 0.5,
		priorMean = 0,
		priorSD = 2.5,
		priorDist = "Cauchy",
		isConstant = isConstant)
	mcmcInits$priorSD[grep("[a-z][a-z]?[0]", parNames)] = 10
	write.csv(mcmcInits, mcmcIFile, row.names=FALSE)
}