# prep_mcmc.r
# March 1 2016
# Matt Talluto
#
# Prepare 4-state data for the MCMC run.
#
# Inputs:
#	calibration_data.rds

# outputs:
#	transitions.txt  (transition dataset for the mcmc)



calib = readRDS("dat/calibration_data.rds")
inits = readRDS("dat/inits.rds")

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

write.csv(mcmcData, "dat/mcmc_calib.txt", row.names=FALSE)


# prep the inits
parNames = colnames(inits)
isConstant = rep(0,length(parNames))
isConstant[grep("[a-z][a-z]?[5-6]", parNames)] = 1
for(i in 1:nrow(inits))
{
	mcmcIFile = paste0("dat/mcmc_inits", i, ".txt")
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
mcmcIFile = "dat/mcmc_inits_int.txt"
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
