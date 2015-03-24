# R CMD BATCH --no-save --no-restore '--args initForFit_name' 5-fit_model.R r.out
# or # R script 5-fit_model.R initForFit_name
#rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)

initForFit <- as.character(args)[1]

#------------------------------
#setwd("/Users/isabelle/Documents/RESEARCH/RECHERCHE/2013-2015 UQAR/QUICCFOR/STModel-Calibration/scripts")
source("3-transition_model.R")
load(paste(initForFit, ".RData", sep=""))
load(paste("../scripts/GenSA_init.RObj", sep = ""))
params = estim.pars$par

#print(getwd())

#------------------------------

# Maximum likelihood estimation
library(GenSA)

#test
cat("starting logLik")
print(model(params, datSel))

estim.pars = GenSA(par = params, fn = model, lower = par_lo, upper= par_hi, control = list(verbose =TRUE, maxit = 7000, smooth=TRUE), dat = datSel)


#save(estim.pars, file="../estimated_params/GenSA_test.rdata")
names(estim.pars$par) = unlist(lapply(names(params), function(x){strsplit(x, split = ".", fixed= TRUE)[[1]][[1]]}))
write.table(estim.pars$par,file=paste("../estimated_params/GenSA_", initForFit, ".txt", sep=""), quote=FALSE, col.names=FALSE)
save(estim.pars, file=paste("../estimated_params/GenSA_", initForFit, ".RData", sep=""))
