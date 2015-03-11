# R CMD BATCH --no-save --no-restore '--args initForFit_name' 5-fit_model.R r.out
# or # R script 5-fit_model.R initForFit_name
#rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)

initForFit <- as.character(args)[1]

#------------------------------
#setwd("/Users/isabelle/Documents/RESEARCH/RECHERCHE/2013-2015 UQAR/QUICCFOR/STModel-Calibration/scripts")
source("3-transition_model.R")
load(paste(initForFit, ".RData", sep=""))
#load(paste("../estimated_params/GenSA_initForFit_", strsplit(fit, "_")[[1]][1], "_0.05.RData", sep = ""))
#params = estim.pars$par

#print(getwd())

#------------------------------

# Maximum likelihood estimation
library(GenSA)

#test
cat("starting logLik")
print(model(params, datSel))

##----
datSel$ENV1.sq = scale(datSel$ENV1^2)
datSel$ENV1.cu = scale(datSel$ENV1^3)
datSel$ENV2.sq = scale(datSel$ENV2^2)
datSel$ENV2.cu = scale(datSel$ENV2^3)

scale_poly = list(means = c(ENV1.sq = mean(datSel$ENV1^2), ENV1.cu = mean(datSel$ENV1^3), ENV2.sq = mean(datSel$ENV2^2), ENV2.cu = mean(datSel$ENV2^3)), vars = c(ENV1.sq = sd(datSel$ENV1^2), ENV1.cu = sd(datSel$ENV1^3), ENV2.sq = sd(datSel$ENV2^2), ENV2.cu = sd(datSel$ENV2^3)))

save(scale_poly, file = paste(initForFit, "_scale_poly.RData", sep = ""))
##----


estim.pars = GenSA(par = params, fn = model, lower = par_lo, upper= par_hi, control = list(verbose =TRUE, maxit = 2000, smooth=FALSE), dat = datSel)


#save(estim.pars, file="../estimated_params/GenSA_test.rdata")
names(estim.pars$par) = unlist(lapply(names(params), function(x){strsplit(x, split = ".", fixed= TRUE)[[1]][[1]]}))
write.table(estim.pars$par,file=paste("../estimated_params/GenSA_", initForFit, ".txt", sep=""))
save(estim.pars, file=paste("../estimated_params/GenSA_", initForFit, ".RData", sep=""))
