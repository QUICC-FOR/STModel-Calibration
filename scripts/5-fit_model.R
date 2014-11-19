rm(list=ls())

source("4-init_params.R")
source("3-transition_model.R")

# Maximum likelihood estimation
library(GenSA)

#test
model(params, data)

estim.pars = GenSA(par = params, fn = model, lower = rep(-30, length(params)), upper= rep(30, length(params)), control = list(verbose =TRUE, max.time = 100, smooth=FALSE), dat = data)


#estim.pars = GenSA(par = params, fn = model, lower = par_lo, upper= par_hi, control = list(verbose =TRUE, max.time = 100, smooth=FALSE), dat = data)


#
#save(coarse, file="../estimated_params/coarse_m3.rdata")
#load("../estimated_params/coarse_m3.rdata")
#write.table(coarse$best_pars,"../estimated_params/par_m3.txt")
#write.table(data.frame(names(coarse$best_pars), unlist(coarse$best_pars)),"../estimated_params/par_m3_forc.txt", sep=" ", row.names=F, col.names=F, quote=F)
#
