rm(list=ls())

source("4-init_params.R")
source("3-transition_model.R")

# Maximum likelihood estimation
library(GenSA)

#test
model(params, data)

estim.pars = GenSA(par = params, fn = model, lower = par_lo, upper= par_hi, control = list(verbose =TRUE, max.time = 2000, smooth=FALSE), dat = data)


estim.pars2 = GenSA(par = estim.pars$par, fn = model, lower = par_lo, upper= par_hi, control = list(verbose =TRUE, max.time = 2000, smooth=FALSE), dat = data)


#save(estim.pars2, file="../estimated_params/GenSA_test.rdata")
names(estim.pars2$par) = unlist(lapply(names(params), function(x){strsplit(x, split = ".", fixed= TRUE)[[1]][[1]]}))
write.table(estim.pars2$par,"../estimated_params/GenSA_test.txt")


#
#save(coarse, file="../estimated_params/coarse_m3.rdata")
#load("../estimated_params/coarse_m3.rdata")
#write.table(coarse$best_pars,"../estimated_params/par_m3.txt")
#write.table(data.frame(names(coarse$best_pars), unlist(coarse$best_pars)),"../estimated_params/par_m3_forc.txt", sep=" ", row.names=F, col.names=F, quote=F)
#
