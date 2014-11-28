rm(list=ls())

source("4-init_params.R")
source("3-transition_model.R")

# Maximum likelihood estimation
library(GenSA)

#test
model(params, data)

estim.pars = GenSA(par = params, fn = model, lower = par_lo, upper= par_hi, control = list(verbose =TRUE, max.time = 2000, smooth=FALSE), dat = data)



#save(estim.pars, file="../estimated_params/GenSA_test.rdata")
names(estim.pars$par) = unlist(lapply(names(params), function(x){strsplit(x, split = ".", fixed= TRUE)[[1]][[1]]}))
write.table(estim.pars$par,file="../estimated_params/GenSA_test.txt")

