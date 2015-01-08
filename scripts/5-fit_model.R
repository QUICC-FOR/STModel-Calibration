#rm(list=ls())

# choice of the SDM
neiborgh == "rf"
neiborgh == "multinom"
neiborgh == "multinom2"

# fit name
fit = 

#------------------------------

source("3-transition_model.R")
source("4-init_params.R")

#------------------------------
# sample data




#------------------------------

# Maximum likelihood estimation
library(GenSA)

#test
cat("starting logLik")
print(model(params, data))

estim.pars = GenSA(par = params, fn = model, lower = par_lo, upper= par_hi, control = list(verbose =TRUE, maxit = 2000, smooth=FALSE), dat = data)


#save(estim.pars, file="../estimated_params/GenSA_test.rdata")
names(estim.pars$par) = unlist(lapply(names(params), function(x){strsplit(x, split = ".", fixed= TRUE)[[1]][[1]]}))
write.table(estim.pars$par,file=paste("../estimated_params/GenSA_", fit, ".txt", sep=""))

