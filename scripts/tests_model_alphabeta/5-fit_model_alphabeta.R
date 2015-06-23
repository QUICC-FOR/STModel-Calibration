# R CMD BATCH --no-save --no-restore '--args initForFit_name ordre step' 5-fit_model.R r.out
# or # R script 5-fit_model.R initForFit_name
#rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)

initForFit <- as.character(args)[1]
#ordre <- as.numeric(args)[2]
step <- as.numeric(args)[2]

#initForFit = "initForFit_rf_0.331"
ordre = 3
#step=1

(name = paste(substr(initForFit, 12, 19),"_", ordre, "_",step, "y_alphabeta",sep=""))
#------------------------------
source("3-transition_model_alphabeta.R")
load(paste(initForFit, ".RData", sep=""))
#------------------------------


par_lo = rep(-50, length(params))
par_hi = rep(50, length(params))

# Maximum likelihood estimation
library(GenSA)

#test
cat("starting logLik")
print(model(params, datSel))

estim.pars = GenSA(par = params, fn = model, lower = par_lo, upper= par_hi, control = list(verbose =TRUE, smooth=FALSE, temperature = 7000, nb.stop.improvement=1000, max.time= 10000, trace.fn = paste("../estimated_params/traceMat_", name,".trMat", sep="")), dat = datSel, step=step)


# max.time = 80000 = 24h

#save(estim.pars, file="../estimated_params/GenSA_test.rdata")
names(estim.pars$par) = unlist(lapply(names(params), function(x){strsplit(x, split = ".", fixed= TRUE)[[1]][[1]]}))
write.table(estim.pars$par,file=paste("../estimated_params/GenSA_", name, ".txt", sep=""), quote=FALSE, col.names=FALSE)
save(estim.pars, par_lo, par_hi, file=paste("../estimated_params/GenSA_", name, ".RData", sep=""))
