# R CMD BATCH --no-save --no-restore '--args initForFit_name option' 5-fit_model.R r.out
# or # R script 5-fit_model.R initForFit_name
#rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)

initForFit <- as.character(args)[1]
option <- as.character(args)[2]
if(is.na(option)) option = ""
#
#initForFit = "initForFit_rf_0.331"
#option = ""

#------------------------------
source(paste("3-transition_model",option,".R", sep =""))
load(paste(initForFit, ".RData", sep=""))
#------------------------------

if(grepl("less", option))
{
params = params[c("ab0", "ab1", "ab2", "ab3","ab4","ab5","ab6", "at0", "at1" , "at2", "at3", "at4", "at5", "at6", "bb0", "bt0", "tt0", "th0", "e0")]
par_lo= par_lo[c("ab0", "ab1", "ab2", "ab3","ab4","ab5","ab6", "at0", "at1" , "at2", "at3", "at4", "at5", "at6", "bb0", "bt0", "tt0", "th0", "e0")]
par_hi = par_hi[c("ab0", "ab1", "ab2", "ab3","ab4","ab5","ab6", "at0", "at1" , "at2", "at3", "at4", "at5", "at6", "bb0", "bt0", "tt0", "th0", "e0")]
}

# Maximum likelihood estimation
library(GenSA)

#test
cat("starting logLik")
print(model(params, datSel))

estim.pars = GenSA(par = params, fn = model, lower = par_lo, upper= par_hi, control = list(verbose =TRUE, smooth=FALSE, max.time = 80000, temperature = 7000,  trace.fn = paste("../estimated_params/traceMat_", initForFit, option, ".trMat", sep="")), dat = datSel)

# max.time = 80000 = 24h

#save(estim.pars, file="../estimated_params/GenSA_test.rdata")
names(estim.pars$par) = unlist(lapply(names(params), function(x){strsplit(x, split = ".", fixed= TRUE)[[1]][[1]]}))
write.table(estim.pars$par,file=paste("../estimated_params/GenSA_", initForFit, option, ".txt", sep=""), quote=FALSE, col.names=FALSE)
save(estim.pars, file=paste("../estimated_params/GenSA_", initForFit, option, ".RData", sep=""))
