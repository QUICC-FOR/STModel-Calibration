# R CMD BATCH --no-save --no-restore '--args initForFit_name ordre step' 5-fit_model.R r.out
# or # R script 5-fit_model.R initForFit_name
#rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)

initForFit <- as.character(args)[1]
ordre <- as.numeric(args)[2]
step <- as.numeric(args)[3]

#initForFit = "initForFit_rf_0.331"
#ordre = 0
#step=1

(name = paste(substr(initForFit, 12, 19),"_", ordre, "_",step, "y",sep=""))
#------------------------------
source(paste("3-transition_model_",ordre,".R", sep =""))
load(paste(initForFit, ".RData", sep=""))
#------------------------------

if(ordre<3)
{
if(ordre==2) params = params[c("ab0", "ab1", "ab2", "ab3","ab4", "at0", "at1" , "at2", "at3", "at4", "bb0", "bb1", "bb2", "bb3", "bb4", "bt0", "bt1", "bt2", "bt3", "bt4", "tt0", "tt1", "tt2", "tt3", "tt4", "th0", "th1", "th2", "th3", "th4", "e0", "e1", "e2", "e3", "e4")]
if(ordre==1) params = params[c("ab0", "ab1", "ab2", "at0", "at1" , "at2", "bb0", "bb1", "bb2",  "bt0", "bt1", "bt2", "tt0", "tt1", "tt2", "th0", "th1", "th2", "e0", "e1", "e2")]
if(ordre==0) params = params[c("ab0", "at0", "bb0", "bt0", "tt0", "th0", "e0")]
}

par_lo = rep(-50, length(params))
par_hi = rep(50, length(params))

# Maximum likelihood estimation
library(GenSA)

#test
cat("starting logLik")
print(model(params, datSel))

estim.pars = GenSA(par = params, fn = model, lower = par_lo, upper= par_hi, control = list(verbose =TRUE, smooth=FALSE, temperature = 7000, nb.stop.improvement=1000, maxit= 10000, trace.fn = paste("../estimated_params/traceMat_", name,".trMat", sep="")), dat = datSel, step=step)


# max.time = 80000 = 24h

#save(estim.pars, file="../estimated_params/GenSA_test.rdata")
names(estim.pars$par) = unlist(lapply(names(params), function(x){strsplit(x, split = ".", fixed= TRUE)[[1]][[1]]}))
write.table(estim.pars$par,file=paste("../estimated_params/GenSA_", name, ".txt", sep=""), quote=FALSE, col.names=FALSE)
save(estim.pars, par_lo, par_hi, file=paste("../estimated_params/GenSA_", name, ".RData", sep=""))
