

args <- commandArgs(trailingOnly = TRUE)

ordre = 2
step= as.numeric(args)[1]
rep =as.numeric(args)[2]
name = paste("rf_all_", ordre, "_", step,"y_rep",rep, sep="")
#------------------------------
source(paste("3-transition_model_",ordre,".R", sep =""))
#------------------------------

##-----------
##
load("datAll.RData")
datSel = dat
veget_pars = read.table(paste("../estimated_params/GenSA_rf_0.339_",ordre,"_",step, "y.txt", sep=""))
params = veget_pars[,2]
names(params) = veget_pars[,1]
##-----------

if(ordre<3)
{
if(ordre==2) params = params[c("ab0", "ab1", "ab2", "ab3","ab4", "at0", "at1" , "at2", "at3", "at4", "bb0", "bb1", "bb2", "bb3", "bb4", "bt0", "bt1", "bt2", "bt3", "bt4", "tt0", "tt1", "tt2", "tt3", "tt4", "th0", "th1", "th2", "th3", "th4", "e0", "e1", "e2", "e3", "e4")]
if(ordre==1) params = params[c("ab0", "ab1", "ab2", "at0", "at1" , "at2", "bb0", "bb1", "bb2",  "bt0", "bt1", "bt2", "tt0", "tt1", "tt2", "th0", "th1", "th2", "e0", "e1", "e2")]
if(ordre==0) params = params[c("ab0", "at0", "bb0", "bt0", "tt0", "th0", "e0")]
}


#test
cat("starting logLik")
print(model(params, datSel, step=step))
#c(y1_334= 13333,  y1_336= 13252, y1_337=13208, y1_339= 13250)
#c(y5_332= 13192, y5_333= 13174, y5_334= 13248, y5_335=13253,  y5_337=13168, y5_339= 13183)
#
#

par_lo = rep(-50, length(params))
par_hi = rep(50, length(params))


# Maximum likelihood estimation
library(GenSA)

estim.pars = GenSA(par = params, fn = model, lower = par_lo, upper= par_hi, control = list(verbose =TRUE, smooth=FALSE, temperature = 7000, nb.stop.improvement=1000, max.time= 80000, trace.fn = paste("../estimated_params/rep_order2_allDat/traceMat_", name,".trMat", sep="")), dat = datSel, step=step)


names(estim.pars$par) = unlist(lapply(names(params), function(x){strsplit(x, split = ".", fixed= TRUE)[[1]][[1]]}))
write.table(estim.pars$par,file=paste("../estimated_params/rep_order2_allDat/GenSA_", name, ".txt", sep=""), quote=FALSE, col.names=FALSE)
save(estim.pars, par_lo, par_hi, file=paste("../estimated_params/rep_order2_allDat/GenSA_", name, ".RData", sep=""))
