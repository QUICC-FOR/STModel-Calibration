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

par_lo = rep(-50, length(params))
par_hi = rep(50, length(params))

# Maximum likelihood estimation
library(GenSA)

#test
cat("starting logLik")
print(model(params, datSel))

estim.pars = GenSA(par = params, fn = model, lower = par_lo, upper= par_hi, control = list(verbose =TRUE, smooth=FALSE, max.time = 80000##-----
#jpeg(paste("../figures/equilibrium_map_", sdm, "_",propData, option, ".jpeg", sep=""), height=3000, width=3000, res = 300)
#
colo = c(AltSS = "pink", 'M wins (coexistence)' = "blue", 'B wins' = "darkgreen", 'B dominates (with M)' ="darkviolet", 'T dominates (with M)' = "violet", 'T wins' = "lightgreen", '??' = 1)

layout(matrix(c(1,2),nr=2,nc=1,byrow=TRUE),heights = c(1,6))

par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
title("",cex=2)
legend("center",legend = levels(coexist),fill = colo[levels(coexist)],bty = "n", cex = 0.8, ncol =2)
par(mar=c(5,5,0,2))

image(x=tpseq, y=ppseq, z = matrix(as.numeric(coexist), ncol = length(ppseq), nrow = length(tpseq)),xlab = "Annual mean temperature (°C)", ylab = "Annual precipitations (mm)", col = colo[levels(coexist)], main = "", xaxt = "n", yaxt="n")
scaled.axis()
#
#dev.off()
, temperature = 7000,  trace.fn = paste("../estimated_params/traceMat_", initForFit, option, ".trMat", sep="")), dat = datSel)

# max.time = 80000 = 24h

#save(estim.pars, file="../estimated_params/GenSA_test.rdata")
names(estim.pars$par) = unlist(lapply(names(params), function(x){strsplit(x, split = ".", fixed= TRUE)[[1]][[1]]}))
write.table(estim.pars$par,file=paste("../estimated_params/GenSA_", initForFit, option, ".txt", sep=""), quote=FALSE, col.names=FALSE)
save(estim.pars, file=paste("../estimated_params/GenSA_", initForFit, option, ".RData", sep=""))
