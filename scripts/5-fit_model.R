# R CMD BATCH --no-save --no-restore '--args choiceSDM subsetProp' 5-fit_model.R r.out
#rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)

## added by MT for using Rscript so we can see console output on the command line
# instead of R CMD BATCH above, use: 
# Rscript 5-fit_model.R choiceSDM subsetProp
# to save to a file instead of printing to the console, use:
# Rscript 5-fit_model.R choiceSDM subsetProp > r.out
if(!interactive()) {
	args = commandArgs(TRUE)
}



# choice of the SDM
#neiborgh == "rf"
#neiborgh == "multinom"
#neiborgh == "multinom2"
neiborgh <- as.character(args)[1]
subsetProp = as.numeric(args[2])


# fit name
fit = paste(neiborgh, subsetProp*100, sep="_")
print(fit)
#------------------------------
#setwd("/Users/isabelle/Documents/RESEARCH/RECHERCHE/2013-2015 UQAR/QUICCFOR/STModel-Calibration/scripts")
source("3-transition_model.R")
source("4-init_params.R")

#print(getwd())

#------------------------------
# sample data (stratified)

library(lhs)
sampl.raw = randomLHS(as.integer(subsetProp*nrow(dat)*1.25), 2)
# transform to the data range
sampl = data.frame(sampl.raw)
sampl[,1] = sampl[,1]*(range(dat$ENV1)[2]-range(dat$ENV1)[1]) + range(dat$ENV1)[1]
sampl[,2] = sampl[,2]*(range(dat$ENV2)[2]-range(dat$ENV2)[1]) + range(dat$ENV2)[1]


# remove from outside the convex hull
library(grDevices)
hull = chull(dat$ENV1, dat$ENV2)
library(splancs)
inPoly <- inout(as.points(sampl[,1], sampl[,2]),as.points(dat$ENV1[hull], dat$ENV2[hull]))
sampl2 = sampl[inPoly,]

cat("sampl asked ", nrow(dat)*subsetProp, "\n")
cat("sampl taken ", nrow(sampl2), "\n")

library(sp)
distance = spDists(as.matrix(dat[,c("ENV1","ENV2")]), sampl2)
#head(distance)

select = rep(NA, nrow(sampl2))
for( i in 1:nrow(sampl2))
{
select[i] = which.min(distance[,i])
distance[select[i],] = 100
}

datSel = dat[select,]

pdf(paste("../figures/climaticSpace_sample_",fit,".pdf", sep=""))
plot(dat$ENV1, dat$ENV2, pch = 19, cex=.5, main = "climatic space", xlab="temperature", ylab = "precipitations")
points(datSel$ENV1, datSel$ENV2, pch = 19, cex=.5, main = "climatic space", xlab="temperature", ylab = "precipitations", col =2)
polygon(dat$ENV1[hull], dat$ENV2[hull],  col = NA)
dev.off()

save(select, file = paste("selectForFit_", fit,sep=""))
#------------------------------

# Maximum likelihood estimation
library(GenSA)

#test
cat("starting logLik")
print(model(params, datSel))

estim.pars = GenSA(par = params, fn = model, lower = par_lo, upper= par_hi, control = list(verbose =TRUE, maxit = 5000, smooth=FALSE), dat = datSel)


#save(estim.pars, file="../estimated_params/GenSA_test.rdata")
names(estim.pars$par) = unlist(lapply(names(params), function(x){strsplit(x, split = ".", fixed= TRUE)[[1]][[1]]}))
write.table(estim.pars$par,file=paste("../estimated_params/GenSA_", fit, ".txt", sep=""))

