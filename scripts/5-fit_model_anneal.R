# R CMD BATCH --no-save --no-restore '--args initForFit_name' 5-fit_model.R r.out
# or # R script 5-fit_model.R initForFit_name
#rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)

initForFit <- as.character(args)[1]

#------------------------------
#setwd("/Users/isabelle/Documents/RESEARCH/RECHERCHE/2013-2015 UQAR/QUICCFOR/STModel-Calibration/scripts")
source("3-transition_model_anneal.R")
load(initForFit)

#print(getwd())

#------------------------------

# Maximum likelihood estimation
var = list(st0 = "st0", st1 = "st1", ENV1 = "ENV1", ENV2 = "ENV2", lik = "predicted", itime = "itime", EB = "EB", ET = "ET", EM = "EM")


# Maximum likelihood estimation
source("fit_fonctions/analyze_function.R")
source("fit_fonctions/likdisplay_ibou.R")
source("fit_fonctions/support_limits_ibou.R")
source("fit_fonctions/likeli.R")
source("fit_fonctions/anneal_simple.R")


#test
cat("starting logLik")
with(as.list(params), {
print(sum(model(ab0, ab1, ab2, ab3,ab4,ab5,ab6,
at0, at1 , at2, at3, at4, at5, at6, 
bb0, bb1, bb2, bb3, bb4, bb5, bb6,
bt0, bt1, bt2, bt3, bt4, bt5, bt6,
tt0, tt1, tt2, tt3, tt4, tt5, tt6, 
t0, t1, t2, t3, t4, t5, t6, 
e0, e1, e2, e3, e4, e5, e6, datSel$st0, datSel$st1, datSel$ET, datSel$EB, datSel$EM, datSel$ENV1, datSel$ENV2, datSel$itime)))
})

estim.pars.anneal= anneal(model = model, par = as.list(params), var = var, source_data = datSel, 
dep_var = "st1", pdf = PDF, par_lo = as.list(par_lo), par_hi = as.list(par_hi), max_iter = 10000, initial_temp = 30, note = "", progress = TRUE, display=FALSE, support = FALSE, min_change = 1, min_drops = 10)


save(estim.pars.anneal, file=paste("../estimated_params/GenSA_anneal_", initForFit, ".RData", sep=""))
