# R CMD BATCH --no-save --no-restore '--args choiceSDM subsetProp' 4-init_params.R r.out
# or # R script 4-init_params.R choiceSDM subsetProp
#rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)

# choice of the SDM
#neiborgh == "rf"
#neiborgh == "multinom"
#neiborgh == "multinom2"
neiborgh <- as.character(args)[1]
subsetProp = as.numeric(args)[2]

# fit name
fit = paste(neiborgh, subsetProp, sep="_")
print(fit)
##-----------
## Choose neigborhood type
##-----------

# neighborhood
if(neiborgh == "rf") pred = read.table("../data/projection_rf_complete.txt", h=T)
if(neiborgh == "multinom") pred = read.table("../data/projection_multimod_complete.txt", h=T)


##-----------
## load dat
##-----------
load("../data/transitions_r1.RData")

dataProj = transitionData
head(dataProj)
dim(dataProj)
#str(dataProj)

# subset 10 degree
select = unique(dataProj$plot[which(dataProj$annual_mean_temp<=10)])
dataProj_subset10 = dataProj[dataProj$plot %in% select,]
pred = pred[pred$plot %in% select,]


# rescale
load("scale_info.Robj")
dat_scale = dataProj_subset10[c("annual_mean_temp", "tot_annual_pp")]
dat_scale = t(apply(dat_scale, 1, function(x) {(x-vars.means[c("annual_mean_temp", "tot_annual_pp")])/vars.sd[c("annual_mean_temp", "tot_annual_pp")]}))
dat_scale = data.frame(dat_scale)
head(dat_scale)
dim(dat_scale)

# remove transitions directes B->T ou T->B
toremove = c(which(dataProj_subset10$state1 == "T" & dataProj_subset10$state2 == "B"), which(dataProj_subset10$state1 == "B" & dataProj_subset10$state2 == "T"))

dat = data.frame(dat_scale[-toremove,])
pred = pred[-toremove,]
dataProj_subset10= dataProj_subset10[-toremove,]

#rename variables
colnames(dat) = c("ENV1", "ENV2")
dat$EB = pred$B
dat$ET = pred$T
dat$EM = pred$M 
dat$st0 = dataProj_subset10$state1
dat$st1 = dataProj_subset10$state2
dat$itime = dataProj_subset10$year2 - dataProj_subset10$year1

head(dat)
dim(dat)

rm(dat_scale, dataProj)


# Evaluate initial parameter values
transitions = paste(dat$st0,dat$st1,sep = "")
sum_transitions = table(transitions)
initState = table(dat$st0)

# estimates for initial parameters
eps_mn = # proba x->R given x is T B or M
(sum_transitions["BR"]/initState["B"]+sum_transitions["MR"]/initState["M"]+sum_transitions["TR"]/initState["T"])/3

theta_mn = # proba M-> B or T
(sum_transitions["MB"] + sum_transitions["MT"])/initState["M"]

thetat_mn = # proba M->T
sum_transitions["MT"]/(sum_transitions["MB"] + sum_transitions["MT"])

betab_mn = # proba T->M
(sum_transitions["TM"]/initState["T"]) * (sum(initState)/(initState["M"] + initState["B"]))

betat_mn = # proba B->M
(sum_transitions["BM"]/initState["B"]) * (sum(initState)/(initState["M"] + initState["T"]))

trRT = sum_transitions["RT"]/initState["R"]
trRB = sum_transitions["RB"]/initState["R"]
trRM = sum_transitions["RB"]/initState["R"]

alphat_mn = trRM/((initState["M"]+initState["T"])*(trRB+trRM))
alphab_mn = (trRB+trRM)/(initState["M"] + initState["B"])

# cf in python:
#alphat, alphab, M, T, B, trRT, trRB, trRM = symbols('alphat alphab M T B trRT trRB trRM')
#eq1 = alphat*(M+T)*alphab*(M+B)
#eq2 = alphab*(M+B)*(1-alphat*(M+T))
#res = solve([Eq(eq1, trRM), Eq(eq2, trRB)], [alphat, alphab])
#simplify(res)

# transform into logits
logit_eps_mn = log(eps_mn/(1-eps_mn))

logit_theta_mn = log(theta_mn/(1-theta_mn))
logit_thetat_mn = log(thetat_mn/(1-thetat_mn))
 
logit_betab_mn = log(betab_mn/(1-betab_mn))
logit_betat_mn = log(betat_mn/(1-betat_mn))

logit_alphab_mn = log(alphab_mn/(1-alphab_mn))
logit_alphat_mn = log(alphat_mn/(1-alphat_mn))



# List initial parameters
#    ab4 = ab5 = ab6 = at2 = at4 = at6 = 0 = tt1 = tt2 = tt3 = tt4 = tt5 = tt6 = t1 = t2 = t3 = t4 = t5 = t6 = e2 = e4 = e6 = e7 = 0
    

params = c(ab0 =as.numeric(logit_alphab_mn), ab1 = 0, ab2 = 0, ab3=0,
at0 = as.numeric(logit_alphat_mn), at1 = 0, at3=0, at5=0,
bb0 = as.numeric(logit_betab_mn), bb1 = 0, bb2 = 0, bb3=0, bb4=0, bb5=0, bb6=0,
bt0 = as.numeric(logit_betat_mn), bt1 = 0, bt2 = 0, bt3=0, bt4=0, bt5=0, bt6=0,
tt0 = as.numeric(logit_thetat_mn), 
t0 = as.numeric(logit_theta_mn), 
e0 = as.numeric(logit_eps_mn), e1 = 0,  e3=0,e5=0)

# bounds
# coeff variation
#cvar = 5
#scaleOfVar = c(rep(abs(logit_alphab_mn*cvar), 7), 
#   rep(abs(logit_alphat_mn*cvar), 7) ,
#   rep(abs(logit_betab_mn*cvar), 7) ,
#   rep(abs(logit_betat_mn*cvar), 7) ,
#   rep(abs(logit_thetat_mn*cvar), 7) ,
#   rep(abs(logit_theta_mn*cvar), 7) ,
#   rep(abs(logit_eps_mn*cvar), 7) )
scaleOfVar = c(rep(200, length(params)))

par_lo = params - scaleOfVar

par_hi = params + scaleOfVar

#-------------------------------------------------------------------
#
# sample data (stratified)
#
#------------------------------
nrow(dataProj_subset10) == nrow(dat)

source("subsample.r")
select = subsample.stratif3D(dataProj_subset10[,c("lon","lat", "annual_mean_temp")], subsetProp, adj = 4.2)


jpeg("../figures/subsample_fit.jpeg", height=5000, width=5000, res=600)
plot(dataProj_subset10[,c("lon","lat")], pch = 20, cex=.2, col = "grey")
points(dataProj_subset10[select,c("lon","lat")], pch = 20, cex=.2, col = 1)
dev.off()


datSel = dat[select,]

#-----------
#coords = cbind(dataProj_subset10$longitude, dataProj_subset10$latitude)
save(datSel, dataProj_subset10, select,params, par_lo, par_hi, file = paste("initForFit_", fit,sep=""))
