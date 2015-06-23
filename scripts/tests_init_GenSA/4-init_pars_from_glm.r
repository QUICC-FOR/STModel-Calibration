rm(list = ls())

# Open data
load("../data/transitions_r1.RData")
# Rename and columns
transitionData$annual_pp = transitionData$tot_annual_pp
transitionData$st0 = transitionData$state1
transitionData$st1 = transitionData$state2

pair.dat <- transitionData

# subset 10 degree
select = unique(pair.dat$plot[which(pair.dat$annual_mean_temp<=10)])
pair.dat_subset10 = pair.dat[pair.dat$plot %in% select,]

# Create transition column
pair.dat_subset10$transition <- paste(pair.dat_subset10$st0,pair.dat_subset10$st1,sep="")

# Datset without filters
pair_dat0 <- pair.dat_subset10

# Graph lim
rg_pp <- range(pair.dat$annual_pp)
rg_tp <- range(pair.dat$annual_mean_temp)



###################################################################
#####    Analyses     GLM                                      #######
###################################################################
library(MASS)
library(ROCR)
library(fmsb)
pal = colorRampPalette(c("lightblue", "yellow", "orange"), space = "rgb")

#####    glm climate                            #######

modelTransition_climate <- function(st0 , st1, pair.dat, name = NULL)
{
print(st0)
print("->")
print(st1)

datst0 = pair.dat[pair.dat$st0%in%st0, ]
datst0$transition = ifelse(datst0$st1 %in% st1, 1, 0)

mod = glm(transition ~ annual_mean_temp + I(scale(annual_mean_temp)^2) + I(scale(annual_mean_temp)^3) + annual_pp + I(scale(annual_pp)^2) + I(scale(annual_pp)^3) + annual_mean_temp:annual_pp, family = "binomial", data =datst0)
stepMod  = stepAIC(mod)
#print(summary(stepMod))


pred = predict(stepMod,new=datst0,"response")
# overall performance
R2 = NagelkerkeR2(stepMod)$R2
#discrimination
perf = performance(prediction(pred, datst0$transition), "auc")
AUC = perf@y.values[[1]]

## selected vars
## selected vars
coeff = summary(stepMod)$coefficients
vars = rownames(coeff)[-1]
effect = coeff[-1,1]
pval = coeff[-1,4]
print(pval)

if(is.null(name)) name = paste(st0, st1, sep = "->")

return(list(mod = stepMod, vars = vars,effect = effect, pval = pval , R2 = R2, AUC = AUC, ranges = apply(datst0[unlist(lapply(1:ncol(datst0), function(x)is.numeric(datst0[,x])))], 2, range), name = name))
}


#------------------------------------------------------------------

pair_dat = pair_dat0
modelTransition = modelTransition_climate
# regeneration
modRT = modelTransition("R", "T", pair.dat = pair_dat)
modRB = modelTransition("R", "B", pair.dat = pair_dat)
modRM = modelTransition("R", "M", pair.dat = pair_dat)
modRR = modelTransition("R", "R", pair.dat = pair_dat)
# exclusion
modMT = modelTransition(c("M"), c("T"), pair.dat = pair_dat)
modMB = modelTransition(c("M"), c("B"), pair.dat = pair_dat)
# colonisation
modTM = modelTransition(c("T"), c("M"), pair.dat = pair_dat)
modBM = modelTransition(c("B"), c("M"), pair.dat = pair_dat)
# disturbance
modR = modelTransition(c("T", "B", "M"), c("R"), pair.dat = pair_dat)
modMR = modelTransition(c("M"), c("R"), pair.dat = pair_dat)
modTR = modelTransition("T", c("R"), pair.dat = pair_dat)
modBR = modelTransition(c("B"), c("R"), pair.dat = pair_dat)

modTT = modelTransition("T", "T", pair.dat = pair_dat)
modBB = modelTransition("B", "B", pair.dat = pair_dat)
modMM = modelTransition("M", "M", pair.dat = pair_dat)


#------------------------------------------------------------------

##-----------
load("../data/transitions_r1.RData")

dataProj = transitionData
head(dataProj)
dim(dataProj)
#str(dataProj)

# subset 10 degree
select = unique(dataProj$plot[which(dataProj$annual_mean_temp<=10)])
dataProj_subset10 = dataProj[dataProj$plot %in% select,]

# rescale
load("scale_info.Robj")
dat_scale = dataProj_subset10[c("annual_mean_temp", "tot_annual_pp")]
dat_scale = t(apply(dat_scale, 1, function(x) {(x-vars.means[c("annual_mean_temp", "tot_annual_pp")])/vars.sd[c("annual_mean_temp", "tot_annual_pp")]}))
dat_scale = data.frame(dat_scale)
#head(dat_scale)
dim(dat_scale)

# remove transitions directes B->T ou T->B
toremove = c(which(dataProj_subset10$state1 == "T" & dataProj_subset10$state2 == "B"), which(dataProj_subset10$state1 == "B" & dataProj_subset10$state2 == "T"))
dataProj_subset10= dataProj_subset10[-toremove,]

#------------------------------------------------------------------
dataSet = cbind(st0 = dataProj_subset10$state1)
#------------------------------------------------------------------
dataProj_subset10$annual_pp = dataProj_subset10$tot_annual_pp
itime = dataProj_subset10[, "year2"] - dataProj_subset10[, "year1"]

#
predT = predict(modRT$mod, newdata = dataProj_subset10[, c("annual_mean_temp", "annual_pp")], type = "response")
predB = predict(modRB$mod, newdata = dataProj_subset10[, c("annual_mean_temp", "annual_pp")], type = "response")
predM = predict(modRM$mod, newdata = dataProj_subset10[, c("annual_mean_temp", "annual_pp")], type = "response")
predR = predict(modRR$mod, newdata = dataProj_subset10[, c("annual_mean_temp", "annual_pp")], type = "response")

datR = cbind(
T=predT[dataProj_subset10$state1=="R"]^itime[dataProj_subset10$state1=="R"], 
B=predB[dataProj_subset10$state1=="R"]^itime[dataProj_subset10$state1=="R"], 
M=predM[dataProj_subset10$state1=="R"]^itime[dataProj_subset10$state1=="R"],
R=predR[dataProj_subset10$state1=="R"]^itime[dataProj_subset10$state1=="R"]
)
head(datR)
Rx = colnames(datR)[apply(datR, 1, function(x){which(rmultinom(1,1,x/(sum(x)))==1)})]
#----------
predT = predict(modTT$mod, newdata = dataProj_subset10[, c("annual_mean_temp", "annual_pp")], type = "response")
predB = rep(0, nrow(dataProj_subset10))
predM = predict(modTM$mod, newdata = dataProj_subset10[, c("annual_mean_temp", "annual_pp")], type = "response")
predR = predict(modTR$mod, newdata = dataProj_subset10[, c("annual_mean_temp", "annual_pp")], type = "response")

datT = cbind(
T=predT[dataProj_subset10$state1=="T"]^itime[dataProj_subset10$state1=="T"], 
B=predB[dataProj_subset10$state1=="T"]^itime[dataProj_subset10$state1=="T"], 
M=predM[dataProj_subset10$state1=="T"]^itime[dataProj_subset10$state1=="T"],
R=predR[dataProj_subset10$state1=="T"]^itime[dataProj_subset10$state1=="T"]
)
head(datT)
Tx = colnames(datT)[apply(datT, 1, function(x){which(rmultinom(1,1,x/(sum(x)))==1)})]
#-----
predT = rep(0, nrow(dataProj_subset10))
predB = predict(modBB$mod, newdata = dataProj_subset10[, c("annual_mean_temp", "annual_pp")], type = "response")
predM = predict(modBM$mod, newdata = dataProj_subset10[, c("annual_mean_temp", "annual_pp")], type = "response")
predR = predict(modBR$mod, newdata = dataProj_subset10[, c("annual_mean_temp", "annual_pp")], type = "response")

datB = cbind(
T=predT[dataProj_subset10$state1=="B"]^itime[dataProj_subset10$state1=="B"], 
B=predB[dataProj_subset10$state1=="B"]^itime[dataProj_subset10$state1=="B"], 
M=predM[dataProj_subset10$state1=="B"]^itime[dataProj_subset10$state1=="B"],
R=predR[dataProj_subset10$state1=="B"]^itime[dataProj_subset10$state1=="B"]
)
head(datB)
Bx = colnames(datB)[apply(datB, 1, function(x){which(rmultinom(1,1,x/(sum(x)))==1)})]
#------
predT = predict(modMT$mod, newdata = dataProj_subset10[, c("annual_mean_temp", "annual_pp")], type = "response")
predB = predict(modMB$mod, newdata = dataProj_subset10[, c("annual_mean_temp", "annual_pp")], type = "response")
predM = predict(modMM$mod, newdata = dataProj_subset10[, c("annual_mean_temp", "annual_pp")], type = "response")
predR = predict(modMR$mod, newdata = dataProj_subset10[, c("annual_mean_temp", "annual_pp")], type = "response")

datM = cbind(
T=predT[dataProj_subset10$state1=="M"]^itime[dataProj_subset10$state1=="M"], 
B=predB[dataProj_subset10$state1=="M"]^itime[dataProj_subset10$state1=="M"], 
M=predM[dataProj_subset10$state1=="M"]^itime[dataProj_subset10$state1=="M"],
R=predR[dataProj_subset10$state1=="M"]^itime[dataProj_subset10$state1=="M"]
)
head(datM)
Mx = colnames(datM)[apply(datM, 1, function(x){which(rmultinom(1,1,x/(sum(x)))==1)})]


#------------------------------------------------------------------
dataSet = data.frame(st0 = dataProj_subset10$state1, itime = itime)
dataSet$st1 = rep(NA, nrow(dataSet))
dataSet[dataProj_subset10$state1=="R", "st1" ] = Rx
dataSet[dataProj_subset10$state1=="T", "st1"] = Tx
dataSet[dataProj_subset10$state1=="B", "st1"] = Bx
dataSet[dataProj_subset10$state1=="M", "st1"] = Mx
head(dataSet)
table(dataSet$st0, dataSet$st1)

#------------------------------------------------------------------
HK <- function (Pred, Obs)
{
	Misc = table(Pred, Obs)

    if (nrow(Misc)!=ncol(Misc)) stop("wrong misclassification table")
    Misc <- unclass(Misc)
    k  <- ncol(Misc)
    Nobs <- apply(Misc, 2, sum)
    Npred <- apply(Misc, 1, sum)
    N <- sum(Nobs)


   HK <- (sum(diag(Misc))/N - sum(as.numeric(Nobs)*as.numeric(Npred))/N/N ) / ( 1 - sum(as.numeric(Nobs)*as.numeric(Nobs))/N/N )

    return(HK)
}

HK(dataSet$st1, dataProj_subset10$state2) #0.8827963

#------------------------------------------------------------------

# neighborhood

pred = read.table("../data/projection_rf_complete.txt", h=T)
#pred = read.table("../data/projection_multimod_complete.txt", h=T)

pred = pred[pred$plot %in% select,]
pred = pred[-toremove,]

dataSet$EB = pred$B
dataSet$ET = pred$T
dataSet$EM = pred$M

load("scale_info.Robj")
dat_scale = dataProj_subset10[c("annual_mean_temp", "tot_annual_pp")]
dat_scale = t(apply(dat_scale, 1, function(x) {(x-vars.means[c("annual_mean_temp", "tot_annual_pp")])/vars.sd[c("annual_mean_temp", "tot_annual_pp")]}))
dat_scale = data.frame(dat_scale)
head(dat_scale)

dataSet$ENV1 = dat_scale$annual_mean_temp
dataSet$ENV2 = dat_scale$tot_annual_pp


head(dataSet)

#------------------------------------------------------------------
# Evaluate initial parameter values
transitions = paste(dataSet$st0,dataSet$st1,sep = "")
sum_transitions = table(transitions)
initState = table(dataSet$st0)

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


# transform into logits
logit_eps_mn = log(eps_mn/(1-eps_mn))

logit_theta_mn = log(theta_mn/(1-theta_mn))
logit_thetat_mn = log(thetat_mn/(1-thetat_mn))
 
logit_betab_mn = log(betab_mn/(1-betab_mn))
logit_betat_mn = log(betat_mn/(1-betat_mn))

logit_alphab_mn = log(alphab_mn/(1-alphab_mn))
logit_alphat_mn = log(alphat_mn/(1-alphat_mn))


   

params = c(ab0 =as.numeric(logit_alphab_mn), ab1 = 0, ab2 = 0, ab3=0, ab4=0, ab5=0, ab6=0, 
at0 = as.numeric(logit_alphat_mn), at1 = 0, at2=0, at3=0, at4=0, at5=0, at6=0,
bb0 = as.numeric(logit_betab_mn), bb1 = 0, bb2 = 0, bb3=0, bb4=0, bb5=0, bb6=0,
bt0 = as.numeric(logit_betat_mn), bt1 = 0, bt2 = 0, bt3=0, bt4=0, bt5=0, bt6=0,
tt0 = as.numeric(logit_thetat_mn), tt1=0, tt2=0, tt3=0, tt4=0, tt5=0, tt6=0,
t0 = as.numeric(logit_theta_mn), t1=0, t2=0, t3=0, t4=0, t5=0, t6=0,
e0 = as.numeric(logit_eps_mn), e1 = 0,  e2=0, e3=0, e4=0, e5=0, e6=0)


scaleOfVar = c(rep(200, 14), rep(50, 35))

par_lo = params - scaleOfVar

par_hi = params + scaleOfVar
#-------------------------------------------------

#source("subsample.r")
#select2 = subsample.stratif3D(dataProj_subset10[,c("lon","lat", "annual_mean_temp")], 0.3, adj = 4.2)
#save(select2, file = "select_fakeDataset_03.RData")
load("select_fakeDataset_03.RData")

datSel = dataSet[select2,]


#-------------------------------------------------
# Maximum likelihood estimation
library(GenSA)
source("3-transition_model.R")

#test
cat("starting logLik")
params[is.na(params)] = 0
par_lo[is.na(par_lo)] = -50
par_hi[is.na(par_hi)] = 50

datSel$EB = rep(0.25, nrow(datSel))
datSel$ET = rep(0.25, nrow(datSel))
datSel$EM = rep(0.25, nrow(datSel))

print(model(params, datSel))



estim.pars = GenSA(par = params, fn = model, lower = par_lo, upper= par_hi, control = list(verbose =TRUE, maxit = 7000, smooth=TRUE), dat = datSel)


#save(estim.pars, file="../estimated_params/GenSA_test.rdata")
names(estim.pars$par) = unlist(lapply(names(params), function(x){strsplit(x, split = ".", fixed= TRUE)[[1]][[1]]}))

write.table(estim.pars$par,file="../estimated_params/GenSA_init_constNeigh.txt")
save(estim.pars, file="../estimated_params/GenSA_init_constNeigh.RData")

write.table(estim.pars$par,file="../estimated_params/GenSA_init.txt")
save(estim.pars, file="../estimated_params/GenSA_init.RData")
