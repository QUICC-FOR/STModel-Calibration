# R CMD BATCH --no-save --no-restore '--args choiceSDM subsetProp' 4-init_params.R r.out
# or # R script 4-initForFit.R subsetProp
#rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)

# choice of the SDM
subsetProp = as.numeric(args)[1]


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

# rescale
load("scale_info.Robj")
dat_scale = dataProj_subset10[c("annual_mean_temp", "tot_annual_pp")]
dat_scale = t(apply(dat_scale, 1, function(x) {(x-vars.means[c("annual_mean_temp", "tot_annual_pp")])/vars.sd[c("annual_mean_temp", "tot_annual_pp")]}))
dat_scale = data.frame(dat_scale)
#head(dat_scale)
dim(dat_scale)


# remove transitions directes B->T ou T->B
trBT = c(which(dataProj_subset10$state1 == "T" & dataProj_subset10$state2 == "B"), which(dataProj_subset10$state1 == "B" & dataProj_subset10$state2 == "T"))

# clean transition time
dataProj_subset10$itime = dataProj_subset10$year2 - dataProj_subset10$year1
rmItime= c(which(dataProj_subset10$itime<5), which(dataProj_subset10$itime>15))

# clean dataset
toremove = unique(c(trBT, rmItime))
dataProj_subset10 = dataProj_subset10[-toremove,]

# create dataset
dat = dat_scale[-toremove,]

#rename variables
colnames(dat) = c("ENV1", "ENV2")

dat$st0 = dataProj_subset10[,"state1"]
dat$st1 = dataProj_subset10[,"state2"]
dat$itime = dataProj_subset10[,"itime"]
dat$plot_id = dataProj_subset10[,"plot"]



#-------------------------------------------------------------------

# neighborhood

pred = read.table("../data/projection_rf_complete.txt", h=T)

pred = pred[pred$plot %in% select,]
pred = pred[-toremove,]

dat$EB = pred$B
dat$ET = pred$T
dat$EM = pred$M 
head(dat)

#
# sample data (stratified)
#
#------------------------------

source("subsample.r")
#select2 = subsample.stratif3D(dataProj_subset10[,c("lon","lat", "annual_mean_temp")], subsetProp, adj = 4.2)

nsampl = as.integer((trunc(subsetProp*100)/100)*nrow(dataProj_subset10))
select2 = subsample.temp.fix(dataProj_subset10[,c("annual_mean_temp")], nsampl)

jpeg(paste("../figures/subsample_",subsetProp,".jpeg", sep=""), height=5000, width=5000, res=600)
plot(dataProj_subset10[,c("lon","lat")], pch = 20, cex=.2, col = "grey")
points(dataProj_subset10[select2,c("lon","lat")], pch = 20, cex=.2, col = 1)
dev.off()


datSel = dat[select2,]
datValid = dat[-select2,]


#-------------------------------------------------------------------

for(neiborgh in c("rf", "cst"))
{

# fit name
fit = paste(neiborgh, subsetProp, sep="_")
print(fit)

if(neiborgh=="cst") 
{
datSel$EB = rep(0.25, nrow(datSel))
datSel$ET = rep(0.25, nrow(datSel))
datSel$EM = rep(0.25, nrow(datSel))
}

# Evaluate initial parameter values
transitions = paste(datSel$st0,datSel$st1,sep = "")
sum_transitions = table(transitions)
initState = table(datSel$st0)

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
    

params = c(ab0 =as.numeric(logit_alphab_mn), ab1 = 0, ab2 = 0, ab3=0, ab4=0, ab5=0, ab6=0, 
at0 = as.numeric(logit_alphat_mn), at1 = 0, at2=0, at3=0, at4=0, at5=0, at6=0,
bb0 = as.numeric(logit_betab_mn), bb1 = 0, bb2 = 0, bb3=0, bb4=0, bb5=0, bb6=0,
bt0 = as.numeric(logit_betat_mn), bt1 = 0, bt2 = 0, bt3=0, bt4=0, bt5=0, bt6=0,
tt0 = as.numeric(logit_thetat_mn), tt1=0, tt2=0, tt3=0, tt4=0, tt5=0, tt6=0,
th0 = as.numeric(logit_theta_mn), th1=0, th2=0, th3=0, th4=0, th5=0, th6=0,
e0 = as.numeric(logit_eps_mn), e1 = 0,  e2=0, e3=0, e4=0, e5=0, e6=0)


scaleOfVar = c(rep(50, 49))

par_lo = 0 - scaleOfVar

par_hi = 0 + scaleOfVar

names(par_lo) = names(par_hi) = names(params)
#-----------

#coords = cbind(dataProj_subset10$longitude, dataProj_subset10$latitude)
save(datSel, dataProj_subset10, select2,params, toremove, select, par_lo, par_hi, fit, file = paste("initForFit_", fit,".RData",sep=""))

}

