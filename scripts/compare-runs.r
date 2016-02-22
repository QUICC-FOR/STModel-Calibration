rm(list = ls())

sdm = "rf"
propData = 0.33
ordre = 3
step = 5
set= 2

#---------------------------------------------------------------------
### alpha beta 
#---------------------------------------------------------------------
(name = paste(sdm,"_", propData, set, "_", ordre, "_",step, "y_alphabeta",sep=""))


folder_ab = "../estimated_params/alphabeta"

load("scale_info.Robj")
## Temp -4 à 10
## PP 700 à 1300
Trange = c(-4, 10)
Tticks = max(Trange)-min(Trange)+1
PPrange = c(700, 1300)
Tbounds = (Trange - vars.means["annual_mean_temp"])/vars.sd["annual_mean_temp"]
PPbounds = (PPrange - vars.means["tot_annual_pp"])/vars.sd["tot_annual_pp"]
PPticks = (max(PPrange)-min(PPrange))/100 +1

tpseq=seq(Tbounds[1],Tbounds[2],l=100)
ppseq=seq(PPbounds[1],PPbounds[2],l=100)

scaled.axis<- function(){
temp = tpseq*vars.sd["annual_mean_temp"]+vars.means["annual_mean_temp"]
precip = ppseq*vars.sd["tot_annual_pp"]+vars.means["tot_annual_pp"]
axis(1, at = seq(Tbounds[1], Tbounds[2], l=Tticks), labels = seq(Trange[1], Trange[2], l = Tticks))
axis(2, at = seq(PPbounds[1],PPbounds[2], l=PPticks), labels = seq(PPrange[1], PPrange[2], l = PPticks))
}

ENV = expand.grid(TP =tpseq , PP = ppseq)
ENV1 = ENV$TP
ENV2 = ENV$PP


###---- alphabeta
veget_pars_ab = read.table(paste(folder_ab, "/GenSA_", name, ".txt", sep=""))
load(paste(folder_ab, "/GenSA_", name, ".RData", sep = ""))
params = veget_pars_ab[,2]
 names(params) = c("ab0", "ab1", "ab2", "ab3","ab4","ab5","ab6",
"at0", "at1" , "at2", "at3", "at4", "at5", "at6", 
"bb0", "bb1", "bb2", "bb3", "bb4", "bb5", "bb6",
"bt0", "bt1", "bt2", "bt3", "bt4", "bt5", "bt6",
"tt0", "tt1", "tt2", "tt3", "tt4", "tt5", "tt6", 
"th0", "th1", "th2", "th3", "th4", "th5", "th6", 
"e0", "e1", "e2", "e3", "e4", "e5", "e6")

    logit_alphab 	= params["ab0"] + params["ab1"]*ENV1 + params["ab2"]*ENV2 + params["ab3"]*ENV1^2 + params["ab4"]*ENV2^2 + params["ab5"]*ENV1^3 + params["ab6"]*ENV2^3
    logit_alphat 	= params["at0"] + params["at1"]*ENV1 + params["at2"]*ENV2 + params["at3"]*ENV1^2 + params["at4"]*ENV2^2 + params["at5"]*ENV1^3 + params["at6"]*ENV2^3
    logit_betab 	= params["bb0"] + params["bb1"]*ENV1 + params["bb2"]*ENV2 + params["bb3"]*ENV1^2 + params["bb4"]*ENV2^2 + params["bb5"]*ENV1^3 + params["bb6"]*ENV2^3
    logit_betat 	= params["bt0"] + params["bt1"]*ENV1 + params["bt2"]*ENV2 + params["bt3"]*ENV1^2 + params["bt4"]*ENV2^2 + params["bt5"]*ENV1^3 + params["bt6"]*ENV2^3
    logit_theta	= params["th0"] + params["th1"]*ENV1 + params["th2"]*ENV2 + params["th3"]*ENV1^2 + params["th4"]*ENV2^2 + params["th5"]*ENV1^3 + params["th6"]*ENV2^3
    logit_thetat	= params["tt0"] + params["tt1"]*ENV1 + params["tt2"]*ENV2 + params["tt3"]*ENV1^2 + params["tt4"]*ENV2^2 + params["tt5"]*ENV1^3 + params["tt6"]*ENV2^3
    logit_eps 	= params["e0"]  + params["e1"]*ENV1 + params["e2"]*ENV2  + params["e3"]*ENV1^2 + params["e4"]*ENV2^2 + params["e5"]*ENV1^3 + params["e6"]*ENV2^3 

load(paste("initForFit_", sdm, "_", propData, set, ".RData", sep=""))
source(paste("3-transition_model_",ordre,".R", sep =""))
logll = model(estim.pars$par, dat, step = step)
logll
estim.pars$count

#------------------------------------------------------------------------------
# long vs short run
#---------------------------------------------------------------------
(name = paste(sdm,"_", propData, set, "_", ordre, "_",step, "y",sep=""))


folder_s = "../estimated_params/avr2015_shortRun"
folder_l = "../estimated_params"
load(paste("initForFit_", sdm, "_", propData, set, ".RData", sep=""))

veget_pars_s = read.table(paste(folder_s, "/GenSA_", name, ".txt", sep=""))
load(paste(folder_s, "/GenSA_", name, ".RData", sep = ""))
source(paste("3-transition_model_",ordre,".R", sep =""))
(logll_s = model(estim.pars$par, dat, step = step))
estim.pars$value

veget_pars_l = read.table(paste(folder_l, "/GenSA_", name, ".txt", sep=""))
load(paste(folder_l, "/GenSA_", name, ".RData", sep = ""))
(logll_l = model(estim.pars$par, dat, step = step))
estim.pars$value

plot(veget_pars_s[,2], veget_pars_l[,2])
cbind(veget_pars_s[,2], veget_pars_l[,2])


