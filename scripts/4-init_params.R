# state transitions
data = as.data.frame(read.table("../data/data_reshaped_RBTM.txt"))
# neighborhood
pred = read.table("../data/data_pred_states_multinom.txt")
#pred = read.table("../data/data_pred_states_randomForest.txt")
#
# remove transitions D->C and C->D
test = numeric(nrow(data))
test[data$st0 == "T" & data$st1 == "B"] = 1
test[data$st0 == "B" & data$st1 == "T"] = 1
data = subset(data, test!=1)
pred = subset(pred, test!=1)

data$ENV1 = scale(data$annual_mean_temp)
data$ENV2 = scale(data$annual_pp)
data$EB = pred$B
data$ET = pred$T
data$EM = pred$M 
#data$Hv = rep(0, nrow(data))
#data$Ha = rep(0, nrow(data))

head(data)

# Evaluate initial parameter values
transitions = paste(data$st0,data$st1,sep = "")
sum_transitions = table(transitions)
tot_transitions = table(data$st0)

# estimates for initial parameters
eps_mn = # proba x->R given x is T B or M
(sum_transitions["BR"]+sum_transitions["MR"]+sum_transitions["TR"])/(tot_transitions["M"] + tot_transitions["B"] + tot_transitions["T"])

thetab_mn = # proba M-> B 
sum_transitions["MB"]/tot_transitions["M"]

thetat_mn = # proba M->T
sum_transitions["MT"]/tot_transitions["M"]

betab_mn = # proba T->M
sum_transitions["TM"]*(tot_transitions["B"]+tot_transitions["M"])/tot_transitions["M"]/sum(tot_transitions)

betat_mn = # proba B->M
sum_transitions["BM"]*(tot_transitions["T"]+tot_transitions["M"])/tot_transitions["M"]/sum(tot_transitions)

trRT = sum_transitions["RT"]/tot_transitions["R"]
trRB = sum_transitions["RB"]/tot_transitions["R"]
trRM = sum_transitions["RB"]/tot_transitions["R"]

alphat_mn = trRM/((tot_transitions["M"]+tot_transitions["T"])*(trRB+trRM))
alphab_mn = (trRB+trRM)/(tot_transitions["M"] + tot_transitions["B"])

# transform into logits
logit_eps_mn = log(eps_mn/(1-eps_mn))

logit_thetab_mn = log(thetab_mn/(1-thetab_mn))
logit_thetat_mn = log(thetat_mn/(1-thetat_mn))
 
logit_betab_mn = log(betab_mn/(1-betab_mn))
logit_betat_mn = log(betat_mn/(1-betat_mn))

logit_alphab_mn = log(alphab_mn/(1-alphab_mn))
logit_alphat_mn = log(alphat_mn/(1-alphat_mn))



# List initial parameters
params = list(ab0 =logit_alphab_mn, ab1 = 0, ab2 = 0, ab3=0, ab4=0, ab5=0, ab6=0,
at0 = logit_alphat_mn, at1 = 0, at2 = 0, at3=0, at4=0, at5=0, at6=0,
bb0 = logit_betab_mn, bb1 = 0, bb2 = 0, bb3=0, bb4=0, bb5=0, bb6=0,
bt0 = logit_betat_mn, bt1 = 0, bt2 = 0, bt3=0, bt4=0, bt5=0, bt6=0,
tt0 = logit_thetat_mn, tt1 = 0, tt2 = 0, tt3 =0, tt4=0, tt5=0, tt6=0,
tb0 = logit_thetab_mn, tb1 = 0, tb2 = 0, tb3=0, tb4=0, tb5=0, tb6=0,
e0 = logit_eps_mn, e1 = 0, e2 = 0, e3=0, e4=0, e5=0, e6=0
)

# coeff variation
cvar = 2

par_lo = list(ab0 = logit_alphab_mn - abs(logit_alphab_mn*cvar),
ab1 = logit_alphab_mn/max(abs(data$ENV1)) - abs(cvar*logit_alphab_mn/max(abs(data$ENV1))),
ab2 = logit_alphab_mn/max(abs(data$ENV2)) - abs(cvar*logit_alphab_mn/max(abs(data$ENV2))), 
ab3 = logit_alphab_mn/max(abs(data$ENV1^2)) - abs(cvar*logit_alphab_mn/max(abs(data$ENV1^2))), 
ab4= logit_alphab_mn/max(abs(data$ENV2^2)) - abs(cvar*logit_alphab_mn/max(abs(data$ENV2^2))),
ab5= logit_alphab_mn/max(abs(data$ENV1^3)) - abs(cvar*logit_alphab_mn/max(abs(data$ENV1^3))),
ab6= logit_alphab_mn/max(abs(data$ENV2^3)) - abs(cvar*logit_alphab_mn/max(abs(data$ENV2^3))),
at0 = logit_alphat_mn - abs(logit_alphat_mn*cvar),
at1 = logit_alphat_mn/max(abs(data$ENV1)) - abs(cvar*logit_alphat_mn/max(abs(data$ENV1))),
at2 = logit_alphat_mn/max(abs(data$ENV2)) - abs(cvar*logit_alphat_mn/max(abs(data$ENV2))),
at3 = logit_alphat_mn/max(abs(data$ENV1^2)) - abs(cvar*logit_alphat_mn/max(abs(data$ENV1^2))),
at4 = logit_alphat_mn/max(abs(data$ENV2^2)) - abs(cvar*logit_alphat_mn/max(abs(data$ENV2^2))),
at5= logit_alphat_mn/max(abs(data$ENV1^3)) - abs(cvar*logit_alphat_mn/max(abs(data$ENV1^3))),
at6= logit_alphat_mn/max(abs(data$ENV2^3)) - abs(cvar*logit_alphat_mn/max(abs(data$ENV2^3))),
bb0 = logit_betab_mn - abs(logit_betab_mn*cvar),
bb1 = logit_betab_mn/max(abs(data$ENV1)) - abs(cvar*logit_betab_mn/max(abs(data$ENV1))),
bb2 = logit_betab_mn/max(abs(data$ENV2)) - abs(cvar*logit_betab_mn/max(abs(data$ENV2))),
bb3 = logit_betab_mn/max(abs(data$ENV1^2)) - abs(cvar*logit_betab_mn/max(abs(data$ENV1^2))),
bb4 = logit_betab_mn/max(abs(data$ENV2^2)) - abs(cvar*logit_betab_mn/max(abs(data$ENV2^2))),
bb5= logit_betab_mn/max(abs(data$ENV1^3)) - abs(cvar*logit_betab_mn/max(abs(data$ENV1^3))),
bb6= logit_betab_mn/max(abs(data$ENV2^3)) - abs(cvar*logit_betab_mn/max(abs(data$ENV2^3))),
bt0 = logit_betat_mn - abs(logit_betat_mn*cvar),
bt1 = logit_betat_mn/max(abs(data$ENV1)) - abs(cvar*logit_betat_mn/max(abs(data$ENV1))),
bt2 =logit_betat_mn/max(abs(data$ENV2)) - abs(cvar*logit_betat_mn/max(abs(data$ENV2))),
bt3 = logit_betat_mn/max(abs(data$ENV1^2)) - abs(cvar*logit_betat_mn/max(abs(data$ENV1^2))),
bt4 = logit_betat_mn/max(abs(data$ENV2^2)) - abs(cvar*logit_betat_mn/max(abs(data$ENV2^2))),
bt5= logit_betat_mn/max(abs(data$ENV1^3)) - abs(cvar*logit_betat_mn/max(abs(data$ENV1^3))),
bt6= logit_betat_mn/max(abs(data$ENV2^3)) - abs(cvar*logit_betat_mn/max(abs(data$ENV2^3))),
tt0 = logit_thetat_mn - abs(logit_thetat_mn*cvar),
tt1 = logit_thetat_mn/max(abs(data$ENV1)) - abs(cvar*logit_thetat_mn/max(abs(data$ENV1))),
tt2 = logit_thetat_mn/max(abs(data$ENV2)) - abs(cvar*logit_thetat_mn/max(abs(data$ENV2))),
tt3 = logit_thetat_mn/max(abs(data$ENV1^2)) - abs(cvar*logit_thetat_mn/max(abs(data$ENV1^2))),
tt4 = logit_thetat_mn/max(abs(data$ENV2^2)) - abs(cvar*logit_thetat_mn/max(abs(data$ENV2^2))),
tt5= logit_thetat_mn/max(abs(data$ENV1^3)) - abs(cvar*logit_thetat_mn/max(abs(data$ENV1^3))),
tt6= logit_thetat_mn/max(abs(data$ENV2^3)) - abs(cvar*logit_thetat_mn/max(abs(data$ENV2^3))),
tb0 = logit_thetab_mn - abs(logit_thetab_mn*cvar),
tb1 = logit_thetab_mn/max(abs(data$ENV1)) - abs(cvar*logit_thetab_mn/max(abs(data$ENV1))),
tb2 = logit_thetab_mn/max(abs(data$ENV2)) - abs(cvar*logit_thetab_mn/max(abs(data$ENV2))),
tb3 = logit_thetab_mn/max(abs(data$ENV1^2)) - abs(cvar*logit_thetab_mn/max(abs(data$ENV1^2))),
tb4 = logit_thetab_mn/max(abs(data$ENV2^2)) - abs(cvar*logit_thetab_mn/max(abs(data$ENV2^2))),
tb5= logit_thetab_mn/max(abs(data$ENV1^3)) - abs(cvar*logit_thetab_mn/max(abs(data$ENV1^3))),
tb6= logit_thetab_mn/max(abs(data$ENV2^3)) - abs(cvar*logit_thetab_mn/max(abs(data$ENV2^3))),
e0 = logit_eps_mn - abs(logit_eps_mn*cvar),
e1 = logit_eps_mn/max(abs(data$ENV1)) - abs(cvar*logit_eps_mn/max(abs(data$ENV1))),
e2 = logit_eps_mn/max(abs(data$ENV2)) - abs(cvar*logit_eps_mn/max(abs(data$ENV2))), 
e3 = logit_eps_mn/max(abs(data$ENV1^2)) - abs(cvar*logit_eps_mn/max(abs(data$ENV1^2))),
e4 = logit_eps_mn/max(abs(data$ENV2^2)) - abs(cvar*logit_eps_mn/max(abs(data$ENV2^2))), 
e5= logit_eps_mn/max(abs(data$ENV1^3)) - abs(cvar*logit_eps_mn/max(abs(data$ENV1^3))),
e6= logit_eps_mn/max(abs(data$ENV2^3)) - abs(cvar*logit_eps_mn/max(abs(data$ENV2^3)))
)

par_hi = list(ab0 = logit_alphab_mn + abs(logit_alphab_mn*cvar),
ab1 = logit_alphab_mn/max(abs(data$ENV1)) + abs(cvar*logit_alphab_mn/max(abs(data$ENV1))),
ab2 = logit_alphab_mn/max(abs(data$ENV2)) + abs(cvar*logit_alphab_mn/max(abs(data$ENV2))), 
ab3 = logit_alphab_mn/max(abs(data$ENV1^2)) + abs(cvar*logit_alphab_mn/max(abs(data$ENV1^2))), 
ab4= logit_alphab_mn/max(abs(data$ENV2^2)) + abs(cvar*logit_alphab_mn/max(abs(data$ENV2^2))),
ab5= logit_alphab_mn/max(abs(data$ENV1^3)) + abs(cvar*logit_alphab_mn/max(abs(data$ENV1^3))),
ab6= logit_alphab_mn/max(abs(data$ENV2^3)) + abs(cvar*logit_alphab_mn/max(abs(data$ENV2^3))),
at0 = logit_alphat_mn + abs(logit_alphat_mn*cvar),
at1 = logit_alphat_mn/max(abs(data$ENV1)) + abs(cvar*logit_alphat_mn/max(abs(data$ENV1))),
at2 = logit_alphat_mn/max(abs(data$ENV2)) + abs(cvar*logit_alphat_mn/max(abs(data$ENV2))),
at3 = logit_alphat_mn/max(abs(data$ENV1^2)) + abs(cvar*logit_alphat_mn/max(abs(data$ENV1^2))),
at4 = logit_alphat_mn/max(abs(data$ENV2^2)) + abs(cvar*logit_alphat_mn/max(abs(data$ENV2^2))),
at5= logit_alphat_mn/max(abs(data$ENV1^3)) + abs(cvar*logit_alphat_mn/max(abs(data$ENV1^3))),
at6= logit_alphat_mn/max(abs(data$ENV2^3)) + abs(cvar*logit_alphat_mn/max(abs(data$ENV2^3))),
bb0 = logit_betab_mn + abs(logit_betab_mn*cvar),
bb1 = logit_betab_mn/max(abs(data$ENV1)) + abs(cvar*logit_betab_mn/max(abs(data$ENV1))),
bb2 = logit_betab_mn/max(abs(data$ENV2)) + abs(cvar*logit_betab_mn/max(abs(data$ENV2))),
bb3 = logit_betab_mn/max(abs(data$ENV1^2)) + abs(cvar*logit_betab_mn/max(abs(data$ENV1^2))),
bb4 = logit_betab_mn/max(abs(data$ENV2^2)) + abs(cvar*logit_betab_mn/max(abs(data$ENV2^2))),
bb5= logit_betab_mn/max(abs(data$ENV1^3)) + abs(cvar*logit_betab_mn/max(abs(data$ENV1^3))),
bb6= logit_betab_mn/max(abs(data$ENV2^3)) + abs(cvar*logit_betab_mn/max(abs(data$ENV2^3))),
bt0 = logit_betat_mn + abs(logit_betat_mn*cvar),
bt1 = logit_betat_mn/max(abs(data$ENV1)) + abs(cvar*logit_betat_mn/max(abs(data$ENV1))),
bt2 =logit_betat_mn/max(abs(data$ENV2)) + abs(cvar*logit_betat_mn/max(abs(data$ENV2))),
bt3 = logit_betat_mn/max(abs(data$ENV1^2)) + abs(cvar*logit_betat_mn/max(abs(data$ENV1^2))),
bt4 = logit_betat_mn/max(abs(data$ENV2^2)) + abs(cvar*logit_betat_mn/max(abs(data$ENV2^2))),
bt5= logit_betat_mn/max(abs(data$ENV1^3)) + abs(cvar*logit_betat_mn/max(abs(data$ENV1^3))),
bt6= logit_betat_mn/max(abs(data$ENV2^3)) + abs(cvar*logit_betat_mn/max(abs(data$ENV2^3))),
tt0 = logit_thetat_mn + abs(logit_thetat_mn*cvar),
tt1 = logit_thetat_mn/max(abs(data$ENV1)) + abs(cvar*logit_thetat_mn/max(abs(data$ENV1))),
tt2 = logit_thetat_mn/max(abs(data$ENV2)) + abs(cvar*logit_thetat_mn/max(abs(data$ENV2))),
tt3 = logit_thetat_mn/max(abs(data$ENV1^2)) + abs(cvar*logit_thetat_mn/max(abs(data$ENV1^2))),
tt4 = logit_thetat_mn/max(abs(data$ENV2^2)) + abs(cvar*logit_thetat_mn/max(abs(data$ENV2^2))),
tt5= logit_thetat_mn/max(abs(data$ENV1^3)) + abs(cvar*logit_thetat_mn/max(abs(data$ENV1^3))),
tt6= logit_thetat_mn/max(abs(data$ENV2^3)) + abs(cvar*logit_thetat_mn/max(abs(data$ENV2^3))),
tb0 = logit_thetab_mn + abs(logit_thetab_mn*cvar),
tb1 = logit_thetab_mn/max(abs(data$ENV1)) + abs(cvar*logit_thetab_mn/max(abs(data$ENV1))),
tb2 = logit_thetab_mn/max(abs(data$ENV2)) + abs(cvar*logit_thetab_mn/max(abs(data$ENV2))),
tb3 = logit_thetab_mn/max(abs(data$ENV1^2)) + abs(cvar*logit_thetab_mn/max(abs(data$ENV1^2))),
tb4 = logit_thetab_mn/max(abs(data$ENV2^2)) + abs(cvar*logit_thetab_mn/max(abs(data$ENV2^2))),
tb5= logit_thetab_mn/max(abs(data$ENV1^3)) + abs(cvar*logit_thetab_mn/max(abs(data$ENV1^3))),
tb6= logit_thetab_mn/max(abs(data$ENV2^3)) + abs(cvar*logit_thetab_mn/max(abs(data$ENV2^3))),
e0 = logit_eps_mn + abs(logit_eps_mn*cvar),
e1 = logit_eps_mn/max(abs(data$ENV1)) + abs(cvar*logit_eps_mn/max(abs(data$ENV1))),
e2 = logit_eps_mn/max(abs(data$ENV2)) + abs(cvar*logit_eps_mn/max(abs(data$ENV2))), 
e3 = logit_eps_mn/max(abs(data$ENV1^2)) + abs(cvar*logit_eps_mn/max(abs(data$ENV1^2))),
e4 = logit_eps_mn/max(abs(data$ENV2^2)) + abs(cvar*logit_eps_mn/max(abs(data$ENV2^2))), 
e5= logit_eps_mn/max(abs(data$ENV1^3)) + abs(cvar*logit_eps_mn/max(abs(data$ENV1^3))),
e6= logit_eps_mn/max(abs(data$ENV2^3)) + abs(cvar*logit_eps_mn/max(abs(data$ENV2^3)))
)

par_initStep =lapply(names(params), function(x){
max(abs(as.numeric(params[x])-as.numeric(par_lo[x])), abs(as.numeric(params[x])-as.numeric(par_hi[x])))
})
names(par_initStep) = names(params)



#ENV1 = data$ENV1
#ENV2 =  data$ENV2
#save(ENV1, ENV2, file ="ENV.RData")

testBounds =
function(par, data)
{
    res = TRUE
    ENV1 = data$ENV1
    ENV2 = data$ENV2
    Hv = data$Hv
    Ha = data$Ha
    ET = data$ET
    EB = data$EB
    EM = data$EM
    

    logit_alphab 	= par$ab0 + par$ab1*ENV1 + par$ab2*ENV2 + par$ab3*ENV1^2 + par$ab4*ENV2^2 + par$ab5*ENV1^3 + par$ab6*ENV2^3
    logit_alphat 	= par$at0 + par$at1*ENV1 + par$at2*ENV2 + par$at3*ENV1^2 + par$at4*ENV2^2 + par$at5*ENV1^3 + par$at6*ENV2^3
    logit_betab 	= par$bb0 + par$bb1*ENV1 + par$bb2*ENV2 + par$bb3*ENV1^2 + par$bb4*ENV2^2 + par$bb5*ENV1^3 + par$bb6*ENV2^3
    logit_betat 	= par$bt0 + par$bt1*ENV1 + par$bt2*ENV2 + par$bt3*ENV1^2 + par$bt4*ENV2^2 + par$bt5*ENV1^3 + par$bt6*ENV2^3
    logit_thetab	= par$tb0 + par$tb1*ENV1 + par$tb2*ENV2 + par$tb3*ENV1^2 + par$tb4*ENV2^2 + par$tb5*ENV1^3 + par$tb6*ENV2^3
    logit_thetat	= par$tt0 + par$tt1*ENV1 + par$tt2*ENV2 + par$tt3*ENV1^2 + par$tt4*ENV2^2 + par$tt5*ENV1^3 +par$ tt6*ENV2^3
    logit_eps 	= par$e0  + par$e1*ENV1 + par$e2*ENV2  + par$e3*ENV1^2 + par$e4*ENV2^2 + par$e5*ENV1^3 + par$e6*ENV2^3
   
    alphab = exp(logit_alphab)/(1+exp(logit_alphab))
    alphat = exp(logit_alphat)/(1+exp(logit_alphat))
    betab = exp(logit_betab)/(1+exp(logit_betab))
    betat = exp(logit_betat)/(1+exp(logit_betat))
    thetab = exp(logit_thetab)/(1+exp(logit_thetab))
    thetat = exp(logit_thetat)/(1+exp(logit_thetat))
    eps = exp(logit_eps)/(1+exp(logit_eps))

    alphab = exp(logit_alphab)/(1+exp(logit_alphab))
    alphat = exp(logit_alphat)/(1+exp(logit_alphat))
    betab = exp(logit_betab)/(1+exp(logit_betab))
    betat = exp(logit_betat)/(1+exp(logit_betat))
    thetab = exp(logit_thetab)/(1+exp(logit_thetab))
    thetat = exp(logit_thetat)/(1+exp(logit_thetat))
    eps = exp(logit_eps)/(1+exp(logit_eps))

    
#    if(sum(alphab==0)>0 | sum(alphab==1)>0 ) res=FALSE
#    if(sum(alphat==0)>0 | sum(alphat==1)>0 ) res=FALSE
#    if(sum(betab==0)>0 | sum(betab==1)>0 ) res=FALSE
#    if(sum(betat==0)>0 | sum(betat==1)>0 ) res=FALSE
#    if(sum(thetab==0)>0 | sum(thetab==1)>0 ) res=FALSE
#    if(sum(thetat==0)>0 | sum(thetat==1)>0 ) res=FALSE
#    if(sum(eps==0)>0 | sum(eps==1)>0 ) res=FALSE   

# conditions for positive probabilities of staying in one state 
# (e.g. (1-eps-thetat)>1)
    condition1 = (eps + betat*(ET+EM)) > 1
    condition2 = (eps + betab*(EB+EM)) > 1
    condition3 = (eps + thetat + thetab) > 1
    condition4 = (alphab*(EB+EM) + alphat*(ET+EM) + alphab*(EB+EM)*alphat*(ET+EM)) > 1
    
    if(sum(condition1, condition2, condition3, condition4)>0) res=FALSE
     
return(res)
}





