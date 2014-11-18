# state transitions
data = as.data.frame(read.table("../data/data_pairs_filter.txt"))
# neighborhood
pred = read.table("../data/pred_states_multinom.txt")
#pred = read.table("../data/data_pred_states_randomForest.txt")


# remove transitions D->C and C->D
test = numeric(nrow(data))
test[data$st0 == "T" & data$st1 == "B"] = 1
test[data$st0 == "B" & data$st1 == "T"] = 1
data = subset(data, test!=1)
pred = subset(pred, test!=1)

#check variables
data$ENV1 = scale(data$annual_mean_temp)
data$ENV2 = scale(data$annual_pp)
data$EB = pred$B
data$ET = pred$T
data$EM = pred$M 
data$st0
data$st1
data$itime = data$yr1 - data$yr0


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
params = c(ab0 =logit_alphab_mn, ab1 = 0, ab2 = 0, ab3=0, ab4=0, ab5=0, ab6=0,
at0 = logit_alphat_mn, at1 = 0, at2 = 0, at3=0, at4=0, at5=0, at6=0,
bb0 = logit_betab_mn, bb1 = 0, bb2 = 0, bb3=0, bb4=0, bb5=0, bb6=0,
bt0 = logit_betat_mn, bt1 = 0, bt2 = 0, bt3=0, bt4=0, bt5=0, bt6=0,
tt0 = logit_thetat_mn, tt1 = 0, tt2 = 0, tt3 =0, tt4=0, tt5=0, tt6=0,
tb0 = logit_thetab_mn, tb1 = 0, tb2 = 0, tb3=0, tb4=0, tb5=0, tb6=0,
e0 = logit_eps_mn, e1 = 0, e2 = 0, e3=0, e4=0, e5=0, e6=0, e7 =0
)

# coeff variation
cvar = 2

par_lo = c(ab0 = logit_alphab_mn - abs(logit_alphab_mn*cvar),
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
e6= logit_eps_mn/max(abs(data$ENV2^3)) - abs(cvar*logit_eps_mn/max(abs(data$ENV2^3))),
e7=0
)

par_hi = c(ab0 = logit_alphab_mn + abs(logit_alphab_mn*cvar),
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
e6= logit_eps_mn/max(abs(data$ENV2^3)) + abs(cvar*logit_eps_mn/max(abs(data$ENV2^3))),
e7=1e-15
)


