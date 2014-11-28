rm(list=ls())
# state transitions
data = read.csv("../data/transitionsFourState.csv")
head(data)
dim(data)
data$annual_mean_temp = apply(data[,c("annual_mean_temp1", "annual_mean_temp2.")], 1, mean)
data$annual_pp = apply(data[,c("annual_pp1", "annual_pp2")], 1, mean)

# neighborhood
pred = read.table("../data/projection_neigbor_rf_temp.txt", h=T)
head(pred)

# remove transitions T->B and B->T
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
data$st0 = data$state1
data$st1 = data$state2
data$itime = data$year2 - data$year1
sum(data$itime<0)

# Evaluate initial parameter values
transitions = paste(data$st0,data$st1,sep = "")
sum_transitions = table(transitions)
initState = table(data$st0)

# estimates for initial parameters
eps_mn = # proba x->R given x is T B or M
(sum_transitions["BR"]/initState["B"]+sum_transitions["MR"]/initState["M"]+sum_transitions["TR"]/initState["T"])/3

theta_mn = # proba M-> B 
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
params = c(ab0 =logit_alphab_mn, ab1 = 0, ab2 = 0, ab3=0, ab4=0, ab5=0, ab6=0,
at0 = logit_alphat_mn, at1 = 0, at2 = 0, at3=0, at4=0, at5=0, at6=0,
bb0 = logit_betab_mn, bb1 = 0, bb2 = 0, bb3=0, bb4=0, bb5=0, bb6=0,
bt0 = logit_betat_mn, bt1 = 0, bt2 = 0, bt3=0, bt4=0, bt5=0, bt6=0,
tt0 = logit_thetat_mn, tt1 = 0, tt2 = 0, tt3 =0, tt4=0, tt5=0, tt6=0,
t0 = logit_theta_mn, t1 = 0, t2 = 0, t3=0, t4=0, t5=0, t6=0,
e0 = logit_eps_mn, e1 = 0, e2 = 0, e3=0, e4=0, e5=0, e6=0, e7 =0
)

# coeff variation
cvar = 5

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
t0 = logit_theta_mn - abs(logit_theta_mn*cvar),
t1 = logit_theta_mn/max(abs(data$ENV1)) - abs(cvar*logit_theta_mn/max(abs(data$ENV1))),
t2 = logit_theta_mn/max(abs(data$ENV2)) - abs(cvar*logit_theta_mn/max(abs(data$ENV2))),
t3 = logit_theta_mn/max(abs(data$ENV1^2)) - abs(cvar*logit_theta_mn/max(abs(data$ENV1^2))),
t4 = logit_theta_mn/max(abs(data$ENV2^2)) - abs(cvar*logit_theta_mn/max(abs(data$ENV2^2))),
t5= logit_theta_mn/max(abs(data$ENV1^3)) - abs(cvar*logit_theta_mn/max(abs(data$ENV1^3))),
t6= logit_theta_mn/max(abs(data$ENV2^3)) - abs(cvar*logit_theta_mn/max(abs(data$ENV2^3))),
e0 = logit_eps_mn - abs(logit_eps_mn*cvar),
e1 = logit_eps_mn/max(abs(data$ENV1)) - abs(cvar*logit_eps_mn/max(abs(data$ENV1))),
e2 = logit_eps_mn/max(abs(data$ENV2)) - abs(cvar*logit_eps_mn/max(abs(data$ENV2))), 
e3 = logit_eps_mn/max(abs(data$ENV1^2)) - abs(cvar*logit_eps_mn/max(abs(data$ENV1^2))),
e4 = logit_eps_mn/max(abs(data$ENV2^2)) - abs(cvar*logit_eps_mn/max(abs(data$ENV2^2))), 
e5= logit_eps_mn/max(abs(data$ENV1^3)) - abs(cvar*logit_eps_mn/max(abs(data$ENV1^3))),
e6= logit_eps_mn/max(abs(data$ENV2^3)) - abs(cvar*logit_eps_mn/max(abs(data$ENV2^3))),
e7=-1e15
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
t0 = logit_theta_mn + abs(logit_theta_mn*cvar),
t1 = logit_theta_mn/max(abs(data$ENV1)) + abs(cvar*logit_theta_mn/max(abs(data$ENV1))),
t2 = logit_theta_mn/max(abs(data$ENV2)) + abs(cvar*logit_theta_mn/max(abs(data$ENV2))),
t3 = logit_theta_mn/max(abs(data$ENV1^2)) + abs(cvar*logit_theta_mn/max(abs(data$ENV1^2))),
t4 = logit_theta_mn/max(abs(data$ENV2^2)) + abs(cvar*logit_theta_mn/max(abs(data$ENV2^2))),
t5= logit_theta_mn/max(abs(data$ENV1^3)) + abs(cvar*logit_theta_mn/max(abs(data$ENV1^3))),
t6= logit_theta_mn/max(abs(data$ENV2^3)) + abs(cvar*logit_theta_mn/max(abs(data$ENV2^3))),
e0 = logit_eps_mn + abs(logit_eps_mn*cvar),
e1 = logit_eps_mn/max(abs(data$ENV1)) + abs(cvar*logit_eps_mn/max(abs(data$ENV1))),
e2 = logit_eps_mn/max(abs(data$ENV2)) + abs(cvar*logit_eps_mn/max(abs(data$ENV2))), 
e3 = logit_eps_mn/max(abs(data$ENV1^2)) + abs(cvar*logit_eps_mn/max(abs(data$ENV1^2))),
e4 = logit_eps_mn/max(abs(data$ENV2^2)) + abs(cvar*logit_eps_mn/max(abs(data$ENV2^2))), 
e5= logit_eps_mn/max(abs(data$ENV1^3)) + abs(cvar*logit_eps_mn/max(abs(data$ENV1^3))),
e6= logit_eps_mn/max(abs(data$ENV2^3)) + abs(cvar*logit_eps_mn/max(abs(data$ENV2^3))),
e7=1e-15
)


