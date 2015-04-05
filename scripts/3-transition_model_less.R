# optimisation
# prerun the model with glm to find some parameter values that are not jointly constrained
# use logits


model = function(params, dat) 

{
    st0 = dat$st0
    st1 = dat$st1
    ET = dat$ET
    EB = dat$EB
    EM = dat$EM
    ENV1 = dat$ENV1
    ENV2 = dat$ENV2
    itime = dat$itime
    
    names(params) = c("ab0", "ab1", "ab2", "ab3","ab4","ab5","ab6",
"at0", "at1" , "at2", "at3", "at4", "at5", "at6", 
"bb0", 
"bt0", 
"tt0", 
"th0", 
"e0")

#bb1 = bb2 = bb3 = bb4 = bb5 = bb6 = 0
#bt1 = bt2 = bt3 = bt4 = bt5 = bt6 = 0
#tt1 = tt2 = tt3 = tt4 = tt5 = tt6 = 0
#t1 = t2 = t3 = t4 = t5 = t6 = 0
#e1 = e2 = e3 = e4 = e5 = e6 = 0
#
#e7 = 0

    
	lik = numeric(length(st0))

    logit_alphab 	= params["ab0"] + params["ab1"]*ENV1 + params["ab2"]*ENV2 + params["ab3"]*ENV1^2 + params["ab4"]*ENV2^2 + params["ab5"]*ENV1^3 + params["ab6"]*ENV2^3
    logit_alphat 	= params["at0"] + params["at1"]*ENV1 + params["at2"]*ENV2 + params["at3"]*ENV1^2 + params["at4"]*ENV2^2 + params["at5"]*ENV1^3 + params["at6"]*ENV2^3
    logit_betab 	= params["bb0"] #+ params["bb1"]*ENV1 + params["bb2"]*ENV2 + params["bb3"]*ENV1^2 + params["bb4"]*ENV2^2 + params["bb5"]*ENV1^3 + params["bb6"]*ENV2^3
    logit_betat 	= params["bt0"] #+ params["bt1"]*ENV1 + params["bt2"]*ENV2 + params["bt3"]*ENV1^2 + params["bt4"]*ENV2^2 + params["bt5"]*ENV1^3 + params["bt6"]*ENV2^3
    logit_theta	= params["th0"] #+ params["t1"]*ENV1 + params["t2"]*ENV2 + params["t3"]*ENV1^2 + params["t4"]*ENV2^2 + params["t5"]*ENV1^3 + params["t6"]*ENV2^3
    logit_thetat	= params["tt0"] #+ params["tt1"]*ENV1 + params["tt2"]*ENV2 + params["tt3"]*ENV1^2 + params["tt4"]*ENV2^2 + params["tt5"]*ENV1^3 + params["tt6"]*ENV2^3
    logit_eps 	= params["e0"]  #+ params["e1"]*ENV1 + params["e2"]*ENV2  + params["e3"]*ENV1^2 + params["e4"]*ENV2^2 + params["e5"]*ENV1^3 + params["e6"]*ENV2^3 + e7*EB

    # compute transitions accounting for interval time and logit transformation
    # ! might be NaN because exp(bigNumber) 
    annualProba <- function(x, itime)
    {
    expx = ifelse(exp(x)==Inf, .Machine$double.xmax, exp(x))
#    1 - (1 - expx/(1+expx))^itime
    proba = expx/(1+expx)
    annualProba = 1 - exp(log(1-proba)/itime)
    return(annualProba)
    }
    alphab = annualProba(logit_alphab, itime)
    alphat = annualProba(logit_alphat, itime)
    betab = annualProba(logit_betab, itime)
    betat = annualProba(logit_betat, itime)
    theta = annualProba(logit_theta, itime)
    thetat = annualProba(logit_thetat, itime)
    eps = annualProba(logit_eps, itime)
       
    
	# Compute the likelihood of observations
	lik[st0 == "B" & st1 == "M"] = (betat*(ET+EM)*(1-eps))[st0 == "B" & st1 == "M"] 
	lik[st0 == "B" & st1 == "R"] = eps[st0 == "B" & st1 == "R"] 	
	lik[st0 == "B" & st1 == "B"] = (1 - eps - betat*(ET+EM)*(1-eps))[st0 == "B" & st1 == "B"]

	lik[st0 == "T" & st1 == "T"] = (1 - eps - betab*(EB+EM)*(1-eps))[st0 == "T" & st1 == "T"] 	
	lik[st0 == "T" & st1 == "M"] = (betab*(EB+EM)*(1-eps))[st0 == "T" & st1 == "M"] 			
	lik[st0 == "T" & st1 == "R"] = eps[st0 == "T" & st1 == "R"] 		
	
	lik[st0 == "M" & st1 == "B"] = (theta*(1-thetat)*(1-eps))[st0 == "M" & st1 == "B"]	
	lik[st0 == "M" & st1 == "T"] = (theta*thetat*(1-eps))[st0 == "M" & st1 == "T"] 	

	lik[st0 == "M" & st1 == "M"] = ((1 - eps)*(1 - theta))[st0 == "M" & st1 == "M"] 			

	lik[st0 == "M" & st1 == "R"] = eps[st0 == "M" & st1 == "R"] 
	
	phib = alphab*(EM + EB)*(1-alphat*(ET+EM))
	phit = alphat*(EM + ET)*(1-alphab*(EB+EM))
	phim = alphab*(EM + EB)*alphat*(EM + ET)
	lik[st0 == "R" & st1 == "B"] = phib[st0 == "R" & st1 == "B"] 	
	lik[st0 == "R" & st1 == "T"] = phit[st0 == "R" & st1 == "T"]	
	lik[st0 == "R" & st1 == "M"] = phim[st0 == "R" & st1 == "M"] 			
	lik[st0 == "R" & st1 == "R"] = (1 - phib - phit - phim)[st0 == "R" & st1 == "R"] 


    # lik might be <0 !!!!
    # for T->T, M->M, B->B and R->R transitions
    # it will give NaN when log
    

    # lik might be equal = 0 -> give -Inf at log   
    # for instance when neighbor (seeds) = 0
    #lik[lik == 0] = .Machine$double.xmin
    
    # calculate sum log likelihood
    # lik might be <0 (in the 1 - other proba)
#    if(sum(lik<0)>0)
#    {
#    sumLL = -.Machine$double.xmax
#    }else{
    sumLL = sum(log(lik))
#    }
    
    if(is.infinite(sumLL)) sumLL = -.Machine$double.xmax
	if(is.nan(sumLL)) sumLL = -.Machine$double.xmax
	if(is.na(sumLL)) print("sumLL: na values!")
#print(sumLL)	
    # return a value to minimize in GenSA
	return(-sumLL)
}


#-----------------------------------------------------------------------------


