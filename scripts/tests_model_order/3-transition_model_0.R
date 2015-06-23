# optimisation
# prerun the model with glm to find some parameter values that are not jointly constrained
# use logits


model = function(params, dat, step= 1) 

{
    st0 = dat$st0
    st1 = dat$st1
    ET = dat$ET
    EB = dat$EB
    EM = dat$EM
    ENV1 = dat$ENV1
    ENV2 = dat$ENV2
    itime = dat$itime
    
    names(params) = c("ab0", 
"at0", 
"bb0",
"bt0", 
"tt0", 
"th0", 
"e0")


    
	lik = numeric(length(st0))

    logit_alphab 	= params["ab0"] 
    logit_alphat 	= params["at0"] 
    logit_betab 	= params["bb0"] 
    logit_betat 	= params["bt0"] 
    logit_theta	= params["th0"] 
    logit_thetat	= params["tt0"] 
    logit_eps 	= params["e0"] 
    
    # compute transitions accounting for interval time and logit transformation
    # ! might be NaN because exp(bigNumber) 
    annualProba <- function(x, itime)
    {
    expx = ifelse(exp(x)==Inf, .Machine$double.xmax, exp(x))
    proba = expx/(1+expx)
#    annualProba = 1 - exp(log(1-proba)/itime)
    proba_itime = 1 - (1 - proba)^(itime/step)
    return(proba_itime)
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
	
#	print(sumLL)
    # return a value to minimize in GenSA
	return(-sumLL)
}


#-----------------------------------------------------------------------------


