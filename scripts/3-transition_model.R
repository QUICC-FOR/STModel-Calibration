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

    ab0 = params["a0"]
    ab1 = params[2]
    ab2 = params[3]
    ab3 = params[4]
    ab4 = params[5]
    ab5 = params[6]
    ab6 = params[7]
    at0 = params[8]
    at1 = params[9]
    at2 = params[10]
    at3 = params[11]
    at4 = params[12]
    at5 = params[13]
    at6 = params[14]
    bb0 = params[15]
    bb1 = params[16]
    bb2 = params[17]
    bb3 = params[18]
    bb4 = params[19]
    bb5 = params[20]
    bb6 = params[21]
    bt0 = params[22]
    bt1 = params[23]
    bt2 = params[24]
    bt3 = params[25]
    bt4 = params[26]
    bt5 = params[27]
    bt6 = params[28]
    tt0 = params[29]
    tt1 = params[30]
    tt2 = params[31]
    tt3 = params[32]
    tt4 = params[33]
    tt5 = params[34]
    tt6 = params[35]
    t0 = params[36]
    t1 = params[37]
    t2 = params[38]
    t3 = params[39]
    t4 = params[40]
    t5 = params[41]
    t6 = params[42]
    e0 = params[43]
    e1 = params[44]
    e2 = params[45]
    e3 = params[46]
    e4 = params[47]
    e5 = params[48]
    e6 = params[49]
#    e7 = params[50]
    
	lik = numeric(length(st0))

    logit_alphab 	= ab0 + ab1*ENV1 + ab2*ENV2 + ab3*ENV1^2 + ab4*ENV2^2 + ab5*ENV1^3 + ab6*ENV2^3
    logit_alphat 	= at0 + at1*ENV1 + at2*ENV2 + at3*ENV1^2 + at4*ENV2^2 + at5*ENV1^3 + at6*ENV2^3
    logit_betab 	= bb0 + bb1*ENV1 + bb2*ENV2 + bb3*ENV1^2 + bb4*ENV2^2 + bb5*ENV1^3 + bb6*ENV2^3
    logit_betat 	= bt0 + bt1*ENV1 + bt2*ENV2 + bt3*ENV1^2 + bt4*ENV2^2 + bt5*ENV1^3 + bt6*ENV2^3
    logit_theta	= t0 + t1*ENV1 + t2*ENV2 + t3*ENV1^2 + t4*ENV2^2 + t5*ENV1^3 + t6*ENV2^3
    logit_thetat	= tt0 + tt1*ENV1 + tt2*ENV2 + tt3*ENV1^2 + tt4*ENV2^2 + tt5*ENV1^3 + tt6*ENV2^3
    logit_eps 	= e0  + e1*ENV1 + e2*ENV2  + e3*ENV1^2 + e4*ENV2^2 + e5*ENV1^3 + e6*ENV2^3  
#    + e7*EB

    # compute transitions accounting for interval time and logit transformation
    # ! might be NaN because exp(bigNumber) 
    annualProba <- function(x, itime)
    {
    expx = ifelse(exp(x)==Inf, .Machine$double.xmax, exp(x))
    1 - (1 - expx/(1+expx))^itime
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
    lik[lik == 0] = .Machine$double.xmin
    
    # calculate sum log likelihood
    # lik might be <0 (in the 1 - other proba)
    if(sum(lik<0)>0)
    {
    sumLL = -1000000
    }else{
    sumLL = sum(log(lik))
    }
    

    # return a value to minimize in GenSA
	return(-sumLL)
}


#-----------------------------------------------------------------------------


