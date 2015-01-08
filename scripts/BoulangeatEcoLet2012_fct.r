



#------------------------------------------------------------------------------------------
# evaluation functions
#------------------------------------------------------------------------------------------
#True Skill Statistic
TSS <-
function (Pred, Obs)
{
 	Misc = unclass(table(Pred, Obs))

    if (dim(Misc)[1] == 1) {
        if (row.names(Misc)[1] == "FALSE")
            Misc <- rbind(Misc, c(0, 0))
        else {
            a <- Misc
            Misc <- c(0, 0)
            Misc <- rbind(Misc, a)
        }
    }
    n <- sum(Misc)
    d <- Misc[1, 1]
    c <- Misc[1, 2]
    b <- Misc[2, 1]
    a <- Misc[2, 2]
    sens <- a/(a + c)
    spec <- d/(b + d)
    K <- (sens + spec) - 1
    return(K)
}

#Hanssen-Kuipers score
#HK <-
#function (Pred, Obs)
#{
#	Misc = table(Pred, Obs)
#
#    if (nrow(Misc)!=ncol(Misc)) stop("wrong misclassification table")
#    Misc <- unclass(Misc)
#    k  <- ncol(Misc)
#    Nobs <- apply(Misc, 2, sum)
#    Npred <- apply(Misc, 1, sum)
#    N <- sum(Nobs)
#
#   HK <- (sum(diag(Misc))/N - sum(Nobs*Npred)/N/N ) / ( 1 - sum(Nobs*Nobs)/N/N )
#
#    return(HK)
#}
#


#
#
# #------------------------------------------------------------------------------------------
# # From probabilities to classes
# #------------------------------------------------------------------------------------------
# #Selection of optimal treshold, adapted from library(BIOMOD)
 CutOff.optim <-
 function (Pred_proba, Obs, ncuts=99)
 {
     if (sum(Obs) == 0)
         stop("\n The observed data only contains 0")

     Quant <- quantile(Pred_proba)
     stat_tab <- as.data.frame(matrix(NA, nrow = ncuts, ncol = 2))

     if (length(unique(Quant)) == 1) {
         CutOff <- Quant[1]

     } else {

      for (j in 0:ncuts) {
             Seuil <- Quant[1] + (j * ((Quant[5] - Quant[1])/(ncuts+1)))
             Pred= ifelse(Pred_proba >= Seuil,1,0)

             Stat <- TSS(Pred, Obs)


     	    if (!is.na(Stat))  {

         	        stat_tab[j+1 , 1] <- Seuil
             	    stat_tab[j+1 , 2] <- Stat
             }

        }

         maxStat <- max(stat_tab[, 2], na.rm = T)
         seuil <- stat_tab[stat_tab[, 2] == maxStat, 1]
         CutOff = seuil[1]
      }
 	return(list(CutOff = CutOff, tab= stat_tab))


 }
#
#
# #Selection of optimal weights for abundance classes
# MultiWeights.optim <-
# function (Obs, Pred_proba){
#
#     nclasses <- ncol(Pred_proba)
#
#     fct=function(x){1/(1+x*x)}     ## to go between 0 and 1
#
#     fct_to_optim <- function(ws){
#
#       ws <- unlist(lapply(ws, fct))
#
#       prediction <- factor(rep(NA, length(Obs)))
#       levels(prediction) <- levels(Obs)
#
#       Pred_proba2 <- t(apply(Pred_proba, 1, function(x){x * ws}))
#       prediction[1:length(prediction)] <- apply(Pred_proba2, 1, function(x){colnames(Pred_proba2)[which.max(x)]} )
#
#
#       Stat <- HK(prediction,Obs)
#
#
#       if(is.na(Stat)) Stat = 0
#
#       return(1-Stat)
#
#     } # end to optim
#
#     Weights <- rep(NA, nclasses)
#     names(Weights) <- colnames(Pred_proba)
#     optima <- optim(rep(0,nclasses), fct_to_optim, method = "SANN", hessian = FALSE)
#     print(paste("convergence",optima$convergence))
#     Weights <-  unlist(lapply(optima$par, fct))
#
#
#     return(Weights)
#
# }
#
#
#
# #
# #------------------------------------------------------------------------------------------
# # Multi model inference
# #------------------------------------------------------------------------------------------
#
# commiteeVote <- function(x){
#   if( (sum(x==0)/length(x))>0.5 ) {
#     res = 0
#   }else{
#     x2=x[x!=0]
#     res = median(x2)
#   }
# return(res)
# }
#
#

