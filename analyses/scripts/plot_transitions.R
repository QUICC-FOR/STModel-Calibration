
attach(read.table("data/par.txt"))
pred = read.table("data/data_pred_gradient.txt")
ENV = pred[,1]
EC = pred[,2]
ED = pred[,3]
EM = pred[,4] 

# Compute the logit
logit_alphac 	= ac0 + ac1*ENV + ac2*ENV^2
logit_alphad 	= ad0 + ad1*ENV + ad2*ENV^2
logit_betac 	= (bc0 + bc1*ENV + bc2*ENV^2)
logit_betad 	= (bd0 + bd1*ENV + bd2*ENV^2)
logit_thetac	= tc0 + tc1*ENV + tc2*ENV^2
logit_thetad	= td0 + td1*ENV + td2*ENV^2
logit_eps 	= e0  + e1*ENV  + e2*ENV^2

# Back transform into probabilities
alphac = exp(logit_alphac)/(1+exp(logit_alphac))
alphad = exp(logit_alphad)/(1+exp(logit_alphad))
betac = exp(logit_betac)/(1+exp(logit_betac))*(EC+EM)
betad = exp(logit_betad)/(1+exp(logit_betad))*(ED+EM)
thetac = exp(logit_thetac)/(1+exp(logit_thetac))
thetad = exp(logit_thetad)/(1+exp(logit_thetad))
eps = exp(logit_eps)/(1+exp(logit_eps))
phic = alphac*(EM + EC)*(1-alphad*(ED+EM))
phid = alphad*(EM + ED)*(1-alphac*(EC+EM))
phim = phic*phid

# Plot the results
x11(height = 6, width = 14)
par(mar=c(5,5,2,1),mfcol = c(1,2))

plot(ENV,phid,type = "l",ylim=c(0,1),cex.axis = 1.25, cex.lab = 1.25, xlab = "Température moyenne annuelle", ylab = "Probabilité",lwd = 2,col = "darkgreen")
lines(ENV,betad,col = "darkred",lwd = 2)
lines(ENV,thetad,col = "darkblue",lwd = 2)
legend("top",bty = "n", col = c("darkgreen","darkred","darkblue"),legend = c("T-->D","C-->M","M-->D"),lty=1,horiz=TRUE,lwd = 3)

plot(ENV,phic,type = "l",ylim=c(0,1),cex.axis = 1.25, cex.lab = 1.25, xlab = "Température moyenne annuelle", ylab = "Probabilité",lwd = 2,col = "darkgreen")
lines(ENV,betac,col = "darkred",lwd = 2)
lines(ENV,thetac,col = "darkblue",lwd = 2)
legend("top",bty = "n", col = c("darkgreen","darkred","darkblue"),legend = c("T-->C","D-->M","M-->C"),lty=1,horiz=TRUE,lwd = 3)

dev.copy2pdf(file = "figures/Transitions.pdf")



