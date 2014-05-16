library(nnet)
# Load data
data = as.data.frame(read.table("data/data_categorical.txt"))
data$E = data$av_annual_mean_tp
data$E2 = data$E^2

# Run the model
SDM = multinom(t1 ~ E + E2, data)

# Get predictions
pred_real = predict(SDM,new=data.frame(E=data$E,E2=data$E2),"probs")
pred_gradient = predict(SDM,new=data.frame(E=seq(-2,6,0.01),E2=seq(-2,6,0.01)^2),"probs")
write.table(pred_real,"data_pred_states.txt")
write.table(cbind(seq(-2,6,0.01),pred_gradient),"data_pred_gradient.txt")

# Plot the results
x11(height = 6, width = 8)
par(mar=c(5,5,2,1))
plot(seq(-2,6,0.01),pred_gradient[,1],type = "l",ylim=c(0,1),cex.axis = 1.25, cex.lab = 1.25, xlab = "Température moyenne annuelle", ylab = "Proportion",lwd = 2)
lines(seq(-2,6,0.01),pred_gradient[,2],col = "darkred",lwd = 2)
lines(seq(-2,6,0.01),pred_gradient[,3],col = "darkblue",lwd = 2)
lines(seq(-2,6,0.01),pred_gradient[,4],col = "darkgreen",lwd = 2)
abline(v = 3.2,lty = 3)
legend("topright",bty = "n", col = c("black","darkred","darkblue","darkgreen"),legend = c("C","D","M","T"),lty=1)

dev.copy2pdf(file = "figures/SDM.pdf")



