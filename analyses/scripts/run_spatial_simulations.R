rm(list=ls())

# Parameters and functions
setwd("/home/DominiqueGravel/Bureau/analyses/data")
par = read.table("par.txt")

setwd("/home/DominiqueGravel/Bureau/analyses/scripts")
source("get_transitions.R")
source("spatial_model.R")
source("plot_map")

# Initial conditions
setwd("/home/DominiqueGravel/Bureau/analyses/data")

BIC = read.table("state_BIC.txt")
XY = BIC[,1:2]
pres = as.character(BIC[,3])
ENV = runif(length(pres),2.8,3.6)
N = nrow(XY)

# Run a short simulation, step by step
plot_map(XY,pres,"2010")
setwd("/home/DominiqueGravel/Bureau/analyses/figures")
dev.copy2pdf(file = "BIC2010.pdf")

for(i in 1:10) {
	shortrun = spatial_model(pres=pres,ENV=ENV,XY=XY,par=par,nstep = 1)
	map = shortrun$map
	plot_map(XY,map,title=2010+i*10)
	dev.copy2pdf(file = paste("BIC",2010+i*10,".pdf",sep=""))
	cat(2010+i*10,'\n')
}

# Run the model over a very large landscape
Xlim = 250
Ylim = 25
N = Xlim*Ylim
X = c(1:Xlim)
Y = c(1:Ylim)
XY = expand.grid(X,Y)
ENV = 6*(XY[,1]-1)/249 
initial_probs = matrix(0.25,nr = 1000, nc = 4)
pres = apply(initial_probs,1,draw)
large_run = spatial_model(pres=pres,ENV=ENV,XY=XY,par=par,nstep = 100)

# Plot the results
X = sort(unique(XY[,1]))
Y = sort(unique(XY[,2]))
pres_int = as.factor(large_run$map)
n = 1
Z = matrix(nr = length(X), nc = length(Y))
for(x in 1:length(X)) 
	for(y in 1:length(Y)) 	
		Z[x,y] = pres_int[XY[,1] == X[x] & XY[,2]==Y[y]]	
x11(width = 25, height = 5)
layout(matrix(c(1,2),nr=2,nc=1,byrow=TRUE),heights = c(1,5))
par(mar=c(0,0,2,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
title(title,cex=2)
legend("center",legend = c("C","M","D","T"),fill = grey(c(0:3)/3),bty = "n",horiz = TRUE,cex = 1.5)
par(mar=c(5,5,0,5))
image(X,Y,Z,cex.axis = 1.5, cex.lab = 1.25, col = grey(c(0:3)/3))
dev.copy2pdf("largeplot.pdf")





