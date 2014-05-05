Xlim = 150
Ylim = 20
N = Xlim*Ylim
X = c(1:Xlim)
Y = c(1:Ylim)
XY = expand.grid(X,Y)
distMat = as.matrix(dist(XY,method = "euclidean", upper = T, diag = T))
ConMat = matrix(0, nr=N,nc=N)
ConMat[distMat<1.5] = 1
diag(ConMat) = 0

S = 2
nsteps = 10000
e1 = 0.1
e2 = 0.1
amax = 1.1
cross = 120
b0 = 1
a = (amax-b0)/cross
a0 = amax - a*XY[,1]
a1 = -0.05

# N: number of cells in the lattice
# presence: Vector of species presence/absence (0: absence, 1: presence, NxS dimension)

# Initiatization of the metaco, distribute trees at random 
pres1 = matrix(0,nrow = N,nc = S)
rand = runif(N,0,1)
pres1[rand<0.5,1] = 1
pres1[rand>0.5,2] = 1
rel = matrix(nr = nsteps,nc = 2)
n=1

quartz(width = 5, height = 5)
# Loop over all time steps
for(time in 1:nsteps) {
	pres0 = pres1

	# Mortality events
	randMort = runif(N,0,1)
	pres1[randMort<e1,1] = 0
	pres1[randMort<e2,2] = 0
	
	# Replacement events
	# Calculate abundance in the neighbourhood
	ConPop = (ConMat%*%pres0)

	# Calculate the recruitment probability
	fitness = matrix(0,nr=N,nc=S)
	fitness[pres0[,1]==1,1] = a0[pres0[,1]==1]
	fitness[pres0[,1]==0,1] = a0[pres0[,1]==0]+a1	 	
	fitness[,2] = b0
	
	recruitProb = numeric(N)
	recruitProb = ConPop[,1]*fitness[,1]/apply(ConPop*fitness,1,sum)
	
	# Pick a species at random and do replacement
	randRecruit = runif(N,0,1)
	pres1[apply(pres1,1,sum)==0 & recruitProb > randRecruit,1] = 1
	pres1[apply(pres1,1,sum)==0 & recruitProb < randRecruit,2] = 1

	# Record each species relative abundance
	rel[time,] = apply(pres1,2,sum)/N

	if(n == 10) {
		z = matrix(nr=Xlim,nc=Ylim)
		for(x in 1:Xlim) for(y in 1:Ylim) z[x,y] = pres1[XY[,1]==x & XY[,2]==y,1]

		t = layout(matrix(c(1,2),nr=2,nc=1),height=c(1,3))
		layout.show(t)
		par(mar=c(0.5,5,1,1))

		sp1 = numeric(Xlim)
		sp2 = numeric(Xlim)
	
		A0 = amax - a*c(1:Xlim)
		sp1[A0>b0] = 1
		sp2[b0>(A0+a1)] = 0.99
	
		plot(c(1:Xlim),sp1,type="l", col = "blue",lwd = 2,ylab="Presence",xaxt = "n",xlab = "")
		lines(c(1:Xlim),sp2,col = "red", lwd = 2)	
		REL = tapply(pres1[,1],INDEX = XY[,c(1,2)],sum)/10	
		for(i in 1:6) points(c(1:Xlim),REL[,i])
		title(main = paste("Time = ",time,sep = ""))

		par(mar=c(5,5,0.5,1))
		image(x=c(1:Xlim),y=c(1:Ylim),z=z,col = c("blue","red"),xlab = "X", ylab = "Y")
		n = 1
		cat(time,'\n')
	}
	else n = n+1
	cat(time,'\n')
}
# nbin = 5
# Ybin = floor(XY[,2]/nbin)
# ind = cbind(XY[,1],Ybin)

# REL = matrix(nr = 150*floor(Ylim/nbin),nc=2)
# n = 1
# for(i in 1:Xlim) 
	# for(j in 0:(floor(Ylim/nbin)-1)) {
		# REL[n,1] =i
		# REL[n,2] = sum(subset(pres1[,1],ind[,1]==i&ind[,2]==j))
		# n = n+1
		# }	
		
# hist(REL[REL[,1]==100,2]/5)


