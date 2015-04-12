subsample.stratif3D <- function(xyz, subsetProp, adj = 1.1)
{
require(lhs)
sampl.raw = randomLHS(as.integer(subsetProp*nrow(xyz)*adj), 3)
# transform to the data range
sampl = data.frame(sampl.raw)
sampl[,1] = sampl[,1]*(range(xyz[,1])[2]-range(xyz[,1])[1]) + range(xyz[,1])[1]
sampl[,2] = sampl[,2]*(range(xyz[,2])[2]-range(xyz[,2])[1]) + range(xyz[,2])[1]
sampl[,3] = sampl[,3]*(range(xyz[,3])[2]-range(xyz[,3])[1]) + range(xyz[,3])[1]


# remove from outside the convex hull (2 plans)
require(grDevices)
hull1 = chull(xyz[,1], xyz[,2])
hull2 = chull(xyz[,3], xyz[,2])

require(splancs)
inPoly1 <- inout(as.points(sampl[,1], sampl[,2]),as.points(xyz[,1][hull1], xyz[,2][hull1]))
inPoly2 <- inout(as.points(sampl[,3], sampl[,2]),as.points(xyz[,3][hull2], xyz[,2][hull2]))

sampl2 = sampl[inPoly1&inPoly2,]

require(sp)

#require(parallel)
select = rep(NA, nrow(sampl2))
for(i in 1:nrow(sampl2)){
#print(i)
distance = spDists(as.matrix(xyz), sampl2[i,])
mini = ifelse(i>1, min(distance[-select[1:(i-1)]]), min(distance) ) 
inds = which(distance==mini)
already_selected = inds%in%select
if(sum(already_selected, na.rm=TRUE)>0) inds = inds[-which(inds%in%select)]
select[i] = ifelse(length(inds==1), inds, sample(inds,1))

}


#require(parallel)
#select = mclapply(1:nrow(sampl2), function(i){
#
#inBuffer = inout(as.points(xy[,1], xy[,2]), as.points(c(rep(sampl2[i,1]-buffer, 2),rep(sampl2[i,1]+buffer, 2)), rep(c(sampl2[i,2]-buffer, sampl2[i,2]+buffer), 2)))
#
#select = ifelse(sum(inBuffer)==0, NA, sample(which(inBuffer),1))
#return(select)
#})

select = unique(unlist(select))
cat("sampl asked ", nrow(xyz)*subsetProp, "\n")
cat("sampl taken ", length(select), "\n")

return(select)
}

subsample.stratif <- function(xy, subsetProp, adj = 1.25)
{
require(lhs)
sampl.raw = randomLHS(as.integer(subsetProp*nrow(xy)*adj), 2)
# transform to the data range
sampl = data.frame(sampl.raw)
sampl[,1] = sampl[,1]*(range(xy[,1])[2]-range(xy[,1])[1]) + range(xy[,1])[1]
sampl[,2] = sampl[,2]*(range(xy[,2])[2]-range(xy[,2])[1]) + range(xy[,2])[1]


# remove from outside the convex hull
require(grDevices)
hull = chull(xy[,1], xy[,2])
require(splancs)
inPoly <- inout(as.points(sampl[,1], sampl[,2]),as.points(xy[,1][hull], xy[,2][hull]))
sampl2 = sampl[inPoly,]

require(sp)

require(parallel)
select = mclapply(1:nrow(sampl2), function(i){

distance = spDists(as.matrix(xy), sampl2[i,])
inds = which(distance==min(distance))
sel = ifelse(length(inds==1), inds, sample(inds,1))

return(sel)
})


#require(parallel)
#select = mclapply(1:nrow(sampl2), function(i){
#
#inBuffer = inout(as.points(xy[,1], xy[,2]), as.points(c(rep(sampl2[i,1]-buffer, 2),rep(sampl2[i,1]+buffer, 2)), rep(c(sampl2[i,2]-buffer, sampl2[i,2]+buffer), 2)))
#
#select = ifelse(sum(inBuffer)==0, NA, sample(which(inBuffer),1))
#return(select)
#})

select = unique(unlist(select))
cat("sampl asked ", nrow(xy)*subsetProp, "\n")
cat("sampl taken ", length(select), "\n")

return(select)
}

subsample.temp <- function(temp, subsetProp)
{
interval = seq(min(temp), max(temp), l= as.integer(subsetProp*1.2*length(temp)/10))

require(parallel)
select = mclapply(2:length(interval), function(i){

inZone = intersect(which(temp>=interval[i-1]), which(temp<interval[i]))

if(length(inZone)==0) select = NA else select = sample(inZone,10, replace=TRUE)
return(select)
})

select = unique(na.omit(unlist(select)))
cat("sampl asked ", length(temp)*subsetProp, "\n")
cat("sampl taken ", length(select), "\n")

return(select)
}

subsample.temp.fix <- function(temp, nsample)
{
interval = seq(min(temp), max(temp), l= (nsample+1)/10)

require(parallel)
select = mclapply(2:length(interval), function(i){

inZone = intersect(which(temp>=interval[i-1]), which(temp<interval[i]))

if(length(inZone)==0) select = NA else select = sample(inZone,10, replace=TRUE)
return(select)
})
bag = 1:length(temp)
select = unlist(select)
select2 = sample(bag[-na.omit(unlist(select))], 10*sum(is.na(select)))
select = c(na.omit(select), select2)
select[duplicated(select)] = sample(bag[-select], sum(duplicated(select)))

select = unique(select)
cat("sampl asked ", nsample, "\n")
cat("sampl taken ", length(select), "\n")

return(select)
}

