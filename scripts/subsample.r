subsample.space <- function(xy, subsetProp, buffer=.5)
{
require(lhs)
sampl.raw = randomLHS(as.integer(subsetProp*nrow(xy)*1.62), 2)
# transform to the data range
sampl = data.frame(sampl.raw)
sampl[,1] = sampl[,1]*(range(xy[,1])[2]-range(xy[,1])[1]) + range(xy[,1])[1]
sampl[,2] = sampl[,2]*(range(xy[,2])[2]-range(xy[,2])[1]) + range(xy[,2])[1]


# remove from outside the convex hull
library(grDevices)
hull = chull(xy[,1], xy[,2])
library(splancs)
inPoly <- inout(as.points(sampl[,1], sampl[,2]),as.points(xy[,1][hull], xy[,2][hull]))
sampl2 = sampl[inPoly,]


require(parallel)
select = mclapply(1:nrow(sampl2), function(i){

inBuffer = inout(as.points(xy[,1], xy[,2]), as.points(c(rep(sampl2[i,1]-buffer, 2),rep(sampl2[i,1]+buffer, 2)), rep(c(sampl2[i,2]-buffer, sampl2[i,2]+buffer), 2)))

select = ifelse(sum(inBuffer)==0, NA, sample(which(inBuffer),1))
return(select)
})

select = unique(na.omit(unlist(select)))
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
