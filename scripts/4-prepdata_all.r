rm(list=ls())
##-----------
## load dat
##-----------
load("../data/transitions_r1.RData")
dataProj = transitionData
head(dataProj)
dim(dataProj)
# subset 10 degree
select = unique(dataProj$plot[which(dataProj$annual_mean_temp<=10)])
dataProj_subset10 = dataProj[dataProj$plot %in% select,]
# rescale
load("scale_info.Robj")
dat_scale = dataProj_subset10[c("annual_mean_temp", "tot_annual_pp")]
dat_scale = t(apply(dat_scale, 1, function(x) {(x-vars.means[c("annual_mean_temp", "tot_annual_pp")])/vars.sd[c("annual_mean_temp", "tot_annual_pp")]}))
dat_scale = data.frame(dat_scale)
#head(dat_scale)
dim(dat_scale)
# remove transitions directes B->T ou T->B
trBT = c(which(dataProj_subset10$state1 == "T" & dataProj_subset10$state2 == "B"), which(dataProj_subset10$state1 == "B" & dataProj_subset10$state2 == "T"))
# clean transition time
dataProj_subset10$itime = dataProj_subset10$year2 - dataProj_subset10$year1
rmItime= c(which(dataProj_subset10$itime<5), which(dataProj_subset10$itime>15))
# clean dataset
toremove = unique(c(trBT, rmItime))
dataProj_subset10 = dataProj_subset10[-toremove,]
# create dataset
dat = dat_scale[-toremove,]
#rename variables
colnames(dat) = c("ENV1", "ENV2")
dat$st0 = dataProj_subset10[,"state1"]
dat$st1 = dataProj_subset10[,"state2"]
dat$itime = dataProj_subset10[,"itime"]
dat$plot_id = dataProj_subset10[,"plot"]
# neighborhood
pred = read.table("../data/projection_rf_complete.txt", h=T)
pred = pred[pred$plot %in% select,]
pred = pred[-toremove,]
dat$EB = pred$B
dat$ET = pred$T
dat$EM = pred$M 
head(dat)

save(dat, file = "datAll.RData")

#
jpeg("../figures/calibration_data_in_climSpace.jpeg", height=3000, width=5000, res=600)

plot(dataProj_subset10$annual_mean_temp, dataProj_subset10$tot_annual_pp, ylab = "precipitations", xlab = "temperature" , cex =.2, pch = 19)
rect(xleft=-4, xright=10, ybottom=700, ytop=1300, col=NA, border = 2, lwd =1.2)


dev.off()