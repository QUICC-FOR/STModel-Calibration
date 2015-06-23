rm(list=ls())
##-----------
## load dat
##-----------

dat = read.table("../data/projection_rf_complete_woLatLon.txt", h=T)
load("scale_info.Robj")
head(dat)
# subset 10 degree
dat = dat[which(dat$annual_mean_temp<=((10-vars.means["annual_mean_temp"])/vars.sd["annual_mean_temp"])),]
#remove transitions directes B->T ou T->B
trBT = c(which(dat$state1 == "T" & dat$state2 == "B"), which(dat$state1 == "B" & dat$state2 == "T"))
# clean transition time
dat$itime = dat$year2 - dat$year1
rmItime= c(which(dat$itime<5), which(dat$itime>15))
# clean dataset
toremove = unique(c(trBT, rmItime))
dat = dat[-toremove,]
dim(dat)
summary(dat)
## attention NA (from soil missing data)

#rename variables
dat$ENV1 = dat[,"annual_mean_temp"]
dat$ENV2 = dat[,"tot_annual_pp"]
dat$st0 = dat[,"state1"]
dat$st1 = dat[,"state2"]
dat$plot_id = dat[,"plot"]

dat = dat[,c("plot_id", "ENV1","ENV2","st0","st1","itime","B","M","R","T","year1","year2","annual_mean_temp",   "tot_annual_pp" ,     "mean_diurnal_range", "ph_2cm"        ,     "slp"            ,    "lat"       ,        "lon")]

head(dat)
dim(dat)

dat = na.omit(dat)
dim(dat)

save(dat, file = "datAll_woLatLon.RData")

#-----
Temp.lim = c((-4-vars.means["annual_mean_temp"])/vars.sd["annual_mean_temp"], (10-vars.means["annual_mean_temp"])/vars.sd["annual_mean_temp"])
Temp.ax = function(x)
{
temp = dat[,"annual_mean_temp"]*vars.sd["annual_mean_temp"]+vars.means["annual_mean_temp"]
axis(1, at = seq((round(min(temp))-vars.means["annual_mean_temp"])/vars.sd["annual_mean_temp"], (round(max(temp))-vars.means["annual_mean_temp"])/vars.sd["annual_mean_temp"], l=16), labels = seq(round(min(temp)), round(max(temp)), l = 16))
}
#-
pp.lim = c((500-vars.means["tot_annual_pp"])/vars.sd["tot_annual_pp"], (2000-vars.means["tot_annual_pp"])/vars.sd["tot_annual_pp"])
pp.ax = function(x)
{
pp = dat[,"tot_annual_pp"]*vars.sd["tot_annual_pp"]+vars.means["tot_annual_pp"]
axis(2, at = seq((round(min(pp))-vars.means["tot_annual_pp"])/vars.sd["tot_annual_pp"], (round(max(pp))-vars.means["tot_annual_pp"])/vars.sd["tot_annual_pp"], l=16), labels = seq(round(min(pp)), round(max(pp)), l = 16))
}

#
jpeg("../figures/calibration_data_in_climSpace.jpeg", height=3000, width=5000, res=600)

plot(dat$annual_mean_temp , dat$tot_annual_pp, ylab = "precipitations", xlab = "temperature" ,cex = .2, pch = 20, xaxt="n", yaxt="n", xlim = Temp.lim, ylim = pp.lim, col = 1)
Temp.ax()
pp.ax()

rect(xleft=(-4-vars.means["annual_mean_temp"])/vars.sd["annual_mean_temp"], 
xright=(10-vars.means["annual_mean_temp"])/vars.sd["annual_mean_temp"], 
ybottom=(700-vars.means["tot_annual_pp"])/vars.sd["tot_annual_pp"], 
ytop=(1300-vars.means["tot_annual_pp"])/vars.sd["tot_annual_pp"],  
col=NA, border = 2, lwd =1.2)


dev.off()