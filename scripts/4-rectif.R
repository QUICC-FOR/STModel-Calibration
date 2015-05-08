rm(list = ls())

sdm = "rf"
propData = 0.33

#-- load dat for datValid
pred = read.table("../data/projection_rf_complete.txt", h=T)

#--
for( i in 1:9)
{
load(paste("initForFit_", sdm, "_", propData, i, ".RData", sep=""))
#-- dat valid --
neig = pred[pred$plot %in% select,][-toremove,]
load("scale_info.Robj")
ENV1 = (dataProj_subset10$annual_mean_temp - vars.means["annual_mean_temp"])/ vars.sd["annual_mean_temp"]
ENV2 = (dataProj_subset10$tot_annual_pp - vars.means["tot_annual_pp"])/ vars.sd["tot_annual_pp"]
dat = data.frame(ENV1 = ENV1 , ENV2 = ENV2, st0 = dataProj_subset10[,"state1"], st1 = dataProj_subset10[,"state2"], itime = dataProj_subset10[,"itime"], plot_id= dataProj_subset10[,"plot"], EB = neig[,"B"], ET = neig[,"T"],EM = neig[,"M"])
datValid = dat[-select2,]
nrow(datValid)+nrow(datSel) == nrow(dataProj_subset10)
#----------------
save(datSel, dataProj_subset10, select2,params, toremove, select, par_lo, par_hi, fit, dat, datValid,file= paste("initForFit_", sdm, "_", propData, i, ".RData", sep=""))
}
