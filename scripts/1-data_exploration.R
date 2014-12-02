
rm(list = ls())
####################################################################
######                Clean data (temporary version)         #######
####################################################################
#
#dat = read.table("../data/data_BA.txt",header = TRUE, sep = ";")
##head(dat)
#dim(dat)
#
####################################################################
######                  Classify plots                       #######
####################################################################
#
## List of interests species
#I_sp  <- c("bop","peb","peg","pet","prp","sal","soa")
#B_sp  <- c("epn","epb","epr", "mel","pig","sab")
#T_sp  <- c("ers","fra","frn","heg","osv","til","cet")
#
##retirer frn
#
## Subset species BA and cover type observed
#BA_dat <- c(which(colnames(dat)=="ame"): which(colnames(dat)=="til"))
#
#### Get BA by hectares (the plot is 400m2)
#dat[,BA_dat] <- dat[,BA_dat]*10000/400
#
#dat$I_tot  <- rowSums(dat[,which(colnames(dat) %in% I_sp)],na.rm=T)
#dat$B_tot  <- rowSums(dat[,which(colnames(dat) %in% B_sp)],na.rm=T)
#dat$T_tot  <- rowSums(dat[,which(colnames(dat) %in% T_sp)],na.rm=T)
#dat$BA_tot <- rowSums(dat[,BA_dat],na.rm=T)
#
## Class into state types
#
#class_fn = function(x) {
#  Itot = x[1]; Btot =x[2]; Ttot=x[3]; BAtot=x[4]
#  classPlot = "Unclass"
#  if(BAtot < 10) classPlot = "R"
#  else if(Btot > 0 & Ttot == 0) classPlot = "B"
#  else if(Ttot > 0 & Btot == 0) classPlot = "T"
#  else if (Btot > 0 & Ttot > 0) classPlot = "M"
#return(classPlot)
#}
#
#dat$class_final  <- as.vector(apply(dat[,c("I_tot", "B_tot", "T_tot", "BA_tot")],1,class_fn))
#
#table(dat$class_final)
##dat[dat$class_final=="Unclass",]
#
#dat = dat[-which(dat$class_final=="Unclass"),]
#dim(dat)
#
#
#
####################################################################
######    Reshape (pairwise) and export data                 #######
####################################################################
#
#liste =split(dat[,c("id_plot", "yr_measured", "lat","lon", "annual_pp", "annual_mean_temp", "class_final")], dat$id_plot)
#nreleves = unlist(lapply(liste, nrow))
#
## remove only one releve in time
#liste = liste[-which(nreleves==1)]
#
### Function to match remeasurements on the same line
#prwise  <- function(x){
#    #climatic columns
#    clim = c(which(colnames(x) == "annual_pp") ,which(colnames(x) == "annual_mean_temp"))
#
#    x = x[order(x$yr_measured),]
#
#    for (i in 1:(nrow(x)-1))
#    {
#    yr0 = x$yr_measured[i]
#    st0 = x$class_final[i]
#    clim0 = x[i, clim ]
#    yr1 = x$yr_measured[i+1]
#    st1 = x$class_final[i+1]
#    clim1 = x[i+1, clim ]
#    if(i ==1) {
#    res = data.frame(idplot = x$id_plot[i], st0, st1, yr0, yr1, lat = x$lat[i], lon = x$lon[i],  t(apply(rbind(clim0, clim1), 2, mean)))
##    res = data.frame(st0, st1, yr0, yr1, lat = x$lat, lon = x$lon,  clim0, clim1)
#    }else{
#    res = rbind(res, data.frame(idplot = x$id_plot[i], st0, st1, yr0, yr1,lat = x$lat[i], lon = x$lon[i], t(apply(rbind(clim0, clim1), 2, mean))))
##    res = data.frame(st0, st1, yr0, yr1, lat = x$lat, lon = x$lon,  clim0, clim1)
#    }
#    }
#
#  return(res)
#}
#
#library(parallel)
#pair_dat = mclapply(liste, prwise)
#pair_dat = (do.call(rbind.data.frame, pair_dat))
#dim(pair_dat)
#head(pair_dat)
#
##check order
#sum(pair_dat$yr0>pair_dat$yr1)
####################################################################
######    filtres                            #######
####################################################################
#
### filtre 0
###-----------
### Rm all plots with no climatic data associated
#
#dat0 = dat[which(!is.na(dat$annual_pp)),]
#liste =split(dat0, dat0$id_plot)
#nreleves = unlist(lapply(liste, nrow))
## remove only one releve in time
#liste = liste[-which(nreleves==1)]
#
#pair_dat0 = mclapply(liste, prwise)
#pair_dat0 = (do.call(rbind.data.frame, pair_dat0))
#
#dim(pair_dat0)
#head(pair_dat0)
#
#
### filtre 1
###-----------
#
### Rm all harvested plots
#dat1  <- dat0[-which(dat0$disturbance%in%c("BRP", "CAM", "CB", "CD", "CDL", "CE", "CJ", "CP", "DLD", "EPC")),]
## BRP = brulis partiel ## ne pas enlever
## CAM = coupe d'amélioration
## CB = coupe par bandes
## CD = coupes en damier
## CDL = coupe à diamètre limité
## CE = coupe partielle et épidémie légère
## CHP = chablis partiel
## CJ = coupe de jardinage
## CP = coupe partielle
## DLD = Coupe à diamètre limite avec dégagement des arbres d'avenir
## DP  = Dépérissement partiel du feuillu
## EL = Épidémie légère
## EPC = Éclaircie précommerciale
## VEP = Verglas partiel
#
#liste =split(dat1, dat1$id_plot)
#nreleves = unlist(lapply(liste, nrow))
## remove only one releve in time
#liste = liste[-which(nreleves==1)]
#
#pair_dat1 = mclapply(liste, prwise)
#pair_dat1 = (do.call(rbind.data.frame, pair_dat1))
#
#head(pair_dat1)
#dim(pair_dat1)
#
#
#
### filtre 2
###-----------
#
### Conserve all plots with drainage 20 to 40
## !!! there are lots of NA, we keep them
#dat2 <- dat1[-which(dat1$drainage<20 | dat1$drainage>=40),]
#
#liste =split(dat2, dat2$id_plot)
#nreleves = unlist(lapply(liste, nrow))
## remove only one releve in time
#liste = liste[-which(nreleves==1)]
#
#pair_dat2 = lapply(liste, prwise)
#pair_dat2 = (do.call(rbind.data.frame, pair_dat2))
#
#dim(pair_dat2)
#head(pair_dat2)
#
##---------------
#
#write.table(dat0,file="../data/data_allyears.txt")
#write.table(pair_dat0,file="../data/data_pairs.txt")
#write.table(dat2,file="../data/data_allyears_filter.txt")
#write.table(pair_dat2,file="../data/data_pairs_filter.txt")
#

#save(pair_dat, pair_dat0, pair_dat1, pair_dat2, file = "pair_dat.RData")
#load("pair_dat.RData")

###################################################################
#####    Reshaping data from STM-Data repo                  #######
###################################################################

# Set wd to STModel-Data folder
setwd("~/Documents/GitHub/STModel-Data")

# Open data
pair.dat <- read.csv("./out_files/transitionsFourState.csv")

# average on climatic data
pair.dat$annual_mean_temp <- rowMeans(subset(pair.dat, select = c(annual_mean_temp1, annual_mean_temp2)), na.rm = TRUE)
pair.dat$annual_pp <- rowMeans(subset(pair.dat, select = c(annual_pp1, annual_pp2)), na.rm = TRUE)

# Rename and clean columns
names(pair.dat)[7:8] <- c("st0","st1")
pair.dat <- pair.dat[,-c(9:12)]

# Create transition column
pair.dat$transition <- paste(pair.dat$st0,pair.dat$st1,sep="")

# Datset without filters
pair_dat0 <- pair.dat

# Graph lim
rg_pp <- range(pair.dat$annual_pp)
rg_tp <- range(pair.dat$annual_mean_temp)

# change wd
setwd("~/Documents/GitHub/STModel-Calibration/scripts")

###################################################################
#####    Analyses    GRAPHS                                 #######
###################################################################

colo = c(rgb(1,165/255,0,0.5), rgb(154/255,205/255,50/255,0.5), rgb(34/255,139/255,34/255,0.5), rgb(148/255,0,211/255,0.5) )
names(colo) = c("R", "T", "B", "M")

figs = function(pair.dat, filename = "../figures/data_explo.pdf" )
{

transitions = paste(pair.dat$st0, pair.dat$st1, sep="")
print(table(transitions))

pdf(file = filename, width = 15 , height=8)

par(mfrow = c(2, 4))
hist(pair.dat$annual_mean_temp[pair.dat$st0=="R"], xlab = "average annual temperature", col = colo["R"], main = "R", freq = TRUE, xlim = c(rg_tp[1],rg_tp[2]))
hist(pair.dat$annual_mean_temp[pair.dat$st0=="T"], xlab = "average annual temperature", col = colo["T"], main = "T", freq = TRUE, xlim = c(rg_tp[1],rg_tp[2]))
hist(pair.dat$annual_mean_temp[pair.dat$st0=="B"], xlab = "average annual temperature", col = colo["B"], main = "B", freq = TRUE, xlim = c(rg_tp[1],rg_tp[2]))
hist(pair.dat$annual_mean_temp[pair.dat$st0=="M"], xlab = "average annual temperature", col = colo["M"], main = "M", freq = TRUE, xlim = c(rg_tp[1],rg_tp[2]))
hist(pair.dat$annual_pp[pair.dat$st0=="R"], xlab = "annual precipitations", col = colo["R"], main = "R", freq = TRUE, xlim = c(rg_pp[1],rg_pp[2]))
hist(pair.dat$annual_pp[pair.dat$st0=="T"], xlab = "annual precipitations", col = colo["T"], main = "T", freq = TRUE, xlim = c(rg_pp[1],rg_pp[2]))
hist(pair.dat$annual_pp[pair.dat$st0=="B"], xlab = "annual precipitations", col = colo["B"], main = "B", freq = TRUE, xlim = c(rg_pp[1],rg_pp[2]))
hist(pair.dat$annual_pp[pair.dat$st0=="M"], xlab = "annual precipitations", col = colo["M"], main = "M", freq = TRUE, xlim = c(rg_pp[1],rg_pp[2]))


par(mfrow = c(2, 3))
hist(pair.dat$annual_mean_temp[transitions=="RT"], xlab = "average annual temperature", col = colo["T"], main = "R->T", freq = TRUE, xlim = c(rg_tp[1],rg_tp[2]))
#lines(density(pair_dat$annual_mean_temp[transitions=="RT"]), col = colo["T"], lwd = 2)
hist(pair.dat$annual_mean_temp[transitions=="RB"], xlab = "average annual temperature", col = colo["B"],  main = "R->B", freq = TRUE, xlim = c(rg_tp[1],rg_tp[2]))
hist(pair.dat$annual_mean_temp[transitions=="RM"], xlab = "average annual temperature", col = colo["M"],  main = "R->M", freq = TRUE, xlim = c(rg_tp[1],rg_tp[2]))

hist(pair.dat$annual_pp[transitions=="RT"], xlab = "annual precipitation", col = colo["T"], main = "R->T", freq = TRUE, xlim = c(rg_pp[1],rg_pp[2]))
hist(pair.dat$annual_pp[transitions=="RB"], xlab = "annual precipitation", col = colo["B"],  main = "R->B", freq = TRUE, xlim = c(rg_pp[1],rg_pp[2]))
hist(pair.dat$annual_pp[transitions=="RM"], xlab = "annual precipitation", col = colo["M"],  main = "R->M", freq = TRUE, xlim =c(rg_pp[1],rg_pp[2]))

#dev.off()
#
#pdf(file = "../figures/data_explo_exclusion.pdf", width = 15 , height=8)

par(mfrow = c(2, 3))
hist(pair.dat$annual_mean_temp[transitions=="MT"], xlab = "average annual temperature", col = colo["T"], main = "M->T", freq = TRUE, xlim = c(rg_tp[1],rg_tp[2]))
hist(pair.dat$annual_mean_temp[transitions=="MB"], xlab = "average annual temperature", col = colo["B"],  main = "M->B", freq = TRUE, xlim = c(rg_tp[1],rg_tp[2]))
hist(pair.dat$annual_mean_temp[transitions=="MM"], xlab = "average annual temperature", col = colo["M"],  main = "M->M", freq = TRUE, xlim = c(rg_tp[1],rg_tp[2]))

hist(pair.dat$annual_pp[transitions=="MT"], xlab = "annual precipitation", col = colo["T"], main = "M->T", freq = TRUE, xlim = c(rg_pp[1],rg_pp[2]))
hist(pair.dat$annual_pp[transitions=="MB"], xlab = "annual precipitation", col = colo["B"],  main = "M->B", freq = TRUE, xlim = c(rg_pp[1],rg_pp[2]))
hist(pair.dat$annual_pp[transitions=="MM"], xlab = "annual precipitation", col = colo["M"],  main = "M->M", freq = TRUE, xlim =c(rg_pp[1],rg_pp[2]))

#dev.off()
#
#pdf(file = "../figures/data_explo_colonisation.pdf", width = 15 , height=8)

par(mfrow = c(2, 4))
hist(pair.dat$annual_mean_temp[transitions=="TT"], xlab = "average annual temperature", col = colo["T"],  main = "T->T", freq = TRUE, xlim = c(rg_tp[1],rg_tp[2]))
hist(pair.dat$annual_mean_temp[transitions=="TM"], xlab = "average annual temperature", col = colo["B"], main = "T->M", freq = TRUE, xlim = c(rg_tp[1],rg_tp[2]))
hist(pair.dat$annual_mean_temp[transitions=="BB"], xlab = "average annual temperature", col = colo["B"],  main = "B->B", freq = TRUE, xlim = c(rg_tp[1],rg_tp[2]))
hist(pair.dat$annual_mean_temp[transitions=="BM"], xlab = "average annual temperature", col = colo["T"],  main = "B->M", freq = TRUE, xlim = c(rg_tp[1],rg_tp[2]))

hist(pair.dat$annual_pp[transitions=="TT"], xlab = "annual precipitation", col = colo["T"], main = "T->T", freq = TRUE, xlim = c(rg_pp[1],rg_pp[2]))
hist(pair.dat$annual_pp[transitions=="TM"], xlab = "annual precipitation", col = colo["B"],  main = "T->M", freq = TRUE, xlim = c(rg_pp[1],rg_pp[2]))
hist(pair.dat$annual_pp[transitions=="BB"], xlab = "annual precipitation", col = colo["B"],  main = "B->B", freq = TRUE, xlim =c(rg_pp[1],rg_pp[2]))
hist(pair.dat$annual_pp[transitions=="BM"], xlab = "annual precipitation", col = colo["T"],  main = "B->M", freq = TRUE, xlim =c(rg_pp[1],rg_pp[2]))

#dev.off()
#
#
#pdf(file = "../figures/data_explo_disturbance.pdf", width = 15 , height=8)

par(mfrow = c(2, 4))
hist(pair.dat$annual_mean_temp[transitions=="RR"], xlab = "average annual temperature", col = colo["R"], main = "R->R", freq = TRUE, xlim = c(rg_tp[1],rg_tp[2]))
hist(pair.dat$annual_mean_temp[transitions=="TR"], xlab = "average annual temperature", col = colo["R"],  main = "T->R", freq = TRUE, xlim = c(rg_tp[1],rg_tp[2]))
hist(pair.dat$annual_mean_temp[transitions=="BR"], xlab = "average annual temperature", col = colo["R"],  main = "B->R", freq = TRUE, xlim = c(rg_tp[1],rg_tp[2]))
hist(pair.dat$annual_mean_temp[transitions=="MR"], xlab = "average annual temperature", col = colo["R"],  main = "M->R", freq = TRUE, xlim = c(rg_tp[1],rg_tp[2]))

hist(pair.dat$annual_pp[transitions=="RR"], xlab = "annual precipitation", col = colo["R"], main = "R->R", freq = TRUE, xlim = c(rg_pp[1],rg_pp[2]))
hist(pair.dat$annual_pp[transitions=="TR"], xlab = "annual precipitation", col = colo["R"],  main = "T->R", freq = TRUE, xlim = c(rg_pp[1],rg_pp[2]))
hist(pair.dat$annual_pp[transitions=="BR"], xlab = "annual precipitation", col = colo["R"],  main = "B->R", freq = TRUE, xlim =c(rg_pp[1],rg_pp[2]))
hist(pair.dat$annual_pp[transitions=="MR"], xlab = "annual precipitation", col = colo["R"],  main = "M->R", freq = TRUE, xlim =c(rg_pp[1],rg_pp[2]))

dev.off()


}


figs(pair_dat0)
figs(pair_dat1, filename = "../figures/data_explo_filterHarvest.pdf")
figs(pair_dat2, filename = "../figures/data_explo_filterHarvest+Drainage.pdf")

###################################################################
#####    Analyses     GLM                                      #######
###################################################################
library(MASS)
library(ROCR)
library(fmsb)
pal = colorRampPalette(c("lightblue", "yellow", "orange"), space = "rgb")

#####    glm climate                            #######

modelTransition_climate <- function(st0 , st1, pair.dat, name = NULL)
{
print(st0)
print("->")
print(st1)

datst0 = pair.dat[pair.dat$st0%in%st0, ]
datst0$transition = ifelse(datst0$st1 %in% st1, 1, 0)

mod = glm(transition ~ annual_mean_temp + I(scale(annual_mean_temp)^2) + I(scale(annual_mean_temp)^3) + annual_pp + I(scale(annual_pp)^2) + I(scale(annual_pp)^3) + annual_mean_temp:annual_pp, family = "binomial", data =datst0)
stepMod  = stepAIC(mod)
#print(summary(stepMod))


pred = predict(stepMod,new=datst0,"response")
# overall performance
R2 = NagelkerkeR2(stepMod)$R2
#discrimination
perf = performance(prediction(pred, datst0$transition), "auc")
AUC = perf@y.values[[1]]

## selected vars
## selected vars
coeff = summary(stepMod)$coefficients
vars = rownames(coeff)[-1]
effect = coeff[-1,1]
pval = coeff[-1,4]
print(pval)

if(is.null(name)) name = paste(st0, st1, sep = "->")

return(list(mod = stepMod, vars = vars,effect = effect, pval = pval , R2 = R2, AUC = AUC, ranges = apply(datst0[unlist(lapply(1:ncol(datst0), function(x)is.numeric(datst0[,x])))], 2, range), name = name))
}

#test
#mod = modelTransition_climate("R", c("T", "M"), pair.dat = pair_dat0)




###################################################################
#####    figures                            #######
###################################################################
pval.star <- function(pval)
{
star=""
if(pval<0.1) star = "."
if(pval<0.05) star = "*"
if(pval<0.01) star= "**"
if(pval<0.001) star = "***"

}


stats <- function(model.list)
{
vars = c("annual_mean_temp", "I(scale(annual_mean_temp)^2)", "I(scale(annual_mean_temp)^3)", "annual_pp" , "I(scale(annual_pp)^2)",  "I(scale(annual_pp)^3)" , "annual_mean_temp:annual_pp")
n = length(model.list)
v = length(vars)

plot(x = c(0,v+8), y  = c(0,n+2), type = "n", main = "Summary", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
text(1, n+1, label = "model")
text(3:(2+v), rep(n+1, v), label = c("T", "T2", "T3","P", "P2", "P3", "TP"))
text(c(v+4, v+6), rep(n+1, 2), label = c("R2", "AUC"))


for (i in 1:n)
{
mod = model.list[[i]]
mod.name = names(model.list)[i]
text(1, n+1-i, label = mod.name)

if(length(mod$vars)!=0)
{
   if(length(mod$vars)==1)
   {
   term = which(vars %in% mod$vars)
    text(2+term, (n+1-i), label = paste("(",ifelse(mod$effect>0, "+", "-"),")", pval.star(mod$pval), sep=""), cex = .6)
   } else {
   for (term in 1:v)
     {

   if(vars[term] %in% mod$vars) text(2+term, (n+1-i), label = paste("(",ifelse(mod$effect[vars[term]]>0, "+", "-"),")", pval.star(mod$pval[vars[term]]), sep=""), cex = .6)

     }
   }
}

text(c(v+4, v+6), rep(n+1-i, 2), label = c(round(mod$R2,2), round(mod$AUC,2)))
}

}


fig_glm <- function(mod)
{

temp = seq(as.numeric(mod$ranges[1, "annual_mean_temp"]), as.numeric(mod$ranges[2, "annual_mean_temp"]), length.out = 50)
pp = seq(as.numeric(mod$ranges[1, "annual_pp"]), as.numeric(mod$ranges[2, "annual_pp"]), length.out = 50)
prob = predict(mod$mod, newdata = data.frame(expand.grid(annual_mean_temp = temp, annual_pp = pp) ), type = "response")

image(x=temp, y=pp, z = matrix(prob, ncol = length(pp), nrow = length(temp)),xlab = "Temperature", ylab = "Precipitations", col = pal(12), main = mod$name)
contour(x=temp, y=pp, z = matrix(prob, ncol = length(pp), nrow = length(temp)), add=TRUE)
}

#test
#fig_glm(mod)

###################################################################
#####    figures  completes                          #######
###################################################################

fig_all_glm <- function(pair_dat, name, modelTransition, stats.fct = stats)
{

#models
# regeneration
modRT = modelTransition("R", "T", pair.dat = pair_dat)
modRB = modelTransition("R", "B", pair.dat = pair_dat)
modRM = modelTransition("R", "M", pair.dat = pair_dat)
# exclusion
modMT = modelTransition(c("M"), c("T"), pair.dat = pair_dat)
modMB = modelTransition(c("M"), c("B"), pair.dat = pair_dat)
# colonisation
modTM = modelTransition(c("T"), c("M"), pair.dat = pair_dat)
modBM = modelTransition(c("B"), c("M"), pair.dat = pair_dat)
# disturbance
modR = modelTransition(c("T", "B", "M"), c("R"), pair.dat = pair_dat)
modMR = modelTransition(c("M"), c("R"), pair.dat = pair_dat)
modTR = modelTransition("T", c("R"), pair.dat = pair_dat)
modBR = modelTransition(c("B"), c("R"), pair.dat = pair_dat)


pdf(paste("../figures/data_explo_glm", name, ".pdf",sep=""), width = 15, height = 8)

layout(matrix(c(1, 2, 3, 4, 5, 9, 6, 7, 10, 12, 8, 11), ncol = 4))
fig_glm(modRT)
fig_glm(modRB)
fig_glm(modRM)

fig_glm(modMT)
fig_glm(modMB)

fig_glm(modTM)
fig_glm(modBM)

fig_glm(modR)
fig_glm(modMR)
fig_glm(modTR)
fig_glm(modBR)

stats.fct(list(RT = modRT, RB = modRB, RM = modRM, MT = modMT, MB = modMB, TM = modTM, BM = modBM, R = modR, TR = modTR, BR = modBR, MR = modMR))

layout(matrix(1, ncol = 1))

stats.fct(list(RT = modRT, RB = modRB, RM = modRM, MT = modMT, MB = modMB, TM = modTM, BM = modBM, R = modR, TR = modTR, BR = modBR, MR = modMR))


dev.off()

}


fig_all_glm(pair_dat0, "", modelTransition = modelTransition_climate)



fig_all_glm(pair_dat1, "_filterHarvest", modelTransition = modelTransition_climate)


fig_all_glm(pair_dat2, "_filterHarvest+Drainage", modelTransition = modelTransition_climate)



#-------------------------------------------------------------------------------#-------------------------------------------------------------------------------

#### Maps

require("ggmap")
datmap = read.csv("~/Documents/GitHub/STModel-Data/out_files/statesFourState.csv")

lon_med <- median(datmap$longitude)
lat_med  <- median(datmap$latitude)

theme_set(theme_grey(base_size=10))
quicc.map = get_map(location = c(lon=lon_med,lat=lat_med), zoom = 4)

ggmap(quicc.map, maprange=FALSE,extent = "normal") %+% datmap + aes(x = longitude, y = latitude) +
#geom_point(size=0.4,aes=0.5) +
facet_wrap(~state) +
stat_density2d(aes(fill = ..level.., alpha = ..level..), bins = 10, geom = 'polygon') +
scale_fill_gradient(low = "red", high = "red4") +
#geom_density2d(color="orange4",size=.1) +
scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))+
xlab("Longitude") + ylab("Latitude") +
theme(legend.position = "none", text = element_text(size = 10))
