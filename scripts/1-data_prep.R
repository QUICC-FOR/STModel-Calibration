###################################################################
#####                        Clean data                     #######
###################################################################

dat = read.table("../data/data_BA.txt",header = TRUE, sep = ";")
#head(dat)
dim(dat)

## Rm all harvested plots 
dat  <- dat[-which(dat$disturbance%in%c("BRP", "CAM", "CB", "CD", "CDL", "CE", "CJ", "CP", "DLD", "EPC")),]
dim(dat)
# BRP = brulis partiel
# CAM = coupe d'amélioration
# CB = coupe par bandes
# CD = coupes en damier
# CDL = coupe à diamètre limité
# CE = coupe partielle et épidémie légère
# CHP = chablis partiel
# CJ = coupe de jardinage
# CP = coupe partielle
# DLD = Coupe à diamètre limite avec dégagement des arbres d'avenir
# DP  = Dépérissement partiel du feuillu
# EL = Épidémie légère
# EPC = Éclaircie précommerciale
# VEP = Verglas partiel


## Rm all plots with no climatic data associated
dat  <- dat[which(!is.na(dat$annual_pp)),]
dim(dat)

## Conserve all plots with drainage 20 to 40
dat <- dat[which(dat$drainage>=20 & dat$drainage<=40),]
dim(dat)

## rm plot with no coordinates constraints
dat <- dat[-which(dat$lon==0),]
dim(dat)

###################################################################
#####                  Classify plots                       #######
###################################################################

# List of interests species
I_sp  <- c("bop","peb","peg","pet","prp","sal","soa")
B_sp  <- c("epn","epb","epr", "mel","pig","sab") 
T_sp  <- c("boj", "chr", "err","ers","fra","frn","heg","osv","til","cet")

# Subset species BA and cover type observed
BA_dat <- c(which(colnames(dat)=="ame"): which(colnames(dat)=="til"))

### Get BA by hectares (the plot is 400m2)
dat[,BA_dat] <- dat[,BA_dat]*10000/400

dat$I_tot  <- rowSums(dat[,which(colnames(dat) %in% I_sp)],na.rm=T)
dat$B_tot  <- rowSums(dat[,which(colnames(dat) %in% B_sp)],na.rm=T)
dat$T_tot  <- rowSums(dat[,which(colnames(dat) %in% T_sp)],na.rm=T)
dat$BA_tot <- rowSums(dat[,BA_dat],na.rm=T)

# Class into state types

class_fn = function(x) { 
  Itot = x[1]; Btot =x[2]; Ttot=x[3]; BAtot=x[4]
  classPlot = "Unclass"
  if(BAtot < 10) classPlot = "R"
  else if(Btot > 0 & Ttot == 0) classPlot = "B"
  else if(Ttot > 0 & Btot == 0) classPlot = "T"
  else if (Btot > 0 & Ttot > 0) classPlot = "M"
return(classPlot)
}

dat$class_final  <- as.vector(apply(dat[,c("I_tot", "B_tot", "T_tot", "BA_tot")],1,class_fn))

table(dat$class_final)
#dat[dat$class_final=="Unclass",]

dat = dat[-which(dat$class_final=="Unclass"),]
dim(dat)

write.table(dat,file="../data/data_allyears_RBTM.txt")


###################################################################
#####    Reshape (pairwise) and export data                 #######
###################################################################

liste =split(dat, dat$id_plot)
nreleves = unlist(lapply(liste, nrow))

# remove only one releve in time
liste = liste[-which(nreleves==1)]

## Function to match remeasurements on the same line 
prwise  <- function(x){
    #climatic columns
    clim = c(which(colnames(x) == "annual_pp") :which(colnames(x) == "annual_max_temp"))
        
    #number of pairs
    npairs = nrow(x)-1

    for (i in 1:npairs)
    {
    yr0 = min(x$yr_measured)
    st0 = x$class_final[which.min(x$yr_measured)]
    clim0 = x[which.min(x$yr_measured), clim ]
    x = x[-which.min(x$yr_measured),]
    yr1 = min(x$yr_measured)
    st1 = x$class_final[which.min(x$yr_measured)]
    clim1 = x[which.min(x$yr_measured), clim ]
    if(i ==1) {
    res = data.frame(st0, st1, yr0, yr1, lat = x$lat, lon = x$lon,  t(apply(rbind(clim0, clim1), 2, mean)))
#    res = data.frame(st0, st1, yr0, yr1, lat = x$lat, lon = x$lon,  clim0, clim1)
    }else{
    res = rbind(res, data.frame(st0, st1, yr0, yr1,lat = x$lat, lon = x$lon, t(apply(rbind(clim0, clim1), 2, mean))))
#    res = data.frame(st0, st1, yr0, yr1, lat = x$lat, lon = x$lon,  clim0, clim1)
    }
    }

  return(res)
}

pair_dat = lapply(liste, prwise)
pair_dat = (do.call(rbind.data.frame, pair_dat))

dim(pair_dat)

# ------------------- write

write.table(pair_dat,file="../data/data_reshaped_RBTM.txt")

# ---------------- herbivores
reshape_dat = read.table("../data/data_reshaped_RBTM.txt")
head(reshape_dat)
dim(reshape_dat)
reshape_dat$id_plot = rownames(reshape_dat)

moose = read.table("../data/densities_moose.txt", h=T, sep="\t")
colnames(moose) = c("NO_ZONE"     ,  "id_plot", "moose")
head(moose)
dim(moose)

deer = read.table("../data/densities_deer.txt", h=T, sep="\t")
colnames(deer) = c("NO_ZONE"     ,  "id_plot", "deer")
head(deer)
dim(deer)

herbivores = merge(moose[,-1], deer[,-1], by = 'id_plot')
dim(herbivores)
head(herbivores)
herbivores[is.na(herbivores$moose),]
herbivores[is.na(herbivores$deer),]

finalTab = merge(reshape_dat, herbivores, by = "id_plot")
head(finalTab)
dim(finalTab)

#library(foreign)
#coords = read.csv("../data/plot_coords.csv")
#head(coords)
#coords$lat[which(coords$lat==0.0)]=NA
#coords = na.omit(coords)
#summary(coords)
#
#fint = merge(finalTab, coords, by = 'id_plot', all.x=TRUE, all.y = FALSE)
#head(unique(fint))

write.table(finalTab,file="../data/data_reshaped_RBTM_herbivores.txt")


###

dat = read.table("../data/data_reshaped_RBTM_herbivores.txt")
head(dat)

transitions = paste(dat$st0, dat$st1, sep="")

par(mfrow = c(2, 4))
boxplot(dat$annual_mean_temp~ifelse(transitions=="RT", 1, 0), xlab = "RT", ylab = "T°", main ="regeneration")
boxplot(dat$annual_mean_temp~ifelse(transitions=="BM", 1, 0), xlab = "BM", ylab = "T°", main ="invasion")
boxplot(dat$annual_mean_temp~ ifelse(transitions=="MT", 1, 0), xlab = "MT", ylab = "T°", main = "succession")
boxplot(dat$annual_mean_temp~ ifelse(transitions=="TT", 1, 0), xlab = "TT", ylab = "T°", main = "stays")

boxplot(dat$annual_mean_temp~ ifelse(transitions=="RB", 1, 0), xlab = "RB", ylab = "T°", main = "regeneration")
boxplot(dat$annual_mean_temp~ ifelse(transitions=="TM", 1, 0), xlab = "TM", ylab = "T°", main ="invasion")
boxplot(dat$annual_mean_temp~ ifelse(transitions=="MB", 1, 0), xlab = "MB", ylab = "T°", main = "succession")
boxplot(dat$annual_mean_temp~ ifelse(transitions=="BB", 1, 0), xlab = "BB", ylab = "T°", main = "stays")

dev.copy2pdf(file="../figures/observation transitions.pdf")
