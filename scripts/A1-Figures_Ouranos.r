setwd("~/Documents/GitHub/STModel-Calibration/")

load("data/RandomForest_temp.rObj")
load("data/Multinom_temp.rObj")

require("randomForest")
require("nnet")
require("reshape2")
require("ggplot2")
require("grid")
require("maptools")
require("RColorBrewer")
require("rgeos")

zone_veg <- readShapePoly("~/Documents/Maitrise/Analyse/dom_ouranos/crop_shapefile/zone_veg.shp")
zone_veg@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
lakes <- readShapePoly("~/Documents/Maitrise/Analyse/dom_ouranos/crop_shapefile/lakes_qc.shp")
lakes@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
region_qc <- readShapePoly("~/Documents/Maitrise/Analyse/dom_ouranos/crop_shapefile/region_qc.shp")
region_qc@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")


#Subset sur l'aire des polygones (prend les lacs avec une superficie supérieur à 0.005 km²)
area <- gArea(lakes, byid=TRUE)
area <- area[area>0.005]
lakes <- lakes[names(area),]

# Transforme pour ggplot2
lakes <- fortify(lakes)
zone_veg <- fortify(zone_veg)
region_qc <- fortify(region_qc)
#rivers <- fortify(rivers)

#################################
## Get calibration spatial zone ##
#################################

data = read.csv("../STModel-Data/out_files/statesFourState.csv")
#data = read.csv("~/Documents/GitHub/STModel-Data/out_files/statesFourState.csv")
head(data)
dim(data)

plots <- unique(data[,2:3])
rg_lat <- range(plots$latitude)
rg_lon <- range(plots$longitude)

#################################
####  Prep data SMD   ###########
#################################

##############
# Carte de proba pour chaque SDM ('Multi', 'Random Forest')

climate_grid <- read.csv("~/Documents/Maitrise/Analyse/dom_ouranos/dom_climate_grid.csv")
climate_grid <- subset(climate_grid,lat>=rg_lat[1] & lat<=rg_lat[2])
climate_grid <- subset(climate_grid,lon>=rg_lon[1] & lon<=rg_lon[2])

# Prob SDM
prob_RF = cbind(climate_grid,as.data.frame(predict(SDM2,new=climate_grid,"prob", OOB=TRUE)))
prob_Multi = cbind(climate_grid,as.data.frame(predict(SDM1,new=climate_grid,"prob", OOB=TRUE)))
prob_RF = data.frame(prob_RF[,-c(3:4)],SDM=rep("Random forest",nrow(prob_RF)))
prob_Multi = data.frame(prob_Multi[,-c(3:4)],SDM=rep("Multinomiale",nrow(prob_Multi)))
prob_clim <- rbind(prob_RF,prob_Multi)

gg_prob_clim <- melt(prob_clim,id=c("lon","lat","SDM"))
gg_prob_RF <- melt(prob_RF[,-7],id=c("lon","lat"))
gg_prob_Multi <- melt(prob_Multi[,-7],id=c("lon","lat"))



#################################
#######  Carto SMD  #############
################################


ggplot_RF = ggplot(gg_prob_RF) +
    geom_polygon(data = subset(region_qc,hole==FALSE), aes(x = long, y = lat, group = group),fill="grey90",colour="grey60",size=0.05) +
    geom_raster(aes(lon,lat,fill=factor(cut(value,11),rev(levels(cut(value,11)))))) +
    facet_wrap(~variable) +
    scale_fill_manual(values=brewer.pal(11,"Spectral"),name="Probability")+
    geom_polygon(data = subset(lakes,hole==FALSE), aes(x = long, y = lat, group = group),fill="lightskyblue",color="deepskyblue2",size=0.1) +
    #geom_polygon(data = zone_veg, aes(x = long, y = lat, group = group),fill=NA,colour="grey90",size=0.2) +
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))+
    coord_equal() +
    xlab("Longitude") + ylab("Latitude")+
    theme(plot.title = element_text(lineheight=.8,face="bold"),
        panel.background = element_rect(fill = "lightskyblue"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey90",size=0.1),
        strip.text.x = element_text(size = 12,face="bold" ,colour = "white"),strip.background = element_rect(colour="black", fill="black"))+
    ggtitle("SDM: Random forest")


ggplot_multi = ggplot(gg_prob_Multi) +
    geom_polygon(data = subset(region_qc,hole==FALSE), aes(x = long, y = lat, group = group),fill="grey90",colour="grey60",size=0.05) +
    geom_raster(aes(lon,lat,fill=factor(cut(value,11),rev(levels(cut(value,11)))))) +
    facet_wrap(~variable) +
    scale_fill_manual(values=brewer.pal(11,"Spectral"),name="Probability")+
    geom_polygon(data = subset(lakes,hole==FALSE), aes(x = long, y = lat, group = group),fill="lightskyblue",color="deepskyblue2",size=0.1) +
    #geom_polygon(data = zone_veg, aes(x = long, y = lat, group = group),fill=NA,colour="grey90",size=0.2) +
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))+
    coord_equal() +
    xlab("Longitude") + ylab("Latitude")+
    theme(plot.title = element_text(lineheight=.8,face="bold"),
        panel.background = element_rect(fill = "lightskyblue"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey90",size=0.1),
        strip.text.x = element_text(size = 12,face="bold" ,colour = "white"),strip.background = element_rect(colour="black", fill="black"))+
    ggtitle("SDM: Multinomial")

ggsave(ggplot_RF,file="../data/proj_RF_qc.pdf",height=8,width=10)
ggsave(ggplot_multi,file="../data/proj_multi_qc.pdf",height=8,width=10)


#################################
#######  Carto State  ###########
#################################

str(prob_clim)
prob <- prob_clim[,c(3:6)]
head(prob)

# Draw states
draw = function(p) {
  states = c("B","M","R","T")
  draw = rmultinom(n=1,size=1,prob=p)
  states[which(draw==1)]
}


climate_grid <- read.csv("~/Documents/Maitrise/Analyse/dom_ouranos/dom_climate_grid.csv")
climate_grid <- subset(climate_grid,lat>=rg_lat[1] & lat<=rg_lat[2])
climate_grid <- subset(climate_grid,lon>=rg_lon[1] & lon<=rg_lon[2])


colors_state = c("darkcyan","palegreen3","black","orange")
States = apply(prob,1,draw)

gg_prob_state <- cbind(prob_clim,States)
gg_prob_state <- subset(gg_prob_state,SDM=='Random forest')
gg_prob_state <- gg_prob_state[,-c(3:7)]

theme_set(theme_grey(base_size=16))
ggplot_state = ggplot(gg_prob_state) +
    geom_polygon(data = subset(region_qc,hole==FALSE), aes(x = long, y = lat, group = group),fill="grey90",colour="grey60",size=0.05) +
    geom_raster(aes(lon,lat,fill=States),alpha=0.9) +
    scale_fill_manual(values=colors_state,name="States")+
    geom_polygon(data = subset(lakes,hole==FALSE), aes(x = long, y = lat, group = group),fill="lightskyblue",color="deepskyblue2",size=0.1) +
    #geom_polygon(data = zone_veg, aes(x = long, y = lat, group = group),fill=NA,colour="grey90",size=0.2) +
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))+
    coord_equal() +
    xlab("Longitude") + ylab("Latitude")+
    theme(plot.title = element_text(lineheight=.8,face="bold"),
        panel.background = element_rect(fill = "lightskyblue"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey90",size=0.1))

ggsave(ggplot_state,file="./data/proj_state_qc.pdf",height=7,width=8)

#################################
####  GLM Transition  ###########
#################################

###################################################################
#####    Reshaping data from STM-Data repo                  #######
###################################################################

# Open data
pair.dat <- read.csv("../STModel-Data/out_files/transitionsFourState.csv")

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

###############################
#####    GLM CLIMATE    #######
###############################

library(MASS)
library(ROCR)
library(fmsb)

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

modRT = modelTransition_climate("R", "T", pair.dat = pair_dat0)
modRB = modelTransition_climate("R", "B", pair.dat = pair_dat0)
modRM = modelTransition_climate("R", "M", pair.dat = pair_dat0)
# exclusion
modMT = modelTransition_climate(c("M"), c("T"), pair.dat = pair_dat0)
modMB = modelTransition_climate(c("M"), c("B"), pair.dat = pair_dat0)
# colonisation
modTM = modelTransition_climate(c("T"), c("M"), pair.dat = pair_dat0)
modBM = modelTransition_climate(c("B"), c("M"), pair.dat = pair_dat0)
# disturbance
modR = modelTransition_climate(c("T", "B", "M"), c("R"), pair.dat = pair_dat0)
modMR = modelTransition_climate(c("M"), c("R"), pair.dat = pair_dat0)
modTR = modelTransition_climate("T", c("R"), pair.dat = pair_dat0)
modBR = modelTransition_climate(c("B"), c("R"), pair.dat = pair_dat0)

models = list(modRT=modRT,modRB=modRB,modRM=modRM,modR=modR)
prob_grid = list()

for (i in 1:length(models)){
    mod = models[[i]]
    temp = seq(-5,25, length.out = 500)
    pp = seq(500,2000, length.out = 500)
    prob = predict(mod$mod, newdata = data.frame(expand.grid(annual_mean_temp = temp, annual_pp = pp) ), type = "response")
    grid = expand.grid(annual_mean_temp = temp, annual_pp = pp)
    prob_grid[[i]] = data.frame(grid,prob=prob,transition=rep(names(models)[i],length(pp)))
}

gg_prob_dat <- do.call(rbind,prob_grid)
gg_prob_dat$class_prob <- cut(gg_prob_dat$prob,11)
gg_prob_dat$transition <- factor(gg_prob_dat$transition,labels=c("R -> T","R -> B","R -> M", "(T,B,M) -> R"))

theme_set(theme_grey(base_size=14))
ggplot(gg_prob_dat, aes(x=annual_mean_temp,y=annual_pp)) + geom_raster(aes(fill=class_prob)) + facet_wrap(~transition) + scale_fill_manual(values = rev(brewer.pal(11,"RdYlBu")),guide = guide_legend(reverse=TRUE),name="Class of probability")+
    xlab("Annual mean temperature (°C)") + ylab("Precipitation (meter)")
ggsave(file="Transition_prob_glm.pdf",width=12,height=8)



