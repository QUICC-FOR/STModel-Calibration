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
#######  Carto SMD  #############
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

ggsave(ggplot_RF,file="./data/proj_RF_qc.pdf",height=8,width=10)
ggsave(ggplot_multi,file="./data/proj_multi_qc.pdf",height=8,width=10)


#################################
#######  Carto State  ###########
#################################

prob <- prob_clim[,c(5,3,6,4)]

# Draw states
draw = function(p) {
  states = c("R","B","T","M")
  draw = rmultinom(n=1,size=1,prob=p)
  states[which(draw==1)]
}


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