###############
# Carte de proba basé sur les deux SDM

climate_grid <- read.csv("~/Documents/Maitrise/Analyse/dom_ouranos/dom_climate_grid.csv")

# Prob SDM
prob_RF = cbind(climate_grid,as.data.frame(predict(SDM2,new=climate_grid,"prob", OOB=TRUE)))
prob_Multi = cbind(climate_grid,as.data.frame(predict(SDM1,new=climate_grid,"prob", OOB=TRUE)))
prob_RF = data.frame(prob_RF[,-c(3:4)],SDM=rep("Random forest",nrow(prob_RF)))
prob_Multi = data.frame(prob_Multi[,-c(3:4)],SDM=rep("Multinomiale",nrow(prob_Multi)))
prob_clim <- rbind(prob_RF,prob_Multi)

require("reshape2")
require("ggplot2")
require("grid")
require("maptools")
require("RColorBrewer")
require("rgeos")

gg_prob_clim <- melt(prob_clim,id=c("lon","lat","SDM"))
gg_prob_RF <- melt(prob_RF[,-7],id=c("lon","lat"))
gg_prob_Multi <- melt(prob_Multi[,-7],id=c("lon","lat"))

zone_veg <- readShapePoly("~/Documents/Maitrise/Analyse/dom_ouranos/crop_shapefile/zone_veg.shp")
zone_veg@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
lakes <- readShapePoly("~/Documents/Maitrise/Analyse/dom_ouranos/crop_shapefile/lakes_qc.shp")
lakes@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

#Subset sur l'aire des polygones
area <- gArea(lakes, byid=TRUE)
area <- area[area>0.005]
lakes <- lakes[names(area),]
# Transforme pour ggplot2
lakes <- fortify(lakes)
zone_veg <- fortify(zone_veg)
#rivers <- fortify(rivers)

#################################
##   Compute convexe Polygon   ##
#################################

require("grDevices")

data = read.csv("~/Documents/GitHub/STModel-Data/out_files/statesFourState.csv")
#data = read.csv("~/Documents/GitHub/STModel-Data/out_files/statesFourState.csv")
head(data)
dim(data)

plots <- unique(data[,2:3])

xy.coords(chull(plots))



#################################
#######     Carto   #############
#################################

ggplot_RF = ggplot(gg_prob_RF) +
    geom_raster(aes(lon,lat,fill=factor(cut(value,9),rev(levels(cut(value,9)))))) +
    facet_wrap(~variable) +
    scale_fill_manual(values=rev(brewer.pal(9,"Reds")),name="Probability")+
    geom_polygon(data = subset(lakes,hole==FALSE), aes(x = long, y = lat, group = group),fill="lightskyblue",colour="dodgerblue2",size=0.05) +
    geom_polygon(data = zone_veg, aes(x = long, y = lat, group = group),fill=NA,colour="grey90",size=0.2) +
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))+
    coord_equal() +
    xlab("Longitude") + ylab("Latitude")+
    theme(plot.title = element_text(lineheight=.8,face="bold"),
        panel.background = element_rect(fill = "lightskyblue"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "deepskyblue3"),
        strip.text.x = element_text(size = 12,face="bold" ,colour = "white"),strip.background = element_rect(colour="black", fill="black"))+
    ggtitle("SDM: Random forest")


ggplot_multi = ggplot(gg_prob_Multi) +
    geom_raster(aes(lon,lat,fill=factor(cut(value,9),rev(levels(cut(value,9)))))) +
    facet_wrap(~variable) +
    scale_fill_manual(values=rev(brewer.pal(9,"Reds")),name="Probability")+
    geom_polygon(data = subset(lakes,hole==FALSE), aes(x = long, y = lat, group = group),fill="lightskyblue",colour="dodgerblue2",size=0.05) +
    geom_polygon(data = zone_veg, aes(x = long, y = lat, group = group),fill=NA,colour="grey90",size=0.2) +
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))+
    coord_equal() +
    xlab("Longitude") + ylab("Latitude")+
    theme(plot.title = element_text(lineheight=.8,face="bold"),
        panel.background = element_rect(fill = "lightskyblue"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "deepskyblue3"),
        strip.text.x = element_text(size = 12,face="bold" ,colour = "white"),strip.background = element_rect(colour="black", fill="black"))+
    ggtitle("SDM: Multinomial")

ggsave(ggplot_RF,file="./Rout_files/proj_RF_qc.pdf",height=8,width=10)
ggsave(ggplot_multi,file="./Rout_files/proj_multi_qc.pdf",height=8,width=10)