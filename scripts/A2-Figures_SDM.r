# Projection of the SDM
# Date: December 12th, 2014

# Two SDM methods:
# - RandomForest
# - Multinomial (nnet package)

# Climatic vars selected for the SDM:
    # 'annual_mean_temp'
    # 'pp_seasonality'
    # 'pp_warmest_quarter'
    # 'mean_diurnal_range'
    # 'tot_annual_pp'
    # 'mean_temperature_wettest_quarter'

####################################
###### Database connection
####################################

setwd('~/Documents/GitHub/STModel-Calibration/scripts/')
source('../../STModel-Data/con_quicc_db_local.r')

#Load librairies
require('reshape2')

# Query
query_SDMClimate_grid  <- " SELECT ST_X(geom) as lon, ST_Y(geom) as lat, val, biovar FROM (
    SELECT biovar, (ST_PixelAsCentroids(ST_Transform(rasters,4326))).* FROM (
    SELECT biovar, ST_Union(ST_Clip(ST_Resample(rast,ref_rast),env_stm.env_plots),'MEAN') as rasters
    FROM
    (SELECT rast, biovar,year_clim FROM clim_rs.clim_allbiovars
     WHERE (year_clim >= 1970 AND year_clim <= 2000)
    AND biovar IN (
    'annual_mean_temp',
    'pp_seasonality',
    'pp_warmest_quarter',
    'mean_diurnal_range',
    'tot_annual_pp',
    'mean_temp_wettest_quarter')) AS rast_noram,
    (SELECT ST_Transform(ST_ConvexHull(ST_Collect(stm_plot_ids.coord_postgis)),4269) as env_plots FROM rdb_quicc.stm_plot_ids) AS env_stm,
    (SELECT rast as ref_rast FROM clim_rs.clim_allbiovars WHERE biovar = 'annual_mean_temp' LIMIT 1) as ref
    WHERE ST_Intersects(rast_noram.rast,env_stm.env_plots)
    GROUP BY biovar) AS union_query
) AS points_query;"

## Send the query to the database
res_SDMClimate_grid <- dbGetQuery(con, query_SDMClimate_grid)
## Time: Approx. 5-15 minutes

# Clean data
res_SDMClimate_grid$biovar <- as.factor(res_SDMClimate_grid$biovar)

# Remove climatic spaces outside of the SDM calibration
stm_dat <- read.csv('../../STModel-Data/out_files/statesFourState.csv')
str(stm_dat)
selectedVars = c("annual_mean_temp", "pp_seasonality", "pp_warmest_quarter", "mean_diurnal_range","annual_pp", "mean_temperatre_wettest_quarter")
stm_dat <- stm_dat[,selectedVars]
rg_amt <- range(stm_dat[,1])
rg_ppseas <- range(stm_dat[,2])
rg_ppwarm <- range(stm_dat[,3])
rg_md <- range(stm_dat[,4])
rg_app <- range(stm_dat[,5])
rg_mtwet <- range(stm_dat[,6])


# # Calib space
spaceClim_SDM <- dcast(res_SDMClimate_grid,lon+lat ~ biovar, value.var="val")
spaceClim_SDM[,5] <- spaceClim_SDM[,5]/10
spaceClim_SDM <- subset(spaceClim_SDM, annual_mean_temp <= rg_amt[2] & annual_mean_temp >= rg_amt[1])
spaceClim_SDM <- subset(spaceClim_SDM, pp_seasonality <= rg_ppseas[2] & pp_seasonality >= rg_ppseas[1])
spaceClim_SDM <- subset(spaceClim_SDM, pp_warmest_quarter <= rg_ppwarm[2] & pp_warmest_quarter >= rg_ppwarm[1])
spaceClim_SDM <- subset(spaceClim_SDM, mean_diurnal_range <= rg_md[2] & mean_diurnal_range >= rg_md[1])
spaceClim_SDM <- subset(spaceClim_SDM, tot_annual_pp <= rg_app[2] & tot_annual_pp >= rg_app[1])
spaceClim_SDM <- subset(spaceClim_SDM, mean_temp_wettest_quarter <= rg_mtwet[2] & mean_temp_wettest_quarter >= rg_mtwet[1])

###### Load shapefiles

require("ggplot2")
require("maptools")
require("raster")
require("reshape2")
require("grid")
require("rgeos")
require("RColorBrewer")

load("../data/shp_lakes_regions_study_area.Robj")

###### Predicted map
names(spaceClim_SDM)[8] <- "annual_pp"
names(spaceClim_SDM)[5] <- "mean_temperatre_wettest_quarter"

require("randomForest")
require("nnet")
require("reshape2")

load('../data/Multinom_6vars_version_b.rObj')
load('../data/RandomForest_7vars.rObj')

pred_multinom <- predict(SDM1.b,new=pred_dat[,-c(1:2)],"prob")
pred_RF <- predict(SDM2,new=pred_dat[,-c(1:2)],"prob")

pred <- pred_multinom

pred <- cbind(pred_dat[,c(1,2)],pred)
pred <- melt(pred,id=c("lon","lat"))
names(pred)<-c("lon","lat","state","prob")
pred$prob <- cut(pred$prob,11)
pred$prob <- factor(pred$prob,levels=rev(levels(pred$prob)))

cols <- brewer.pal(11,"Spectral")

map_pred <- ggplot(pred) +
geom_polygon(data = regions_usa, aes(x = long, y = lat, group = group),fill="grey80",colour="grey50",size=0.1) +
geom_polygon(data = regions_can, aes(x = long, y = lat, group = group),fill="grey80",colour="grey50",size=0.1) +
geom_raster(aes(lon,lat,fill=prob)) +
facet_wrap(~state)+
scale_fill_manual(values=cols,name="Class") +
geom_polygon(data = subset(lakes_usa,hole==FALSE), aes(x = long, y = lat, group = group),fill="lightskyblue",colour="dodgerblue4",size=0.1) +
geom_polygon(data = subset(lakes_can,hole==FALSE), aes(x = long, y = lat, group = group),fill="lightskyblue",colour="dodgerblue4",size=0.1) +
scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))+
coord_equal() +
xlab("Longitude") + ylab("Latitude")+
theme(plot.title = element_text(face="bold"),
    panel.background = element_rect(fill = "lightskyblue"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90",size=0.1),
    strip.text.x = element_text(size = 12,face="bold" ,colour = "white"),strip.background = element_rect(colour="black", fill="black"))

ggsave(map_pred,file="../figures/multinom_proj_SDM.jpg",width=15,height=12)
#«
###### Climatic map

map_var <- function(dat,var,pal,n){

    datvar <- subset(dat,biovar == var)
    datvar$val <- cut(datvar$val,n, dig.lab = 2)
    datvar$val <- factor(datvar$val,levels=rev(levels(datvar$val)))

    cols <- brewer.pal(n,pal)

    ggplot(datvar) +
    geom_polygon(data = regions_usa, aes(x = long, y = lat, group = group),fill="grey80",colour="grey50",size=0.1) +
    geom_polygon(data = regions_can, aes(x = long, y = lat, group = group),fill="grey80",colour="grey50",size=0.1) +
    geom_raster(aes(lon,lat,fill=val, order = val)) +
    scale_fill_manual(values=cols,name="Class") +
    geom_polygon(data = subset(lakes_usa,hole==FALSE), aes(x = long, y = lat, group = group),fill="lightskyblue",colour="dodgerblue4",size=0.1) +
    geom_polygon(data = subset(lakes_can,hole==FALSE), aes(x = long, y = lat, group = group),fill="lightskyblue",colour="dodgerblue4",size=0.1) +
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))+
    coord_equal() +
    xlab("Longitude") + ylab("Latitude")+
    theme(plot.title = element_text(face="bold"),
        panel.background = element_rect(fill = "lightskyblue"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey90",size=0.1),
        strip.text.x = element_text(size = 12,face="bold" ,colour = "white"),strip.background = element_rect(colour="black", fill="black"))+
    ggtitle(paste(var,"(Horizon 1970-2000)"))
}

res_SDMClimate_grid$biovar <- factor(res_SDMClimate_grid$biovar,labels=c(
                    "Annual mean temperature (°C)",
                    "Mean diurnal range (JD)",
                    "Mean temperature of the wettest period (°C)",
                    "Precipitation seasonality (mm)",
                    "Precipitation of the warmest quarter (mm)",
                    "Annual precipitation (mm)"))


 for(i in 1:nlevels(res_SDMClimate_grid$biovar)){
    map_var(res_SDMClimate_grid,levels(res_SDMClimate_grid$biovar)[i],"Spectral",11)
    ggsave(paste("../figures/var",i,"_1970-2000.jpg",sep=""),height=5,width=7)
 }
