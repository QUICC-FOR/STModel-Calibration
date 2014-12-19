rm(list=ls())
# Load data
#dataProj = as.data.frame(read.table("../data/data_pairs_filter.txt"))
#dataProj$E = scale(dataProj$annual_mean_temp)
#dataProj$P = scale(dataProj$annual_pp)

#data = read.csv("../data/statesFourState.csv")
data = read.csv("~/Documents/GitHub/STModel-Data/out_files/statesFourState.csv")
head(data)
dim(data)

setwd("~/Documents/GitHub/STModel-Calibration/scripts")

# ----------------------
# Clean data
# ---------------------

# Clean Undefined state
dat_wo_U <- subset(data, state != "U")
dat_wo_U$state <- droplevels(dat_wo_U$state)

# ----------------------
### choice of variables
# ----------------------

varCor = cor(dat_wo_U[,-c(1:5)]) # check correlation between variables

# ----------------------
### ACP & Correlation
# ----------------------

library("ade4")
var.pca = dudi.pca(dat_wo_U[,-c(1:5)], scannf=FALSE, nf = 5)
var.pca$eig / sum(var.pca$eig)

# ----------------------
### Analyze de coinertie
# ----------------------

com <- acm.disjonctif(dat_wo_U[,c("plot","state")])

com <- dat_wo_U[,"state"]
com <- data.frame(rec=seq(1,length(com),1),state=com,val=rep(1,length(com)))
com <- dcast(com,rec~state, value.var="val")
com[is.na(com)]<-0
com <- com[,-1]

coi <- niche(var.pca,com,scannf=F,nf=4)
s.class(coi$ls,fac=dat_wo_U$state,cpoint=0,cstar=0,col=colors_state)
s.arrow(30*coi$co[selectedVars,],clab=0.8,add.plot=TRUE)
s.label(coi$li,add.plot=TRUE)

coi$eig/sum(coi$eig)
s.corcircle(coi$co, clab = 0.5)


# s.corcircle(var.pca$co, clab = 0.5)
# contrib = inertia.dudi(var.pca, row = FALSE, col = TRUE)$col.abs

# nonCorVars = intersect(names(varCor[which(abs(varCor[,"annual_mean_temp"])<0.7),"annual_mean_temp"]), names(varCor[which(abs(varCor[,"pp_seasonality"])<0.7),"pp_seasonality"]))

# contrib[nonCorVars,]

# nonCorVars = intersect(nonCorVars, names(varCor[which(abs(varCor[,"pp_warmest_quarter"])<0.7),"pp_warmest_quarter"]))

# nonCorVars = intersect(nonCorVars, names(varCor[which(abs(varCor[,"mean_diurnal_range"])<0.7),"mean_diurnal_range"]))

# nonCorVars = intersect(nonCorVars, names(varCor[which(abs(varCor[,"annual_pp"])<0.7),"annual_pp"]))

selectedVars = c("annual_mean_temp", "pp_seasonality", "pp_warmest_quarter", "mean_diurnal_range","annual_pp", "mean_temperatre_wettest_quarter")

varCor2 = cor(dat_wo_U[, selectedVars])
varCor2

datSel = dat_wo_U[,c("state",selectedVars)]


colors_state = c("darkcyan","palegreen3","black","orange","pink")

s.class(var.pca$li,fac=datSel$state,cpoint=0,cstar=0,col=colors_state)
s.arrow(10*var.pca$co[selectedVars,],clab=0.8,,add.plot=TRUE)


###########################################################################################################

# ----------------------
### Explo Data et PCA Steve
# ----------------------

require(reshape2)
require(ggplot2)

ggdata <- melt(datSel_wo_U,id=c("state"))
ggdata_all <- data[,5:34]

# Histograme
hists = ggplot(ggdata) + geom_histogram(aes(x=value,fill=state)) + facet_grid(state~variable,scales="free") + scale_fill_brewer(palette="Accent","State")
ggsave(hists,file="../figures/explo_hist_selVars_by_state.pdf",width=15,height=7)

#install.packages("GGally")
require(GGally)

ggplot <- function(...) ggplot2::ggplot(...) + scale_color_brewer(palette="Accent")

ggpairs(ggdata_all,
        columns=2:ncol(ggdata_all),
        lower = list(continuous = "density"), # data.frame with variables
        title="Climate exploration by state", # title of the plot
        colour = "state") # aesthetics, ggplot2 style

dev.copy2pdf(height=16,width=18,out.type="pdf")
dev.off()