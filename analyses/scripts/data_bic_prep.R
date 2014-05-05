

setwd("/home/DominiqueGravel/Bureau/analyses/data")
bic = read.table("tree_coord_bic.txt",header=T,sep = ";")
bic = subset(bic, bic$DHP>100)
ID = paste(bic$BORX,bic$BORY,sep="-")
BA = 25*((bic$DHP/2000)^2)
ID_list = unique(ID)
XY = matrix(as.numeric(unlist(strsplit(ID_list,"-"))),nr=length(ID_list),nc=2,byrow=T)

T_sp  <- c("BEPA","QURU","POBA","POGR","POTR","PRPE","SOAM","SOAU","SODE","ACSP","ACPE")
C_sp  <- c("ABBA","PIAB","PIGL","PIMA","PIRE","PIST","THOC")
D_sp  <- c("ACRU","ACSA")

state_sp = numeric(nrow(bic))

state_sp[bic$Esp%in%C_sp]="C"
state_sp[bic$Esp%in%D_sp]="D"
state_sp[bic$Esp%in%T_sp]="T"

table_bic = as.data.frame(tapply(X=BA,INDEX = list(ID,state_sp),FUN = sum))

sum.rm = function(x) sum(x,na.rm=T)
BA_quad = apply(table_bic,1,sum.rm)
table_bic_rel = table_bic/BA_quad

state_quad = numeric(nrow(table_bic_rel))
state_quad[table_bic_rel[,1]>1/3]="Unclass"
state_quad[table_bic_rel[,2]>2/3]="C"
state_quad[table_bic_rel[,3]>2/3]="D"
state_quad[table_bic_rel[,4]>2/3]="T"
state_quad[table_bic_rel[,2]>2/3]="C"
state_quad[state_quad==0]="M"

write.table(cbind(XY,state_quad),"state_BIC.txt")

####################################################

setwd("/home/DominiqueGravel/Bureau/analyses/data")
sut = read.table("tree_coord_sut.txt",header=T,sep = ";")
sut = subset(sut, sut$DHP>100)
ID = paste(sut$BORX,sut$BORY,sep="-")
BA = 25*((sut$DHP/2000)^2)
ID_list = unique(ID)
XY = matrix(as.numeric(unlist(strsplit(ID_list,"-"))),nr=length(ID_list),nc=2,byrow=T)

T_sp  <- c("BEPA","BEAL","beal","QURU","POBA","POGR","POTR","PRPE","SOAM","SOAU","SODE","ACSP","ACPE","bepa","sode","acpe","AMSP")
C_sp  <- c("ABBA","PIAB","PIGL","PIMA","PIRE","PIST","THOC","abba","piru","tsca","PIRU","TSCA")
D_sp  <- c("ACRU","ACSA","FAGR","fagr")

state_sp = numeric(nrow(sut))

state_sp[sut$Esp%in%C_sp]="C"
state_sp[sut$Esp%in%D_sp]="D"
state_sp[sut$Esp%in%T_sp]="T"

table_sut = as.data.frame(tapply(X=BA,INDEX = list(ID,state_sp),FUN = sum))

sum.rm = function(x) sum(x,na.rm=T)
BA_quad = apply(table_sut,1,sum.rm)
table_sut_rel = table_sut/BA_quad

state_quad = numeric(nrow(table_sut_rel))
state_quad[table_sut_rel[,1]>1/3]="Unclass"
state_quad[table_sut_rel[,2]>2/3]="C"
state_quad[table_sut_rel[,3]>2/3]="D"
state_quad[table_sut_rel[,4]>2/3]="T"
state_quad[state_quad==0]="M"

XY = rbind(XY,c(40,560))
XY = rbind(XY,c(100,740))
state_quad = c(state_quad,c("M","M"))

state_quad = state_quad[order(XY[,1],XY[,2])]
XY = XY[order(XY[,1],XY[,2]),]

write.table(cbind(XY,state_quad),"state_SUT.txt")



