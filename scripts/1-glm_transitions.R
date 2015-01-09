
rm(list = ls())

# Open data
pair.dat <- read.csv("../data/transitionsFourState.csv")

# subset 10 degree
select = unique(dat_wo_U$plot[which(dat_wo_U$annual_mean_temp<=10)])
dat_subset10 = dat_wo_U[dat_wo_U$plot %in% select,]

# Rename and clean columns
names(pair.dat)[7:8] <- c("st0","st1")
#pair.dat <- pair.dat[,-c(9:12)]

# Create transition column
pair.dat$transition <- paste(pair.dat$st0,pair.dat$st1,sep="")

# Datset without filters
pair_dat0 <- pair.dat

# Graph lim
rg_pp <- range(pair.dat$annual_pp)
rg_tp <- range(pair.dat$annual_mean_temp)



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


colo = c(R = rgb(.5,.5,.5,.5), T = rgb(1,0,0,.5), B = rgb(0.2,.8,.2,.5), M = rgb(0,0,1,.5))


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


pdf(paste("../figures/glm_transitions", name, ".pdf",sep=""), width = 15, height = 8)

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




#-------------------------------------------------------------------------------#-------------------------------------------------------------------------------
