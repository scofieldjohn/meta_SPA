####load packages
library(pwr)
library(grid)
library(weightr)
library(meta)
library(metafor)
library(ez)
library(BayesMed)
library(BayesFactor)
library(Hmisc)
library(devtools)
devtools::install_github("RobbievanAert/puniform")
library(puniform)
devtools::install_github("doomlab/MOTE")
library(MOTE)


####load data
metaALL <- read.csv("~/John M.S. Experimental Psychology/Research and Articles/Evolutionary Psychology/SPA Meta-Analysis/Data analysis/data analysis round 3/metaALL.csv")

####misc exploratory tests
#simple regression with ES and N
Regmodel = lm(Recalculated.Eta~N, data=metaALL)
summary(Regmodel) #no sig differences
#Bayes regression
BFRegmodel = regressionBF(formula = Recalculated.Eta~N,
                        data = metaALL)
BFRegmodel
plot(BFRegmodel)

#t-test for research design
ttestmodel = t.test(metaALL$Recalculated.Eta~metaALL$Design)
ttestmodel
tapply(metaALL$Recalculated.Eta, list(metaALL$Design), mean)
tapply(metaALL$Recalculated.Eta, list(metaALL$Design), sd)
tapply(metaALL$Recalculated.Eta, list(metaALL$Design), length)
CohensD = (0.2019733-0.1131946)/(sqrt(((48*0.14018604^2)+(40*0.09381229^2))/(49+41-2)))
CohensD 
#Bayes t-test
ttestBF(formula = Recalculated.Eta~Design, data = metaALL)

ttestmodel2 = t.test(metaALL$Individual.Power~metaALL$Design)
ttestmodel2
tapply(metaALL$Individual.Power, list(metaALL$Design), mean)
tapply(metaALL$Individual.Power, list(metaALL$Design), sd)
tapply(metaALL$Individual.Power, list(metaALL$Design), length)
CohensD2 = (0.8125557-0.7480047)/(sqrt(((37*0.2188377^2)+(36*0.2465742^2))/(37+38-2)))
CohensD2
#BF
ttestBF(formula = Individual.Power~Design, data = metaALL)

#correlations
options(scipen=999)
cor.test(metaALL$Individual.Power , metaALL$Recalculated.Eta, method = "pearson")
cor.test(metaALL$Individual.Power , metaALL$N, method = "pearson")
cor.test(metaALL$Fisher.s.z , metaALL$SE.of.Z, method = "pearson")
jzs_cor(metaALL$Individual.Power, metaALL$Recalculated.Eta)
jzs_cor(metaALL$Individual.Power, metaALL$N)
jzs_cor(metaALL$Fisher.s.z, metaALL$SE.of.Z)

cor.test(metaALL$Recalculated.Eta , metaALL$N, method = "pearson")
jzs_cor(metaALL$Recalculated.Eta , metaALL$N)



####function to convert z back to eta
z2e = function(z){
  corrFE=(exp(z/0.5)-1)/(exp(z/0.5)+1)
  etaFE=corrFE^2
  etaFE
}

####initial full meta to estimate heterogeneity
#only interpreting heterogeneity
#not meta-analytic effect size
meta = metagen(Fisher.s.z, SE.of.Z, data = metaALL, method.tau="PM")
meta #plus the CI's for individual studies are calculated from normal approximations


#translating CI's back to eta to put into spreadsheet
IndividualCIlo = z2e(meta$lower)
IndividualCIhi = z2e(meta$upper)
IndividualCIlo
IndividualCIhi

#don't report meta results like these below until you exclude outliers
#FE ES
z2e(meta$TE.fixed)
#FE CI lo
z2e(meta$lower.fixed)
#FE CI hi
z2e(meta$upper.fixed)
#RE ES
z2e(meta$TE.random)
#RE CI lo
z2e(meta$lower.random)
#RE CI hi
z2e(meta$upper.random)

####Individually calculated all non-central CI's from test statistic
#calculate individual CI's for all experiments for non-central distribution
#t = t-test value
#n = sample size
#a = alpha type 1 error rate, p < .05 example
#k = number of decimals
#dfm = df model (between subjects, numerator)
#dfe = df error (within subjects, denominator)
#f = F ratio statistic
##
##dependent t from t
##the same as single t this formula is d differences
d.single.tt(t = 4.66, n = 90, a = .05, k = 3)
##independent t from t
d.ind.tt(t = 1.71, n1 = 27, n2 = 27, a = .05, k = 2)
##eta, r squared, ICCs
Eta.anova(dfm = 4, dfe = 95, f = 3.232
, a = .05, k = 3)



####run metafor first to test for influence
#only testing for influence here
metaALLmodelFE = rma(yi=Fisher.s.z, vi=Var.of.Z, data = metaALL, method="FE")
metaALLmodelRE = rma(yi=Fisher.s.z, vi=Var.of.Z, data = metaALL, method="PM")

#testing for influence with full dataset
ALLinflueFE =influence(metaALLmodelFE)
ALLinflueRE =influence(metaALLmodelRE)
ALLinflueFE #22,28,33,34,38,48,50,51,52
ALLinflueRE #28,50,52
plot(ALLinflueFE)
plot(ALLinflueRE)

nooutall = metaALL[-c(28,50,52), ] 

####breaking up into more homogenous subgroups to deal with heterogeneity
nooutbetween = subset(nooutall, Design == "Between")
nooutwithin = subset(nooutall, Design == "Within")
##using nooutall nooutbetween and nooutwithin as data with no outliers


####meta for all studies with no outliers
ALLmetamodel = metagen(Fisher.s.z, SE.of.Z, data = nooutall, method.tau="PM")
ALLmetamodel
ALLmetamodel$w.fixed
ALLmetamodel$w.random
#FE ES
z2e(ALLmetamodel$TE.fixed)
#FE CI lo
z2e(ALLmetamodel$lower.fixed)
#FE CI hi
z2e(ALLmetamodel$upper.fixed)
#RE ES
z2e(ALLmetamodel$TE.random)
#RE CI lo
z2e(ALLmetamodel$lower.random)
#RE CI hi
z2e(ALLmetamodel$upper.random)

####meta for between studies with no outliers
BWmetamodel = metagen(Fisher.s.z, SE.of.Z, data = nooutbetween, method.tau="PM")
BWmetamodel
BWmetamodel$w.fixed
BWmetamodel$w.random
#FE ES
z2e(BWmetamodel$TE.fixed)
#FE CI lo
z2e(BWmetamodel$lower.fixed)
#FE CI hi
z2e(BWmetamodel$upper.fixed)
#RE ES
z2e(BWmetamodel$TE.random)
#RE CI lo
z2e(BWmetamodel$lower.random)
#RE CI hi
z2e(BWmetamodel$upper.random)

####meta for within studies with no outliers
WImetamodel = metagen(Fisher.s.z, SE.of.Z, data = nooutwithin, method.tau="PM")
WImetamodel
WImetamodel$w.fixed
WImetamodel$w.random
#FE ES
z2e(WImetamodel$TE.fixed)
#FE CI lo
z2e(WImetamodel$lower.fixed)
#FE CI hi
z2e(WImetamodel$upper.fixed)
#RE ES
z2e(WImetamodel$TE.random)
#RE CI lo
z2e(WImetamodel$lower.random)
#RE CI hi
z2e(WImetamodel$upper.random)


####Deriving individual CIs (same as meta model) 
####and p for each experiment and the number of #studies where p>.05
#CI's based on normal approximation
options(scipen=999)
ci.all<-data.frame(ci(TE = metaALL$Fisher.s.z, seTE = metaALL$SE.of.Z))
ci.all$p
NumNonSig.all<-with(ci.all, c(sum(p>.05)))
NumSig.all<-86-NumNonSig.all
NumSig.all
(NumSig.all/86)*100 ##Percent of studies with significant findings




####p-curve done from p-curve.com  R script is provided on site.
#test statistics to put into application for between studies
F(3, 80) = 4.84
t(48) = -1.4
F(1, 62)=0.21
F(2, 57) = 6.64
F(2, 54) = 2.98
t(71) = 2.09
F(1, 80) = 8.32
F(1, 34) = 6.69
F(1,32)=4.26
F(3,132) = 12.01
F(2,99) = 11.14
F(3,103) = 8.92
F(2, 104) = 4.55
F(1,165) = 13.15
F(2,78) = 3.76
F(1,99)=4.35
F(1,66) = 4.03
F(1, 78) = 3.70
F(1, 108) = 11.68
F(1, 170) = 10.59
F(2, 147) = 6.89
F( 5, 294) = 4.41
F(2,144) = 6.05
F(1, 96) = 9.04
F(1, 48) = 5.22
F(1,68) = 5.84
F(2, 66)= 7.51
F(2,69) = 15.47
F(1, 37) = 9.17
F(2, 72) = 9.57
F(1, 58) = 17.97
F(2, 215) = 4.26
F(1, 98) = 0.13
F(1,76) = 15.96
F(1,70) = 2.83
F(1, 87) = 15.38
F(3, 216) = 4.08
F(3, 76) = 1.11
F(3, 116) = 1.27
F(3, 116) = 7.80
F(1,125) = 5.85
F(1,40) = 0.21
F(2,159) = 3.52
t(52) = 1.71
F(4, 95) = 11.024
F(2, 117) = 1.89
F(2, 177) = 3.045
F(4, 95) = 3.232



#test statistics to put into application for within subjects
F(1, 30) = 17.42
F(2, 138)= 17.2
F(2, 112) = 22.95
F(1, 88) = 5.05
F(2, 58) = 3.17
F(2,70) = 11.69
t(27) = 2.51
t(38) = 2.30
F(2, 92) = 28.4
F(2, 96) = 19.9
t(75)=4.66
F(2,94) = 9.219
t(47) = 2.14
F(2,94) = 6.433
F(1,54) = 34.98
F(2, 60) = 16.18
F(2, 146) = 11.50
F(2, 126) = 13.68
F(1, 79) = 11.85
F(1,118)=45.61
F(1,70) = 37.49
F(1,43) = 12.47
F(1,31) = 4.48
F(1,27)=7.36
F(1, 113) = 7.13
F(1, 27)= 4.72
F(1, 37) = 16.34
F(1, 39) = 8.04
F(1, 49) = 29.88
F(1, 23)= 5.70
t(29) = 1.71
t(29) = 0.47
F(1, 78) = 9.13
F(1,37) = 8.08
F(2, 100) = 4.50
F(2, 142) = 3.62
F(1, 50) = 8.48
F(1, 62) = 4.47
F(3,62) = 15.20

#for all studies, just combine the two above

####TES done in excel #BW .095, WI .171
pwr.f2.test(u = 2, v = 58, f2 = .171, sig.level = 0.05)
pwr.f2.test(u = 2, v = 57, f2 = .095, sig.level = 0.05)
pwr.f2.test(u = 2, v = 60, f2 = .171, sig.level = 0.05)
pwr.f2.test(u = 2, v = 54, f2 = .095, sig.level = 0.05)

pwr.f2.test(u = 2, v = 60, f2 = .171, sig.level = 0.05)
pwr.f2.test(u = 2, v = 104, f2 = .095, sig.level = 0.05)
pwr.f2.test(u = 2, v = 146, f2 = .171, sig.level = 0.05)
pwr.f2.test(u = 2, v = 126, f2 = .171, sig.level = 0.05)

pwr.f2.test(u = 1, v = 27, f2 = .171, sig.level = 0.05)
pwr.f2.test(u = 1, v = 108, f2 = .095, sig.level = 0.05)
pwr.f2.test(u = 1, v = 113, f2 = .171, sig.level = 0.05)
pwr.f2.test(u = 1, v = 170, f2 = .095, sig.level = 0.05)
pwr.f2.test(u = 1, v = 27, f2 = .171, sig.level = 0.05)

pwr.f2.test(u = 2, v = 147, f2 = .095, sig.level = 0.05)
pwr.f2.test(u = 1, v = 37, f2 = .171, sig.level = 0.05)
pwr.f2.test(u = 1, v = 39, f2 = .171, sig.level = 0.05)
pwr.f2.test(u = 1, v = 49, f2 = .171, sig.level = 0.05)

pwr.f2.test(u = 2, v = 66, f2 = .095, sig.level = 0.05)
pwr.f2.test(u = 2, v = 95, f2 = .095, sig.level = 0.05)
pwr.f2.test(u = 2, v = 69, f2 = .095, sig.level = 0.05)
pwr.f2.test(u = 1, v = 37, f2 = .095, sig.level = 0.05)


#####Funnel plot
#This will save as a .png file in the documents folder
png(file="funnel.png", width = 1600, height = 1400, res = 200)
layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = TRUE), respect = TRUE)
funnel(ALLmetamodel, lty.fixed = 1, lty.random = 5, xlim = c(-0.6, 1.1), 
       xlab = "Effect size", ylab = "Standard Error (SE)", cex = .75, 
       col = 1, bg = 1, contour=c(0.00001, 0.95), col.contour=c("snow3", "white"), 
       family = "serif", cex.lab = 1.5, main = "Full Sample", cex.main = 1.8) 
funnel(BWmetamodel, lty.fixed = 1, lty.random = 5, xlim = c(-0.6, 1.1),
       xlab = "Effect size",  ylab = "Standard Error (SE)", cex = .75, 
       col = 1, bg = 1, contour=c(0.00001, 0.95), col.contour=c("snow3", "white"),
       family = "serif", cex.lab = 1.5, main = "Between-Subjects Design", cex.main = 1.8) 
funnel(WImetamodel, lty.fixed = 1, lty.random = 5, xlim = c(-0.6, 1.1),
       xlab = "Effect size",  ylab = "Standard Error (SE)", cex = .75, col = 1, bg = 1, 
       contour=c(0.00001, 0.95), col.contour=c("snow3", "white"), family = "serif", cex.lab = 1.5,
       main = "Within-Subjects Design", cex.main = 1.8) 



#p-uniform ##will have to report raw correlation coefficients and sample size
#no option to put eta.
#ri = correlation coefficients
#ni = sample size
#side = "Right"
#method = "P"
#alpha = 0.05

puniALL = puniform(ri=nooutall$r, ni=nooutall$N,
                         side = "right", method = "P",
                         alpha = 0.05, plot = TRUE)
puniALL
(puniALL$est)^2
(puniALL$ci.lb)^2
(puniALL$ci.ub)^2

puniBW = puniform(ri=nooutbetween$r, ni=nooutbetween$N,
                  side = "right", method = "P",
                  alpha = 0.05, plot = TRUE)
puniBW
(puniBW$est)^2
(puniBW$ci.lb)^2
(puniBW$ci.ub)^2

puniWI = puniform(ri=nooutwithin$r, ni=nooutwithin$N,
                  side = "right", method = "P",
                  alpha = 0.05, plot = TRUE)
puniWI
(puniWI$est)^2
(puniWI$ci.lb)^2
(puniWI$ci.ub)^2


#Trim and fill  
tf.all<-trimfill(ALLmetamodel)
tf.bw<-trimfill(BWmetamodel)
tf.wi<-trimfill(WImetamodel)

tf.all
z2e(tf.all$TE.random)
z2e(tf.all$lower.random)
z2e(tf.all$upper.random)
z2e(tf.all$TE.fixed)
z2e(tf.all$lower.fixed)
z2e(tf.all$upper.fixed)

tf.bw
z2e(tf.bw$TE.random)
z2e(tf.bw$lower.random)
z2e(tf.bw$upper.random)
z2e(tf.bw$TE.fixed)
z2e(tf.bw$lower.fixed)
z2e(tf.bw$upper.fixed)

tf.wi
z2e(tf.wi$TE.random)
z2e(tf.wi$lower.random)
z2e(tf.wi$upper.random)
z2e(tf.wi$TE.fixed)
z2e(tf.wi$lower.fixed)
z2e(tf.wi$upper.fixed)


#Funnel plot after trim and fill
#This will save as a .png file in the documents folder
png(file="funneltrimfill.png", width = 1600, height = 1400, res = 200)
layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = TRUE), respect = TRUE)
funnel(tf.all, lty.fixed = 1, lty.random = 5, xlim = c(-0.6, 1.1), 
       xlab = "Effect size", ylab = "Standard Error (SE)", cex = .75, 
       col = 1, bg = 1, contour=c(0.00001, 0.95), col.contour=c("snow3", "white"), 
       family = "serif", cex.lab = 1.5, main = "Full Sample", cex.main = 1.8) 
funnel(tf.bw, lty.fixed = 1, lty.random = 5, xlim = c(-0.6, 1.1),
       xlab = "Effect size",  ylab = "Standard Error (SE)", cex = .75, 
       col = 1, bg = 1, contour=c(0.00001, 0.95), col.contour=c("snow3", "white"),
       family = "serif", cex.lab = 1.5, main = "Between-Subjects Design", cex.main = 1.8) 
funnel(tf.wi, lty.fixed = 1, lty.random = 5, xlim = c(-0.6, 1.1),
       xlab = "Effect size",  ylab = "Standard Error (SE)", cex = .75, col = 1, bg = 1, 
       contour=c(0.00001, 0.95), col.contour=c("snow3", "white"), family = "serif", cex.lab = 1.5,
       main = "Within-Subjects Design", cex.main = 1.8) 

#PET-PEESE 
#all
PET.all<-lm(Fisher.s.z~SE.of.Z, weights = 1/Var.of.Z, data = nooutall)
summary(PET.all)
confint(PET.all)

PEESE.all<-lm(Fisher.s.z~Var.of.Z, weights = 1/Var.of.Z, data = nooutall)
summary(PEESE.all)
confint(PEESE.all)
z2e(0.27562) #ES
z2e(0.2130469) #lo
z2e(0.3381838) #hi

#bw
PET.BW<-lm(Fisher.s.z~SE.of.Z, weights = 1/Var.of.Z, data = nooutbetween)
summary(PET.BW)
confint(PET.BW)

PEESE.BW<-lm(Fisher.s.z~Var.of.Z, weights = 1/Var.of.Z, data = nooutbetween)
summary(PEESE.BW)
confint(PEESE.BW)
z2e(0.24741) #ES
z2e(0.1700323) #lo
z2e(0.3247855) #hi

#wi
PET.WI<-lm(Fisher.s.z~SE.of.Z, weights = 1/Var.of.Z, data = nooutwithin)
summary(PET.WI)
confint(PET.WI)

PEESE.WI<-lm(Fisher.s.z~Var.of.Z, weights = 1/Var.of.Z, data = nooutwithin)
summary(PEESE.WI)
confint(PEESE.WI)
z2e(0.45592) #ES
z2e(0.32388) #lo
z2e(0.5879555) #hi


####selection models 

#all no out
selectionmodelallFE = weightfunct(effect = nooutall$Fisher.s.z, v = nooutall$Var.of.Z, fe=TRUE)
z2e(selectionmodelallFE$adj_est) #first number
z2e(selectionmodelallFE$ci.lb_adj) #first number
z2e(selectionmodelallFE$ci.ub_adj) #first number
selectionmodelallRE = weightfunct(effect = nooutall$Fisher.s.z, v = nooutall$Var.of.Z, fe=FALSE)
z2e(selectionmodelallRE$adj_est) #second number
z2e(selectionmodelallRE$ci.lb_adj) #second number
z2e(selectionmodelallRE$ci.ub_adj) #second number

#BW no out
selectionmodelbwFE = weightfunct(effect = nooutbetween$Fisher.s.z, v = nooutbetween$Var.of.Z, fe=TRUE)
z2e(selectionmodelbwFE$adj_est) #first number
z2e(selectionmodelbwFE$ci.lb_adj) #first number
z2e(selectionmodelbwFE$ci.ub_adj) #first number
selectionmodelbwRE = weightfunct(effect = nooutbetween$Fisher.s.z, v = nooutbetween$Var.of.Z, fe=FALSE)
z2e(selectionmodelbwRE$adj_est) #second number
z2e(selectionmodelbwRE$ci.lb_adj) #second number
z2e(selectionmodelbwRE$ci.ub_adj) #second number

#WI no out
selectionmodelwiFE = weightfunct(effect = nooutwithin$Fisher.s.z, v = nooutwithin$Var.of.Z, fe=TRUE)
z2e(selectionmodelwiFE$adj_est) #first number
z2e(selectionmodelwiFE$ci.lb_adj) #first number
z2e(selectionmodelwiFE$ci.ub_adj) #first number
selectionmodelwiRE = weightfunct(effect = nooutwithin$Fisher.s.z, v = nooutwithin$Var.of.Z, fe=FALSE)
z2e(selectionmodelwiRE$adj_est) #second number
z2e(selectionmodelwiRE$ci.lb_adj) #second number
z2e(selectionmodelwiRE$ci.ub_adj) #second number


####forest plots. However, these are presented in Fisher's Z, not eta


#BW
bwauthor = as.vector(nooutbetween$Author)
forest(BWmetamodel, studlab = bwauthor, comb.fixed = F, overall = TRUE,
       xlim=c(0.0, 1.00), hetstat = F, fontsize = 11)


#wi
wiauthor = as.vector(nooutwithin$Author)
forest(WImetamodel, studlab = wiauthor, comb.fixed = F, overall = TRUE,
       xlim=c(0.0, 1.00), hetstat = F)



