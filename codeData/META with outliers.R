
library(meta)
library(metafor)
####load data
metaALL <- read.csv("~/John M.S. Experimental Psychology/Research and Articles/Evolutionary Psychology/SPA Meta-Analysis/Data analysis/data analysis round 3/metaALL.csv")

####function to convert z back to eta
z2e = function(z){
  corrFE=(exp(z/0.5)-1)/(exp(z/0.5)+1)
  etaFE=corrFE^2
  etaFE
}

yesout = metaALL

####breaking up into more homogenous subgroups to deal with heterogeneity
yesoutbetween = subset(metaALL, Design == "Between")
yesoutwithin = subset(metaALL, Design == "Within")



####meta for all studies
ALLmetamodel.x = metagen(Fisher.s.z, SE.of.Z, data = metaALL, method.tau="PM")
ALLmetamodel.x
ALLmetamodel.x$w.fixed
ALLmetamodel.x$w.random
#FE ES
z2e(ALLmetamodel.x$TE.fixed)
#FE CI lo
z2e(ALLmetamodel.x$lower.fixed)
#FE CI hi
z2e(ALLmetamodel.x$upper.fixed)
#RE ES
z2e(ALLmetamodel.x$TE.random)
#RE CI lo
z2e(ALLmetamodel.x$lower.random)
#RE CI hi
z2e(ALLmetamodel.x$upper.random)

####meta for between studies
BWmetamodel.x = metagen(Fisher.s.z, SE.of.Z, data = yesoutbetween, method.tau="PM")
BWmetamodel.x
BWmetamodel.x$w.fixed
BWmetamodel.x$w.random
#FE ES
z2e(BWmetamodel.x$TE.fixed)
#FE CI lo
z2e(BWmetamodel.x$lower.fixed)
#FE CI hi
z2e(BWmetamodel.x$upper.fixed)
#RE ES
z2e(BWmetamodel.x$TE.random)
#RE CI lo
z2e(BWmetamodel.x$lower.random)
#RE CI hi
z2e(BWmetamodel.x$upper.random)

####meta for within studies
WImetamodel.x = metagen(Fisher.s.z, SE.of.Z, data = yesoutwithin, method.tau="PM")
WImetamodel.x
WImetamodel.x$w.fixed
WImetamodel.x$w.random
#FE ES
z2e(WImetamodel.x$TE.fixed)
#FE CI lo
z2e(WImetamodel.x$lower.fixed)
#FE CI hi
z2e(WImetamodel.x$upper.fixed)
#RE ES
z2e(WImetamodel.x$TE.random)
#RE CI lo
z2e(WImetamodel.x$lower.random)
#RE CI hi
z2e(WImetamodel.x$upper.random)


#Trim and fill  
tf.all.x<-trimfill(ALLmetamodel.x)
tf.bw.x<-trimfill(BWmetamodel.x)
tf.wi.x<-trimfill(WImetamodel.x)

tf.all.x
z2e(tf.all.x$TE.random)
z2e(tf.all.x$lower.random)
z2e(tf.all.x$upper.random)
z2e(tf.all.x$TE.fixed)
z2e(tf.all.x$lower.fixed)
z2e(tf.all.x$upper.fixed)

tf.bw.x
z2e(tf.bw.x$TE.random)
z2e(tf.bw.x$lower.random)
z2e(tf.bw.x$upper.random)
z2e(tf.bw.x$TE.fixed)
z2e(tf.bw.x$lower.fixed)
z2e(tf.bw.x$upper.fixed)

tf.wi.x
z2e(tf.wi.x$TE.random)
z2e(tf.wi.x$lower.random)
z2e(tf.wi.x$upper.random)
z2e(tf.wi.x$TE.fixed)
z2e(tf.wi.x$lower.fixed)
z2e(tf.wi.x$upper.fixed)


#PET-PEESE 
#all
PET.all.x<-lm(Fisher.s.z~SE.of.Z, weights = 1/Var.of.Z, data = metaALL)
summary(PET.all.x)
confint(PET.all.x)

PEESE.all.x<-lm(Fisher.s.z~Var.of.Z, weights = 1/Var.of.Z, data = metaALL)
summary(PEESE.all.x)
confint(PEESE.all.x)
z2e(0.29969) #ES
z2e(0.2248061) #lo
z2e(0.3745742) #hi

#bw
PET.BW.x<-lm(Fisher.s.z~SE.of.Z, weights = 1/Var.of.Z, data = yesoutbetween)
summary(PET.BW.x)
confint(PET.BW.x)

PEESE.BW.x<-lm(Fisher.s.z~Var.of.Z, weights = 1/Var.of.Z, data = yesoutbetween)
summary(PEESE.BW.x)
confint(PEESE.BW.x)
z2e(0.25618) #ES
z2e(0.1654194) #lo
z2e(0.3469461) #hi

#wi
PET.WI.x<-lm(Fisher.s.z~SE.of.Z, weights = 1/Var.of.Z, data = yesoutwithin)
summary(PET.WI.x)
confint(PET.WI.x)

PEESE.WI.x<-lm(Fisher.s.z~Var.of.Z, weights = 1/Var.of.Z, data = yesoutwithin)
summary(PEESE.WI.x)
confint(PEESE.WI.x)
z2e(0.53187) #ES
z2e(0.3794364) #lo
z2e(0.6843123) #hi
