# calculate counterfactuals for simple models

rm(heardata_pred_center)
rm(heardata_pred)

heardata_pred <- data.frame(
  PTA_pred,
  PTAz_pred,
  rep(age_pred, each=NPTApred),
  rep(agez_pred, each=NPTApred),
  c(rep(1,NPTApred*Nagepred),rep(2,NPTApred*Nagepred),rep(3,NPTApred*Nagepred),rep(4,NPTApred*Nagepred),rep(5,NPTApred*Nagepred),rep(6,NPTApred*Nagepred),rep(7,NPTApred*Nagepred)),
  c(rep(0,NPTApred*Nagepred*7),rep(1,NPTApred*Nagepred*7)))
Ntotalpred = dim(heardata_pred)[1]
colnames(heardata_pred) <- c("PTA","PTAz","age", "agez","center","gender")
heardata_pred$agefactor <- cut.default(heardata_pred$age, seq(from=min_age-1,to=max_age+1,length.out=4))

ps_pred = matrix(nrow = Nsamples, ncol = Ntotalpred)
ps_pred_center = matrix(nrow = Nsamples, ncol = Ntotalpred)
pn_pred = matrix(nrow = Nsamples, ncol = Ntotalpred)
pn_pred_center = matrix(nrow = Nsamples, ncol = Ntotalpred)
logit_ps_pred = matrix(nrow = Nsamples, ncol = Ntotalpred)
logit_ps_pred_center = matrix(nrow = Nsamples, ncol = Ntotalpred)
logit_pn_pred = matrix(nrow = Nsamples, ncol = Ntotalpred)
logit_pn_pred_center = matrix(nrow = Nsamples, ncol = Ntotalpred)

for (i in 1:Ntotalpred)
{logit_ps_pred[,i] = parsfit$beta_0[] +
  parsfit$beta_age[]*heardata_pred$agez[i] + 
  parsfit$beta_gender[]*heardata_pred$gender[i] + 
  parsfit$beta_PTA[]*heardata_pred$PTAz[i]
logit_pn_pred[,i] = logit_ps_pred[,i] + 
  parsfit$beta_cond[]
logit_ps_pred_center[,i] = 
  parsfit$gamma_0[,heardata_pred$center[i]] +
  parsfit$beta_gender[]*heardata_pred$gender[i] + 
  parsfit$beta_age[]*heardata_pred$agez[i] + 
  parsfit$beta_PTA[]*heardata_pred$PTAz[i]
logit_pn_pred_center[,i] = logit_ps_pred_center[,i] +
  parsfit$beta_cond[]

ps_pred[,i] = 1/16+(1-1/16)*inv.logit(logit_ps_pred[,i])
pn_pred[,i] = 1/16+(1-apply(parsfit$plapse[,],1,mean)-1/16)*inv.logit(logit_pn_pred[,i])
ps_pred_center[,i] = 1/16+(1-1/16)*inv.logit(logit_ps_pred_center[,i])
pn_pred_center[,i] = 1/16+(1-1/16-parsfit$plapse[,heardata_pred$center[i]])*inv.logit(logit_pn_pred_center[,i])
}  

heardata_pred$agefactor <- cut.default(heardata_pred$age, seq(from=min_age-1,to=max_age+1,length.out=4))
heardata_pred$center <- as.factor(heardata_pred$center)
levels(heardata_pred$center)[levels(heardata_pred$center)=="1"] <- "Marseille"
levels(heardata_pred$center)[levels(heardata_pred$center)=="2"] <- "Lille"
levels(heardata_pred$center)[levels(heardata_pred$center)=="3"] <- "Paris"
levels(heardata_pred$center)[levels(heardata_pred$center)=="4"] <- "Clermont"
levels(heardata_pred$center)[levels(heardata_pred$center)=="5"] <- "Lyon"
levels(heardata_pred$center)[levels(heardata_pred$center)=="6"] <- "Bordeaux"
levels(heardata_pred$center)[levels(heardata_pred$center)=="7"] <- "Toulouse"

heardata_pred$PC_silence<-t(apply(100*ps_pred[],2,quantile,probs=c(0.025,0.5,0.975),na.rm = TRUE)) #the median line with 95% credible intervals
colnames(heardata_pred$PC_silence)<-c("PC_silence_Cinf","PC_silence_med","PC_silence_Csup")
heardata_pred$PC_noise<-t(apply(100*pn_pred[],2,quantile,probs=c(0.025,0.5,0.975),na.rm = TRUE)) #the median line with 95% credible intervals
colnames(heardata_pred$PC_noise)<-c("PC_noise_Cinf","PC_noise_med","PC_noise_Csup")

heardata_pred_center = heardata_pred
heardata_pred_center$PC_silence<-t(apply(100*ps_pred_center[],2,quantile,probs=c(0.025,0.5,0.975),na.rm = TRUE)) #the median line with 95% credible intervals
colnames(heardata_pred_center$PC_silence)<-c("PC_silence_Cinf","PC_silence_med","PC_silence_Csup")
heardata_pred$PC_silence<-t(apply(100*ps_pred[],2,quantile,probs=c(0.025,0.5,0.975),na.rm = TRUE)) #the median line with 95% credible intervals
colnames(heardata_pred$PC_silence)<-c("PC_silence_Cinf","PC_silence_med","PC_silence_Csup")
heardata_pred_center$PC_noise<-t(apply(100*pn_pred_center[],2,quantile,probs=c(0.025,0.5,0.975),na.rm = TRUE)) #the median line with 95% credible intervals
colnames(heardata_pred_center$PC_noise)<-c("PC_noise_Cinf","PC_noise_med","PC_noise_Csup")
heardata_pred$PC_noise<-t(apply(100*pn_pred[],2,quantile,probs=c(0.025,0.5,0.975),na.rm = TRUE)) #the median line with 95% credible intervals
colnames(heardata_pred$PC_noise)<-c("PC_noise_Cinf","PC_noise_med","PC_noise_Csup")

heardata_pred_center=aggregate(.~agefactor*PTA*center,data=heardata_pred_center,mean)
#heardata_pred=aggregate(.~agefactor*PTA*center,data=heardata_pred,mean)
heardata_pred=aggregate(.~agefactor*PTA,data=heardata_pred,mean)
