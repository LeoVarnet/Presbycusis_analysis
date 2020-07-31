# calculate counterfactuals for simple models

if (exists("heardata_counterf_center")){rm(heardata_counterf_center)}
if (exists("heardata_counterf")){rm(heardata_counterf)}

heardata_counterf <- data.frame(
  PTA_pred,
  PTAz_pred,
  rep(age_pred, each=NPTApred),
  rep(agez_pred, each=NPTApred),
  c(rep(1,NPTApred*Nagepred),rep(2,NPTApred*Nagepred),rep(3,NPTApred*Nagepred),rep(4,NPTApred*Nagepred),rep(5,NPTApred*Nagepred),rep(6,NPTApred*Nagepred),rep(7,NPTApred*Nagepred)),
  c(rep(0,NPTApred*Nagepred*7),rep(1,NPTApred*Nagepred*7)))
Ntotalpred = dim(heardata_counterf)[1]
colnames(heardata_counterf) <- c("PTA","PTAz","age", "agez","center","gender")
# heardata_counterf$agefactor <- cut.default(heardata_counterf$age, seq(from=min_age-1,to=max_age+1,length.out=4))
# levels(heardata_counterf$agefactor) <- c(paste("age = ", floor(agefactor_cutoff[1]), "-", floor(agefactor_cutoff[2]),"y"),
#                                      paste("age = ", floor(agefactor_cutoff[2]), "-", floor(agefactor_cutoff[3]),"y"),
#                                      paste("age = ", floor(agefactor_cutoff[3]), "-", floor(agefactor_cutoff[4]),"y"))#paste("age (y) =", levels(heardata$agefactor))

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
  parsfit$beta_age[]*heardata_counterf$agez[i] + 
  parsfit$beta_gender[]*heardata_counterf$gender[i] + 
  parsfit$beta_PTA[]*heardata_counterf$PTAz[i]
logit_pn_pred[,i] = logit_ps_pred[,i] + 
  parsfit$beta_cond[] + 
  parsfit$beta_agecond[]*heardata_counterf$agez[i] + 
  parsfit$beta_condPTA[]*heardata_counterf$PTAz[i] +
  parsfit$beta_agecondPTA[]*heardata_counterf$agez[i]*heardata_counterf$PTAz[i] 
logit_ps_pred_center[,i] = 
  parsfit$gamma_0[,heardata_counterf$center[i]] +
  parsfit$beta_gender[]*heardata_counterf$gender[i] + 
  parsfit$beta_age[]*heardata_counterf$agez[i] + 
  parsfit$beta_PTA[]*heardata_counterf$PTAz[i]
logit_pn_pred_center[,i] = logit_ps_pred_center[,i] +
  parsfit$beta_cond[] + 
  parsfit$beta_agecond[]*heardata_counterf$agez[i] + 
  parsfit$beta_condPTA[]*heardata_counterf$PTAz[i] +
  parsfit$beta_agecondPTA[]*heardata_counterf$agez[i]*heardata_counterf$PTAz[i] 

ps_pred[,i] = 1/16+(1-1/16)*inv.logit(logit_ps_pred[,i])
#pn_pred[,i] = 1/16+(1-apply(parsfit$plapse[,],1,mean)-1/16)*inv.logit(logit_pn_pred[,i])
#pn_pred[,i] = 1/16+(1-1/16-parsfit$plapse[,1])*inv.logit(logit_pn_pred[,i])
pn_pred[,i] = 1/16+(1-1/16-parsfit$plapse)*inv.logit(logit_pn_pred[,i])
ps_pred_center[,i] = 1/16+(1-1/16)*inv.logit(logit_ps_pred_center[,i])
#pn_pred_center[,i] = 1/16+(1-1/16-parsfit$plapse[,heardata_counterf$center[i]])*inv.logit(logit_pn_pred_center[,i])
#pn_pred_center[,i] = 1/16+(1-1/16-parsfit$plapse[,1])*inv.logit(logit_pn_pred_center[,i])
pn_pred_center[,i] = 1/16+(1-1/16-parsfit$plapse)*inv.logit(logit_pn_pred_center[,i])
}  

heardata_counterf$agefactor <- cut.default(heardata_counterf$age, agefactor_cutoff)
levels(heardata_counterf$agefactor) <- c(paste("age = ", floor(agefactor_cutoff[1]), "-", floor(agefactor_cutoff[2]),"y"),
                                     paste("age = ", floor(agefactor_cutoff[2]), "-", floor(agefactor_cutoff[3]),"y"),
                                     paste("age = ", floor(agefactor_cutoff[3]), "-", floor(agefactor_cutoff[4]),"y"))#paste("age (y) =", levels(heardata$agefactor))
heardata_counterf$center <- as.factor(heardata_counterf$center)
levels(heardata_counterf$center)[levels(heardata_counterf$center)=="1"] <- "Marseille"
levels(heardata_counterf$center)[levels(heardata_counterf$center)=="2"] <- "Lille"
levels(heardata_counterf$center)[levels(heardata_counterf$center)=="3"] <- "Paris"
levels(heardata_counterf$center)[levels(heardata_counterf$center)=="4"] <- "Clermont"
levels(heardata_counterf$center)[levels(heardata_counterf$center)=="5"] <- "Lyon"
levels(heardata_counterf$center)[levels(heardata_counterf$center)=="6"] <- "Bordeaux"
levels(heardata_counterf$center)[levels(heardata_counterf$center)=="7"] <- "Toulouse"

heardata_counterf$PC_silence<-t(apply(100*ps_pred[],2,quantile,probs=c(0.025,0.5,0.975),na.rm = TRUE)) #the median line with 95% credible intervals
colnames(heardata_counterf$PC_silence)<-c("PC_silence_Cinf","PC_silence_med","PC_silence_Csup")
heardata_counterf$PC_noise<-t(apply(100*pn_pred[],2,quantile,probs=c(0.025,0.5,0.975),na.rm = TRUE)) #the median line with 95% credible intervals
colnames(heardata_counterf$PC_noise)<-c("PC_noise_Cinf","PC_noise_med","PC_noise_Csup")

heardata_counterf_center = heardata_counterf
heardata_counterf_center$PC_silence<-t(apply(100*ps_pred_center[],2,quantile,probs=c(0.025,0.5,0.975),na.rm = TRUE)) #the median line with 95% credible intervals
colnames(heardata_counterf_center$PC_silence)<-c("PC_silence_Cinf","PC_silence_med","PC_silence_Csup")
heardata_counterf$PC_silence<-t(apply(100*ps_pred[],2,quantile,probs=c(0.025,0.5,0.975),na.rm = TRUE)) #the median line with 95% credible intervals
colnames(heardata_counterf$PC_silence)<-c("PC_silence_Cinf","PC_silence_med","PC_silence_Csup")
heardata_counterf_center$PC_noise<-t(apply(100*pn_pred_center[],2,quantile,probs=c(0.025,0.5,0.975),na.rm = TRUE)) #the median line with 95% credible intervals
colnames(heardata_counterf_center$PC_noise)<-c("PC_noise_Cinf","PC_noise_med","PC_noise_Csup")
heardata_counterf$PC_noise<-t(apply(100*pn_pred[],2,quantile,probs=c(0.025,0.5,0.975),na.rm = TRUE)) #the median line with 95% credible intervals
colnames(heardata_counterf$PC_noise)<-c("PC_noise_Cinf","PC_noise_med","PC_noise_Csup")

heardata_counterf_center=aggregate(.~agefactor*PTA*center,data=heardata_counterf_center,mean)
#heardata_counterf=aggregate(.~agefactor*PTA*center,data=heardata_counterf,mean)
heardata_counterf=aggregate(.~agefactor*PTA,data=heardata_counterf,mean)
