# Setup -------------------------------------------------------------------

library(ggplot2) 
library(ggExtra)
library(rstan)
library(shinystan)
#library(tidyverse)
library(cowplot)
library(dagitty)
library(loo)
library(gtools)

rm(list=ls())
setwd("C:/Users/Léo/ownCloud/Professionnel/Projet Genopath/011_Presbyacousie_Lorenzi/R_analyses/Presbycusis_analysis")
theme_set(theme_bw())
heardata = read.table('mat_full.txt')
colnames(heardata) <- c("group", "age", "PTA", "subject", "center", "gender", "NC_silence", "NC_noise")

Ntrials = max(heardata$NC_silence); # 48 trials per condition
Ncenter = max((heardata$center)); # 7 centers
heardata$group = factor(heardata$group)
heardata$subject = factor(heardata$subject)
heardata$PTA = as.numeric(heardata$PTA)
heardata$age = as.numeric(heardata$age)
heardata$center = factor(heardata$center)
heardata$gender = factor(heardata$gender)
heardata$NC_silence = as.numeric(heardata$NC_silence)
heardata$NC_noise = as.numeric(heardata$NC_noise)
heardata$PC_silence = 100*heardata$NC_silence/Ntrials
heardata$PC_noise = 100*heardata$NC_noise/Ntrials
heardata$lPC_silence <- qlogis(heardata$PC_silence/100)
heardata$lPC_noise <- qlogis(heardata$PC_noise/100)
min_age = min(heardata$age) # subjects 42- to 92-yo
max_age = max(heardata$age)
levels(heardata$group)[levels(heardata$group)=="0"] <- "NH"
levels(heardata$group)[levels(heardata$group)=="1"] <- "HI"
levels(heardata$gender)[levels(heardata$gender)=="1"] <- "Man"
levels(heardata$gender)[levels(heardata$gender)=="2"] <- "Woman"
levels(heardata$center)[levels(heardata$center)=="1"] <- "Marseille"
levels(heardata$center)[levels(heardata$center)=="2"] <- "Lille"
levels(heardata$center)[levels(heardata$center)=="3"] <- "Paris"
levels(heardata$center)[levels(heardata$center)=="4"] <- "Clermont"
levels(heardata$center)[levels(heardata$center)=="5"] <- "Lyon"
levels(heardata$center)[levels(heardata$center)=="6"] <- "Bordeaux"
levels(heardata$center)[levels(heardata$center)=="7"] <- "Toulouse"


# transform continuous variable age into a 3-level factor for representation purpose
heardata$agefactor <- cut.default(heardata$age, seq(from=min_age-1,to=max_age+1,length.out=4))

heardata2 <- data.frame(  unlist(list(heardata$group,heardata$group)),
                          c(heardata$age,heardata$age),
                          c(heardata$PTA,heardata$PTA),
                          unlist(list(heardata$subject,heardata$subject)),
                          unlist(list(heardata$center,heardata$center)),
                          unlist(list(heardata$gender,heardata$gender)),
                          c(heardata$NC_silence,heardata$NC_noise),
                          c(heardata$PC_silence,heardata$PC_noise),
                          c(rep(1,length(heardata$NC_silence)),rep(2,length(heardata$NC_noise))))
colnames(heardata2) <- c("group", "age", "PTA", "subject", "center", "gender", "NC", "PC", 'condition')
heardata2$condition <- factor(heardata2$condition)
levels(heardata2$condition)[levels(heardata2$condition)=="1"] <- "Silence"
levels(heardata2$condition)[levels(heardata2$condition)=="2"] <- "Noise"
heardata2$agefactor <- cut.default(heardata2$age, seq(from=min_age-1,to=max_age+1,length.out=4))
heardata2$lPC <- qlogis(heardata2$PC/100)

# Zscored data and counterfactual predictors

PTA_pred = seq(from = 0, to = 70, by = 7)
NPTApred = length(PTA_pred)
heardata$PTAz = (heardata$PTA-mean(heardata$PTA))/sd(heardata$PTA)
heardata2$PTAz = (heardata2$PTA-mean(heardata2$PTA))/sd(heardata2$PTA)
PTAz_pred = (PTA_pred-mean(heardata2$PTA))/sd(heardata2$PTA)

age_pred = seq(from = 40, to = 100, by = 6)
Nagepred = length(age_pred)
heardata$agez = (heardata$age-mean(heardata$age))/sd(heardata$age)
heardata2$agez = (heardata2$age-mean(heardata2$age))/sd(heardata2$age)
agez_pred = (age_pred-mean(heardata2$age))/sd(heardata2$age)

group_pred = c(0,1)
cond_pred = c(0,1)

# Data as a list for feeding Rstan models

data <- list( agez = heardata$agez[heardata$group=="HI"],
              PTAz = heardata$PTAz[heardata$group=="HI"],
              agez = heardata$agez[heardata$group=="HI"],
              PC_silence = heardata$PC_silence[heardata$group=="HI"],
              PC_noise = heardata$PC_noise[heardata$group=="HI"],
              lPC_silence = heardata$lPC_silence[heardata$group=="HI"],
              lPC_noise = heardata$lPC_noise[heardata$group=="HI"],
              gender = as.numeric(heardata$gender[heardata$group=="HI"])-1,
              N = length(heardata$PTAz[heardata$group=="HI"]),
              NCs = heardata$NC_silence[heardata$group=="HI"],
              NCn = heardata$NC_noise[heardata$group=="HI"],
              N_NH = length(heardata$NC_noise[heardata$group=="NH"]),
              NCn_NH = heardata$NC_noise[heardata$group=="NH"],
              center = as.numeric(heardata$center[heardata$group=="HI"]),
              Ntrials = Ntrials,
              Ncenter = Ncenter,
              prior_only = 0)

# Functions

Rsquared <- function(data,parsfit){
  SSerr = ((data$lPC_silence-apply(parsfit$eta_s,2,mean))^2 + (data$lPC_noise-apply(parsfit$eta_n,2,mean))^2)/2
  SSerr = sum(SSerr[SSerr!=Inf])
  meanlPC = (mean(data$lPC_silence[data$lPC_silence!=Inf])+mean(data$lPC_noise[data$lPC_noise!=Inf]))/2
  SStot = ((data$lPC_silence-meanlPC)^2 + (data$lPC_noise-meanlPC)^2)/2
  SStot = sum(SStot[SStot!=Inf])
  Rsquared = 1-SSerr/SStot
}


# 1. Plot PTA vs. age -----------------------------------------------------

p1 <- ggplot() +
  geom_count(data = heardata, aes(x=age, y=PTA, color=group)) +
  scale_size_area(max_size = 2,name = "N",n.breaks = 3) +
  labs(y="PTA (dB HL)", x = "age (years)") +
  xlim(40, 95) + ylim(0, 70)

p1xhist <- axis_canvas(p1, axis = "x") + 
  geom_histogram(data = heardata,
                 aes(x = age, color = group, fill = group, alpha=0.1),
                 binwidth = 2, position="identity")+
  xlim(40, 95)
p1xSD <- axis_canvas(p1xhist, axis = "x") + 
  geom_boxplot(data = heardata, aes(x = age,color = group), outlier.shape = NA)+
  xlim(40, 95)
p1yhist <- axis_canvas(p1, axis = "y", coord_flip = TRUE) + 
  geom_histogram(data = heardata,
                 aes(x = PTA, color = group, fill = group, alpha=0.1),
                 binwidth = 2, position="identity") +
  coord_flip() + xlim(0, 70)
p1ySD <- axis_canvas(p1yhist, axis = "y", coord_flip = TRUE) + 
  geom_boxplot(data = heardata, aes(x = PTA, color = group), outlier.shape = NA) +
  coord_flip() + xlim(0, 70)
p1xhist <- insert_xaxis_grob(p1xhist, p1xSD, grid::unit(.2, "null"), position = "top")
p1 <- insert_xaxis_grob(p1, p1xhist, grid::unit(.2, "null"), position = "top")
p1yhist <- insert_yaxis_grob(p1yhist, p1ySD, grid::unit(.2, "null"), position = "right")
p1 <- insert_yaxis_grob(p1, p1yhist, grid::unit(.2, "null"), position = "right")
ggdraw(p1)

# 2. Plot PTA vs. PC --------------------

bin2d = c(25,25)
p2.1 <- ggplot(data = heardata2,aes(x=PTA, y=PC)) +
  geom_bin2d(bins = bin2d) +
  theme_bw() + 
  #xlim(0,70) + ylim(0,100) + 
  xlab('PTA (dB HL)') + ylab('% correct') +
  scale_fill_continuous(type = "viridis",limits=NULL) + facet_grid(condition~group)
p2.1
p2.1a <- ggplot(data = heardata2,aes(x=age, y=PC)) +
  geom_bin2d(bins = bin2d) +
  theme_bw() + 
  #xlim(0,70) + ylim(0,100) + 
  xlab('age (y)') + ylab('% correct') +
  scale_fill_continuous(type = "viridis",limits=NULL) + facet_grid(condition~group)
p2.1a

p2.2 <- ggplot(data = heardata2,aes(x=PTA, y=PC, color=age)) +
  geom_point() +
  theme_bw() + 
  #xlim(0,70) + ylim(0,100) + 
  xlab('PTA (dB HL)') + ylab('% correct') +
  scale_fill_continuous(type = "viridis",limits=NULL) + facet_grid(condition~group)
p2.2
p2.2a <- ggplot(data = heardata2,aes(x=age, y=PC, color=PTA)) +
  geom_point() +
  theme_bw() + 
  #xlim(0,70) + ylim(0,100) + 
  xlab('age (y)') + ylab('% correct') +
  scale_fill_continuous(type = "viridis",limits=NULL) + facet_grid(condition~group)
p2.2a

p2.3 <- ggplot(data = heardata2,aes(x=PTA, y=PC, color=center)) +
  geom_point() +
  theme_bw() + 
  #xlim(0,70) +  + 
  xlab('PTA (dB HL)') + ylab('% correct') +
  scale_fill_continuous(type = "viridis",limits=NULL) + facet_grid(condition~group)
p2.3
p2.3a <- ggplot(data = heardata2,aes(x=age, y=PC, color=center)) +
  geom_point() +
  theme_bw() + 
  #xlim(0,70) + ylim(0,100) + 
  xlab('age (y)') + ylab('% correct') +
  scale_fill_continuous(type = "viridis",limits=NULL) + facet_grid(condition~group)
p2.3a

p2.4 <- ggplot(data = heardata2,aes(x=PTA, y=age, z=PC)) + 
  stat_summary_2d() + 
  facet_grid(~condition)
p2.4

p2.5 <- ggplot(data = subset(heardata2,group=="HI"),aes(x=PTA, y=PC, color=group)) + 
  geom_count() +
  #scale_size(name = "N")+
  scale_size_area(max_size = 2,name = "N",n.breaks = 5)+
#geom_point() +
  facet_grid(condition~agefactor)
p2.5

p2.5a <- ggplot(data = subset(heardata2,group=="HI"),aes(x=PTA, y=PC, color=center)) + 
  geom_count() +
  #scale_size(name = "N")+
  scale_size_area(max_size = 2,name = "N",n.breaks = 5)+
  #geom_point() +
  facet_grid(condition~agefactor)
p2.5a

p2.6 <- ggplot(data = heardata2,aes(x=PTA, y=lPC, color=group)) + 
  geom_count() +
  #scale_size(name = "N")+
  scale_size_area(max_size = 2,name = "N",n.breaks = 5)+
  #geom_point() +
  facet_grid(condition~agefactor)
p2.6

p2.6a <- ggplot(data = heardata2,aes(x=PTA, y=lPC, color=center)) + 
  geom_count() +
  #scale_size(name = "N")+
  scale_size_area(max_size = 2,name = "N",n.breaks = 5)+
  #geom_point() +
  facet_grid(condition~agefactor)
p2.6a

## 3. Plot masking effect --------------

p3.0a <- ggplot(data = subset(heardata, group=="HI"),aes(x=PC_silence, y=PC_noise, color=center)) + 
  geom_point() +
  scale_size_area(max_size = 2,name = "N",n.breaks = 5) +
  geom_abline(intercept = 0, slope = 1,linetype="dashed", color = "black")+
  xlim(0, 100) + ylim(0, 100)
p3.0a

p3.0b <- ggplot(data = subset(heardata, group=="HI"),aes(x=lPC_silence, y=lPC_noise, color=center)) + 
  geom_count() +
  scale_size_area(max_size = 2,name = "N",n.breaks = 5) +
  geom_abline(intercept = 0, slope = 1,linetype="dashed", color = "black")
p3.0b

p3.1a <- ggplot(data = subset(heardata,group=="HI")) + 
  geom_point(aes(x=PTA, y=lPC_silence, color=center)) +
  geom_point(aes(x=PTA, y=lPC_noise, color=center), shape=1) +
  #geom_segment(aes(x=PTA, xend=PTA, y=lPC_noise, yend=lPC_silence, color=center)) +
  #scale_size(name = "N")+
  ylab('logit(PC)')+
  scale_size_area(max_size = 2,name = "N",n.breaks = 5) +
  facet_grid(~agefactor)
p3.1a

p3.1b <- ggplot(data = subset(heardata,group=="HI")) + 
  geom_point(aes(x=PTA, y=PC_silence, color=center)) +
  geom_point(aes(x=PTA, y=PC_noise, color=center), shape=1) +
  #geom_segment(aes(x=PTA, xend=PTA, y=lPC_noise, yend=lPC_silence, color=center)) +
  #scale_size(name = "N")+
  ylab('PC')+
  scale_size_area(max_size = 2,name = "N",n.breaks = 5) +
  facet_grid(~agefactor)
p3.1b

p3.1c <- ggplot(data = subset(heardata,group=="HI")) + 
  geom_point(aes(x=PTA, y=PC_silence, color=center)) +
  geom_point(aes(x=PTA, y=PC_noise, color=center), shape=1) +
  #geom_segment(aes(x=PTA, xend=PTA, y=lPC_noise, yend=lPC_silence, color=center)) +
  #scale_size(name = "N")+
  ylab('PC')+
  scale_size_area(max_size = 2,name = "N",n.breaks = 5) +
  facet_grid(center~agefactor)+ 
  theme(legend.position = "none")
p3.1c

## 4. Simple model (w only main effects) -----------------

m1.1 <- stan_model(file = 'm1.1.stan')

fit.m1.1 <- sampling(m1.1,
                       data = data,
                       chains = 4,             # number of Markov chains
                       warmup = 3000,          # number of warmup iterations per chain
                       iter = 7000,            # total number of iterations per chain
                       refresh = 1000)
# diagnosis

parameters = c("beta_0","beta_PTA","beta_age","beta_cond","beta_gender","plapse")#,"beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA","plapse")#c("beta_0","beta_PTA","beta_age","beta_cond","beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA")
print(fit.m1.1, pars = parameters)
#launch_shinystan(fit.m1.1)
plot(fit.m1.1, show_density = TRUE, show_outer_line = FALSE, ci_level= 0.95,outer_level= 0.99, pars = parameters)+ ggtitle("m1.1") #+ coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlim(-3,1)

# counterfactual predictions

parsfit<-extract(fit.m1.1,pars=rev(fit.m1.1@model_pars))#c("beta_0","beta_PTA","beta_age","beta_cond","beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA","gamma_0","gamma_PTA","gamma_age","gamma_cond","gamma_condPTA","gamma_agePTA","gamma_agecond","gamma_agecondPTA"))
Nsamples = length(parsfit$beta_0)

parsfit$gamma_0 <- matrix(0, Nsamples, 7)
parsfit$plapse <- matrix(parsfit$plapse,nrow=Nsamples,ncol=7,byrow=TRUE)

source("counterfactuals_simple.R")
 
p3.1b +
  geom_line(data=heardata_pred,aes(x=PTA, y=PC_silence_med), inherit.aes = FALSE, color='black', 'size'=1)+
  geom_line(data=heardata_pred,aes(x=PTA, y=PC_noise_med), inherit.aes = FALSE, color='black', 'size'=1, linetype='longdash')+
  geom_ribbon(data=heardata_pred,aes(x=PTA, ymin=PC_silence_Cinf, ymax=PC_silence_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='black')+
  geom_ribbon(data=heardata_pred,aes(x=PTA, ymin=PC_noise_Cinf, ymax=PC_noise_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='black')
ylim(0, 100)

# Compute goodness of fit measures

Rsquared_m1.1 <- Rsquared(data,parsfit)
log_lik_1.1 <- extract_log_lik(fit.m1.1)
loo_1.1 <- loo(log_lik_1.1)
print(loo_1.1)

## 5. Hier. GLM w only main effects -----------------

m1.1hi <- stan_model(file = 'm1.1hi_bis.stan')

fit.m1.1hi <- sampling(m1.1hi,
                     data = data,
                     chains = 4,             # number of Markov chains
                     warmup = 3000,          # number of warmup iterations per chain
                     iter = 7000,            # total number of iterations per chain
                     refresh = 1000)
# diagnosis

parameters = c("beta_0","gamma_0","beta_PTA","beta_age","beta_cond","beta_gender","plapse")#,"beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA","plapse")#c("beta_0","beta_PTA","beta_age","beta_cond","beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA")
print(fit.m1.1hi, pars = parameters)
launch_shinystan(fit.m1.1hi)
plot(fit.m1.1hi, show_density = TRUE, show_outer_line = FALSE, ci_level= 0.95,outer_level= 0.99, pars = parameters) + ggtitle("m1.1hi") #+ coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1))#+xlim(-3,1)

# counterfactual predictions

parsfit<-extract(fit.m1.1hi,pars=rev(fit.m1.1hi@model_pars))#c("beta_0","beta_PTA","beta_age","beta_cond","beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA","gamma_0","gamma_PTA","gamma_age","gamma_cond","gamma_condPTA","gamma_agePTA","gamma_agecond","gamma_agecondPTA"))
Nsamples = length(parsfit$beta_0)

source("counterfactuals_simple.R")

p3.1c +
  geom_line(data=heardata_pred_center,aes(x=PTA, y=PC_silence_med, color=center), inherit.aes = FALSE,  'size'=1)+
  geom_line(data=heardata_pred_center,aes(x=PTA, y=PC_noise_med, color=center), inherit.aes = FALSE, linetype='longdash', 'size'=1)+
  geom_ribbon(data=heardata_pred_center,aes(x=PTA, ymin=PC_silence_Cinf, ymax=PC_silence_Csup, fill=center), inherit.aes = FALSE, color=NA, alpha=0.1)+
  geom_ribbon(data=heardata_pred_center,aes(x=PTA, ymin=PC_noise_Cinf, ymax=PC_noise_Csup,  fill=center), inherit.aes = FALSE, color=NA, alpha=0.1)+
  ylim(0, 100)

p3.1b +
  geom_line(data=heardata_pred_center,aes(x=PTA, y=PC_silence_med, color=center), inherit.aes = FALSE, 'size'=1)+
  geom_line(data=heardata_pred_center,aes(x=PTA, y=PC_noise_med, color=center), inherit.aes = FALSE, 'size'=1, linetype='longdash')

p3.1b +
  geom_line(data=heardata_pred,aes(x=PTA, y=PC_silence_med), inherit.aes = FALSE, color='blue', 'size'=1)+
  geom_line(data=heardata_pred,aes(x=PTA, y=PC_noise_med), inherit.aes = FALSE, color='blue', linetype='dotted', 'size'=1)+
  geom_ribbon(data=heardata_pred,aes(x=PTA, ymin=PC_silence_Cinf, ymax=PC_silence_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='blue')+
  geom_ribbon(data=heardata_pred,aes(x=PTA, ymin=PC_noise_Cinf, ymax=PC_noise_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='blue')

p3.1a +
  geom_line(data=heardata_pred,aes(x=PTA, y=logit(PC_silence_med/100)), inherit.aes = FALSE, color='blue', 'size'=1)+
  geom_line(data=heardata_pred,aes(x=PTA, y=logit(PC_noise_med/100)), inherit.aes = FALSE, color='blue', linetype='dotted', 'size'=1)+
  geom_ribbon(data=heardata_pred,aes(x=PTA, ymin=logit(PC_silence_Cinf/100), ymax=logit(PC_silence_Csup/100)), inherit.aes = FALSE, color=NA, alpha=0.1, fill='blue')+
  geom_ribbon(data=heardata_pred,aes(x=PTA, ymin=logit(PC_noise_Cinf/100), ymax=logit(PC_noise_Csup/100)), inherit.aes = FALSE, color=NA, alpha=0.1, fill='blue')

# Compute goodness of fit measures

Rsquared_m1.1hi <- Rsquared(data,parsfit)
log_lik_1.1hi <- extract_log_lik(fit.m1.1hi)
pwll_1.1hi <- apply(log_lik_1.1hi, 2, mean)
waic_1.1hi <- waic(log_lik_1.1hi)
pwwaic_1.1hi <- waic_1.1hi$pointwise[,"elpd_waic"]
loo_1.1hi <- loo(log_lik_1.1hi)
print(loo_1.1hi)

## 6. Hier. GLM w only main effects (without PTA) -----------------

m1.2hi <- stan_model(file = 'm1.2hi_bis.stan')

fit.m1.2hi <- sampling(m1.2hi,
                       data = data,
                       chains = 4,             # number of Markov chains
                       warmup = 3000,          # number of warmup iterations per chain
                       iter = 7000,            # total number of iterations per chain
                       refresh = 1000)
# diagnosis

parameters = c("beta_0","gamma_0","beta_age","beta_cond","beta_gender","plapse")#,"beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA","plapse")#c("beta_0","beta_PTA","beta_age","beta_cond","beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA")
print(fit.m1.2hi, pars = parameters)
#launch_shinystan(fit.m1.2hi)
#plot(fit.m1.2hi, pars = parameters)+ coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1))#+xlim(-3,1)
plot(fit.m1.2hi, show_density = TRUE, show_outer_line = FALSE, ci_level= 0.95,outer_level= 0.99, pars = parameters) + ggtitle("m1.2hi") #+ coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1))#+xlim(-3,1)

# counterfactual predictions

parsfit<-extract(fit.m1.2hi,pars=rev(fit.m1.2hi@model_pars))#c("beta_0","beta_PTA","beta_age","beta_cond","beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA","gamma_0","gamma_PTA","gamma_age","gamma_cond","gamma_condPTA","gamma_agePTA","gamma_agecond","gamma_agecondPTA"))
Nsamples = length(parsfit$beta_0)

source("counterfactuals_simple.R")

p3.1c +
  geom_line(data=heardata_pred_center,aes(x=PTA, y=PC_silence_med, color=center), inherit.aes = FALSE,  'size'=1)+
  geom_line(data=heardata_pred_center,aes(x=PTA, y=PC_noise_med, color=center), inherit.aes = FALSE, linetype='longdash', 'size'=1)+
  geom_ribbon(data=heardata_pred_center,aes(x=PTA, ymin=PC_silence_Cinf, ymax=PC_silence_Csup, fill=center), inherit.aes = FALSE, color=NA, alpha=0.1)+
  geom_ribbon(data=heardata_pred_center,aes(x=PTA, ymin=PC_noise_Cinf, ymax=PC_noise_Csup,  fill=center), inherit.aes = FALSE, color=NA, alpha=0.1)+
  ylim(0, 100)

p3.1b +
  geom_line(data=heardata_pred_center,aes(x=PTA, y=PC_silence_med, color=center), inherit.aes = FALSE, 'size'=1)+
  geom_line(data=heardata_pred_center,aes(x=PTA, y=PC_noise_med, color=center), inherit.aes = FALSE, 'size'=1, linetype='longdash')

p3.1b +
  geom_line(data=heardata_pred,aes(x=PTA, y=PC_silence_med), inherit.aes = FALSE, color='blue', 'size'=1)+
  geom_line(data=heardata_pred,aes(x=PTA, y=PC_noise_med), inherit.aes = FALSE, color='blue', linetype='dotted', 'size'=1)+
  geom_ribbon(data=heardata_pred,aes(x=PTA, ymin=PC_silence_Cinf, ymax=PC_silence_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='blue')+
  geom_ribbon(data=heardata_pred,aes(x=PTA, ymin=PC_noise_Cinf, ymax=PC_noise_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='blue')

# Compute goodness of fit measures

Rsquared_m1.2hi <- Rsquared(data,parsfit)
log_lik_1.2hi <- extract_log_lik(fit.m1.2hi)
waic_1.2hi <- waic(log_lik_1.2hi)
loo_1.2hi <- loo(log_lik_1.2hi)
print(loo_1.2hi)

## 7. Hier. GLM w only main effects (without age) -----------------

m1.3hi <- stan_model(file = 'm1.3hi_bis.stan')

fit.m1.3hi <- sampling(m1.3hi,
                       data = data,
                       chains = 4,             # number of Markov chains
                       warmup = 3000,          # number of warmup iterations per chain
                       iter = 7000,            # total number of iterations per chain
                       refresh = 1000)
# diagnosis

parameters = c("beta_PTA","beta_cond","beta_gender")#,"beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA","plapse")#c("beta_0","beta_PTA","beta_age","beta_cond","beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA")
parameters = c("beta_0","gamma_0","beta_PTA","beta_cond","beta_gender","plapse")#,"beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA","plapse")#c("beta_0","beta_PTA","beta_age","beta_cond","beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA")
print(fit.m1.3hi, pars = parameters)
#launch_shinystan(fit.m1.3hi)
plot(fit.m1.3hi, pars = parameters)+ coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1))#+xlim(-3,1)
plot(fit.m1.3hi, show_density = TRUE, show_outer_line = FALSE, ci_level= 0.95,outer_level= 0.99, pars = parameters) + ggtitle("m1.3hi") #+ coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1))#+xlim(-3,1)

# counterfactual predictions

parsfit<-extract(fit.m1.3hi,pars=rev(fit.m1.3hi@model_pars))#c("beta_0","beta_PTA","beta_age","beta_cond","beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA","gamma_0","gamma_PTA","gamma_age","gamma_cond","gamma_condPTA","gamma_agePTA","gamma_agecond","gamma_agecondPTA"))
Nsamples = length(parsfit$beta_0)

source("counterfactuals_simple.R")

p3.1c +
  geom_line(data=heardata_pred_center,aes(x=PTA, y=PC_silence_med, color=center), inherit.aes = FALSE,  'size'=1)+
  geom_line(data=heardata_pred_center,aes(x=PTA, y=PC_noise_med, color=center), inherit.aes = FALSE, linetype='longdash', 'size'=1)+
  geom_ribbon(data=heardata_pred_center,aes(x=PTA, ymin=PC_silence_Cinf, ymax=PC_silence_Csup, fill=center), inherit.aes = FALSE, color=NA, alpha=0.1)+
  geom_ribbon(data=heardata_pred_center,aes(x=PTA, ymin=PC_noise_Cinf, ymax=PC_noise_Csup,  fill=center), inherit.aes = FALSE, color=NA, alpha=0.1)+
  ylim(0, 100)

p3.1b +
  geom_line(data=heardata_pred_center,aes(x=PTA, y=PC_silence_med, color=center), inherit.aes = FALSE, 'size'=1)+
  geom_line(data=heardata_pred_center,aes(x=PTA, y=PC_noise_med, color=center), inherit.aes = FALSE, 'size'=1, linetype='longdash')

p3.1b +
  geom_line(data=heardata_pred,aes(x=PTA, y=PC_silence_med), inherit.aes = FALSE, color='blue', 'size'=1)+
  geom_line(data=heardata_pred,aes(x=PTA, y=PC_noise_med), inherit.aes = FALSE, color='blue', linetype='dotted', 'size'=1)+
  geom_ribbon(data=heardata_pred,aes(x=PTA, ymin=PC_silence_Cinf, ymax=PC_silence_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='blue')+
  geom_ribbon(data=heardata_pred,aes(x=PTA, ymin=PC_noise_Cinf, ymax=PC_noise_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='blue')

# Compute goodness of fit measures

Rsquared_m1.3hi <- Rsquared(data,parsfit)
log_lik_1.3hi <- extract_log_lik(fit.m1.3hi)
waic_1.3hi <- waic(log_lik_1.3hi)
loo_1.3hi <- loo(log_lik_1.3hi)
print(loo_1.3hi)

## 8. Full hier. GLM on intercept -----------------

m2.1hi <- stan_model(file = 'm2.1hi_bis.stan')

fit.m2.1hi <- sampling(m2.1hi,
                       data = data,
                       chains = 4,             # number of Markov chains
                       warmup = 3000,          # number of warmup iterations per chain
                       iter = 7000,            # total number of iterations per chain
                       refresh = 1000)
# diagnosis

parameters = c("beta_0","gamma_0","beta_PTA","beta_age","beta_agePTA","beta_cond","beta_condPTA","beta_agecond","beta_agecondPTA","beta_gender","plapse")#,"beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA","plapse")#c("beta_0","beta_PTA","beta_age","beta_cond","beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA")
print(fit.m2.1hi, pars = parameters)
#launch_shinystan(fit.m2.1hi)
plot(fit.m2.1hi, show_density = TRUE, show_outer_line = FALSE, ci_level= 0.95,outer_level= 0.99, pars = parameters) + ggtitle("m2.1hi") #+ coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1))#+xlim(-3,1)

# counterfactual predictions

parsfit<-extract(fit.m2.1hi,pars=rev(fit.m2.1hi@model_pars))#c("beta_0","beta_PTA","beta_age","beta_cond","beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA","gamma_0","gamma_PTA","gamma_age","gamma_cond","gamma_condPTA","gamma_agePTA","gamma_agecond","gamma_agecondPTA"))
Nsamples = length(parsfit$beta_0)

source("counterfactuals_full.R")

p3.1c +
  geom_line(data=heardata_pred_center,aes(x=PTA, y=PC_silence_med, color=center), inherit.aes = FALSE,  'size'=1)+
  geom_line(data=heardata_pred_center,aes(x=PTA, y=PC_noise_med, color=center), inherit.aes = FALSE, linetype='longdash', 'size'=1)+
  geom_ribbon(data=heardata_pred_center,aes(x=PTA, ymin=PC_silence_Cinf, ymax=PC_silence_Csup, fill=center), inherit.aes = FALSE, color=NA, alpha=0.1)+
  geom_ribbon(data=heardata_pred_center,aes(x=PTA, ymin=PC_noise_Cinf, ymax=PC_noise_Csup,  fill=center), inherit.aes = FALSE, color=NA, alpha=0.1)+
  ylim(0, 100)

p3.1b +
  geom_line(data=heardata_pred_center,aes(x=PTA, y=PC_silence_med, color=center), inherit.aes = FALSE, 'size'=1)+
  geom_line(data=heardata_pred_center,aes(x=PTA, y=PC_noise_med, color=center), inherit.aes = FALSE, 'size'=1, linetype='longdash')

p3.1b +
  geom_line(data=heardata_pred,aes(x=PTA, y=PC_silence_med), inherit.aes = FALSE, color='blue', 'size'=1)+
  geom_line(data=heardata_pred,aes(x=PTA, y=PC_noise_med), inherit.aes = FALSE, color='blue', linetype='dotted', 'size'=1)+
  geom_ribbon(data=heardata_pred,aes(x=PTA, ymin=PC_silence_Cinf, ymax=PC_silence_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='blue')+
  geom_ribbon(data=heardata_pred,aes(x=PTA, ymin=PC_noise_Cinf, ymax=PC_noise_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='blue')

# Compute goodness of fit measures

Rsquared_m2.1hi <- Rsquared(data,parsfit)
log_lik_2.1hi <- extract_log_lik(fit.m2.1hi)
pwll_2.1hi <- apply(log_lik_2.1hi, 2, mean)
waic_2.1hi <- waic(log_lik_2.1hi)
pwwaic_2.1hi <- waic_2.1hi$pointwise[,"elpd_waic"]
loo_2.1hi <- loo(log_lik_2.1hi)
print(loo_2.1hi)
