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
#library(psych)

rm(list=ls())
setwd("C:/Users/Léo/ownCloud/Professionnel/Projet Genopath/011_Presbyacousie_Lorenzi/R_analyses/Presbycusis_analysis")
theme_set(theme_bw())
heardata = read.table('mat_full.txt')
colnames(heardata) <- c("group", "age", "PTA", "subject", "center", "gender", "NC_silence", "NC_noise")
ESII = read.table('mat_ESII.txt')


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
heardata$groupcenter = with(heardata, interaction(group,  center, sep = " - "))
heardata$ESII_s <- ESII$V3[1:(nrow(ESII)/2)]
heardata$ESII_n <- ESII$V3[((nrow(ESII)/2)+1):nrow(ESII)]

# transform continuous variable age into a 3-level factor for representation purpose
agefactor_cutoff = seq(from=min_age-1,to=max_age+1,length.out=4)
heardata$agefactor <- cut.default(heardata$age, agefactor_cutoff)
levels(heardata$agefactor) <- c(paste("age = ", floor(agefactor_cutoff[1]), "-", floor(agefactor_cutoff[2]),"y"),
                                paste("age = ", floor(agefactor_cutoff[2]), "-", floor(agefactor_cutoff[3]),"y"),
                                paste("age = ", floor(agefactor_cutoff[3]), "-", floor(agefactor_cutoff[4]),"y"))#paste("age (y) =", levels(heardata$agefactor))

heardata2 <- data.frame(  unlist(list(heardata$group,heardata$group)),
                          c(heardata$age,heardata$age),
                          c(heardata$PTA,heardata$PTA),
                          unlist(list(heardata$subject,heardata$subject)),
                          unlist(list(heardata$center,heardata$center)),
                          unlist(list(heardata$gender,heardata$gender)),
                          c(heardata$NC_silence,heardata$NC_noise),
                          c(heardata$PC_silence,heardata$PC_noise),
                          c(rep(1,length(heardata$NC_silence)),rep(2,length(heardata$NC_noise))),
                          ESII$V3)
colnames(heardata2) <- c("group", "age", "PTA", "subject", "center", "gender", "NC", "PC", 'condition', 'ESII')
heardata2$condition <- factor(heardata2$condition)
levels(heardata2$condition)[levels(heardata2$condition)=="1"] <- "Silence"
levels(heardata2$condition)[levels(heardata2$condition)=="2"] <- "Noise"
heardata2$agefactor <- cut.default(heardata2$age, agefactor_cutoff)
levels(heardata2$agefactor) <- c(paste("age = ", floor(agefactor_cutoff[1]), "-", floor(agefactor_cutoff[2]),"y"),
                                paste("age = ", floor(agefactor_cutoff[2]), "-", floor(agefactor_cutoff[3]),"y"),
                                paste("age = ", floor(agefactor_cutoff[3]), "-", floor(agefactor_cutoff[4]),"y"))#paste("age (y) =", levels(heardata$agefactor))

heardata2$lPC <- qlogis(heardata2$PC/100)

# Zscored data and counterfactual predictors

PTA_pred = seq(from = 0, to = 70, by = 7)
NPTApred = length(PTA_pred)
heardata$PTAz = (heardata$PTA-mean(heardata$PTA))/sd(heardata$PTA)
heardata2$PTAz = (heardata2$PTA-mean(heardata2$PTA))/sd(heardata2$PTA)
PTAz_pred = (PTA_pred-mean(heardata2$PTA))/sd(heardata2$PTA)

age_pred = seq(from = min_age, to = max_age, by = 6)
Nagepred = length(age_pred)
heardata$agez = (heardata$age-mean(heardata$age))/sd(heardata$age)
heardata2$agez = (heardata2$age-mean(heardata2$age))/sd(heardata2$age)
agez_pred = (age_pred-mean(heardata2$age))/sd(heardata2$age)

ESII_pred = seq(from = 0, to = 1, by = 0.1)
NESIIpred = length(ESII_pred)
heardata$ESII_nz = (heardata$ESII_n-mean(heardata2$ESII))/sd(heardata2$ESII)
heardata$ESII_sz = (heardata$ESII_s-mean(heardata2$ESII))/sd(heardata2$ESII)
heardata2$ESIIz = (heardata2$ESII-mean(heardata2$ESII))/sd(heardata2$ESII)
ESIIz_pred = (ESII_pred-mean(heardata2$ESII))/sd(heardata2$ESII)

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
              ESII_nz = heardata$ESII_nz[heardata$group=="HI"],
              ESII_sz = heardata$ESII_sz[heardata$group=="HI"],
              prior_only = 0)
#What to do with missing values?
# data$NCn[data$NCn==3] <- 100 # Ignore the missing data istead of replacing it with 3%
data$Nmissing = sum(data$NCn==100)

# Inclusion thresholds

InclusionHI = read.table('mat_inclusionHI.txt')
colnames(InclusionHI) <- c("age", "PTA_female", "PTA_male")
InclusionHI$PTA_male = as.numeric(InclusionHI$PTA_male)
InclusionHI$PTA_female = as.numeric(InclusionHI$PTA_female)
InclusionHI$age = as.numeric(InclusionHI$age)
InclusionNH = read.table('mat_inclusionNH.txt')
colnames(InclusionNH) <- c("age", "PTA_female", "PTA_male")
InclusionNH$PTA_male = as.numeric(InclusionNH$PTA_male)
InclusionNH$PTA_female = as.numeric(InclusionNH$PTA_female)
InclusionNH$age = as.numeric(InclusionNH$age)

# Functions

Rsquared <- function(data,parsfit){
  SSerr = c((data$lPC_silence-apply(parsfit$eta_s,2,mean))^2,(data$lPC_noise-apply(parsfit$eta_n,2,mean))^2)
  SSerr = mean(SSerr[SSerr!=Inf])
  meanlPC = mean(c(data$lPC_silence[data$lPC_silence!=Inf],data$lPC_noise[data$lPC_noise!=Inf]))
  SStot = c((data$lPC_silence-mean(c(data$lPC_silence[data$lPC_silence!=Inf])))^2, (data$lPC_noise-mean(c(data$lPC_noise[data$lPC_noise!=Inf])))^2)
  SStot = mean(SStot[SStot!=Inf])
  Rsquared = 1-SSerr/SStot
}

Rsquared_noise <- function(data,parsfit){
  SSerr = (data$lPC_noise-apply(parsfit$eta_n,2,mean))^2
  SSerr = mean(SSerr[SSerr!=Inf])
  meanlPC = mean(data$lPC_noise[data$lPC_noise!=Inf])
  SStot = (data$lPC_noise-mean(data$lPC_noise[data$lPC_noise!=Inf]))^2
  SStot = mean(SStot[SStot!=Inf])
  Rsquared_noise = 1-SSerr/SStot
}

# Plot PTA vs. age -----------------------------------------------------

p1 <- ggplot() + geom_ribbon(aes(x=InclusionHI$age,ymin=InclusionHI$PTA_male,ymax=rep(70,1,100)), fill="gray",alpha=0.3)+ 
  geom_ribbon(aes(x=InclusionHI$age,ymin=InclusionHI$PTA_female,ymax=rep(70,1,100)), fill="gray",alpha=0.2)
# p1 <- p1 + geom_ribbon(aes(x=InclusionNH$age,ymin=rep(0,1,100),ymax=InclusionNH$PTA_male), fill="#FF0000",alpha=0.1)+ 
#   geom_ribbon(aes(x=InclusionNH$age,ymin=rep(0,1,100),ymax=InclusionNH$PTA_female), fill="#FF0000",alpha=0.2)
#   #scale_alpha(range=c(0,0.5))
p1 <- p1 + geom_count(data = heardata, aes(x=age, y=PTA, color=group, alpha=0.5)) +
  scale_size_area(max_size = 2,name = "N",n.breaks = 3) +
  labs(y="PTA (dB HL)", x = "age (years)") + guides(alpha = FALSE, size = FALSE)

p1 <- p1 + xlim(40, 95) + ylim(0, 70)

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

ggsave("FigureInclusion.pdf", width = 6, height = 5)

# Plot ESII vs. PTA --------------------

p0 <- ggplot(data=heardata) + 
  geom_count(aes(x=PTA, y=ESII_n, color=group, alpha=0.5)) +
  geom_count(aes(x=PTA, y=ESII_s, color=group, alpha=0.5)) +
  #scale_size_area(max_size = 2,name = "N",n.breaks = 3) +
  labs(x = "PTA", y = "ESII") + guides(alpha = FALSE, size = FALSE)

p0 <- p0 + xlim(0, 70) + ylim(0, 1)
p0

# Plot PTA vs. PC --------------------

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

## Summarize effects of group and center on PTA, age and scores ----

aggregate( formula = cbind(PTA, age) ~ group, data = heardata, FUN = function(x) c(mean = mean(x), sd = sd(x), med= median(x)))
aggregate( formula = cbind(PC_silence, PC_noise) ~ group, data = heardata, FUN = function(x) c(mean = mean(x), sd = sd(x), med= median(x)))

aggregate( formula = cbind(PTA, age) ~ group + center, data = heardata, FUN = function(x) c(mean = mean(x), sd = sd(x), med= median(x)))
aggregate( formula = cbind(PC_silence, PC_noise) ~ group + center, data = heardata, FUN = function(x) c(mean = mean(x), sd = sd(x), med= median(x)))

## Plot PTA vs masking effect --------------

p3.0a <- ggplot(data = heardata, aes(x=PC_silence, y=PC_noise, color=group,alpha=0.5)) + 
  geom_count() +
  scale_size_area(max_size = 2,name = "N",n.breaks = 5) +
  geom_abline(intercept = 0, slope = 1,linetype="dashed", color = "black")+
  xlim(0, 100) + ylim(0, 100) + guides(alpha = FALSE, size = FALSE)+
  labs(y="Score in noise (% correct)", x = "Score in quiet (% correct)")
  
p3.0axhist <- axis_canvas(p3.0a, axis = "x") + 
  geom_histogram(data = heardata,
                 aes(x = PC_silence, color = group, fill = group, alpha = 0.1),
                 binwidth = 5, position="identity")+
  xlim(0, 100)
p3.0axSD <- axis_canvas(p3.0axhist, axis = "x") + 
  geom_boxplot(data = heardata, aes(x = PC_silence,color = group), outlier.shape = NA)+
  xlim(0, 100)
p3.0ayhist <- axis_canvas(p3.0a, axis = "y", coord_flip = TRUE) + 
  geom_histogram(data = heardata,
                 aes(x = PC_noise, color = group, fill = group, alpha=0.1),
                 binwidth = 5, position="identity") +
  coord_flip() + xlim(0, 100)
p3.0aySD <- axis_canvas(p3.0ayhist, axis = "y", coord_flip = TRUE) + 
  geom_boxplot(data = heardata, aes(x = PC_noise, color = group), outlier.shape = NA) +
  coord_flip() + xlim(0, 100)
p3.0axhist <- insert_xaxis_grob(p3.0axhist, p3.0axSD, grid::unit(.2, "null"), position = "top")
p3.0a <- insert_xaxis_grob(p3.0a, p3.0axhist, grid::unit(.2, "null"), position = "top")
p3.0ayhist <- insert_yaxis_grob(p3.0ayhist, p3.0aySD, grid::unit(.2, "null"), position = "right")
p3.0a <- insert_yaxis_grob(p3.0a, p3.0ayhist, grid::unit(.2, "null"), position = "right")
ggdraw(p3.0a)

p3.0b <- ggplot(data = subset(heardata, group=="HI"),aes(x=lPC_silence, y=lPC_noise, color=center)) + 
  geom_count() +
  scale_size_area(max_size = 2,name = "N",n.breaks = 5) +
  geom_abline(intercept = 0, slope = 1,linetype="dashed", color = "black")

p3.1a <- ggplot(data = subset(heardata,group=="HI")) + 
  geom_point(aes(x=PTA, y=lPC_silence, color=center)) +
  geom_point(aes(x=PTA, y=lPC_noise, color=center), shape=1) +
  #geom_segment(aes(x=PTA, xend=PTA, y=lPC_noise, yend=lPC_silence, color=center)) +
  #scale_size(name = "N")+
  ylab('logit(PC)')+
  scale_size_area(max_size = 2,name = "N",n.breaks = 5) +
  facet_grid(~agefactor)

p3.1b <- ggplot(data = subset(heardata,group=="HI")) + 
  geom_count(aes(x=PTA, y=PC_silence), alpha = 1) +#, color=center
  geom_count(aes(x=PTA, y=PC_noise), shape=1, alpha = 1) +#, color=center
  #geom_segment(aes(x=PTA, xend=PTA, y=lPC_noise, yend=lPC_silence, color=center)) +
  #scale_size(name = "N")
  ylim(0, 100) + guides(alpha = FALSE, size = FALSE)+
  scale_alpha(range = c(0, 1))+
  ylab('Intellitest scores (% correct)')+
  xlab('PTA (dB HL)')+
  scale_size_area(max_size = 1.5,name = "N",n.breaks = 5) +
  facet_grid(~agefactor) 

p3.1c <- ggplot(data = subset(heardata,group=="HI")) + 
  geom_count(aes(x=PTA, y=PC_silence, color=center)) +
  geom_count(aes(x=PTA, y=PC_noise, color=center), shape=1) +
  #geom_segment(aes(x=PTA, xend=PTA, y=lPC_noise, yend=lPC_silence, color=center)) +
  #scale_size(name = "N")+
  ylim(0, 100) + guides(alpha = FALSE, size = FALSE)+
  scale_alpha(range = c(0, 1))+
  ylab('Intellitest scores (% correct)')+
  xlab('PTA (dB HL)')+
  scale_size_area(max_size = 1.5,name = "N",n.breaks = 5) +
  facet_grid(center~agefactor)+ 
  theme(legend.position = "none")

## Plot ESII vs masking effect --------------

p4.1b <- ggplot(data = subset(heardata,group=="HI")) + 
  geom_count(aes(x=ESII_s, y=PC_silence), alpha = 1) +#, color=center
  geom_count(aes(x=ESII_n, y=PC_noise), shape=1, alpha = 1) +#, color=center
  #geom_segment(aes(x=PTA, xend=PTA, y=lPC_noise, yend=lPC_silence, color=center)) +
  #scale_size(name = "N")
  ylim(0, 100) + guides(alpha = FALSE, size = FALSE)+
  scale_alpha(range = c(0, 1))+
  ylab('Intellitest scores (% correct)')+
  xlab('ESII')+
  scale_size_area(max_size = 1.5,name = "N",n.breaks = 5) +
  facet_grid(~agefactor)

## 1.1 Simple non-hier model (w only main effects) -----------------

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

parsfit$gamma_0 <- cbind(parsfit$beta_0,parsfit$beta_0,parsfit$beta_0,parsfit$beta_0,parsfit$beta_0,parsfit$beta_0,parsfit$beta_0)

source("counterfactuals_simple.R")
 
p3.1b +
  geom_line(data=heardata_counterf,aes(x=PTA, y=PC_silence_med), inherit.aes = FALSE, color='black', 'size'=1)+
  geom_line(data=heardata_counterf,aes(x=PTA, y=PC_noise_med), inherit.aes = FALSE, color='black', 'size'=1, linetype='longdash')+
  geom_ribbon(data=heardata_counterf,aes(x=PTA, ymin=PC_silence_Cinf, ymax=PC_silence_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='black')+
  geom_ribbon(data=heardata_counterf,aes(x=PTA, ymin=PC_noise_Cinf, ymax=PC_noise_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='black')
ylim(0, 100)

p3.1b +
  geom_point(data=subset(heardata,group=="HI"),aes(x=PTA, y=PC_noise_pred), inherit.aes = FALSE, color='blue')
# 
#   geom_line(data=heardata_counterf,aes(x=PTA, y=PC_silence_med), inherit.aes = FALSE, color='black', 'size'=1)+
#   geom_line(data=heardata_counterf,aes(x=PTA, y=PC_noise_med), inherit.aes = FALSE, color='black', 'size'=1, linetype='longdash')+
#   geom_ribbon(data=heardata_counterf,aes(x=PTA, ymin=PC_silence_Cinf, ymax=PC_silence_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='black')+
#   geom_ribbon(data=heardata_counterf,aes(x=PTA, ymin=PC_noise_Cinf, ymax=PC_noise_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='black')
ylim(0, 100)


# Compute goodness of fit measures

Rsquared_m1.1 <- Rsquared(data,parsfit)
log_lik_1.1 <- extract_log_lik(fit.m1.1)
loo_1.1 <- loo(log_lik_1.1)
print(loo_1.1)

## 2.1 Full non-hier. GLM on intercept -----------------

m2.1 <- stan_model(file = 'm2.1.stan')

fit.m2.1 <- sampling(m2.1,
                       data = data,
                       chains = 4,             # number of Markov chains
                       warmup = 3000,          # number of warmup iterations per chain
                       iter = 7000,            # total number of iterations per chain
                       refresh = 1000)
# diagnosis

parameters = c("beta_0","beta_PTA","beta_age","beta_agePTA","beta_cond","beta_condPTA","beta_agecond","beta_agecondPTA","beta_gender","plapse")#,"beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA","plapse")#c("beta_0","beta_PTA","beta_age","beta_cond","beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA")
print(fit.m2.1, pars = parameters)
#launch_shinystan(fit.m2.1hi)
plot(fit.m2.1, show_density = TRUE, show_outer_line = FALSE, ci_level= 0.95,outer_level= 0.99, pars = parameters) + ggtitle("m2.1") #+ coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1))#+xlim(-3,1)

# counterfactual predictions

parsfit<-extract(fit.m2.1,pars=rev(fit.m2.1@model_pars))#c("beta_0","beta_PTA","beta_age","beta_cond","beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA","gamma_0","gamma_PTA","gamma_age","gamma_cond","gamma_condPTA","gamma_agePTA","gamma_agecond","gamma_agecondPTA"))
Nsamples = length(parsfit$beta_0)

parsfit$gamma_0 <- cbind(parsfit$beta_0,parsfit$beta_0,parsfit$beta_0,parsfit$beta_0,parsfit$beta_0,parsfit$beta_0,parsfit$beta_0)

source("counterfactuals_full.R")

p3.1c +
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_silence_med, color=center), inherit.aes = FALSE,  'size'=1)+
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_noise_med, color=center), inherit.aes = FALSE, linetype='longdash', 'size'=1)+
  geom_ribbon(data=heardata_counterf_center,aes(x=PTA, ymin=PC_silence_Cinf, ymax=PC_silence_Csup, fill=center), inherit.aes = FALSE, color=NA, alpha=0.1)+
  geom_ribbon(data=heardata_counterf_center,aes(x=PTA, ymin=PC_noise_Cinf, ymax=PC_noise_Csup,  fill=center), inherit.aes = FALSE, color=NA, alpha=0.1)+
  ylim(0, 100)

p3.1b +
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_silence_med, color=center), inherit.aes = FALSE, 'size'=1)+
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_noise_med, color=center), inherit.aes = FALSE, 'size'=1, linetype='longdash')

p3.1b +
  geom_line(data=heardata_counterf,aes(x=PTA, y=PC_silence_med), inherit.aes = FALSE, color='blue', 'size'=1)+
  geom_line(data=heardata_counterf,aes(x=PTA, y=PC_noise_med), inherit.aes = FALSE, color='blue', linetype='dotted', 'size'=1)+
  geom_ribbon(data=heardata_counterf,aes(x=PTA, ymin=PC_silence_Cinf, ymax=PC_silence_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='blue')+
  geom_ribbon(data=heardata_counterf,aes(x=PTA, ymin=PC_noise_Cinf, ymax=PC_noise_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='blue')

# Compute goodness of fit measures

Rsquared_m2.1 <- Rsquared(data,parsfit)
log_lik_2.1 <- extract_log_lik(fit.m2.1)
loo_2.1 <- loo(log_lik_2.1)
print(loo_2.1)


## 1.1hi Simple hier. GLM w only main effects -----------------

m1.1hi <- stan_model(file = 'm1.1hi_quater.stan')

fit.m1.1hi <- sampling(m1.1hi,
                     data = data,
                     chains = 4,             # number of Markov chains
                     warmup = 3000,          # number of warmup iterations per chain
                     iter = 7000,            # total number of iterations per chain
                     refresh = 1000)
# diagnosis

parameters = c("beta_0","gamma_0","beta_PTA","beta_age","beta_cond","beta_gender","plapse","sigma_0")#,"beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA","plapse")#c("beta_0","beta_PTA","beta_age","beta_cond","beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA")
print(fit.m1.1hi, pars = parameters)
#launch_shinystan(fit.m1.1hi)
plot(fit.m1.1hi, show_density = TRUE, show_outer_line = FALSE, ci_level= 0.95,outer_level= 0.99, pars = parameters) + ggtitle("m1.1hi") + xlim(-3,3) #+ coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1))#+xlim(-3,1)

# counterfactual predictions & posterior predictions

parsfit<-extract(fit.m1.1hi,pars=rev(fit.m1.1hi@model_pars))#c("beta_0","beta_PTA","beta_age","beta_cond","beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA","gamma_0","gamma_PTA","gamma_age","gamma_cond","gamma_condPTA","gamma_agePTA","gamma_agecond","gamma_agecondPTA"))
Nsamples = length(parsfit$beta_0)

source("posterior_predictions.R")
source("counterfactuals_simple.R")

p3.1c_fit <- p3.1c +
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_silence_med, color=center), inherit.aes = FALSE,  'size'=1)+
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_noise_med, color=center), inherit.aes = FALSE, linetype='longdash', 'size'=1)+
  geom_ribbon(data=heardata_counterf_center,aes(x=PTA, ymin=PC_silence_Cinf, ymax=PC_silence_Csup, fill=center), inherit.aes = FALSE, color=NA, alpha=0.1)+
  geom_ribbon(data=heardata_counterf_center,aes(x=PTA, ymin=PC_noise_Cinf, ymax=PC_noise_Csup,  fill=center), inherit.aes = FALSE, color=NA, alpha=0.1)+
  ylim(0, 100)
p3.1c_fit 

p3.1b +
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_silence_med, color=center), inherit.aes = FALSE, 'size'=1)+
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_noise_med, color=center), inherit.aes = FALSE, 'size'=1, linetype='longdash')

p3.1b_fit <- p3.1b +
  geom_line(data=heardata_counterf,aes(x=PTA, y=PC_silence_med), inherit.aes = FALSE, color='blue', 'size'=1)+
  geom_line(data=heardata_counterf,aes(x=PTA, y=PC_noise_med), inherit.aes = FALSE, color='blue', linetype='dotted', 'size'=1)+
  geom_ribbon(data=heardata_counterf,aes(x=PTA, ymin=PC_silence_Cinf, ymax=PC_silence_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='blue')+
  geom_ribbon(data=heardata_counterf,aes(x=PTA, ymin=PC_noise_Cinf, ymax=PC_noise_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='blue')
p3.1b_fit 

p3.1b +
  geom_count(data=subset(heardata,group=="HI"),aes(x=PTA, y=PC_noise_pred), inherit.aes = FALSE, color='blue', shape=1, alpha = 1)+
  geom_count(data=subset(heardata,group=="HI"),aes(x=PTA, y=PC_silence_pred), inherit.aes = FALSE, color='blue', alpha = 1)+
  ylim(0, 100) + xlim(PTA_pred[1], PTA_pred[length(PTA_pred)])

# Compute goodness of fit measures

Rsquared_m1.1hi <- Rsquared(data,parsfit)
Rsquaredn_m1.1hi <- Rsquared_noise(data,parsfit)
log_lik_1.1hi <- extract_log_lik(fit.m1.1hi)
pwll_1.1hi <- apply(log_lik_1.1hi, 2, mean)
waic_1.1hi <- waic(log_lik_1.1hi)
pwwaic_1.1hi <- waic_1.1hi$pointwise[,"elpd_waic"]
loo_1.1hi <- loo(log_lik_1.1hi)
print(loo_1.1hi)

## extract slope ----

# calculate PTA slope for each iteration

PTAslope_s = (1-1/16)*parsfit$beta_PTA
PTAslope_n = (1-1/16-parsfit$plapse)*parsfit$beta_PTA
PTAslope_diff = PTAslope_s - PTAslope_n

C = apply(PTAslope_s,2,quantile,probs=c(0.025,0.5,0.975),na.rm = TRUE) #the median line with 95% credible intervals

p <- ggplot(data.frame(PTAslope=c(PTAslope_s,PTAslope_n),cond=c(rep("Silence",Nsamples),rep("Noise",Nsamples))), aes(x=PTAslope,color=cond)) + 
  geom_density()

## 1.2hi Simple hier. GLM w only main effects (without PTA) -----------------

m1.2hi <- stan_model(file = 'm1.2hi_quater.stan')

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
plot(fit.m1.2hi, show_density = TRUE, show_outer_line = FALSE, ci_level= 0.95,outer_level= 0.99, pars = parameters) + ggtitle("m1.2hi") #+ coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1))#+xlim(-3,1)

# counterfactual predictions & posterior predictions

parsfit<-extract(fit.m1.2hi,pars=rev(fit.m1.2hi@model_pars))#c("beta_0","beta_PTA","beta_age","beta_cond","beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA","gamma_0","gamma_PTA","gamma_age","gamma_cond","gamma_condPTA","gamma_agePTA","gamma_agecond","gamma_agecondPTA"))
Nsamples = length(parsfit$beta_0)

source("posterior_predictions.R")
source("counterfactuals_simple.R")

p3.1c +
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_silence_med, color=center), inherit.aes = FALSE,  'size'=1)+
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_noise_med, color=center), inherit.aes = FALSE, linetype='longdash', 'size'=1)+
  geom_ribbon(data=heardata_counterf_center,aes(x=PTA, ymin=PC_silence_Cinf, ymax=PC_silence_Csup, fill=center), inherit.aes = FALSE, color=NA, alpha=0.1)+
  geom_ribbon(data=heardata_counterf_center,aes(x=PTA, ymin=PC_noise_Cinf, ymax=PC_noise_Csup,  fill=center), inherit.aes = FALSE, color=NA, alpha=0.1)+
  ylim(0, 100)

p3.1b +
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_silence_med, color=center), inherit.aes = FALSE, 'size'=1)+
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_noise_med, color=center), inherit.aes = FALSE, 'size'=1, linetype='longdash')

p3.1b +
  geom_line(data=heardata_counterf,aes(x=PTA, y=PC_silence_med), inherit.aes = FALSE, color='blue', 'size'=1)+
  geom_line(data=heardata_counterf,aes(x=PTA, y=PC_noise_med), inherit.aes = FALSE, color='blue', linetype='dotted', 'size'=1)+
  geom_ribbon(data=heardata_counterf,aes(x=PTA, ymin=PC_silence_Cinf, ymax=PC_silence_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='blue')+
  geom_ribbon(data=heardata_counterf,aes(x=PTA, ymin=PC_noise_Cinf, ymax=PC_noise_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='blue')

p3.1b +
  geom_count(data=subset(heardata,group=="HI"),aes(x=PTA, y=PC_noise_pred), inherit.aes = FALSE, color='blue', shape=1, alpha = 1)+
  geom_count(data=subset(heardata,group=="HI"),aes(x=PTA, y=PC_silence_pred), inherit.aes = FALSE, color='blue', alpha = 1)+
  ylim(0, 100) + xlim(PTA_pred[1], PTA_pred[length(PTA_pred)])

# Compute goodness of fit measures

Rsquared_m1.2hi <- Rsquared(data,parsfit)
Rsquaredn_m1.2hi <- Rsquared_noise(data,parsfit)
log_lik_1.2hi <- extract_log_lik(fit.m1.2hi)
pwll_1.2hi <- apply(log_lik_1.2hi, 2, mean)
waic_1.2hi <- waic(log_lik_1.2hi)
loo_1.2hi <- loo(log_lik_1.2hi)
print(loo_1.2hi)

## 1.3hi Simple hier. GLM w only main effects (without age) -----------------

m1.3hi <- stan_model(file = 'm1.3hi_quater.stan')

fit.m1.3hi <- sampling(m1.3hi,
                       data = data,
                       chains = 4,             # number of Markov chains
                       warmup = 3000,          # number of warmup iterations per chain
                       iter = 7000,            # total number of iterations per chain
                       refresh = 1000)
# diagnosis

parameters = c("beta_0","gamma_0","beta_PTA","beta_cond","plapse","sigma_0")#,"beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA","plapse")#c("beta_0","beta_PTA","beta_age","beta_cond","beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA")
print(fit.m1.3hi, pars = parameters)
#launch_shinystan(fit.m1.3hi)
plot(fit.m1.3hi, show_density = TRUE, show_outer_line = FALSE, ci_level= 0.95,outer_level= 0.99, pars = parameters) + ggtitle("m1.3hi") + xlim(-3,3) #+ coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1))#

# counterfactual predictions & posterior predictions

parsfit<-extract(fit.m1.3hi,pars=rev(fit.m1.3hi@model_pars))#c("beta_0","beta_PTA","beta_age","beta_cond","beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA","gamma_0","gamma_PTA","gamma_age","gamma_cond","gamma_condPTA","gamma_agePTA","gamma_agecond","gamma_agecondPTA"))
Nsamples = length(parsfit$beta_0)

source("posterior_predictions.R")
source("counterfactuals_simple.R")

p3.1c +
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_silence_med, color=center), inherit.aes = FALSE,  'size'=1)+
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_noise_med, color=center), inherit.aes = FALSE, linetype='longdash', 'size'=1)+
  geom_ribbon(data=heardata_counterf_center,aes(x=PTA, ymin=PC_silence_Cinf, ymax=PC_silence_Csup, fill=center), inherit.aes = FALSE, color=NA, alpha=0.1)+
  geom_ribbon(data=heardata_counterf_center,aes(x=PTA, ymin=PC_noise_Cinf, ymax=PC_noise_Csup,  fill=center), inherit.aes = FALSE, color=NA, alpha=0.1)+
  ylim(0, 100)

p3.1b +
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_silence_med, color=center), inherit.aes = FALSE, 'size'=1)+
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_noise_med, color=center), inherit.aes = FALSE, 'size'=1, linetype='longdash')

p3.1b +
  geom_line(data=heardata_counterf,aes(x=PTA, y=PC_silence_med), inherit.aes = FALSE, color='blue', 'size'=1)+
  geom_line(data=heardata_counterf,aes(x=PTA, y=PC_noise_med), inherit.aes = FALSE, color='blue', linetype='dotted', 'size'=1)+
  geom_ribbon(data=heardata_counterf,aes(x=PTA, ymin=PC_silence_Cinf, ymax=PC_silence_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='blue')+
  geom_ribbon(data=heardata_counterf,aes(x=PTA, ymin=PC_noise_Cinf, ymax=PC_noise_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='blue')

p3.1b +
  geom_count(data=subset(heardata,group=="HI"),aes(x=PTA, y=PC_noise_pred), inherit.aes = FALSE, color='blue', shape=1, alpha = 1)+
  geom_count(data=subset(heardata,group=="HI"),aes(x=PTA, y=PC_silence_pred), inherit.aes = FALSE, color='blue', alpha = 1)+
  ylim(0, 100) + xlim(PTA_pred[1], PTA_pred[length(PTA_pred)])

# Compute goodness of fit measures

Rsquared_m1.3hi <- Rsquared(data,parsfit)
Rsquaredn_m1.3hi <- Rsquared_noise(data,parsfit)
log_lik_1.3hi <- extract_log_lik(fit.m1.3hi)
pwll_1.3hi <- apply(log_lik_1.3hi, 2, mean)
waic_1.3hi <- waic(log_lik_1.3hi)
loo_1.3hi <- loo(log_lik_1.3hi)
print(loo_1.3hi)

## 2.1hi Full hier. GLM on intercept -----------------

m2.1hi <- stan_model(file = 'm2.1hi_quater.stan')

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
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_silence_med, color=center), inherit.aes = FALSE,  'size'=1)+
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_noise_med, color=center), inherit.aes = FALSE, linetype='longdash', 'size'=1)+
  geom_ribbon(data=heardata_counterf_center,aes(x=PTA, ymin=PC_silence_Cinf, ymax=PC_silence_Csup, fill=center), inherit.aes = FALSE, color=NA, alpha=0.1)+
  geom_ribbon(data=heardata_counterf_center,aes(x=PTA, ymin=PC_noise_Cinf, ymax=PC_noise_Csup,  fill=center), inherit.aes = FALSE, color=NA, alpha=0.1)+
  ylim(0, 100)

p3.1b +
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_silence_med, color=center), inherit.aes = FALSE, 'size'=1)+
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_noise_med, color=center), inherit.aes = FALSE, 'size'=1, linetype='longdash')

p3.1b +
  geom_line(data=heardata_counterf,aes(x=PTA, y=PC_silence_med), inherit.aes = FALSE, color='blue', 'size'=1)+
  geom_line(data=heardata_counterf,aes(x=PTA, y=PC_noise_med), inherit.aes = FALSE, color='blue', linetype='dotted', 'size'=1)+
  geom_ribbon(data=heardata_counterf,aes(x=PTA, ymin=PC_silence_Cinf, ymax=PC_silence_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='blue')+
  geom_ribbon(data=heardata_counterf,aes(x=PTA, ymin=PC_noise_Cinf, ymax=PC_noise_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='blue')

# Compute goodness of fit measures

Rsquared_m2.1hi <- Rsquared(data,parsfit)
Rsquaredn_m2.1hi <- Rsquared_noise(data,parsfit)
log_lik_2.1hi <- extract_log_lik(fit.m2.1hi)
pwll_2.1hi <- apply(log_lik_2.1hi, 2, mean)
waic_2.1hi <- waic(log_lik_2.1hi)
pwwaic_2.1hi <- waic_2.1hi$pointwise[,"elpd_waic"]
loo_2.1hi <- loo(log_lik_2.1hi)
pwelpdloo_2.1hi <- loo_2.1hi$pointwise[,"elpd_loo"]
print(loo_2.1hi)

## 2.2hi Full hier. GLM on intercept (without cond*PTA and age*cond*PTA) -----------------

m2.2hi <- stan_model(file = 'm2.2hi_quater.stan')

fit.m2.2hi <- sampling(m2.2hi,
                       data = data,
                       chains = 4,             # number of Markov chains
                       warmup = 3000,          # number of warmup iterations per chain
                       iter = 7000,            # total number of iterations per chain
                       refresh = 1000)
# diagnosis

parameters = c("beta_0","gamma_0","beta_PTA","beta_age","beta_agePTA","beta_cond","beta_condPTA","beta_agecond","beta_agecondPTA","beta_gender","plapse")#,"beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA","plapse")#c("beta_0","beta_PTA","beta_age","beta_cond","beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA")
print(fit.m2.2hi, pars = parameters)
#launch_shinystan(fit.m2.1hi)
plot(fit.m2.2hi, show_density = TRUE, show_outer_line = FALSE, ci_level= 0.95,outer_level= 0.99, pars = parameters) + ggtitle("m2.2hi") #+ coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1))#+xlim(-3,1)

# counterfactual predictions

parsfit<-extract(fit.m2.2hi,pars=rev(fit.m2.2hi@model_pars))#c("beta_0","beta_PTA","beta_age","beta_cond","beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA","gamma_0","gamma_PTA","gamma_age","gamma_cond","gamma_condPTA","gamma_agePTA","gamma_agecond","gamma_agecondPTA"))
Nsamples = length(parsfit$beta_0)

source("counterfactuals_full.R")

p3.1c +
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_silence_med, color=center), inherit.aes = FALSE,  'size'=1)+
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_noise_med, color=center), inherit.aes = FALSE, linetype='longdash', 'size'=1)+
  geom_ribbon(data=heardata_counterf_center,aes(x=PTA, ymin=PC_silence_Cinf, ymax=PC_silence_Csup, fill=center), inherit.aes = FALSE, color=NA, alpha=0.1)+
  geom_ribbon(data=heardata_counterf_center,aes(x=PTA, ymin=PC_noise_Cinf, ymax=PC_noise_Csup,  fill=center), inherit.aes = FALSE, color=NA, alpha=0.1)+
  ylim(0, 100)

p3.1b +
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_silence_med, color=center), inherit.aes = FALSE, 'size'=1)+
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_noise_med, color=center), inherit.aes = FALSE, 'size'=1, linetype='longdash')

p3.1b +
  geom_line(data=heardata_counterf,aes(x=PTA, y=PC_silence_med), inherit.aes = FALSE, color='blue', 'size'=1)+
  geom_line(data=heardata_counterf,aes(x=PTA, y=PC_noise_med), inherit.aes = FALSE, color='blue', linetype='dotted', 'size'=1)+
  geom_ribbon(data=heardata_counterf,aes(x=PTA, ymin=PC_silence_Cinf, ymax=PC_silence_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='blue')+
  geom_ribbon(data=heardata_counterf,aes(x=PTA, ymin=PC_noise_Cinf, ymax=PC_noise_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='blue')

# Compute goodness of fit measures

Rsquared_m2.2hi <- Rsquared(data,parsfit)
log_lik_2.2hi <- extract_log_lik(fit.m2.2hi)
pwll_2.2hi <- apply(log_lik_2.2hi, 2, mean)
waic_2.2hi <- waic(log_lik_2.2hi)
pwwaic_2.2hi <- waic_2.2hi$pointwise[,"elpd_waic"]
loo_2.2hi <- loo(log_lik_2.2hi)
pwelpdloo_2.2hi <- loo_2.2hi$pointwise[,"elpd_loo"]
print(loo_2.2hi)

## 2.3hi Hier. GLM on intercept (with main effect + cond*PTA) -----------------

m2.3hi <- stan_model(file = 'm2.3hi_quater.stan')

fit.m2.3hi <- sampling(m2.3hi,
                       data = data,
                       chains = 4,             # number of Markov chains
                       warmup = 3000,          # number of warmup iterations per chain
                       iter = 7000,            # total number of iterations per chain
                       refresh = 1000)
# diagnosis

parameters = c("beta_0","gamma_0","beta_PTA","beta_age","beta_agePTA","beta_cond","beta_condPTA","beta_agecond","beta_agecondPTA","beta_gender","plapse")#,"beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA","plapse")#c("beta_0","beta_PTA","beta_age","beta_cond","beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA")
print(fit.m2.3hi, pars = parameters)
#launch_shinystan(fit.m2.1hi)
plot(fit.m2.3hi, show_density = TRUE, show_outer_line = FALSE, ci_level= 0.95,outer_level= 0.99, pars = parameters) + ggtitle("m2.3hi") #+ coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1))#+xlim(-3,1)

# counterfactual predictions

parsfit<-extract(fit.m2.3hi,pars=rev(fit.m2.3hi@model_pars))#c("beta_0","beta_PTA","beta_age","beta_cond","beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA","gamma_0","gamma_PTA","gamma_age","gamma_cond","gamma_condPTA","gamma_agePTA","gamma_agecond","gamma_agecondPTA"))
Nsamples = length(parsfit$beta_0)

source("counterfactuals_full.R")

p3.1c +
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_silence_med, color=center), inherit.aes = FALSE,  'size'=1)+
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_noise_med, color=center), inherit.aes = FALSE, linetype='longdash', 'size'=1)+
  geom_ribbon(data=heardata_counterf_center,aes(x=PTA, ymin=PC_silence_Cinf, ymax=PC_silence_Csup, fill=center), inherit.aes = FALSE, color=NA, alpha=0.1)+
  geom_ribbon(data=heardata_counterf_center,aes(x=PTA, ymin=PC_noise_Cinf, ymax=PC_noise_Csup,  fill=center), inherit.aes = FALSE, color=NA, alpha=0.1)+
  ylim(0, 100)

p3.1b +
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_silence_med, color=center), inherit.aes = FALSE, 'size'=1)+
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_noise_med, color=center), inherit.aes = FALSE, 'size'=1, linetype='longdash')

p3.1b +
  geom_line(data=heardata_counterf,aes(x=PTA, y=PC_silence_med), inherit.aes = FALSE, color='blue', 'size'=1)+
  geom_line(data=heardata_counterf,aes(x=PTA, y=PC_noise_med), inherit.aes = FALSE, color='blue', linetype='dotted', 'size'=1)+
  geom_ribbon(data=heardata_counterf,aes(x=PTA, ymin=PC_silence_Cinf, ymax=PC_silence_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='blue')+
  geom_ribbon(data=heardata_counterf,aes(x=PTA, ymin=PC_noise_Cinf, ymax=PC_noise_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='blue')

# Compute goodness of fit measures

Rsquared_m2.3hi <- Rsquared(data,parsfit)
log_lik_2.3hi <- extract_log_lik(fit.m2.3hi)
pwll_2.3hi <- apply(log_lik_2.3hi, 2, mean)
waic_2.3hi <- waic(log_lik_2.3hi)
pwwaic_2.3hi <- waic_2.3hi$pointwise[,"elpd_waic"]
loo_2.3hi <- loo(log_lik_2.3hi)
pwelpdloo_2.3hi <- loo_2.3hi$pointwise[,"elpd_loo"]
print(loo_2.3hi)

## weight ratio -------

mean(parsfit$beta_PTA/parsfit$beta_age)

## compare predictions ----------

heardata2_fit = subset(heardata2, group=="HI")
heardata2_fit$pwll_2.1hi = pwll_2.1hi
heardata2_fit$pwelpdloo_2.1hi = pwelpdloo_2.1hi
heardata2_fit$gamma2 <- rep(apply(parsfit$gamma2_0,2,mean),2)
ggplot(data = heardata2_fit, aes(x=PTA, y=gamma2, color=center))+
  geom_point(aes(shape=condition))+
  scale_shape_manual(values=c(16, 1))+
  facet_grid(center~agefactor)
ggplot(data = heardata2_fit, aes(x=pwelpdloo_2.1hi, y=gamma2))+  geom_point(aes(color=condition))

# error of fit
heardataHI=subset(heardata,heardata$group=="HI")
heardataHI$PC_silence_pred <- apply(parsfit$p_s,2,mean)*100
heardataHI$PC_noise_pred <- apply(parsfit$p_n,2,mean)*100
heardataHI$lPC_silence_pred <- apply(parsfit$eta_s,2,mean)
heardataHI$lPC_noise_pred <- apply(parsfit$eta_n,2,mean)

p3.1b + 
  geom_point(aes(x=PTA, y=heardataHI$PC_silence_pred, color=center,alpha=0.5)) +
  geom_point(aes(x=PTA, y=heardataHI$PC_noise_pred, color=center,alpha=0.5), shape=1)
  
ggplot(data = subset(heardataHI,NC_noise>3)) +
  geom_point(aes(x=PTA, y=PC_silence-PC_silence_pred, color=center))+
  geom_point(aes(x=PTA, y=PC_noise-PC_noise_pred, color=center), shape=1)+  
  facet_grid(~agefactor)

## -- 3.1hi Simple hier. GLM ESII -------

m3.1hi <- stan_model(file = 'm3.1hi.stan')

fit.m3.1hi <- sampling(m3.1hi,
                       data = data,
                       chains = 4,             # number of Markov chains
                       warmup = 3000,          # number of warmup iterations per chain
                       iter = 7000,            # total number of iterations per chain
                       refresh = 1000)
# diagnosis

parameters = c("beta_0","gamma_0","beta_ESII","beta_cond","beta_age","beta_gender","plapse","sigma_0")#,"beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA","plapse")#c("beta_0","beta_PTA","beta_age","beta_cond","beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA")
print(fit.m3.1hi, pars = parameters)
#launch_shinystan(fit.m1.1hi)
plot(fit.m3.1hi, show_density = TRUE, show_outer_line = FALSE, ci_level= 0.95,outer_level= 0.99, pars = parameters) + ggtitle("m3.1hi") + xlim(-3,3) #+ coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1))#+xlim(-3,1)

# counterfactual predictions

parsfit<-extract(fit.m3.1hi,pars=rev(fit.m3.1hi@model_pars))#c("beta_0","beta_PTA","beta_age","beta_cond","beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA","gamma_0","gamma_PTA","gamma_age","gamma_cond","gamma_condPTA","gamma_agePTA","gamma_agecond","gamma_agecondPTA"))
Nsamples = length(parsfit$beta_0)

source("counterfactuals_ESII.R")

p4.1b_fit <- p4.1b +
  geom_line(data=heardata_counterf,aes(x=ESII, y=PC_silence_med), inherit.aes = FALSE, color='blue', 'size'=1)+
  geom_line(data=heardata_counterf,aes(x=ESII, y=PC_noise_med), inherit.aes = FALSE, color='blue', linetype='dotted', 'size'=1)+
  geom_ribbon(data=heardata_counterf,aes(x=ESII, ymin=PC_silence_Cinf, ymax=PC_silence_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='blue')+
  geom_ribbon(data=heardata_counterf,aes(x=ESII, ymin=PC_noise_Cinf, ymax=PC_noise_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='blue')
p4.1b_fit 

p4.1b_fit <- p4.1b +
  #geom_point(data=subset(heardata,group=="HI"),aes(x=ESII_s, y=PC_silence_pred), inherit.aes = FALSE, color='blue')+
  geom_point(data=subset(heardata,group=="HI"),aes(x=ESII_n, y=PC_noise_pred), inherit.aes = FALSE, color='blue')
p4.1b_fit 

## 4.1hi Full hier. GLM on intercept -----------------

m4.1hi <- stan_model(file = 'm4.1hi.stan')

fit.m4.1hi <- sampling(m4.1hi,
                       data = data,
                       chains = 4,             # number of Markov chains
                       warmup = 3000,          # number of warmup iterations per chain
                       iter = 7000,            # total number of iterations per chain
                       refresh = 1000)
# diagnosis

parameters = c("beta_0","gamma_0","beta_ESII","beta_age","beta_ageESII","beta_cond","beta_condESII","beta_agecond","beta_agecondESII","beta_gender","plapse")#,"beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA","plapse")#c("beta_0","beta_PTA","beta_age","beta_cond","beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA")
print(fit.m4.1hi, pars = parameters)
#launch_shinystan(fit.m2.1hi)
plot(fit.m4.1hi, show_density = TRUE, show_outer_line = FALSE, ci_level= 0.95,outer_level= 0.99, pars = parameters) + ggtitle("m4.1hi") #+ coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1))#+xlim(-3,1)

# counterfactual predictions

parsfit<-extract(fit.m4.1hi,pars=rev(fit.m4.1hi@model_pars))#c("beta_0","beta_PTA","beta_age","beta_cond","beta_condPTA","beta_agePTA","beta_agecond","beta_agecondPTA","gamma_0","gamma_PTA","gamma_age","gamma_cond","gamma_condPTA","gamma_agePTA","gamma_agecond","gamma_agecondPTA"))
Nsamples = length(parsfit$beta_0)

source("counterfactuals_full_ESII.R")

p4.1b_fit <- p4.1b +
  #geom_point(data=subset(heardata,group=="HI"),aes(x=ESII_s, y=PC_silence_pred), inherit.aes = FALSE, color='blue')+
  geom_point(data=subset(heardata,group=="HI"),aes(x=ESII_n, y=PC_noise_pred), inherit.aes = FALSE, color='blue')
p4.1b_fit 

p4.1c +
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_silence_med, color=center), inherit.aes = FALSE,  'size'=1)+
  geom_line(data=heardata_counterf_center,aes(x=PTA, y=PC_noise_med, color=center), inherit.aes = FALSE, linetype='longdash', 'size'=1)+
  geom_ribbon(data=heardata_counterf_center,aes(x=PTA, ymin=PC_silence_Cinf, ymax=PC_silence_Csup, fill=center), inherit.aes = FALSE, color=NA, alpha=0.1)+
  geom_ribbon(data=heardata_counterf_center,aes(x=PTA, ymin=PC_noise_Cinf, ymax=PC_noise_Csup,  fill=center), inherit.aes = FALSE, color=NA, alpha=0.1)+
  ylim(0, 100)

p4.1b +
  geom_line(data=heardata_counterf_center,aes(x=ESII, y=PC_silence_med, color=center), inherit.aes = FALSE, 'size'=1)+
  geom_line(data=heardata_counterf_center,aes(x=ESII, y=PC_noise_med, color=center), inherit.aes = FALSE, 'size'=1, linetype='longdash')

p4.1b +
  geom_line(data=heardata_counterf,aes(x=ESII, y=PC_silence_med), inherit.aes = FALSE, color='blue', 'size'=1)+
  geom_line(data=heardata_counterf,aes(x=ESII, y=PC_noise_med), inherit.aes = FALSE, color='blue', linetype='dotted', 'size'=1)+
  geom_ribbon(data=heardata_counterf,aes(x=ESII, ymin=PC_silence_Cinf, ymax=PC_silence_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='blue')+
  geom_ribbon(data=heardata_counterf,aes(x=ESII, ymin=PC_noise_Cinf, ymax=PC_noise_Csup), inherit.aes = FALSE, color=NA, alpha=0.1, fill='blue')

# Compute goodness of fit measures

Rsquared_m4.1hi <- Rsquared(data,parsfit)
log_lik_4.1hi <- extract_log_lik(fit.m4.1hi)
pwll_4.1hi <- apply(log_lik_4.1hi, 2, mean)
waic_4.1hi <- waic(log_lik_4.1hi)
pwwaic_4.1hi <- waic_4.1hi$pointwise[,"elpd_waic"]
loo_4.1hi <- loo(log_lik_4.1hi)
pwelpdloo_4.1hi <- loo_4.1hi$pointwise[,"elpd_loo"]
print(loo_4.1hi)

## Generate Figures for manuscript -----------

plot_grid(p1, p3.0a, labels=c("A", "B"), ncol = 2, nrow = 1)

ggsave("FigureResults.pdf", width = 11, height = 5)

p3.1b_fit <- p3.1b_fit + facet_grid(group~agefactor)+theme(strip.text.y = element_text(color ="lightgray"))
plot_grid(p3.1b_fit, p3.1c_fit, labels=c("A", "B"), ncol = 1, nrow = 2, rel_heights = c(1, 3))

ggsave("FigureModelFit.pdf", width = 6, height = 10)
