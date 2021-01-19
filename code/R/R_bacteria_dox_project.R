
### Bacteria-Dox Interactions

# Ryan Blaustein
# ryan.blaustein1@gmail.com
# last edited: 01-18-21

### call packages
library(DescTools)
library(emmeans)
library(ggplot2)
library(gplots)
library(gtools)
library(gridExtra)
library(lme4)
library(lmerTest)
library(pbkrtest)
library(reshape2)
library(stringr)

### set working directory to 'bacteria_doxo_interactions/data_for_R' (https://github.com/hartmann-lab/bacteria_doxo_interactions)


########################################
### Individual Strains & Doxorubicin ###
########################################

### load data

# dox standard curve
dox_curve = read.table("individual_strains/1_dox_st_curve.txt", sep='\t',h=T,row.names=1,check=F,comment='')

# bacteria growth curves
mono_data = read.table("individual_strains/2_growth_curves.txt", sep='\t',h=T,row.names=1,check=F,comment='')

# aerobic growth
aerobic_growth = read.table("individual_strains/3_aerobic_growth.txt", sep='\t',h=T,row.names=1,check=F,comment='')

# efflux pump inhibitor (PABN)
pabn = read.table("individual_strains/4_pabn.txt",  sep='\t',h=T,row.names=1,check=F,comment='')

# reporter strains in spent media
sm_data = read.table("individual_strains/5_spent_media.txt", sep='\t',h=T,row.names=1,check=F,comment='')


###
### Dox in broth media
###

## concentration stability over time (negative control)
neg_ct = data.frame(trt = c("Control", "Low", "Medium", "High", "Control", "Low", "Medium", "High"),
                    Time = c("0 h", "0 h", "0 h", "0 h", "24 h", "24 h", "24 h", "24 h"),
                    dox = c(0, 14.32, 103.62, 238.35, 0, 16.17, 103.93, 241.43),
                    dox_se = c(0.92, 1.67, 4.81, 4.72, 0.15, 0.77, 2.22, 21.68))
neg_ct$trt = factor(neg_ct$trt, levels = c("Control", "Low", "Medium", "High"))

# plot (Fig. S1)
ggplot(neg_ct, aes(x = trt, y = dox, fill = Time)) +
  geom_bar(stat="identity", position=position_dodge2(preserve="single")) +
  geom_errorbar(stat="identity", position=position_dodge2(preserve="single"),
                aes(ymax=dox+dox_se, ymin=dox-dox_se)) +
  xlab("Dox Treatment") +
  ylab("Dox µM") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = c(0.15, 0.8))

## standard curve (Fig. S2 A)
ggplot(dox_curve,
       aes(x=Conc_uM, y=OD480/1000)) +
  geom_point(size=3) +
  geom_smooth(method="lm") +
  theme_bw() +
  xlab("Dox µM") +
  ylab("λ480") +
  theme(axis.title=element_text(size = 20),
        axis.text=element_text(size = 18))

# regression parameters
summary(lm(OD480~Conc_uM, data=dox_curve))


###
### Bacterial growth & drug transformation in model system (anaerobic)
###

## prep table
mono_data_avgs = mono_data[which(is.na(mono_data$OD600_corr_avg) == F),]
mono_data_avgs$Treatment = factor(mono_data_avgs$Treatment,
                                  levels = c("Ctrl", "Dox_10µM", "Dox_100µM", "Dox_250µM"))

### strains with no drug (control)
ggplot(mono_data_avgs[grep("Ctrl", mono_data_avgs$Treatment),],
       aes(x=Time_Point, y=OD600_corr_avg/1000, col=Strain)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymax=OD600_corr_avg/1000 + OD600_corr_se/1000,
                    ymin=OD600_corr_avg/1000 - OD600_corr_se/1000), width=1.1) +
  geom_line(lty=1) +
  theme_bw() +
  xlab("Time (hr)") +
  ylab("OD600") +
  theme(axis.title=element_text(size = 18),
        axis.text=element_text(size = 16),
        legend.text=element_text(size = 16),
        legend.title=element_text(size = 16)) +
  theme(legend.position = "none")


### growth curves with drug

# C. innocuum
CI = ggplot(mono_data_avgs[grep("innocuum", mono_data_avgs$Strain),],
            aes(x=Time_Point, y=OD600_corr_avg/1000, col=Treatment)) +
  geom_point(size=2.7) +
  geom_errorbar(aes(ymax=OD600_corr_avg/1000 + OD600_corr_se/1000,
                    ymin=OD600_corr_avg/1000 - OD600_corr_se/1000), width=1.2) +
  geom_line(lty=1, aes(col=Treatment)) +
  geom_point(aes(y=Dox_conc_avg/450), size=2.7, shape=17) +
  geom_line(aes(y=Dox_conc_avg/450), lty=2) +
  geom_errorbar(aes(ymax=Dox_conc_avg/450 + Dox_conc_se/450, 
                    ymin=Dox_conc_avg/450 - Dox_conc_se/450), width=1.2) +
  theme_bw() +
  xlab("Time (hr)") +
  ylab("OD600") +
  scale_y_continuous(limits = c(-0.02, 0.58),
                     sec.axis=sec_axis(~.*450, name="Dox uM")) +
  scale_color_manual(values=c("black", "blue", "orange", "red"),
                     labels=c("Control", "Dox 10µM", "Dox 100µM", "Dox 250µM")) +
  theme(axis.title=element_text(size = 13),
        axis.text=element_text(size = 12),
        legend.text=element_text(size = 12),
        legend.title=element_text(size = 12)) +
  theme(legend.position = "none") +
  theme(axis.title = element_text(size=0))

# E. coli
EC = ggplot(mono_data_avgs[grep("coli", mono_data_avgs$Strain),],
            aes(x=Time_Point, y=OD600_corr_avg/1000, col=Treatment)) +
  geom_point(size=2.7) +
  geom_errorbar(aes(ymax=OD600_corr_avg/1000 + OD600_corr_se/1000,
                    ymin=OD600_corr_avg/1000 - OD600_corr_se/1000), width=1.2) +
  geom_line(lty=1, aes(col=Treatment)) +
  geom_point(aes(y=Dox_conc_avg/450), size=2.7, shape=17) +
  geom_line(aes(y=Dox_conc_avg/450), lty=2) +
  geom_errorbar(aes(ymax=Dox_conc_avg/450 + Dox_conc_se/450, 
                    ymin=Dox_conc_avg/450 - Dox_conc_se/450), width=1.2) +
  theme_bw() +
  xlab("Time (hr)") +
  ylab("OD600") +
  scale_y_continuous(limits = c(-0.02, 0.58),
                     sec.axis=sec_axis(~.*450, name="Dox uM")) +
  scale_color_manual(values=c("black", "blue", "orange", "red")) +
  theme(axis.title=element_text(size = 13),
        axis.text=element_text(size = 12),
        legend.text=element_text(size = 12),
        legend.title=element_text(size = 12)) +
  theme(legend.position = "none") +
  theme(axis.title = element_text(size=0))

# E. faecium
EF = ggplot(mono_data_avgs[grep("faecium", mono_data_avgs$Strain),],
            aes(x=Time_Point, y=OD600_corr_avg/1000, col=Treatment)) +
  geom_point(size=2.7) +
  geom_errorbar(aes(ymax=OD600_corr_avg/1000 + OD600_corr_se/1000,
                    ymin=OD600_corr_avg/1000 - OD600_corr_se/1000), width=1.2) +
  geom_line(lty=1, aes(col=Treatment)) +
  geom_point(aes(y=Dox_conc_avg/450), size=2.7, shape=17) +
  geom_line(aes(y=Dox_conc_avg/450), lty=2) +
  geom_errorbar(aes(ymax=Dox_conc_avg/450 + Dox_conc_se/450, 
                    ymin=Dox_conc_avg/450 - Dox_conc_se/450), width=1.2) +
  theme_bw() +
  xlab("Time (hr)") +
  ylab("OD600") +
  scale_y_continuous(limits = c(-0.02, 0.58),
                     sec.axis=sec_axis(~.*450, name="Dox uM")) +
  scale_color_manual(values=c("black", "blue", "orange", "red")) +
  theme(axis.title=element_text(size = 13),
        axis.text=element_text(size = 12),
        legend.text=element_text(size = 12),
        legend.title=element_text(size = 12)) +
  theme(legend.position = "none") +
  theme(axis.title = element_text(size=0))

# K. pneumoniae
KP = ggplot(mono_data_avgs[grep("pneumoniae", mono_data_avgs$Strain),],
            aes(x=Time_Point, y=OD600_corr_avg/1000, col=Treatment)) +
  geom_point(size=2.7) +
  geom_errorbar(aes(ymax=OD600_corr_avg/1000 + OD600_corr_se/1000,
                    ymin=OD600_corr_avg/1000 - OD600_corr_se/1000), width=1.2) +
  geom_line(lty=1, aes(col=Treatment)) +
  geom_point(aes(y=Dox_conc_avg/450), size=2.7, shape=17) +
  geom_line(aes(y=Dox_conc_avg/450), lty=2) +
  geom_errorbar(aes(ymax=Dox_conc_avg/450 + Dox_conc_se/450, 
                    ymin=Dox_conc_avg/450 - Dox_conc_se/450), width=1.2) +
  theme_bw() +
  xlab("Time (hr)") +
  ylab("OD600") +
  scale_y_continuous(limits = c(-0.02, 0.58),
                     sec.axis=sec_axis(~.*450, name="Dox uM")) +
  scale_color_manual(values=c("black", "blue", "orange", "red")) +
  theme(axis.title=element_text(size = 13),
        axis.text=element_text(size = 12),
        legend.text=element_text(size = 12),
        legend.title=element_text(size = 12)) +
  theme(legend.position = "none") +
  theme(axis.title = element_text(size=0))

# L. rhamnosus
Lac = ggplot(mono_data_avgs[grep("Lacto", mono_data_avgs$Strain),],
             aes(x=Time_Point, y=OD600_corr_avg/1000, col=Treatment)) +
  geom_point(size=2.7) +
  geom_errorbar(aes(ymax=OD600_corr_avg/1000 + OD600_corr_se/1000,
                    ymin=OD600_corr_avg/1000 - OD600_corr_se/1000), width=1.2) +
  geom_line(lty=1, aes(col=Treatment)) +
  geom_point(aes(y=Dox_conc_avg/450), size=2.7, shape=17) +
  geom_line(aes(y=Dox_conc_avg/450), lty=2) +
  geom_errorbar(aes(ymax=Dox_conc_avg/450 + Dox_conc_se/450, 
                    ymin=Dox_conc_avg/450 - Dox_conc_se/450), width=1.2) +
  theme_bw() +
  xlab("Time (hr)") +
  ylab("OD600") +
  scale_y_continuous(limits = c(-0.02, 0.58),
                     sec.axis=sec_axis(~.*450, name="Dox uM")) +
  scale_color_manual(values=c("black", "blue", "orange", "red")) +
  theme(axis.title=element_text(size = 13),
        axis.text=element_text(size = 12),
        legend.text=element_text(size = 12),
        legend.title=element_text(size = 12)) +
  theme(legend.position = "none") +
  theme(axis.title = element_text(size=0))

## plot (Fig. 1)
grid.arrange(KP, EC, EF, CI, Lac, nrow=2)


### stats: bacterial growth

## linear mixed effects (LME) model 

# K. pneumoniae
summary(lmer(OD600_corr ~ Time_Point*Treatment + (1|Treatment),
             data = mono_data_stats[grep("K. pneum", mono_data_stats$Strain),]))
# E. coli
summary(lmer(OD600_corr ~ Time_Point*Treatment + (1|Treatment),
             data = mono_data_stats[grep("E. coli", mono_data_stats$Strain),]))
# E. faecium
summary(lmer(OD600_corr ~ Time_Point*Treatment + (1|Treatment),
             data = mono_data_stats[grep("E. fae", mono_data_stats$Strain),]))
# C. innocuum
summary(lmer(OD600_corr ~ Time_Point*Treatment + (1|Treatment),
             data = mono_data_stats[grep("C. inn", mono_data_stats$Strain),]))
# L. rhamnosus
summary(lmer(OD600_corr ~ Time_Point*Treatment + (1|Treatment),
             data = mono_data_stats[grep("Lac", mono_data_stats$Strain),]))


### stats: drug transformation

## prep table: area under the curve (AUC)
mono_data_stats = mono_data
mono_data_stats$Time_Point = factor(as.character(mono_data_stats$Time_Point),
                                    levels = c("0", "3", "6", "12", "18", "24"))
mono_data_stats$Treatment = factor(mono_data_stats$Treatment,
                                   levels = c("Ctrl", "Dox_10µM", "Dox_100µM", "Dox_250µM"))
mono_data_stats_AUC = mono_data_stats

## normalize to controls

# KP
mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]
mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("Ctrl", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Dox_conc
mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Dox_conc = 
  mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Dox_conc - mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("Ctrl", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Dox_conc
mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Dox_conc =
  mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Dox_conc - mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("Ctrl", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Dox_conc
mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Dox_conc =
  mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Dox_conc - mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("Ctrl", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Dox_conc

# EC
mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]
mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("Ctrl", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Dox_conc
mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Dox_conc = 
  mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Dox_conc - mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("Ctrl", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Dox_conc
mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Dox_conc =
  mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Dox_conc - mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("Ctrl", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Dox_conc
mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Dox_conc =
  mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Dox_conc - mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("Ctrl", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Dox_conc

# AUC:KP
AUC(x=mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("Ctrl", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Time_hr,
    y=mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("Ctrl", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Dox_conc)
AUC(x=mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Time_hr,
    y=mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Dox_conc)
AUC(x=mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Time_hr,
    y=mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Dox_conc)
AUC(x=mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Time_hr,
    y=mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Dox_conc)

# AUC:EC
AUC(x=mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("Ctrl", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Time_hr,
    y=mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("Ctrl", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Dox_conc)
AUC(x=mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Time_hr,
    y=mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Dox_conc)
AUC(x=mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Time_hr,
    y=mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Dox_conc)
AUC(x=mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Time_hr,
    y=mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Dox_conc)

# AUC by replicate
AUC_stats = data.frame(strain = c(rep("EC", 9), rep("KP", 9)), 
                       dox = c(rep("10", 3), rep("100", 3), rep("250", 3), rep("10", 3), rep("100", 3), rep("250", 3)),
                       AUC = c(AUC(x=mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),][grep("1", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Time_hr,
                                   y=mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),][grep("1", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Dox_conc),
                               AUC(x=mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),][grep("2", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Time_hr,
                                   y=mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),][grep("2", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Dox_conc),
                               AUC(x=mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),][grep("3", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Time_hr,
                                   y=mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),][grep("3", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Dox_conc),
                               AUC(x=mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),][grep("1", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Time_hr,
                                   y=mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),][grep("1", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Dox_conc),
                               AUC(x=mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),][grep("2", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Time_hr,
                                   y=mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),][grep("2", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Dox_conc),
                               AUC(x=mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),][grep("3", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Time_hr,
                                   y=mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),][grep("3", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Dox_conc),
                               AUC(x=mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),][grep("1", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Time_hr,
                                   y=mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),][grep("1", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Dox_conc),
                               AUC(x=mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),][grep("2", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Time_hr,
                                   y=mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),][grep("2", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Dox_conc),
                               AUC(x=mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),][grep("3", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Time_hr,
                                   y=mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),][grep("3", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("E. coli", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Dox_conc),
                               AUC(x=mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),][grep("1", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Time_hr,
                                   y=mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),][grep("1", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Dox_conc),
                               AUC(x=mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),][grep("2", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Time_hr,
                                   y=mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),][grep("2", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Dox_conc),
                               AUC(x=mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),][grep("3", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Time_hr,
                                   y=mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),][grep("3", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("10µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Dox_conc),
                               AUC(x=mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),][grep("1", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Time_hr,
                                   y=mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),][grep("1", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Dox_conc),
                               AUC(x=mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),][grep("2", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Time_hr,
                                   y=mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),][grep("2", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Dox_conc),
                               AUC(x=mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),][grep("3", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Time_hr,
                                   y=mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),][grep("3", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("100µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Dox_conc),
                               AUC(x=mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),][grep("1", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Time_hr,
                                   y=mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),][grep("1", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Dox_conc),
                               AUC(x=mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),][grep("2", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Time_hr,
                                   y=mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),][grep("2", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Dox_conc),
                               AUC(x=mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),][grep("3", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Time_hr,
                                   y=mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),][grep("3", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),][grep("250µM", mono_data_stats_AUC[grep("K. pneum", mono_data_stats_AUC$Strain),]$Treatment),]$Rep),]$Dox_conc)))

## AUC KP vs. EC: t-test & ANOVA
t.test(AUC_stats[grep("EC", AUC_stats$strain),][grep("^10$", AUC_stats[grep("EC", AUC_stats$strain),]$dox),]$AUC,
       AUC_stats[grep("KP", AUC_stats$strain),][grep("^10$", AUC_stats[grep("KP", AUC_stats$strain),]$dox),]$AUC)
t.test(AUC_stats[grep("EC", AUC_stats$strain),][grep("100", AUC_stats[grep("EC", AUC_stats$strain),]$dox),]$AUC,
       AUC_stats[grep("KP", AUC_stats$strain),][grep("100", AUC_stats[grep("KP", AUC_stats$strain),]$dox),]$AUC)
t.test(AUC_stats[grep("EC", AUC_stats$strain),][grep("250", AUC_stats[grep("EC", AUC_stats$strain),]$dox),]$AUC,
       AUC_stats[grep("KP", AUC_stats$strain),][grep("250", AUC_stats[grep("KP", AUC_stats$strain),]$dox),]$AUC)
summary(aov(AUC ~ strain, AUC_stats[grep("^10$", AUC_stats$dox),]))
summary(aov(AUC ~ strain, AUC_stats[grep("100", AUC_stats$dox),]))
summary(aov(AUC ~ strain, AUC_stats[grep("250", AUC_stats$dox),]))


### LME model for drug transformation by strain

# K. pneumoniae
summary(lmer(Dox_conc ~ Time_Point*Treatment + (1|Treatment),
             data = mono_data_stats[grep("K. pneum", mono_data_stats$Strain),]))
# E. coli
summary(lmer(Dox_conc ~ Time_Point*Treatment + (1|Treatment),
             data = mono_data_stats[grep("E. coli", mono_data_stats$Strain),]))

# remove time 0 since values were approximate
mono_data_stats = mono_data_stats[-c(grep("0", mono_data_stats$Time_Point)),]
mono_data_stats$Time_Point = factor(as.character(mono_data_stats$Time_Point),
                                    levels = c("3", "6", "12", "18", "24"))

# E. faecium
summary(lmer(Dox_conc ~ Time_Point*Treatment + (1|Treatment),
             data = mono_data_stats[grep("E. fae", mono_data_stats$Strain),]))
# E. faecium: dox reduction by treatment
summary(lm(Dox_conc ~ Time_Point,
           data = mono_data_stats[grep("E. fae", mono_data_stats$Strain),][grep("100", mono_data_stats[grep("E. fae", mono_data_stats$Strain),]$Treatment),]))
summary(lm(Dox_conc ~ Time_Point, 
           data = mono_data_stats[grep("E. fae", mono_data_stats$Strain),][grep("250", mono_data_stats[grep("E. fae", mono_data_stats$Strain),]$Treatment),]))
# C. innocuum
summary(lmer(Dox_conc ~ Time_Point*Treatment + (1|Treatment),
             data = mono_data_stats[grep("C. inn", mono_data_stats$Strain),]))
# L. rhamnosus
summary(lmer(Dox_conc ~ Time_Point*Treatment + (1|Treatment),
             data = mono_data_stats[grep("Lac", mono_data_stats$Strain),]))


###
### Aerobic growth conditions
###

## prep table
aerobe_stats = data.frame(strain = aerobic_growth$Strain, 
                          trt = aerobic_growth$Trt_char,
                          time = aerobic_growth$Time_Char, 
                          val = aerobic_growth$Dox_conc)
aerobe_stats$strain = factor(aerobe_stats$strain,
                             levels = c("neg_control", "EF", "EC", "KP"))

# anova: low trt
summary(aov(val ~ strain,
            aerobe_stats[grep("A", aerobe_stats$time),][grep("Low", aerobe_stats[grep("A", aerobe_stats$time),]$trt),]))
summary(aov(val ~ strain,
            aerobe_stats[grep("B", aerobe_stats$time),][grep("Low", aerobe_stats[grep("B", aerobe_stats$time),]$trt),]))

# anova: med trt
summary(aov(val ~ strain,
            aerobe_stats[grep("A", aerobe_stats$time),][grep("Med", aerobe_stats[grep("A", aerobe_stats$time),]$trt),]))
summary(aov(val ~ strain,
            aerobe_stats[grep("B", aerobe_stats$time),][grep("Med", aerobe_stats[grep("B", aerobe_stats$time),]$trt),]))
TukeyHSD(aov(val ~ strain,
             aerobe_stats[grep("B", aerobe_stats$time),][grep("Med", aerobe_stats[grep("B", aerobe_stats$time),]$trt),]))

# anova: high trt
summary(aov(val ~ strain,
            aerobe_stats[grep("A", aerobe_stats$time),][grep("Hi", aerobe_stats[grep("A", aerobe_stats$time),]$trt),]))
summary(aov(val ~ strain,
            aerobe_stats[grep("B", aerobe_stats$time),][grep("Hi", aerobe_stats[grep("B", aerobe_stats$time),]$trt),]))
TukeyHSD(aov(val ~ strain,
             aerobe_stats[grep("B", aerobe_stats$time),][grep("Hi", aerobe_stats[grep("B", aerobe_stats$time),]$trt),]))

## clean data for plotting
aerobic_growth = aerobic_growth[-c(which(is.na(aerobic_growth$Dox_avg_corr))),]
aerobic_growth$Strain = factor(aerobic_growth$Strain,
                               levels = c("neg_control", "EF", "EC", "KP"))

# low dox
aerobe_low = ggplot(aerobic_growth[grep("Low", aerobic_growth$Trt_char),], aes(x = Time_Char, y = Dox_avg_corr, fill = Strain)) +
  geom_bar(stat="identity", position=position_dodge2(preserve="single"),
           alpha = 0.8, col = "black") +
  geom_errorbar(stat="identity", position=position_dodge2(preserve="single"),
                aes(ymax=Dox_avg_corr+Dox_se, ymin=Dox_avg_corr-Dox_se), alpha = 0.8) +
  xlab("Incubation time (h)") +
  ylab("Dox µM") +
  theme_bw() +
  scale_fill_manual(values = c("gray", "orange", "skyblue", "darkblue")) +
  scale_x_discrete(labels = c("12", "24")) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) +
  theme(legend.position = "none")

# med dox
aerobe_med = ggplot(aerobic_growth[grep("Med", aerobic_growth$Trt_char),], aes(x = Time_Char, y = Dox_avg_corr, fill = Strain)) +
  geom_bar(stat="identity", position=position_dodge2(preserve="single"),
           alpha = 0.8, col = "black") +
  geom_errorbar(stat="identity", position=position_dodge2(preserve="single"),
                aes(ymax=Dox_avg_corr+Dox_se, ymin=Dox_avg_corr-Dox_se), alpha = 0.8) +
  xlab("Incubation time (h)") +
  ylab("Dox µM") +
  theme_bw() +
  scale_fill_manual(values = c("gray", "orange", "skyblue", "darkblue")) +
  scale_x_discrete(labels = c("12", "24")) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) +
  theme(legend.position = "none")

# high dox
aerobe_high = ggplot(aerobic_growth[grep("High", aerobic_growth$Trt_char),], aes(x = Time_Char, y = Dox_avg_corr, fill = Strain)) +
  geom_bar(stat="identity", position=position_dodge2(preserve="single"),
           alpha = 0.8, col = "black") +
  geom_errorbar(stat="identity", position=position_dodge2(preserve="single"),
                aes(ymax=Dox_avg_corr+Dox_se, ymin=Dox_avg_corr-Dox_se), alpha = 0.8) +
  xlab("Incubation time (h)") +
  ylab("Dox µM") +
  theme_bw() +
  scale_fill_manual(values = c("gray", "orange", "skyblue", "darkblue")) +
  scale_x_discrete(labels = c("12", "24")) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) +
  theme(legend.position = "none")

# plot (Fig. S3)
grid.arrange(aerobe_low, aerobe_med, aerobe_high, nrow=1)


###
### PABN
###

## treatment effects on each strain
summary(aov(OD600_corr ~ Group, data = pabn[grep("EC", pabn$Strain),]))
TukeyHSD(aov(OD600_corr ~ Group, data = pabn[grep("EC", pabn$Strain),]))

summary(aov(OD600_corr ~ Group, data = pabn[grep("EF", pabn$Strain),]))
TukeyHSD(aov(OD600_corr ~ Group, data = pabn[grep("EF", pabn$Strain),]))

summary(aov(OD600_corr ~ Group, data = pabn[grep("KP", pabn$Strain),]))
TukeyHSD(aov(OD600_corr ~ Group, data = pabn[grep("KP", pabn$Strain),]))

## clean data for plotting
pabn = pabn[-c(which(is.na(pabn$OD600_avg))),]
pabn = pabn[-c(grep("PABN2",pabn$Group)),]
pabn$Group = factor(pabn$Group, levels = c("Control", "Dox", "PABN", "Dox+PABN"))
pabn$Strain = factor(pabn$Strain, levels = c("KP", "EC", "EF"))

# plot (Fig. S4)
ggplot(pabn, aes(x = Strain, y = OD600_avg, fill = Group)) +
  geom_bar(stat="identity", position=position_dodge2(preserve="single"),
           alpha = 0.5, aes(col = Group)) +
  geom_errorbar(stat="identity", position=position_dodge2(preserve="single"),
                aes(ymax=OD600_avg+OD600_se, ymin=OD600_avg-OD600_se, col = Group), alpha = 0.8) +
  xlab("Taxon") +
  ylab("OD600") +
  theme_bw() +
  ylim(0, 580) +
  scale_x_discrete(labels = c("K. pneumoniae", "E. coli", "E. faecium")) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        axis.text.x = element_text(face = "italic"))


###
### Reporter strain growth in spent media of drug-resistant strains (Fig. 2)
###

## prep table
sm_data_avgs = sm_data[which(is.na(sm_data$OD600_corr_avg) == F),]

# set factor levels for strain and spent media sample
sm_data_avgs$Strain = factor(sm_data_avgs$Strain,
                             levels = c("C. innocuum", "L. rhamnosus", "E. faecium"))
sm_data_avgs$Spent_Media_Sample = factor(sm_data_avgs$Spent_Media_Sample,
                                         levels = c("KP", "KP+Dox", "EC", "EC+Dox", "EF", "EF+Dox"))

## plot (Fig. 2)
ggplot(sm_data_avgs[-c(grep("faecium",sm_data_avgs$Strain)),],
       aes(x=Strain, y=OD600_corr_avg/1000, color=Spent_Media_Sample, fill=Spent_Media_Sample)) +
  geom_bar(stat="identity", position=position_dodge2(preserve="single"),
           alpha = 0.5) + 
  geom_errorbar(stat="identity", position=position_dodge2(preserve="single"),
                aes(ymax=OD600_corr_avg/1000 + OD600_corr_se/1000,
                    ymin=OD600_corr_avg/1000 - OD600_corr_se/1000,
                    color = Spent_Media_Sample)) +
  theme_bw() +  
  ylab("OD600 at 24 h in 50% spent media") +
  xlab("Reporter strain") +
  ylim(0, 0.37) +
  scale_color_manual(values=c("darkblue", "darkblue", "skyblue", "skyblue", "orange", "orange")) +
  scale_fill_manual(values=c("darkblue", "pink", "skyblue", "pink", "orange", "pink")) +
  theme(axis.title=element_text(size = 14), 
        axis.text=element_text(size = 14),
        axis.text.x=element_text(face="italic"),
        legend.text=element_text(size = 13),
        legend.title=element_text(size = 13)) +
  theme(legend.position = "none")


### reporter strain stats for restored growth across different spent media

## C. innocuum

# ANOVA
summary(aov(OD600_corr ~ Spent_Media_Sample, data = sm_data[c(grep("innoc",sm_data$Strain)),]))
TukeyHSD(aov(OD600_corr ~ Spent_Media_Sample, data = sm_data[c(grep("innoc",sm_data$Strain)),]))
# Mann-Whitney: EF
wilcox.test(OD600_corr ~ Spent_Media_Sample, data = sm_data[grep("innoc",sm_data$Strain),][grep("EF", sm_data[c(grep("innoc",sm_data$Strain)),]$Spent_Media_Sample),])
# Mann-Whitney: EC
wilcox.test(OD600_corr ~ Spent_Media_Sample, data = sm_data[grep("innoc",sm_data$Strain),][grep("EC", sm_data[c(grep("innoc",sm_data$Strain)),]$Spent_Media_Sample),])
# Mann-Whitney: KP
wilcox.test(OD600_corr ~ Spent_Media_Sample, data = sm_data[grep("innoc",sm_data$Strain),][grep("KP", sm_data[c(grep("innoc",sm_data$Strain)),]$Spent_Media_Sample),])

## L. rhamnosus

# ANOVA
summary(aov(OD600_corr ~ Spent_Media_Sample, data = sm_data[c(grep("rham",sm_data$Strain)),]))
TukeyHSD(aov(OD600_corr ~ Spent_Media_Sample, data = sm_data[c(grep("rham",sm_data$Strain)),]))
# Mann-Whitney: EF
wilcox.test(OD600_corr ~ Spent_Media_Sample, data = sm_data[grep("rham",sm_data$Strain),][grep("EF", sm_data[c(grep("rham",sm_data$Strain)),]$Spent_Media_Sample),])
# Mann-Whitney: EC
wilcox.test(OD600_corr ~ Spent_Media_Sample, data = sm_data[grep("rham",sm_data$Strain),][grep("EC", sm_data[c(grep("rham",sm_data$Strain)),]$Spent_Media_Sample),])
# Mann-Whitney: KP
wilcox.test(OD600_corr ~ Spent_Media_Sample, data = sm_data[grep("rham",sm_data$Strain),][grep("KP", sm_data[c(grep("rham",sm_data$Strain)),]$Spent_Media_Sample),])



###############################
### Mixed-microbial communities
###############################

### load data

# inoculation concentrations
inoc <- read.table("mixed_community/incoculation_conc.txt", sep='\t',h=T,row.names=1,check=F,comment='')

## OD600 data
OD_raw = read.table('mixed_community/OD_data.txt', header = T, sep = '\t', row.names = 1, comment.char = '')

# 16S taxa table: raw data
taxa_raw = read.table('mixed_community/taxon_assignments_R1.txt', header = T, sep = '\t', row.names = 1, comment.char = '')


###
### Bacterial starting concentrations
###

## set factor 
inoc$Strain = factor(inoc$Strain,
                     levels = c("C_innocuum", "Lactobacillus", "E_faecalis", "E_coli",  "K_pneumoniae"))

## plot (Fig. S5 A)
ggplot(inoc[grep("GAM", inoc$Media),], aes(x = Strain, y =Avg, fill = Strain)) + 
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_errorbar(aes(ymax = Avg+SE,
                    ymin = Avg-SE),
                width = 0.2, position = position_dodge(width = 0.8)) +
  #ylab(bquote('log(CFU'~ml^-1)) +
  ylab(bquote('log( CFU'~ml^-1~')')) +
  scale_x_discrete(labels = c("C. innocuum", "L. rhamnosus", "E. faecalis", "E. coli", "K. pneumoniae")) +
  scale_fill_manual(values = c("darkred", "pink", "orange", "skyblue", "darkblue")) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, face="italic"),
        axis.ticks.x = element_blank(), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 0),
        legend.position = "none")

## ANOVA for different strain concentrations
summary(aov(Log_CFU_ml ~ Strain, inoc[grep("GAM", inoc$Media),]))


###
### 16S amplicon raw data pre-processing
###

## clean table

# separate taxon assignments from the dataset
taxa_sep = taxa_raw[,grep("D_", colnames(taxa_raw))]

# taxon summary stats
sum(apply(taxa_sep, 1, sum))
dim(taxa_sep)

# filter only target bacteria strains and combine counts by group
taxa_target = taxa_raw[,grep("Clostridium|Enterococcus|Escherichia|Enterobacteriaceae.__|Klebsiella|Lactobacillus", colnames(taxa_raw))]
colnames(taxa_target)
taxa_target_counts = data.frame(Clostridium = c(taxa_target[,3] + taxa_target[,4] + taxa_target[,5] + taxa_target[,6]),
                                Lactobacillus = taxa_target[,2],
                                Enterococcus = taxa_target[,1],
                                Escherichia = taxa_target[,7],
                                Klebsiella = taxa_target[,8],
                                counts_target = apply(taxa_sep[,c(grep("Clostridium|Enterococcus|Escherichia|Enterobacteriaceae.__|Klebsiella|Lactobacillus", colnames(taxa_raw)))], 1, sum),
                                counts_other = apply(taxa_sep[,-c(grep("Clostridium|Enterococcus|Escherichia|Enterobacteriaceae.__|Klebsiella|Lactobacillus", colnames(taxa_raw)))], 1, sum))
taxa_target_counts$total_counts = taxa_target_counts$counts_target+taxa_target_counts$counts_other
sum(taxa_target_counts$total_counts)
head(taxa_target_counts)

# add metadata
taxa_target_counts_meta = data.frame(taxa_target_counts,
                                     taxa_raw[,c(214:220)])
taxa_target_counts_meta = taxa_target_counts_meta[mixedorder(taxa_target_counts_meta$sample.code),]
colnames(taxa_target_counts_meta)
head(taxa_target_counts_meta)

# clean counts by "control-pool" samples (i.e., positives for R3 -- pool 2 on 3/27)
taxa_target_counts_meta[grep("ctrl", taxa_target_counts_meta$sample.code),]
taxa_target_counts_meta[grep("pool", taxa_target_counts_meta$sample.code),]

taxa_counts_clean = taxa_target_counts_meta
taxa_counts_clean$Enterococcus[grep("R3", taxa_target_counts_meta$sample.code)] = taxa_counts_clean$Enterococcus[grep("R3", taxa_target_counts_meta$sample.code)]-26
taxa_counts_clean$Enterococcus[grep("R3", taxa_target_counts_meta$sample.code)][which(taxa_counts_clean$Enterococcus[grep("R3", taxa_target_counts_meta$sample.code)]<0)] =
  rep(0, length(which(taxa_counts_clean$Enterococcus[grep("R3", taxa_target_counts_meta$sample.code)]<0)))
taxa_counts_clean$Escherichia[grep("R3", taxa_target_counts_meta$sample.code)] = taxa_counts_clean$Escherichia[grep("R3", taxa_target_counts_meta$sample.code)]-2
taxa_counts_clean$Escherichia[grep("R3", taxa_target_counts_meta$sample.code)][which(taxa_counts_clean$Escherichia[grep("R3", taxa_target_counts_meta$sample.code)]<0)] =
  rep(0, length(which(taxa_counts_clean$Escherichia[grep("R3", taxa_target_counts_meta$sample.code)]<0)))


### compute relative abundances of target strains by community

# Community 1 (C1)
C1_RA = matrix(ncol = 5, nrow = dim(taxa_counts_clean[grep("C1",taxa_counts_clean$microbial.community),])[1])
for (i in 1:dim(taxa_counts_clean[grep("C1",taxa_counts_clean$microbial.community),])[1]) {
  C1_RA[i,1] = taxa_counts_clean[grep("C1",taxa_counts_clean$microbial.community),][i,1]/sum(taxa_counts_clean[grep("C1",taxa_counts_clean$microbial.community),][i,1:5])
  C1_RA[i,2] = taxa_counts_clean[grep("C1",taxa_counts_clean$microbial.community),][i,2]/sum(taxa_counts_clean[grep("C1",taxa_counts_clean$microbial.community),][i,1:5])
  C1_RA[i,3] = taxa_counts_clean[grep("C1",taxa_counts_clean$microbial.community),][i,3]/sum(taxa_counts_clean[grep("C1",taxa_counts_clean$microbial.community),][i,1:5])
  C1_RA[i,4] = taxa_counts_clean[grep("C1",taxa_counts_clean$microbial.community),][i,4]/sum(taxa_counts_clean[grep("C1",taxa_counts_clean$microbial.community),][i,1:5])
  C1_RA[i,5] = taxa_counts_clean[grep("C1",taxa_counts_clean$microbial.community),][i,5]/sum(taxa_counts_clean[grep("C1",taxa_counts_clean$microbial.community),][i,1:5])
} 
C1_RA

# C2
C2_RA = matrix(ncol = 5, nrow = dim(taxa_counts_clean[grep("C2",taxa_counts_clean$microbial.community),])[1])
for (i in 1:dim(taxa_counts_clean[grep("C2",taxa_counts_clean$microbial.community),])[1]) {
  C2_RA[i,1] = taxa_counts_clean[grep("C2",taxa_counts_clean$microbial.community),][i,1]/sum(taxa_counts_clean[grep("C2",taxa_counts_clean$microbial.community),][i,c(1,2,4,5)])
  C2_RA[i,2] = taxa_counts_clean[grep("C2",taxa_counts_clean$microbial.community),][i,2]/sum(taxa_counts_clean[grep("C2",taxa_counts_clean$microbial.community),][i,c(1,2,4,5)])
  #3
  C2_RA[i,4] = taxa_counts_clean[grep("C2",taxa_counts_clean$microbial.community),][i,4]/sum(taxa_counts_clean[grep("C2",taxa_counts_clean$microbial.community),][i,c(1,2,4,5)])
  C2_RA[i,5] = taxa_counts_clean[grep("C2",taxa_counts_clean$microbial.community),][i,5]/sum(taxa_counts_clean[grep("C2",taxa_counts_clean$microbial.community),][i,c(1,2,4,5)])
} 
C2_RA[is.na(C2_RA)] <- 0
C2_RA

# C3
C3_RA = matrix(ncol = 5, nrow = dim(taxa_counts_clean[grep("C3",taxa_counts_clean$microbial.community),])[1])
for (i in 1:dim(taxa_counts_clean[grep("C3",taxa_counts_clean$microbial.community),])[1]) {
  C3_RA[i,1] = taxa_counts_clean[grep("C3",taxa_counts_clean$microbial.community),][i,1]/sum(taxa_counts_clean[grep("C3",taxa_counts_clean$microbial.community),][i,c(1,2,3,5)])
  C3_RA[i,2] = taxa_counts_clean[grep("C3",taxa_counts_clean$microbial.community),][i,2]/sum(taxa_counts_clean[grep("C3",taxa_counts_clean$microbial.community),][i,c(1,2,3,5)])
  C3_RA[i,3] = taxa_counts_clean[grep("C3",taxa_counts_clean$microbial.community),][i,3]/sum(taxa_counts_clean[grep("C3",taxa_counts_clean$microbial.community),][i,c(1,2,3,5)])
  #4
  C3_RA[i,5] = taxa_counts_clean[grep("C3",taxa_counts_clean$microbial.community),][i,5]/sum(taxa_counts_clean[grep("C3",taxa_counts_clean$microbial.community),][i,c(1,2,3,5)])
} 
C3_RA[is.na(C3_RA)] <- 0
C3_RA

# C4
C4_RA = matrix(ncol = 5, nrow = dim(taxa_counts_clean[grep("C4",taxa_counts_clean$microbial.community),])[1])
for (i in 1:dim(taxa_counts_clean[grep("C4",taxa_counts_clean$microbial.community),])[1]) {
  C4_RA[i,1] = taxa_counts_clean[grep("C4",taxa_counts_clean$microbial.community),][i,1]/sum(taxa_counts_clean[grep("C4",taxa_counts_clean$microbial.community),][i,c(1,2,3,4)])
  C4_RA[i,2] = taxa_counts_clean[grep("C4",taxa_counts_clean$microbial.community),][i,2]/sum(taxa_counts_clean[grep("C4",taxa_counts_clean$microbial.community),][i,c(1,2,3,4)])
  C4_RA[i,3] = taxa_counts_clean[grep("C4",taxa_counts_clean$microbial.community),][i,3]/sum(taxa_counts_clean[grep("C4",taxa_counts_clean$microbial.community),][i,c(1,2,3,4)])
  C4_RA[i,4] = taxa_counts_clean[grep("C4",taxa_counts_clean$microbial.community),][i,4]/sum(taxa_counts_clean[grep("C4",taxa_counts_clean$microbial.community),][i,c(1,2,3,4)])
  #5
} 
C4_RA[is.na(C4_RA)] <- 0
C4_RA

# C5
C5_RA = matrix(ncol = 5, nrow = dim(taxa_counts_clean[grep("C5",taxa_counts_clean$microbial.community),])[1])
for (i in 1:dim(taxa_counts_clean[grep("C5",taxa_counts_clean$microbial.community),])[1]) {
  C5_RA[i,1] = taxa_counts_clean[grep("C5",taxa_counts_clean$microbial.community),][i,1]/sum(taxa_counts_clean[grep("C5",taxa_counts_clean$microbial.community),][i,c(1,2,3)])
  C5_RA[i,2] = taxa_counts_clean[grep("C5",taxa_counts_clean$microbial.community),][i,2]/sum(taxa_counts_clean[grep("C5",taxa_counts_clean$microbial.community),][i,c(1,2,3)])
  C5_RA[i,3] = taxa_counts_clean[grep("C5",taxa_counts_clean$microbial.community),][i,3]/sum(taxa_counts_clean[grep("C5",taxa_counts_clean$microbial.community),][i,c(1,2,3)])
  #4
  #5
} 
C5_RA[is.na(C5_RA)] <- 0
C5_RA

## combine relative abundance matrices
rel_abund = rbind(C1_RA, C2_RA, C3_RA, C4_RA, C5_RA)
rel_abund = as.data.frame(rel_abund)
colnames(rel_abund) = c("RA_CI", "RA_Lac", "RA_EF", "RA_EC", "RA_KP")

# add to 'taxa_counts_clean'
taxa_RA = data.frame(taxa_counts_clean[grep("C1|C2|C3|C4|C5", taxa_counts_clean$microbial.community),],
                     rel_abund)
taxa_RA$sample.code = as.factor(as.character(taxa_RA$sample.code))
taxa_RA$microbial.community = as.character(taxa_RA$microbial.community)
taxa_RA$microbial.community = as.factor(taxa_RA$microbial.community)
taxa_RA$microbial.community

## summary stats

# total counts by sample
sum(taxa_RA$total_counts)
mean(taxa_RA$total_counts)
sd(taxa_RA$total_counts)/sqrt(595)

# target strains counts by sample
sum(taxa_RA$counts_target)
mean(taxa_RA$counts_target)
sd(taxa_RA$counts_target)/sqrt(595)

# proportion counts assigned to target strains
sum(taxa_RA$counts_target)/sum(taxa_RA$total_counts)

## re-compute relative abundances with adjustment for 'below detection limit' (BDL) when counts are 0
taxa_RA_new = taxa_RA
taxa_RA_new$RA_CI[which(taxa_RA_new$RA_CI == 0)] = (taxa_RA_new[which(taxa_RA_new$RA_CI == 0),1]+0.5)/taxa_RA_new[which(taxa_RA_new$RA_CI == 0),8]
taxa_RA_new$RA_Lac[which(taxa_RA_new$RA_Lac == 0)] = (taxa_RA_new[which(taxa_RA_new$RA_Lac == 0),2]+0.5)/taxa_RA_new[which(taxa_RA_new$RA_Lac == 0),8]
taxa_RA_new$RA_EF[which(taxa_RA_new$RA_EF == 0)] = (taxa_RA_new[which(taxa_RA_new$RA_EF == 0),3]+0.5)/taxa_RA_new[which(taxa_RA_new$RA_EF == 0),8]
taxa_RA_new$RA_EC[which(taxa_RA_new$RA_EC == 0)] = (taxa_RA_new[which(taxa_RA_new$RA_EC == 0),4]+0.5)/taxa_RA_new[which(taxa_RA_new$RA_EC == 0),8]
taxa_RA_new$RA_KP[which(taxa_RA_new$RA_KP == 0)] = (taxa_RA_new[which(taxa_RA_new$RA_KP == 0),5]+0.5)/taxa_RA_new[which(taxa_RA_new$RA_KP == 0),8]

# clean BDL data for target strains by community
taxa_RA_new[grep("C2",taxa_RA_new$microbial.community),]$RA_EF = rep(0, length(taxa_RA_new[grep("C2",taxa_RA_new$microbial.community),]$RA_EF))
taxa_RA_new[grep("C3",taxa_RA_new$microbial.community),]$RA_EC = rep(0, length(taxa_RA_new[grep("C3",taxa_RA_new$microbial.community),]$RA_EC))
taxa_RA_new[grep("C4",taxa_RA_new$microbial.community),]$RA_KP = rep(0, length(taxa_RA_new[grep("C4",taxa_RA_new$microbial.community),]$RA_KP))
taxa_RA_new[grep("C5",taxa_RA_new$microbial.community),]$RA_EC = rep(0, length(taxa_RA_new[grep("C5",taxa_RA_new$microbial.community),]$RA_EC))
taxa_RA_new[grep("C5",taxa_RA_new$microbial.community),]$RA_KP = rep(0, length(taxa_RA_new[grep("C5",taxa_RA_new$microbial.community),]$RA_KP))

# adjust to 100% RA by sample
taxa_RA_new = data.frame(taxa_RA_new[,1:15],
                         taxa_RA_new[,16:20]/apply(taxa_RA_new[,16:20], 1, sum))
apply(taxa_RA_new[,16:20], 1, sum)

## check exclusion info for duplicate amplicons by sample ("X")
taxa_RA_new[c(which(taxa_RA_new$exclusion.info == "X"), which(taxa_RA_new$exclusion.info == "X")-1),]

# combine technical replicates (where applicable) -- average amplicon data by sample code
taxa_RA_new$sample.code = as.character(taxa_RA_new$sample.code)
taxa_RA_new$sample.code[grep("X", taxa_RA_new$sample.code)] = substr(taxa_RA_new$sample.code[grep("X", taxa_RA_new$sample.code)], start = 1, stop = nchar(taxa_RA_new$sample.code[grep("X", taxa_RA_new$sample.code)])-1)
taxa_RA_new$sample.code = as.factor(taxa_RA_new$sample.code)

taxa_merge_code = data.frame(Clostridium = tapply(taxa_RA_new$RA_CI, taxa_RA_new$sample.code, mean),
                             Lactobacillus = tapply(taxa_RA_new$RA_Lac, taxa_RA_new$sample.code, mean),
                             Enterococcus = tapply(taxa_RA_new$RA_EF, taxa_RA_new$sample.code, mean),
                             Escherichia = tapply(taxa_RA_new$RA_EC, taxa_RA_new$sample.code, mean),
                             Klebsiella = tapply(taxa_RA_new$RA_KP, taxa_RA_new$sample.code, mean))

## re-compute RA for target strains
taxa_merge_rep = data.frame(CI_mean = tapply(taxa_merge_code[,1], substr(rownames(taxa_merge_code), start = 1, stop = nchar(rownames(taxa_merge_code)) - 3), mean),
                            Lac_mean = tapply(taxa_merge_code[,2], substr(rownames(taxa_merge_code), start = 1, stop = nchar(rownames(taxa_merge_code)) - 3), mean),
                            EF_mean = tapply(taxa_merge_code[,3], substr(rownames(taxa_merge_code), start = 1, stop = nchar(rownames(taxa_merge_code)) - 3), mean),
                            EC_mean = tapply(taxa_merge_code[,4], substr(rownames(taxa_merge_code), start = 1, stop = nchar(rownames(taxa_merge_code)) - 3), mean),
                            KP_mean = tapply(taxa_merge_code[,5], substr(rownames(taxa_merge_code), start = 1, stop = nchar(rownames(taxa_merge_code)) - 3), mean))

# re-order data
taxa_merge_rep = taxa_merge_rep[mixedorder(rownames(taxa_merge_rep)),]


###
### OD600 and 16S data analysis
###

## merge data
OD_16S = merge.data.frame(taxa_merge_code, OD_raw, by=0, all.x=TRUE)

# re-set negative OD values
OD_16S$OD600_corr[which(OD_16S$OD600_corr < 0)] = rep(0, length(which(OD_16S$OD600_corr < 0)))

## add values for "relative growth" (16S proportion * OD600)
OD_16S$CI_rel = OD_16S$Clostridium*OD_16S$OD600_corr
OD_16S$Lac_rel = OD_16S$Lactobacillus*OD_16S$OD600_corr
OD_16S$EF_rel = OD_16S$Enterococcus*OD_16S$OD600_corr
OD_16S$EC_rel = OD_16S$Escherichia*OD_16S$OD600_corr
OD_16S$KP_rel = OD_16S$Klebsiella*OD_16S$OD600_corr

## log transform "relative growth"
OD_16S$log_CI_rel = log(OD_16S$CI_rel+1)
OD_16S$log_Lac_rel = log(OD_16S$Lac_rel+1)
OD_16S$log_EF_rel = log(OD_16S$EF_rel+1)
OD_16S$log_EC_rel = log(OD_16S$EC_rel+1)
OD_16S$log_KP_rel = log(OD_16S$KP_rel+1)

## export OD_16S for editing based on qpcr vals (i.e., discount CI in G3 samples where it may have been diluted out)
# re-order
OD_16S = OD_16S[mixedorder(OD_16S$Row.names),]
head(OD_16S)
# write
#write.table(OD_16S, 'mixed_community/OD_16S_updated.txt', sep="\t", col.names=NA, quote = FALSE)

## read-in corrected table; set as OD_16S
OD_16S_corr = read.table('mixed_community/OD_16S_updated_edit.txt', header = T, sep = '\t', row.names = 1, comment.char = '')

# re-name
OD_16S = OD_16S_corr
OD_16S$Row.names = as.character(OD_16S$Row.names)

## compute summary stats
OD_16S_rep = data.frame(CI_mean = tapply(OD_16S$Clostridium, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), function(x) {mean(x, na.rm = T)}),
                        Lac_mean  = tapply(OD_16S$Lactobacillus, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), mean),
                        EF_mean  = tapply(OD_16S$Enterococcus, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), mean),
                        EC_mean  = tapply(OD_16S$Escherichia, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), mean),
                        KP_mean  = tapply(OD_16S$Klebsiella, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), mean),
                        
                        CI_rel_mean  = tapply(OD_16S$CI_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), function(x) {mean(x, na.rm = T)}),
                        Lac_rel_mean  = tapply(OD_16S$Lac_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), mean),
                        EF_rel_mean  = tapply(OD_16S$EF_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), mean),
                        EC_rel_mean  = tapply(OD_16S$EC_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), mean),
                        KP_rel_mean  = tapply(OD_16S$KP_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), mean),
                        
                        log_CI_rel_mean  = tapply(OD_16S$log_CI_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), function(x) {mean(x, na.rm = T)}),
                        log_Lac_rel_mean  = tapply(OD_16S$log_Lac_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), mean),
                        log_EF_rel_mean  = tapply(OD_16S$log_EF_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), mean),
                        log_EC_rel_mean  = tapply(OD_16S$log_EC_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), mean),
                        log_KP_rel_mean  = tapply(OD_16S$log_KP_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), mean),
                        
                        OD_600_mean  = tapply(OD_16S$OD600_corr, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), mean),
                        Dox_mean = tapply(OD_16S$Dox_R1_corr, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), mean),
                        Time_mean = tapply(OD_16S$Time, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), mean),
                        
                        CI_sd = tapply(OD_16S$Clostridium, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), function(x) {mean(x, na.rm = T)}),
                        Lac_sd  = tapply(OD_16S$Lactobacillus, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), sd),
                        EF_sd  = tapply(OD_16S$Enterococcus, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), sd),
                        EC_sd  = tapply(OD_16S$Escherichia, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), sd),
                        KP_sd  = tapply(OD_16S$Klebsiella, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), sd),
                        
                        CI_rel_sd  = tapply(OD_16S$CI_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), function(x) {mean(x, na.rm = T)}),
                        Lac_rel_sd  = tapply(OD_16S$Lac_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), sd),
                        EF_rel_sd  = tapply(OD_16S$EF_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), sd),
                        EC_rel_sd  = tapply(OD_16S$EC_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), sd),
                        KP_rel_sd  = tapply(OD_16S$KP_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), sd),
                        
                        log_CI_rel_sd  = tapply(OD_16S$log_CI_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), function(x) {mean(x, na.rm = T)}),
                        log_Lac_rel_sd  = tapply(OD_16S$log_Lac_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), sd),
                        log_EF_rel_sd  = tapply(OD_16S$log_EF_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), sd),
                        log_EC_rel_sd  = tapply(OD_16S$log_EC_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), sd),
                        log_KP_rel_sd  = tapply(OD_16S$log_KP_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), sd),
                        
                        log_CI_rel_se  = tapply(OD_16S$log_CI_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), function(x) {sd(x, na.rm = T)})/
                          sqrt(tapply(OD_16S$Clostridium, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), function(x) {length(which(x > 0))})),
                        log_Lac_rel_se  = tapply(OD_16S$log_Lac_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), sd)/
                          sqrt(tapply(OD_16S$Lac_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), length)),
                        log_EF_rel_se  = tapply(OD_16S$log_EF_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), sd)/
                          sqrt(tapply(OD_16S$EF_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), length)),
                        log_EC_rel_se  = tapply(OD_16S$log_EC_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), sd)/
                          sqrt(tapply(OD_16S$EC_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), length)),
                        log_KP_rel_se  = tapply(OD_16S$log_KP_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), sd)/
                          sqrt(tapply(OD_16S$KP_rel, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), length)),
                        
                        OD_600_sd  = tapply(OD_16S$OD600_corr, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), sd),
                        Dox_sd = tapply(OD_16S$Dox_R1_corr, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), sd),
                        
                        OD_600_se  = tapply(OD_16S$OD600_corr, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), sd)/
                          sqrt(tapply(OD_16S$OD600_corr, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), length)),
                        Dox_se = tapply(OD_16S$Dox_R1_corr, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), sd)/
                          sqrt(tapply(OD_16S$Dox_R1_corr, substr(OD_16S$Row.names, start = 1, stop = nchar(OD_16S$Row.names) - 3), length)))

## re-order
OD_16S_rep = OD_16S_rep[mixedorder(rownames(OD_16S_rep)),]
head(OD_16S_rep)
dim(OD_16S_rep)

## write-out table
#write.table(OD_16S_rep, 'mixed_community/OD_16S_updated_edit_rep.txt', sep="\t", col.names=NA, quote = FALSE)

## "melt" for plotting
OD_16S_melt = melt(t(OD_16S_rep[,grep("log",colnames(OD_16S_rep))]))
OD_16S_melt = data.frame(OD_16S_melt[grep("mean", OD_16S_melt$Var1),],
                         OD_16S_melt[grep("se", OD_16S_melt$Var1),])
OD_16S_melt$Var2 == OD_16S_melt$Var2.1

## clean table
head(OD_16S_melt)
OD_16S_melt = OD_16S_melt[,-c(4:5)]
colnames(OD_16S_melt) = c("taxa", "sample", "mean", "se")

## split metadata
OD_16S_melt = data.frame(OD_16S_melt,
                         as.data.frame(str_split_fixed(as.character(OD_16S_melt$sample), "-", 4)),
                         Time = rep(OD_16S_rep$Time, each = 5))
#OD_16S_melt = OD_16S_melt[-c(which(is.na(OD_16S_melt$mean) == TRUE)),]
head(OD_16S_melt)
dim(OD_16S_melt)

## call G1 (no dox) to control treatment
OD_16S_melt$V3[which(OD_16S_melt$V3 == "NA")] = rep("T1", length(which(OD_16S_melt$V3 == "NA")))

## write-out table for manual edit before plotting
#write.table(OD_16S_melt, 'mixed_community/OD_16S_updated_edit_rep_melt.txt', sep="\t", col.names=NA, quote = FALSE)


###
### Plot relative bacterial growth within communities by treatment group
###

## read-in corrected table
OD_16S_melt_clean = read.table('mixed_community/OD_16S_updated_edit_rep_melt_clean.txt', header = T, sep = '\t', row.names = 1, comment.char = '')

# reset factor for pollting
OD_16S_melt_clean$taxa = as.character(OD_16S_melt_clean$taxa)
OD_16S_melt_clean$taxa = factor(OD_16S_melt_clean$taxa,
                                levels = c("log_CI_rel_mean", "log_Lac_rel_mean", "log_EF_rel_mean", "log_EC_rel_mean", "log_KP_rel_mean")) 

## prep plot panels

C1 = ggplot(OD_16S_melt_clean[grep("C1", OD_16S_melt_clean$sample),],  
            aes(x = Time, y = mean, fill = taxa, col = taxa, shape = V3, linetype = V3)) +
  geom_point(size = 2.5, col = "black", alpha = 0.7) +
  geom_line(size = 0.5, alpha = 0.6) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se,
                    colour = taxa), width = 1, alpha = 0.5) +
  scale_shape_manual(values = c(21, 24, 23, 22),
                     labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  scale_fill_manual(values = c("darkred", "pink", "orange", "skyblue", "darkblue", "darkgray"),
                    labels = c("Clostridium (S)", "Lactobacillus (S)", "Enterococcus (R)", "Escherichia (R*)", "Klebsiella (R*)", "Doxo uM")) +
  scale_colour_manual(values = c("darkred", "pink", "orange", "skyblue", "darkblue", "darkgray"),
                      labels = c("Clostridium (S)", "Lactobacillus (S)", "Enterococcus (R)", "Escherichia (R*)", "Klebsiella (R*)", "Doxo uM")) +
  scale_linetype_manual(values = c(1, 2, 3, 4),
                        labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  theme_bw() +
  xlab("Cumulative time (h)") +
  ylab("log(RA*OD600)") +
  scale_x_continuous(breaks = c(0, 12, 24, 36, 48, 60, 72)) +
  #facet_grid(.~V2, scale="free_x", space="free") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = "none")

C2 = ggplot(OD_16S_melt_clean[intersect(grep("C2", OD_16S_melt_clean$sample), grep("CI|Lac|EC|KP|dox", OD_16S_melt_clean$taxa)),],  
            aes(x = Time, y = mean, fill = taxa, col = taxa, shape = V3, linetype = V3)) +
  geom_point(size = 2.5, col = "black", alpha = 0.7) +
  geom_line(size = 0.5, alpha = 0.6) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se,
                    colour = taxa), width = 1, alpha = 0.5) +
  scale_shape_manual(values = c(21, 24, 23, 22),
                     labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  scale_fill_manual(values = c("darkred", "pink", "skyblue", "darkblue", "darkgray"),
                    labels = c("Clostridium (S)", "Lactobacillus (S)", "Escherichia (R*)", "Klebsiella (R*)", "Doxo uM")) +
  scale_colour_manual(values = c("darkred", "pink", "skyblue", "darkblue", "darkgray"),
                      labels = c("Clostridium (S)", "Lactobacillus (S)", "Escherichia (R*)", "Klebsiella (R*)", "Doxo uM")) +
  scale_linetype_manual(values = c(1, 2, 3, 4),
                        labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  theme_bw() +
  xlab("Cumulative time (h)") +
  ylab("log(RA*OD600)") +
  scale_x_continuous(breaks = c(0, 12, 24, 36, 48, 60, 72)) +
  #facet_grid(.~V2, scale="free_x", space="free") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = "none")

C3 = ggplot(OD_16S_melt_clean[intersect(grep("C3", OD_16S_melt_clean$sample), grep("CI|Lac|EF|KP|dox", OD_16S_melt_clean$taxa)),],  
            aes(x = Time, y = mean, fill = taxa, col = taxa, shape = V3, linetype = V3)) +
  geom_point(size = 2.5, col = "black", alpha = 0.7) +
  geom_line(size = 0.5, alpha = 0.6) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se,
                    colour = taxa), width = 1, alpha = 0.5) +
  scale_shape_manual(values = c(21, 24, 23, 22),
                     labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  scale_fill_manual(values = c("darkred", "pink", "orange", "darkblue", "darkgray"),
                    labels = c("Clostridium (S)", "Lactobacillus (S)", "Enterococcus (R)", "Klebsiella (R*)", "Doxo uM")) +
  scale_colour_manual(values = c("darkred", "pink", "orange", "darkblue", "darkgray"),
                      labels = c("Clostridium (S)", "Lactobacillus (S)", "Enterococcus (R)", "Klebsiella (R*)", "Doxo uM")) +
  scale_linetype_manual(values = c(1, 2, 3, 4),
                        labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  theme_bw() +
  xlab("Cumulative time (h)") +
  ylab("log(RA*OD600)") +
  scale_x_continuous(breaks = c(0, 12, 24, 36, 48, 60, 72)) +
  #facet_grid(.~V2, scale="free_x", space="free") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = "none")

C4 = ggplot(OD_16S_melt_clean[intersect(grep("C4", OD_16S_melt_clean$sample), grep("CI|Lac|EF|EC|dox", OD_16S_melt_clean$taxa)),],  
            aes(x = Time, y = mean, fill = taxa, col = taxa, shape = V3, linetype = V3)) +
  geom_point(size = 2.5, col = "black", alpha = 0.7) +
  geom_line(size = 0.5, alpha = 0.6) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se,
                    colour = taxa), width = 1, alpha = 0.5) +
  scale_shape_manual(values = c(21, 24, 23, 22),
                     labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  scale_fill_manual(values = c("darkred", "pink", "orange", "darkblue","darkgray"),
                    labels = c("Clostridium (S)", "Lactobacillus (S)", "Enterococcus (R)", "Escherichia (R*)", "Doxo uM")) +
  scale_colour_manual(values = c("darkred", "pink", "orange", "skyblue", "darkgray"),
                      labels = c("Clostridium (S)", "Lactobacillus (S)", "Enterococcus (R)", "Escherichia (R*)", "Doxo uM")) +
  scale_linetype_manual(values = c(1, 2, 3, 4),
                        labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  theme_bw() +
  xlab("Cumulative time (h)") +
  ylab("log(RA*OD600)") +
  scale_x_continuous(breaks = c(0, 12, 24, 36, 48, 60, 72)) +
  #facet_grid(.~V2, scale="free_x", space="free") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = "none")

C5 = ggplot(OD_16S_melt_clean[intersect(grep("C5", OD_16S_melt_clean$sample), grep("CI|Lac|EF|dox", OD_16S_melt_clean$taxa)),],  
            aes(x = Time, y = mean, fill = taxa, col = taxa, shape = V3, linetype = V3)) +
  geom_point(size = 2.5, col = "black", alpha = 0.7) +
  geom_line(size = 0.5, alpha = 0.6) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se,
                    colour = taxa), width = 1, alpha = 0.5) +
  scale_shape_manual(values = c(21, 24, 23, 22),
                     labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  scale_fill_manual(values = c("darkred", "pink", "orange", "darkgray"),
                    labels = c("Clostridium (S)", "Lactobacillus (S)", "Enterococcus (R)", "Doxo uM")) +
  scale_colour_manual(values = c("darkred", "pink", "orange", "darkgray"),
                      labels = c("Clostridium (S)", "Lactobacillus (S)", "Enterococcus (R)", "Doxo uM")) +
  scale_linetype_manual(values = c(1, 2, 3, 4),
                        labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  theme_bw() +
  xlab("Cumulative time (h)") +
  ylab("log(RA*OD600)") +
  scale_x_continuous(breaks = c(0, 12, 24, 36, 48, 60, 72)) +
  #facet_grid(.~V2, scale="free_x", space="free") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = "none")

## plot all panels (Fig. S6)
grid.arrange(C1, C2, C3, C4, C5, nrow = 2)


###
### Plot without low treatment (clean figure)
###

C1 = ggplot(OD_16S_melt_clean[grep("C1", OD_16S_melt_clean$sample),][-c(grep("T2", OD_16S_melt_clean[grep("C1", OD_16S_melt_clean$sample),]$V3)),],  
            aes(x = Time, y = mean, fill = taxa, col = taxa, shape = V3, linetype = V3)) +
  geom_point(size = 2.5, col = "black", alpha = 0.7) +
  geom_line(size = 0.5, alpha = 0.6) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se,
                    colour = taxa), width = 1, alpha = 0.5) +
  scale_shape_manual(values = c(21, 23, 22, 24),
                     labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  scale_fill_manual(values = c("darkred", "pink", "orange", "skyblue", "darkblue", "darkgray"),
                    labels = c("Clostridium (S)", "Lactobacillus (S)", "Enterococcus (R)", "Escherichia (R*)", "Klebsiella (R*)", "Doxo uM")) +
  scale_colour_manual(values = c("darkred", "pink", "orange", "skyblue", "darkblue", "darkgray"),
                      labels = c("Clostridium (S)", "Lactobacillus (S)", "Enterococcus (R)", "Escherichia (R*)", "Klebsiella (R*)", "Doxo uM")) +
  scale_linetype_manual(values = c(1, 3, 4, 2),
                        labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  theme_bw() +
  xlab("Cumulative time (h)") +
  ylab("log(RA*OD600)") +
  scale_x_continuous(breaks = c(0, 12, 24, 36, 48, 60, 72)) +
  #facet_grid(.~V2, scale="free_x", space="free") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = "none")

C2 = ggplot(OD_16S_melt_clean[intersect(grep("C2", OD_16S_melt_clean$sample), grep("CI|Lac|EC|KP|dox", OD_16S_melt_clean$taxa)),][-c(grep("T2", OD_16S_melt_clean[intersect(grep("C2", OD_16S_melt_clean$sample), grep("CI|Lac|EC|KP|dox", OD_16S_melt_clean$taxa)),]$V3)),],
            aes(x = Time, y = mean, fill = taxa, col = taxa, shape = V3, linetype = V3)) +
  geom_point(size = 2.5, col = "black", alpha = 0.7) +
  geom_line(size = 0.5, alpha = 0.6) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se,
                    colour = taxa), width = 1, alpha = 0.5) +
  scale_shape_manual(values = c(21, 23, 22, 24),
                     labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  scale_fill_manual(values = c("darkred", "pink", "skyblue", "darkblue", "darkgray"),
                    labels = c("Clostridium (S)", "Lactobacillus (S)", "Escherichia (R*)", "Klebsiella (R*)", "Doxo uM")) +
  scale_colour_manual(values = c("darkred", "pink", "skyblue", "darkblue", "darkgray"),
                      labels = c("Clostridium (S)", "Lactobacillus (S)", "Escherichia (R*)", "Klebsiella (R*)", "Doxo uM")) +
  scale_linetype_manual(values = c(1, 3, 4, 2),
                        labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  theme_bw() +
  xlab("Cumulative time (h)") +
  ylab("log(RA*OD600)") +
  scale_x_continuous(breaks = c(0, 12, 24, 36, 48, 60, 72)) +
  #facet_grid(.~V2, scale="free_x", space="free") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = "none")

C3 = ggplot(OD_16S_melt_clean[intersect(grep("C3", OD_16S_melt_clean$sample), grep("CI|Lac|EF|KP|dox", OD_16S_melt_clean$taxa)),][-c(grep("T2", OD_16S_melt_clean[intersect(grep("C3", OD_16S_melt_clean$sample), grep("CI|Lac|EF|KP|dox", OD_16S_melt_clean$taxa)),]$V3)),],  
            aes(x = Time, y = mean, fill = taxa, col = taxa, shape = V3, linetype = V3)) +
  geom_point(size = 2.5, col = "black", alpha = 0.7) +
  geom_line(size = 0.5, alpha = 0.6) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se,
                    colour = taxa), width = 1, alpha = 0.5) +
  scale_shape_manual(values = c(21, 23, 22, 24),
                     labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  scale_fill_manual(values = c("darkred", "pink", "orange", "darkblue", "darkgray"),
                    labels = c("Clostridium (S)", "Lactobacillus (S)", "Enterococcus (R)", "Klebsiella (R*)", "Doxo uM")) +
  scale_colour_manual(values = c("darkred", "pink", "orange", "darkblue", "darkgray"),
                      labels = c("Clostridium (S)", "Lactobacillus (S)", "Enterococcus (R)", "Klebsiella (R*)", "Doxo uM")) +
  scale_linetype_manual(values = c(1, 3, 4, 2),
                        labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  theme_bw() +
  xlab("Cumulative time (h)") +
  ylab("log(RA*OD600)") +
  scale_x_continuous(breaks = c(0, 12, 24, 36, 48, 60, 72)) +
  #facet_grid(.~V2, scale="free_x", space="free") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = "none")

C4 = ggplot(OD_16S_melt_clean[intersect(grep("C4", OD_16S_melt_clean$sample), grep("CI|Lac|EF|EC|dox", OD_16S_melt_clean$taxa)),][-c(grep("T2", OD_16S_melt_clean[intersect(grep("C4", OD_16S_melt_clean$sample), grep("CI|Lac|EF|EC|dox", OD_16S_melt_clean$taxa)),]$V3)),],  
            aes(x = Time, y = mean, fill = taxa, col = taxa, shape = V3, linetype = V3)) +
  geom_point(size = 2.5, col = "black", alpha = 0.7) +
  geom_line(size = 0.5, alpha = 0.6) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se,
                    colour = taxa), width = 1, alpha = 0.5) +
  scale_shape_manual(values = c(21, 23, 22, 24),
                     labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  scale_fill_manual(values = c("darkred", "pink", "orange", "darkblue","darkgray"),
                    labels = c("Clostridium (S)", "Lactobacillus (S)", "Enterococcus (R)", "Escherichia (R*)", "Doxo uM")) +
  scale_colour_manual(values = c("darkred", "pink", "orange", "skyblue", "darkgray"),
                      labels = c("Clostridium (S)", "Lactobacillus (S)", "Enterococcus (R)", "Escherichia (R*)", "Doxo uM")) +
  scale_linetype_manual(values = c(1, 3, 4, 2),
                        labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  theme_bw() +
  xlab("Cumulative time (h)") +
  ylab("log(RA*OD600)") +
  scale_x_continuous(breaks = c(0, 12, 24, 36, 48, 60, 72)) +
  #facet_grid(.~V2, scale="free_x", space="free") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = "none")

C5 = ggplot(OD_16S_melt_clean[intersect(grep("C5", OD_16S_melt_clean$sample), grep("CI|Lac|EF|dox", OD_16S_melt_clean$taxa)),][-c(grep("T2", OD_16S_melt_clean[intersect(grep("C5", OD_16S_melt_clean$sample), grep("CI|Lac|EF|dox", OD_16S_melt_clean$taxa)),]$V3)),],  
            aes(x = Time, y = mean, fill = taxa, col = taxa, shape = V3, linetype = V3)) +
  geom_point(size = 2.5, col = "black", alpha = 0.7) +
  geom_line(size = 0.5, alpha = 0.6) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se,
                    colour = taxa), width = 1, alpha = 0.5) +
  scale_shape_manual(values = c(21, 23, 22, 24),
                     labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  scale_fill_manual(values = c("darkred", "pink", "orange", "darkgray"),
                    labels = c("Clostridium (S)", "Lactobacillus (S)", "Enterococcus (R)", "Doxo uM")) +
  scale_colour_manual(values = c("darkred", "pink", "orange", "darkgray"),
                      labels = c("Clostridium (S)", "Lactobacillus (S)", "Enterococcus (R)", "Doxo uM")) +
  scale_linetype_manual(values = c(1, 3, 4, 2),
                        labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  theme_bw() +
  xlab("Cumulative time (h)") +
  ylab("log(RA*OD600)") +
  scale_x_continuous(breaks = c(0, 12, 24, 36, 48, 60, 72)) +
  #facet_grid(.~V2, scale="free_x", space="free") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = "none")

## plot all panels (Fig. 5)
grid.arrange(C1, C2, C3, C4, C5, nrow = 2)


###
### Taxa stats: relative growth
###

taxa_stats = data.frame(OD_16S_corr,
                        as.data.frame(str_split_fixed(as.character(OD_16S_corr$Row.names), "-", 5)))
taxa_stats$V4 = factor(taxa_stats$V4,
                       levels = c("TP0", "TP3", "TP6", "TP12", "TP24"))

## add factor to differentiate controls from all treatment groups (for CI and LR stats)
test = taxa_stats
test$V6 = rep("Low", length(test$V5))
test$V6[grep("T2|T3|T4", test$V3)] = rep("High", length(test$V6[grep("T2|T3|T4", test$V3)]))

## Effects of doxorubicin (any level) on relative growth of CI/LR by community: LME

# CI-G2
summary(lmer(log_CI_rel ~ V6*V4 + (V4/V5) + (1|V6), data = test[grep("C1-G2", test$Row.names),]))
summary(lmer(log_CI_rel ~ V6*V4 + (V4/V5) + (1|V6), data = test[grep("C2-G2", test$Row.names),]))
summary(lmer(log_CI_rel ~ V6*V4 + (V4/V5) + (1|V6), data = test[grep("C3-G2", test$Row.names),]))
summary(lmer(log_CI_rel ~ V6*V4 + (V4/V5) + (1|V6), data = test[grep("C4-G2", test$Row.names),]))
summary(lmer(log_CI_rel ~ V6*V4 + (V4/V5) + (1|V6), data = test[grep("C5-G2", test$Row.names),]))

# CI-G3
summary(lmer(log_CI_rel ~ V6*V4 + (V4/V5) + (1|V6), data = test[grep("C1-G3", test$Row.names),]))
summary(lmer(log_CI_rel ~ V6*V4 + (V4/V5) + (1|V6), data = test[grep("C2-G3", test$Row.names),]))
summary(lmer(log_CI_rel ~ V6*V4 + (V4/V5) + (1|V6), data = test[grep("C3-G3", test$Row.names),]))
summary(lmer(log_CI_rel ~ V6*V4 + (V4/V5) + (1|V6), data = test[grep("C4-G3", test$Row.names),]))
summary(lmer(log_CI_rel ~ V6*V4 + (V4/V5) + (1|V6), data = test[grep("C5-G3", test$Row.names),]))

## LR-G2
summary(lmer(log_Lac_rel ~ V6*V4 + (V4/V5) + (1|V6), data = test[grep("C1-G2", test$Row.names),]))
summary(lmer(log_Lac_rel ~ V6*V4 + (V4/V5) + (1|V6), data = test[grep("C2-G2", test$Row.names),]))
summary(lmer(log_Lac_rel ~ V6*V4 + (V4/V5) + (1|V6), data = test[grep("C3-G2", test$Row.names),]))
summary(lmer(log_Lac_rel ~ V6*V4 + (V4/V5) + (1|V6), data = test[grep("C4-G2", test$Row.names),]))
summary(lmer(log_Lac_rel ~ V6*V4 + (V4/V5) + (1|V6), data = test[grep("C5-G2", test$Row.names),]))

# LR-G3
summary(lmer(log_Lac_rel ~ V6*V4 + (V4/V5) + (1|V6), data = test[grep("C1-G3", test$Row.names),]))
summary(lmer(log_Lac_rel ~ V6*V4 + (V4/V5) + (1|V6), data = test[grep("C2-G3", test$Row.names),]))
summary(lmer(log_Lac_rel ~ V6*V4 + (V4/V5) + (1|V6), data = test[grep("C3-G3", test$Row.names),]))
summary(lmer(log_Lac_rel ~ V6*V4 + (V4/V5) + (1|V6), data = test[grep("C4-G3", test$Row.names),]))
summary(lmer(log_Lac_rel ~ V6*V4 + (V4/V5) + (1|V6), data = test[grep("C5-G3", test$Row.names),]))


## Effects of treatment (levels differentiated) on relative growth of individual strains by community: LME

# EF-G3
summary(lmer(log_EF_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C1-G3", taxa_stats$Row.names),]))
summary(lmer(log_EF_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C3-G3", taxa_stats$Row.names),]))
summary(lmer(log_EF_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C4-G3", taxa_stats$Row.names),]))
summary(lmer(log_EF_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C5-G3", taxa_stats$Row.names),]))

# CI-G2
summary(lmer(log_CI_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C1-G2", taxa_stats$Row.names),]))
summary(lmer(log_CI_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C2-G2", taxa_stats$Row.names),]))
summary(lmer(log_CI_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C3-G2", taxa_stats$Row.names),]))
summary(lmer(log_CI_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C4-G2", taxa_stats$Row.names),]))
summary(lmer(log_CI_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C5-G2", taxa_stats$Row.names),]))

# CI-G3
summary(lmer(log_CI_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C1-G3", taxa_stats$Row.names),]))
summary(lmer(log_CI_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C2-G3", taxa_stats$Row.names),]))
summary(lmer(log_CI_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C3-G3", taxa_stats$Row.names),]))
summary(lmer(log_CI_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C4-G3", taxa_stats$Row.names),]))
summary(lmer(log_CI_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C5-G3", taxa_stats$Row.names),]))

# Lac-G2
summary(lmer(log_Lac_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C1-G2", taxa_stats$Row.names),]))
summary(lmer(log_Lac_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C2-G2", taxa_stats$Row.names),]))
summary(lmer(log_Lac_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C3-G2", taxa_stats$Row.names),]))
summary(lmer(log_Lac_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C4-G2", taxa_stats$Row.names),]))
summary(lmer(log_Lac_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C5-G2", taxa_stats$Row.names),]))

# Lac-G3
summary(lmer(log_Lac_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C1-G3", taxa_stats$Row.names),]))
summary(lmer(log_Lac_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C2-G3", taxa_stats$Row.names),]))
summary(lmer(log_Lac_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C3-G3", taxa_stats$Row.names),]))
summary(lmer(log_Lac_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C4-G3", taxa_stats$Row.names),]))
summary(lmer(log_Lac_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C5-G3", taxa_stats$Row.names),]))


## Effects of treatment on relative growth of members in each community (C1-C5): LME

#C1
summary(lmer(log_KP_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C1-G2", taxa_stats$Row.names),]))
summary(lmer(log_EC_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C1-G2", taxa_stats$Row.names),]))
summary(lmer(log_EF_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C1-G2", taxa_stats$Row.names),]))
summary(lmer(log_Lac_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C1-G2", taxa_stats$Row.names),]))
summary(lmer(log_CI_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C1-G2", taxa_stats$Row.names),]))

#C2
summary(lmer(log_KP_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C2-G2", taxa_stats$Row.names),]))
summary(lmer(log_EC_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C2-G2", taxa_stats$Row.names),]))
summary(lmer(log_Lac_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C2-G2", taxa_stats$Row.names),]))
summary(lmer(log_CI_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C2-G2", taxa_stats$Row.names),]))

#C3
summary(lmer(log_KP_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C3-G2", taxa_stats$Row.names),]))
summary(lmer(log_EF_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C3-G2", taxa_stats$Row.names),]))
summary(lmer(log_Lac_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C3-G2", taxa_stats$Row.names),]))
summary(lmer(log_CI_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C3-G2", taxa_stats$Row.names),]))

#C4
summary(lmer(log_EC_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C4-G2", taxa_stats$Row.names),]))
summary(lmer(log_EF_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C4-G2", taxa_stats$Row.names),]))
summary(lmer(log_Lac_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C4-G2", taxa_stats$Row.names),]))
summary(lmer(log_CI_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C4-G2", taxa_stats$Row.names),]))

#C5
summary(lmer(log_EF_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C5-G2", taxa_stats$Row.names),]))
summary(lmer(log_Lac_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C5-G2", taxa_stats$Row.names),]))
summary(lmer(log_CI_rel ~ V3*V4 + (V4/V5) + (1|V3), data =  taxa_stats[grep("C5-G2", taxa_stats$Row.names),]))


###
### OD600 and Dox data prep
###

## prep data
colnames(OD_raw)[1] = c("Row.names")
OD_raw$Row.names = as.character(OD_raw$Row.names)

# set df with summary stats
OD_rep = data.frame(OD_600_mean  = tapply(OD_raw$OD600_corr, OD_raw$Row.names, function(x) {mean(x, na.rm = T)}),
                    OD_600_sd  = tapply(OD_raw$OD600_corr, OD_raw$Row.names, function(x) {sd(x, na.rm = T)}),
                    OD_600_se  = tapply(OD_raw$OD600_corr, OD_raw$Row.names, function(x) {sd(x, na.rm = T)})/sqrt(3),
                    
                    Dox_mean = tapply(OD_raw$Dox_R1_corr, OD_raw$Row.names, function(x) {mean(x, na.rm = T)}),
                    Dox_sd = tapply(OD_raw$Dox_R1_corr, OD_raw$Row.names, function(x) {sd(x, na.rm = T)}),
                    Dox_se = tapply(OD_raw$Dox_R1_corr, OD_raw$Row.names, function(x) {sd(x, na.rm = T)})/sqrt(3),
                    
                    Time_mean = tapply(OD_raw$Time, OD_raw$Row.names, mean))

## write-out df for manual edit before plotting
#write.table(OD_rep, 'mixed_community/OD_data_rep.txt', sep="\t", col.names=NA, quote = FALSE)


###
### OD600 stats
###

## OD stats
OD_stats = data.frame(OD_raw, 
                      as.data.frame(str_split_fixed(as.character(OD_raw$Row.names), "-", 4)))
head(OD_stats) 

# set time as factor
OD_stats$V4 = factor(OD_stats$V4,
                     levels = c("TP0", "TP3", "TP6", "TP12", "TP24")) 

# LME for time*community with random effect for community (Generation 1)
summary(lmer(OD600_corr ~ V1*V4 + (V4/Rep) + (1|V1), data = OD_stats[grep("G1", OD_stats$Row.names),]))
# LME for time*community with random effect for community (Generation 2)
summary(lmer(OD600_corr ~ V1*V4 + (V4/Rep) + (1|V1), data = OD_stats[grep("G2-T1", OD_stats$Row.names),]))
# LME for time*community with random effect for community (Generation 3)
summary(lmer(OD600_corr ~ V1*V4 + (V4/Rep) + (1|V1), data = OD_stats[grep("G3-T1", OD_stats$Row.names),]))

# LME for fixed effect of time*community*treatment with random effects for community and treatment
summary(lmer(OD600_corr ~ V3*V1*V4 + (V4/Rep) + (1|V3) + (1|V1), data = OD_stats[grep("G2", OD_stats$Row.names),]))
summary(lmer(OD600_corr ~ V3*V1*V4 + (V4/Rep) + (1|V3) + (1|V1), data = OD_stats[grep("G3", OD_stats$Row.names),]))

# LME for fixed effect of time*treatment with random effect for treatment (Generation 2)
summary(lmer(OD600_corr ~ V3*V4 + (V4/Rep) + (1|V3), data = OD_stats[grep("C1-G2", OD_stats$Row.names),]))
summary(lmer(OD600_corr ~ V3*V4 + (V4/Rep) + (1|V3), data = OD_stats[grep("C2-G2", OD_stats$Row.names),]))
summary(lmer(OD600_corr ~ V3*V4 + (V4/Rep) + (1|V3), data = OD_stats[grep("C3-G2", OD_stats$Row.names),]))
summary(lmer(OD600_corr ~ V3*V4 + (V4/Rep) + (1|V3), data = OD_stats[grep("C4-G2", OD_stats$Row.names),]))
summary(lmer(OD600_corr ~ V3*V4 + (V4/Rep) + (1|V3), data = OD_stats[grep("C5-G2", OD_stats$Row.names),]))

# LME for fixed effect of time*treatment with random effect for treatment (Generation 3)
summary(lmer(OD600_corr ~ V3*V4 + (V4/Rep) + (1|V3), data = OD_stats[grep("C1-G3", OD_stats$Row.names),]))
summary(lmer(OD600_corr ~ V3*V4 + (V4/Rep) + (1|V3), data = OD_stats[grep("C2-G3", OD_stats$Row.names),]))
summary(lmer(OD600_corr ~ V3*V4 + (V4/Rep) + (1|V3), data = OD_stats[grep("C3-G3", OD_stats$Row.names),]))
summary(lmer(OD600_corr ~ V3*V4 + (V4/Rep) + (1|V3), data = OD_stats[grep("C4-G3", OD_stats$Row.names),]))
summary(lmer(OD600_corr ~ V3*V4 + (V4/Rep) + (1|V3), data = OD_stats[grep("C5-G3", OD_stats$Row.names),]))

# LME for time*community with random effect for community (Generation 2)
summary(lmer(OD600_corr ~ V1*V4 + (V4/Rep) + (1|V1), data = OD_stats[grep("G2-T2", OD_stats$Row.names),]))
summary(lmer(OD600_corr ~ V1*V4 + (V4/Rep) + (1|V1), data = OD_stats[grep("G2-T3", OD_stats$Row.names),]))
summary(lmer(OD600_corr ~ V1*V4 + (V4/Rep) + (1|V1), data = OD_stats[grep("G2-T4", OD_stats$Row.names),]))

# LME for time*community with random effect for community (Generation 3)
summary(lmer(OD600_corr ~ V1*V4 + (V4/Rep) + (1|V1), data = OD_stats[grep("G3-T1", OD_stats$Row.names),]))
summary(lmer(OD600_corr ~ V1*V4 + (V4/Rep) + (1|V1), data = OD_stats[grep("G3-T2", OD_stats$Row.names),]))
summary(lmer(OD600_corr ~ V1*V4 + (V4/Rep) + (1|V1), data = OD_stats[grep("G3-T3", OD_stats$Row.names),]))
summary(lmer(OD600_corr ~ V1*V4 + (V4/Rep) + (1|V1), data = OD_stats[grep("G3-T4", OD_stats$Row.names),]))

## read-in corrected table
OD_rep_clean = read.table('mixed_community/OD_data_rep_clean.txt', header = T, sep = '\t', row.names = 1, comment.char = '')
head(OD_rep_clean)
dim(OD_rep_clean)

## prep data for plot

# C1
C1 = ggplot(OD_rep_clean[grep("C1", OD_rep_clean$Sample_ID),], 
            aes(x = Time, y = OD_600_mean/1000, col = Trt, fill = Trt)) +
  geom_point(size = 2.25, alpha = 0.8, pch = 21, alpha = 0.7) +
  geom_line(alpha = 0.7, lty = 1) +
  geom_errorbar(aes(ymax = OD_600_mean/1000 + OD_600_se/1000, 
                    ymin = OD_600_mean/1000 - OD_600_se/1000), width = 1, alpha = 0.7) +
  ylim(0, 0.55) +
  theme_bw() +
  xlab("Time (h)") +
  ylab("OD600") +
  scale_fill_manual(values = c("Black", "Blue", "Orange", "Red"),
                    labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  scale_color_manual(values = c("Black", "Blue", "Orange", "Red"),
                     labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.position = "none")

# C2
C2 = ggplot(OD_rep_clean[grep("C2", OD_rep_clean$Sample_ID),], 
            aes(x = Time, y = OD_600_mean/1000, col = Trt, fill = Trt)) +
  geom_point(size = 2.25, alpha = 0.8, pch = 21, alpha = 0.7) +
  geom_line(alpha = 0.7, lty = 1) +
  geom_errorbar(aes(ymax = OD_600_mean/1000 + OD_600_se/1000, 
                    ymin = OD_600_mean/1000 - OD_600_se/1000), width = 1, alpha = 0.7) +
  ylim(0, 0.55) +
  theme_bw() +
  xlab("Time (h)") +
  ylab("OD600") +
  scale_fill_manual(values = c("Black", "Blue", "Orange", "Red"),
                    labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  scale_color_manual(values = c("Black", "Blue", "Orange", "Red"),
                     labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.position = "none")

# C3
C3 = ggplot(OD_rep_clean[grep("C3", OD_rep_clean$Sample_ID),], 
            aes(x = Time, y = OD_600_mean/1000, col = Trt, fill = Trt)) +
  geom_point(size = 2.25, alpha = 0.8, pch = 21, alpha = 0.7) +
  geom_line(alpha = 0.7, lty = 1) +
  geom_errorbar(aes(ymax = OD_600_mean/1000 + OD_600_se/1000, 
                    ymin = OD_600_mean/1000 - OD_600_se/1000), width = 1, alpha = 0.7) +
  ylim(0, 0.55) +
  theme_bw() +
  xlab("Time (h)") +
  ylab("OD600") +
  scale_fill_manual(values = c("Black", "Blue", "Orange", "Red"),
                    labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  scale_color_manual(values = c("Black", "Blue", "Orange", "Red"),
                     labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.position = "none")

# C4
C4 = ggplot(OD_rep_clean[grep("C4", OD_rep_clean$Sample_ID),], 
            aes(x = Time, y = OD_600_mean/1000, col = Trt, fill = Trt)) +
  geom_point(size = 2.25, alpha = 0.8, pch = 21, alpha = 0.7) +
  geom_line(alpha = 0.7, lty = 1) +
  geom_errorbar(aes(ymax = OD_600_mean/1000 + OD_600_se/1000, 
                    ymin = OD_600_mean/1000 - OD_600_se/1000), width = 1, alpha = 0.7) +
  ylim(0, 0.55) +
  theme_bw() +
  xlab("Time (h)") +
  ylab("OD600") +
  scale_fill_manual(values = c("Black", "Blue", "Orange", "Red"),
                    labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  scale_color_manual(values = c("Black", "Blue", "Orange", "Red"),
                     labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.position = "none")

# C5
C5 = ggplot(OD_rep_clean[grep("C5", OD_rep_clean$Sample_ID),], 
            aes(x = Time, y = OD_600_mean/1000, col = Trt, fill = Trt)) +
  geom_point(size = 2.25, alpha = 0.8, pch = 21, alpha = 0.7) +
  geom_line(alpha = 0.7, lty = 1) +
  geom_errorbar(aes(ymax = OD_600_mean/1000 + OD_600_se/1000, 
                    ymin = OD_600_mean/1000 - OD_600_se/1000), width = 1, alpha = 0.7) +
  ylim(0, 0.55) +
  theme_bw() +
  xlab("Time (h)") +
  ylab("OD600") +
  scale_fill_manual(values = c("Black", "Blue", "Orange", "Red"),
                    labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  scale_color_manual(values = c("Black", "Blue", "Orange", "Red"),
                     labels = c("Control", "Low (10uM)", "Med (100uM)", "High (250uM)")) +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.position = "none")

## plot all (Fig. 4)
grid.arrange(C1, C2, C3, C4, C5, nrow = 2)


###
### Dox stats
###

### Area under the curve (AUC)

## prep table
AUC_dox = data.frame(comm = rep(c("1", "1", "1", "2", "2", "2", "3", "3", "3", "4", "4", "4", "5", "5", "5"), 3),
                     dox = c(rep("low", 15), rep("med", 15), rep("high", 15)),
                     AUC = c(AUC(OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T2", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time[c(1,2,3,5)],
                                 OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr[c(1,2,3,5)]),
                             AUC(OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T3", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C1-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C1-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C2-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C2-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C3-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C3-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C4-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time,
                                 OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),][grep("1", OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr),
                             AUC(OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time[c(1,2,3,5)],
                                 OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C4-G2", OD_stats$Row.names),]$V3),][grep("2", OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr[c(1,2,3,5)]),
                             AUC(OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),]$Rep),]$Time[c(1,3,4,5)],
                                 OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),][grep("3", OD_stats[grep("C5-G2", OD_stats$Row.names),][grep("T4", OD_stats[grep("C5-G2", OD_stats$Row.names),]$V3),]$Rep),]$Dox_R1_corr[c(1,3,4,5)])))

# ANOVA for effects of community
summary(aov(AUC~comm, AUC_dox[grep("low", AUC_dox$dox),]))
TukeyHSD(aov(AUC~comm, AUC_dox[grep("low", AUC_dox$dox),]))

summary(aov(AUC~comm, AUC_dox[grep("med", AUC_dox$dox),]))
TukeyHSD(aov(AUC~comm, AUC_dox[grep("med", AUC_dox$dox),]))

summary(aov(AUC~comm, AUC_dox[grep("high", AUC_dox$dox),]))
TukeyHSD(aov(AUC~comm, AUC_dox[grep("high", AUC_dox$dox),]))


### Drug transformation (changes in concentration) over time

# set time as factor
OD_stats$V4 = factor(OD_stats$V4, levels = c("TP0", "TP3", "TP6", "TP12", "TP24"))

# LME for fixed effect of time*treatment with random effect for treatment
summary(lmer(Dox_R1_corr ~ V3*V4 + (V4/Rep) + (1|V3), data = OD_stats[grep("C1-G2", OD_stats$Row.names),]))
summary(lmer(Dox_R1_corr ~ V3*V4 + (V4/Rep) + (1|V3), data = OD_stats[grep("C2-G2", OD_stats$Row.names),]))
summary(lmer(Dox_R1_corr ~ V3*V4 + (V4/Rep) + (1|V3), data = OD_stats[grep("C3-G2", OD_stats$Row.names),]))
summary(lmer(Dox_R1_corr ~ V3*V4 + (V4/Rep) + (1|V3), data = OD_stats[grep("C4-G2", OD_stats$Row.names),]))
summary(lmer(Dox_R1_corr ~ V3*V4 + (V4/Rep) + (1|V3), data = OD_stats[grep("C5-G2", OD_stats$Row.names),]))

# log transformed
summary(lmer(log(Dox_R1_corr+1) ~ V3*V4 + (V4/Rep) + (1|V3), data = OD_stats[grep("C1-G2", OD_stats$Row.names),]))
summary(lmer(log(Dox_R1_corr+1) ~ V3*V4 + (V4/Rep) + (1|V3), data = OD_stats[grep("C2-G2", OD_stats$Row.names),]))
summary(lmer(log(Dox_R1_corr+1) ~ V3*V4 + (V4/Rep) + (1|V3), data = OD_stats[grep("C3-G2", OD_stats$Row.names),]))
summary(lmer(log(Dox_R1_corr+1) ~ V3*V4 + (V4/Rep) + (1|V3), data = OD_stats[grep("C4-G2", OD_stats$Row.names),]))
summary(lmer(log(Dox_R1_corr+1) ~ V3*V4 + (V4/Rep) + (1|V3), data = OD_stats[grep("C5-G2", OD_stats$Row.names),]))

# LME for fixed effect of time*community and time*treatment with random effects for community and treatment
summary(lmer(Dox_R1_corr ~ V3*V1*V4 + (V4/Rep) + (1|V3) + (1|V1), data = OD_stats[grep("G2", OD_stats$Row.names),]))

# LME for fixed effect of time*community with random effects for community
summary(lmer(Dox_R1_corr ~ V1*V4 + (V4/Rep) + (1|V1), data = OD_stats[grep("G2-T2", OD_stats$Row.names),]))
summary(lmer(Dox_R1_corr ~ V1*V4 + (V4/Rep) + (1|V1), data = OD_stats[grep("G2-T3", OD_stats$Row.names),]))
summary(lmer(Dox_R1_corr ~ V1*V4 + (V4/Rep) + (1|V1), data = OD_stats[grep("G2-T4", OD_stats$Row.names),]))

# post-hoc pairwise comparisons
test = lmer(Dox_R1_corr ~ V1*V4 + (V4/Rep) + (1|V1), data = OD_stats[grep("G2-T3", OD_stats$Row.names),])
emmeans(test, list(pairwise ~ V1*V4), adjust = "tukey")
test2 = lmer(Dox_R1_corr ~ V1*V4 + (V4/Rep) + (1|V1), data = OD_stats[grep("G2-T4", OD_stats$Row.names),])
emmeans(test2, list(pairwise ~ V1*V4), adjust = "tukey")

# log transformed
summary(lmer(log(Dox_R1_corr+1) ~ V1*V4 + (V4/Rep) + (1|V1), data = OD_stats[grep("G2-T2", OD_stats$Row.names),]))
summary(lmer(log(Dox_R1_corr+1) ~ V1*V4 + (V4/Rep) + (1|V1), data = OD_stats[grep("G2-T3", OD_stats$Row.names),]))
summary(lmer(log(Dox_R1_corr+1) ~ V1*V4 + (V4/Rep) + (1|V1), data = OD_stats[grep("G2-T4", OD_stats$Row.names),]))

## prep data for plot

# low trt
T2 = ggplot(OD_rep_clean[grep("G2-T2", OD_rep_clean$Sample_ID),], 
            aes(x = TP, y = Dox_mean, col = Com, fill = Com)) +
  geom_point(size = 2.5, alpha = 0.7, pch = 21) +
  geom_line(aes(y = Dox_mean), alpha = 0.6, lty = 1) +
  geom_errorbar(aes(ymax = Dox_mean + Dox_se, ymin = Dox_mean - Dox_se), 
                width=0.8, alpha = 0.5, lty = 1) +
  theme_bw() +
  xlab("Time (h)") +
  ylab("Dox µM") +
  ylim(0, 20) +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.position = "none")

# med trt
T3 = ggplot(OD_rep_clean[grep("G2-T3", OD_rep_clean$Sample_ID),], 
            aes(x = TP, y = Dox_mean, col = Com, fill = Com)) +
  geom_point(size = 2.5, alpha = 0.7, pch = 21) +
  geom_line(aes(y = Dox_mean), alpha = 0.6, lty = 1) +
  geom_errorbar(aes(ymax = Dox_mean + Dox_se, ymin = Dox_mean - Dox_se), 
                width=0.8, alpha = 0.5, lty = 1) +
  theme_bw() +
  xlab("Time (h)") +
  ylab("Dox µM") +
  ylim(0, 115) +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.position = "none")

# high trt
T4 = ggplot(OD_rep_clean[grep("G2-T4", OD_rep_clean$Sample_ID),], 
            aes(x = TP, y = Dox_mean, col = Com, fill = Com)) +
  geom_point(size = 2.5, alpha = 0.7, pch = 21) +
  geom_line(aes(y = Dox_mean), alpha = 0.6, lty = 1) +
  geom_errorbar(aes(ymax = Dox_mean + Dox_se, ymin = Dox_mean - Dox_se), 
                width=0.8, alpha = 0.5, lty = 1) +
  theme_bw() +
  xlab("Time (h)") +
  ylab("Dox µM") +
  ylim(0, 290) +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.position = "none")

## plot all (Fig. 6)
grid.arrange(T2, T3, T4, nrow = 1)


###
### Dox assay sensitivity
###

## monoculture data points
mono_calc = mono_data[grep("Ctrl", mono_data$Treatment),]$Dox_conc
mono_calc = mono_calc[-c(which(is.na(mono_calc)))]
#mono_calc = mono_calc[-c(which(mono_calc == 0))]
mean(mono_calc)
sd(mono_calc)
sd(mono_calc)/sqrt(length(mono_calc))
range(mono_calc)

## mixed community data points
mixed_calc = OD_stats[grep("G2-T1", rownames(OD_stats)), grep("Dox", colnames(OD_stats))]$Dox_conc
mean(mixed_calc)
sd(mixed_calc)
sd(mixed_calc)/sqrt(length(mixed_calc))
range(mixed_calc)

## combine all data for assay on samples with no dox
mean(c(mono_calc, mixed_calc))
sd(c(mono_calc, mixed_calc))
sd(c(mono_calc, mixed_calc))/sqrt(length(c(mono_calc, mixed_calc)))
range(c(mono_calc, mixed_calc))

# plot (Fig. S2 B)
ggplot(data.frame(val = c(mono_calc, mixed_calc),
                  dataset = c(rep("mono", length(mono_calc)), rep("mixed", length(mixed_calc)))),
       aes(x=val)) +
  geom_histogram(#aes(y=..density..),
    binwidth = 1, color = "black", fill = "gray", alpha = 0.6) +
  ylab("Sample count") +
  xlab("Detected dox µM") +
  geom_vline(xintercept = -1.091452, lty = 2, size = 1.5, color = "blue") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 16))


###
### 16S amplicon data: stats
###

## prep table
amplicon_stats = data.frame(OD_16S_corr[,2:6])
rownames(amplicon_stats) = OD_16S_corr$Row.names
#amplicon_stats = taxa_merge_code
#amplicon_stats = melt(t(taxa_merge_rep))
#amplicon_stats$Var2 = as.character(amplicon_stats$Var2)
amplicon_stats = data.frame(amplicon_stats,
                            do.call(rbind, strsplit(as.character(rownames(amplicon_stats)), "-", fixed=TRUE)))
amplicon_stats$sample = substr(rownames(amplicon_stats), start = 1, stop = nchar(rownames(amplicon_stats))-3)

## subset for times 0, 12, 24, 36, 48, 60, 72; set time values
amplicon_stats = amplicon_stats[grep("G1-NA-TP0|TP12|TP24", rownames(amplicon_stats)),]
amplicon_stats$time = rep("0", length(amplicon_stats$sample))
amplicon_stats$time[grep("G1-NA-TP12", amplicon_stats$sample)] = rep("12", length(amplicon_stats$time[grep("G1-NA-TP12", amplicon_stats$sample)]))
amplicon_stats$time[grep("G1-NA-TP24", amplicon_stats$sample)] = rep("24", length(amplicon_stats$time[grep("G1-NA-TP24", amplicon_stats$sample)]))
amplicon_stats$time[grep("G2-T1-TP12|G2-T2-TP12|G2-T3-TP12|G2-T4-TP12", amplicon_stats$sample)] = rep("36", length(amplicon_stats$time[grep("G2-T1-TP12|G2-T2-TP12|G2-T3-TP12|G2-T4-TP12", amplicon_stats$sample)]))
amplicon_stats$time[grep("G2-T1-TP24|G2-T2-TP24|G2-T3-TP24|G2-T4-TP24", amplicon_stats$sample)] = rep("48", length(amplicon_stats$time[grep("G2-T1-TP24|G2-T2-TP24|G2-T3-TP24|G2-T4-TP24", amplicon_stats$sample)]))
amplicon_stats$time[grep("G3-T1-TP12|G3-T2-TP12|G3-T3-TP12|G3-T4-TP12", amplicon_stats$sample)] = rep("60", length(amplicon_stats$time[grep("G3-T1-TP12|G3-T2-TP12|G3-T3-TP12|G3-T4-TP12", amplicon_stats$sample)]))
amplicon_stats$time[grep("G3-T1-TP24|G3-T2-TP24|G3-T3-TP24|G3-T4-TP24", amplicon_stats$sample)] = rep("72", length(amplicon_stats$time[grep("G3-T1-TP24|G3-T2-TP24|G3-T3-TP24|G3-T4-TP24", amplicon_stats$sample)]))


### LME for fixed effect of time*treatment with random effect for treatment (Gen 2 and Gen 3)

## By strain

# KP
summary(lmer(Klebsiella ~ X1*X3*time + (time/X5) + (1|X3) + (1|X1), data = amplicon_stats[grep("C1-G2|C1-G3|C2-G2|C2-G3|C3-G2|C3-G3", amplicon_stats$sample),]))
# EC
summary(lmer(Escherichia ~ X1*X3*time + (time/X5) + (1|X3) + (1|X1), data = amplicon_stats[grep("C1-G2|C1-G3|C2-G2|C2-G3|C4-G2|C4-G3", amplicon_stats$sample),]))
# KP
summary(lmer(Enterococcus ~ X1*X3*time + (time/X5) + (1|X3) + (1|X1), data = amplicon_stats[grep("C1-G2|C1-G3|C4-G2|C4-G3|C3-G2|C3-G3|C5-G2|C5-G3", amplicon_stats$sample),]))

## By community

# C1
summary(lmer(Klebsiella ~ X3*time + (time/X5) + (1|X3) + (1|X2), data = amplicon_stats[grep("C1-G2|C1-G3", amplicon_stats$sample),]))
summary(lmer(Escherichia ~ X3*time + (time/X5) + (1|X3) + (1|X2), data = amplicon_stats[grep("C1-G2|C1-G3", amplicon_stats$sample),]))
summary(lmer(Enterococcus ~ X3*time + (time/X5) + (1|X3) + (1|X2), data = amplicon_stats[grep("C1-G2|C1-G3", amplicon_stats$sample),]))
summary(lmer(Lactobacillus ~ X3*time + (time/X5) + (1|X3) + (1|X2), data = amplicon_stats[grep("C1-G2|C1-G3", amplicon_stats$sample),]))
summary(lmer(Clostridium ~ X3*time + (time/X5) + (1|X3) + (1|X2), data = amplicon_stats[grep("C1-G2|C1-G3", amplicon_stats$sample),]))
# C2
summary(lmer(Klebsiella ~ X3*time + (time/X5) + (1|X3) + (1|X2), data = amplicon_stats[grep("C2-G2|C2-G3", amplicon_stats$sample),]))
summary(lmer(Escherichia ~ X3*time + (time/X5) + (1|X3) + (1|X2), data = amplicon_stats[grep("C2-G2|C2-G3", amplicon_stats$sample),]))
summary(lmer(Lactobacillus ~ X3*time + (time/X5) + (1|X3) + (1|X2), data = amplicon_stats[grep("C2-G2|C2-G3", amplicon_stats$sample),]))
summary(lmer(Clostridium ~ X3*time + (time/X5) + (1|X3) + (1|X2), data = amplicon_stats[grep("C2-G2|C2-G3", amplicon_stats$sample),]))
# C3
summary(lmer(Klebsiella ~  X3*time + (time/X5) + (1|X3) + (1|X2), data = amplicon_stats[grep("C3-G2|C3-G3", amplicon_stats$sample),]))
summary(lmer(Enterococcus ~ X3*time + (time/X5) + (1|X3) + (1|X2), data = amplicon_stats[grep("C3-G2|C3-G3", amplicon_stats$sample),]))
summary(lmer(Lactobacillus ~ X3*time + (time/X5) + (1|X3) + (1|X2), data = amplicon_stats[grep("C3-G2|C3-G3", amplicon_stats$sample),]))
summary(lmer(Clostridium ~ X3*time + (time/X5) + (1|X3) + (1|X2), data = amplicon_stats[grep("C3-G2|C3-G3", amplicon_stats$sample),]))
# C4
summary(lmer(Escherichia ~ X3*time + (time/X5) + (1|X3) + (1|X2), data = amplicon_stats[grep("C4-G2|C4-G3", amplicon_stats$sample),]))
summary(lmer(Enterococcus ~ X3*time + (time/X5) + (1|X3) + (1|X2), data = amplicon_stats[grep("C4-G2|C4-G3", amplicon_stats$sample),]))
summary(lmer(Lactobacillus ~ X3*time + (time/X5) + (1|X3) + (1|X2), data = amplicon_stats[grep("C4-G2|C4-G3", amplicon_stats$sample),]))
summary(lmer(Clostridium ~ X3*time + (time/X5) + (1|X3) + (1|X2), data = amplicon_stats[grep("C4-G2|C4-G3", amplicon_stats$sample),]))
# C5
summary(lmer(Enterococcus ~ X3*time + (time/X5) + (1|X3), data = amplicon_stats[grep("C5-G2|C5-G3", amplicon_stats$sample),]))
summary(lmer(Lactobacillus ~ X3*time + (time/X5) + (1|X3) + (1|X2), data = amplicon_stats[grep("C5-G2|C5-G3", amplicon_stats$sample),]))
summary(lmer(Clostridium ~ X3*time + (time/X5) + (1|X3) + (1|X2), data = amplicon_stats[grep("C5-G2|C5-G3", amplicon_stats$sample),]))


### relative abundance data visualization

## re-compute for target strains using corrected 16S data

## prep data frame
taxa_merge_code = data.frame(OD_16S_corr[,2:6])
rownames(taxa_merge_code) = OD_16S_corr$Row.names

taxa_merge_rep = data.frame(CI_mean = tapply(taxa_merge_code[,1], substr(rownames(taxa_merge_code), start = 1, stop = nchar(rownames(taxa_merge_code)) - 3), mean),
                            Lac_mean = tapply(taxa_merge_code[,2], substr(rownames(taxa_merge_code), start = 1, stop = nchar(rownames(taxa_merge_code)) - 3), mean),
                            EF_mean = tapply(taxa_merge_code[,3], substr(rownames(taxa_merge_code), start = 1, stop = nchar(rownames(taxa_merge_code)) - 3), mean),
                            EC_mean = tapply(taxa_merge_code[,4], substr(rownames(taxa_merge_code), start = 1, stop = nchar(rownames(taxa_merge_code)) - 3), mean),
                            KP_mean = tapply(taxa_merge_code[,5], substr(rownames(taxa_merge_code), start = 1, stop = nchar(rownames(taxa_merge_code)) - 3), mean))

# re-order data
taxa_merge_rep = taxa_merge_rep[mixedorder(rownames(taxa_merge_rep)),]

## prepare data frame (skip above for plot)
amplicon_df = taxa_merge_rep
amplicon_df = melt(t(taxa_merge_rep))
amplicon_df$Var2 = as.character(amplicon_df$Var2)
amplicon_df = data.frame(amplicon_df, do.call(rbind, strsplit(as.character(amplicon_df$Var2), "-", fixed=TRUE)))

## subset for times 0, 12, 24, 36, 48, 60, 72
amplicon_df = amplicon_df[grep("G1-NA-TP0|TP12|TP24", amplicon_df$Var2),]

## add categories from sample names
amplicon_df$X1 = as.character(amplicon_df$X1)
amplicon_df$X2 = as.character(amplicon_df$X2)
amplicon_df$X3 = as.character(amplicon_df$X3)
amplicon_df$X4 = as.character(amplicon_df$X4)

paste_df = c(1:dim(amplicon_df)[1])
for (i in 1:dim(amplicon_df)[1]) {
  paste_df[i] = paste(amplicon_df[i,c(4,6,5,7)], collapse = "-")
}
amplicon_df$paste_df = paste_df

amplicon_df$X3[grep("NA", amplicon_df$X3)] = rep("Gen. 1", length(amplicon_df$X3[grep("NA", amplicon_df$X3)]))
amplicon_df$X3[grep("T1", amplicon_df$X3)] = rep("0 µM", length(amplicon_df$X3[grep("T1", amplicon_df$X3)]))
amplicon_df$X3[grep("T2", amplicon_df$X3)] = rep("10 µM", length(amplicon_df$X3[grep("T2", amplicon_df$X3)]))
amplicon_df$X3[grep("T3", amplicon_df$X3)] = rep("100 µM", length(amplicon_df$X3[grep("T3", amplicon_df$X3)]))
amplicon_df$X3[grep("T4", amplicon_df$X3)] = rep("250 µM", length(amplicon_df$X3[grep("T4", amplicon_df$X3)]))

# set factor for plot order
amplicon_df$X3 = as.factor(amplicon_df$X3)
amplicon_df$X3 = factor(amplicon_df$X3,
                        levels = c("Gen. 1", "0 µM", "10 µM", "100 µM", "250 µM"))

## plot starting content (Fig. S5 B)
ggplot(amplicon_df[grep("G1-TP0", amplicon_df$paste_df),], 
       aes(x = X1, y = value, fill = Var1)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  xlab("Time (h)") +
  ylab("Relative abundance") +
  theme_bw() +
  scale_fill_manual(values = c("darkred", "pink", "orange", "skyblue", "darkblue")) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        axis.title.x = element_blank())

## plot changes over time

# C1
C1 = ggplot(amplicon_df[grep("C1", amplicon_df$X1),], 
            aes(x = paste_df, y = value, fill = Var1)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  xlab("Time (h)") +
  ylab("Rel. abundance") +
  facet_grid(.~X3, scale="free_x", space="free") +
  scale_x_discrete(labels = c(36, 48, 60, 72)) +
  theme_bw() +
  scale_fill_manual(values = c("darkred", "pink", "orange", "skyblue", "darkblue")) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        strip.text.x = element_text(size = 9, margin = margin(.05, 0, .05, 0, "cm")),
        legend.position = "none",
        axis.title.x = element_blank())

# C2
C2 = ggplot(amplicon_df[grep("C2", amplicon_df$X1),], 
            aes(x = paste_df, y = value, fill = Var1)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  xlab("Time (h)") +
  ylab("Rel. abundance") +
  facet_grid(.~X3, scale="free_x", space="free") +
  scale_x_discrete(labels = c(36, 48, 60, 72)) +
  theme_bw() +
  scale_fill_manual(values = c("darkred", "pink", "orange", "skyblue", "darkblue")) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        strip.text.x = element_text(size = 9, margin = margin(.05, 0, .05, 0, "cm")),
        legend.position = "none",
        axis.title.x = element_blank())

# C3
C3 = ggplot(amplicon_df[grep("C3", amplicon_df$X1),], 
            aes(x = paste_df, y = value, fill = Var1)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  xlab("Time (h)") +
  ylab("Rel. abundance") +
  facet_grid(.~X3, scale="free_x", space="free") +
  scale_x_discrete(labels = c(36, 48, 60, 72)) +
  theme_bw() +
  scale_fill_manual(values = c("darkred", "pink", "orange", "skyblue", "darkblue")) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        strip.text.x = element_text(size = 9, margin = margin(.05, 0, .05, 0, "cm")),
        legend.position = "none",
        axis.title.x = element_blank())

# C4
C4 = ggplot(amplicon_df[grep("C4", amplicon_df$X1),], 
            aes(x = paste_df, y = value, fill = Var1)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  xlab("Time (h)") +
  ylab("Rel. abundance") +
  facet_grid(.~X3, scale="free_x", space="free") +
  scale_x_discrete(labels = c(36, 48, 60, 72)) +
  theme_bw() +
  scale_fill_manual(values = c("darkred", "pink", "orange", "skyblue", "darkblue")) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        strip.text.x = element_text(size = 9, margin = margin(.05, 0, .05, 0, "cm")),
        legend.position = "none",
        axis.title.x = element_blank())

# C5
C5 = ggplot(amplicon_df[grep("C5", amplicon_df$X1),], 
            aes(x = paste_df, y = value, fill = Var1)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  xlab("Time (h)") +
  ylab("Rel. abundance") +
  facet_grid(.~X3, scale="free_x", space="free") +
  scale_x_discrete(labels = c(36, 48, 60, 72)) +
  theme_bw() +
  scale_fill_manual(values = c("darkred", "pink", "orange", "skyblue", "darkblue")) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        strip.text.x = element_text(size = 9, margin = margin(.05, 0, .05, 0, "cm")),
        legend.position = "none",
        axis.title.x = element_blank())

## plot all
grid.arrange(C1, C2, C3, C4, C5, ncol = 1)


###
### Ratios of drug-resistant strains
###

## prep table
test = taxa_merge_code
test$KP_EC = log(test$Klebsiella/test$Escherichia, 10)
test$KP_EF = log(test$Klebsiella/test$Enterococcus, 10)
test$EC_EF = log(test$Escherichia/test$Enterococcus, 10)

# clean for stats
test_stats = data.frame(test[,6:8],
                        do.call(rbind, strsplit(as.character(rownames(test)), "-", fixed=TRUE)))
test_stats = test_stats[grep("G2|G3", rownames(test_stats)),]
test_stats = test_stats[grep("TP6|TP12|TP24", rownames(test_stats)),]
#test_stats = test_stats[grep("TP12|TP24", rownames(test_stats)),]
head(test_stats)

# set factor for time
test_stats$time = rep("T36", length(test_stats$X4))
test_stats$time[grep("G2-T1-TP24|G2-T2-TP24|G2-T3-TP24|G2-T4-TP24", rownames(test_stats))] = rep("T48", length(test_stats$time[grep("G2-T1-TP24|G2-T2-TP24|G2-T3-TP24|G2-T4-TP24", rownames(test_stats))]))
test_stats$time[grep("G3-T1-TP12|G3-T2-TP12|G3-T3-TP12|G3-T4-TP12", rownames(test_stats))] = rep("T60", length(test_stats$time[grep("G3-T1-TP12|G3-T2-TP12|G3-T3-TP12|G3-T4-TP12", rownames(test_stats))]))
test_stats$time[grep("G3-T1-TP24|G3-T2-TP24|G3-T3-TP24|G3-T4-TP24", rownames(test_stats))] = rep("T72", length(test_stats$time[grep("G3-T1-TP24|G3-T2-TP24|G3-T3-TP24|G3-T4-TP24", rownames(test_stats))]))
test_stats$time = factor(test_stats$time)
head(test_stats)

### Changes in ratios by treatment group: LME

# KP_EC
summary(lmer(KP_EC ~ X3*X1*time + (time/X5) + (1|X3), 
             data = test_stats[grep("C1|C2", test_stats$X1), c(1,4:9)]))
summary(lmer(KP_EC ~ X3*time + (time/X5) + (1|X3), 
             data = test_stats[grep("C1", test_stats$X1), c(1,4:9)]))
summary(lmer(KP_EC ~ X3*time + (time/X5) + (1|X3), 
             data = test_stats[grep("C2", test_stats$X1), c(1,4:9)]))

# KP_EF
summary(lmer(KP_EF ~ X3*X1*time + (time/X5) + (1|X3), 
             data = test_stats[grep("C1|C3", test_stats$X1), c(2,4:9)]))
summary(lmer(KP_EF ~ X3*time + (time/X5) + (1|X3), 
             data = test_stats[grep("C1", test_stats$X1), c(2,4:9)]))
summary(lmer(KP_EF ~ X3*time + (time/X5) + (1|X3), 
             data = test_stats[grep("C3", test_stats$X1), c(2,4:9)]))

# EC_EF
summary(lmer(EC_EF ~ X3*X1*time + (time/X5) + (1|X3), 
             data = test_stats[grep("C1|C4", test_stats$X1), c(3,4:9)]))
summary(lmer(EC_EF ~ X3*time + (time/X5) + (1|X3), 
             data = test_stats[grep("C1", test_stats$X1), c(3,4:9)]))
summary(lmer(EC_EF ~ X3*time + (time/X5) + (1|X3), 
             data = test_stats[grep("C4", test_stats$X1), c(3,4:9)]))

### prep for plot
taxa_merge_code = test
test_merge_rep = data.frame(CI_mean = tapply(taxa_merge_code[,1], substr(rownames(taxa_merge_code), start = 1, stop = nchar(rownames(taxa_merge_code)) - 3), mean),
                            Lac_mean = tapply(taxa_merge_code[,2], substr(rownames(taxa_merge_code), start = 1, stop = nchar(rownames(taxa_merge_code)) - 3), mean),
                            EF_mean = tapply(taxa_merge_code[,3], substr(rownames(taxa_merge_code), start = 1, stop = nchar(rownames(taxa_merge_code)) - 3), mean),
                            EC_mean = tapply(taxa_merge_code[,4], substr(rownames(taxa_merge_code), start = 1, stop = nchar(rownames(taxa_merge_code)) - 3), mean),
                            KP_mean = tapply(taxa_merge_code[,5], substr(rownames(taxa_merge_code), start = 1, stop = nchar(rownames(taxa_merge_code)) - 3), mean),
                            KP_EC_mean = tapply(taxa_merge_code[,6], substr(rownames(taxa_merge_code), start = 1, stop = nchar(rownames(taxa_merge_code)) - 3), mean),
                            KP_EF_mean = tapply(taxa_merge_code[,7], substr(rownames(taxa_merge_code), start = 1, stop = nchar(rownames(taxa_merge_code)) - 3), mean),
                            EC_EF_mean = tapply(taxa_merge_code[,8], substr(rownames(taxa_merge_code), start = 1, stop = nchar(rownames(taxa_merge_code)) - 3), mean))
test = test_merge_rep
test = data.frame(test,
                  do.call(rbind, strsplit(as.character(rownames(test)), "-", fixed=TRUE)))
test = test[grep("NA-TP0|TP6|TP12|TP24", rownames(test)),]

# set factor for time point
test$X4 = factor(test$X4,
                 levels = c("TP0", "TP3", "TP6", "TP12", "TP24"))
test = test[order(test$X4),]
test = test[order(test$X2),]
test = test[order(test$X3),]

# new table for heatmap
test2 = rbind(test[grep("C1", test$X1),6],
              test[grep("C2", test$X1),6],
              test[grep("C1", test$X1),7],
              test[grep("C3", test$X1),7],
              test[grep("C1", test$X1),8],
              test[grep("C4", test$X1),8])
rownames(test2) = c("C1: KP/EC", "C2: KP/EC",
                    "C1: KP/EF", "C3: KP/EF",
                    "C1: EC/EF", "C4: EC/EF")
colnames(test2) = rownames(test[grep("C1", test$X1),10:11])

# set color palette
mypalette <- colorRampPalette(c("blue","white","red"))(n = 100)
mycolors = seq(-1.5, 1.5, length.out = 101)

# plot key for heatmap color
hist(seq(-1.5, 1.5, length.out = 1000), breaks = mycolors, 
     col = mypalette, 
     border = mypalette, 
     ylim = c(0, 0.05),
     xlim = c(-1.5, 1.5),
     xlab = NULL, ylab = NULL, 
     axes = FALSE,
     main = NA)
axis(side = 1, labels = c("-1.5", NA, "-0.5", NA, "0.5", NA, "1.5"),
     at = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5),
     las = 1, cex = 4)

# plot heatmap (Fig. 7)
heatmap.2(test2[,-c(grep("TP6", colnames(test2)))],
          Rowv = F,
          Colv = F,
          keysize = 1,
          key.title = NA,
          key.xlab = NA,
          symkey=FALSE, 
          cexRow = 1.5,
          cexCol = 2,
          colsep = c(3,7,11,15),
          srtCol = 90,
          sepcolor = "black",
          labCol = c("0", "12", "24",
                     rep(c("36", "48", "60", "72"), 4)),
          density.info='none',
          trace = "none",
          col = mypalette,
          breaks = mycolors,
          RowSideColors = rep(c("#F8766D", "darkgray"),3),
          labRow = NA,
          lwid = c(0.5,7),
          offsetCol = -0.25,
          offsetRow = -38)

