#loading library-----------------
library(tidyverse)
library(asreml)
source('./script/outlier.R')
library(data.table)
library(ASRgenomics)
library(flextable)
library(dplyr)
library(SmartEDA)
library(explore)
library(patchwork)
library(janitor)
library(cowplot)
library(RColorBrewer)
#reading field design and wavelength data and filtering year 16 ("wavelength data available for year 16)

design<- read.csv('./data/design.csv') %>% filter(year!= 16)
averagedspectra <- read.csv('./data/Averaged_Spectra.csv')

averagedspectra <- averagedspectra %>%
  separate_wider_delim(
    Spectra,
    delim = "_",
    names = c("Spectra", "plotnum"),
    cols_remove = T
  ) %>%
  select(-plotnum) %>%
  group_by(Spectra) %>%
  summarise_all(mean) %>%
  ungroup() %>%
  rename(PlotID = Spectra)

designnew<- left_join(design,averagedspectra, by = "PlotID") # Apply left_join dplyr function

#finding 7 different check varieties that are replicated in all block inside each set in location "EF"

is_check_ef <- table(design[design$loc == 'EF', ]$name2) > 1
checks_ef <- names(is_check_ef[is_check_ef == TRUE])


##finding 7 different check varieties that are replicated in all block inside each set in location "MW"

is_check_mw<- table(design[design$loc == 'MW', ]$name2) > 1
checks_mw <- names(is_check_mw[is_check_mw == TRUE])


is_check_efsingle <- table(design[design$loc == 'EF',]$name2) >16
check_efsingle<- names(is_check_efsingle[is_check_efsingle == TRUE])
#so we can see check var "spx" is replicated 17 times while other 6 check varieties
#are replicated 16th times one in each block inside four sets in each location "EF" and " MW".

#write.csv(designnew, './data/designnew.csv')

designnew<- fread('./data/designnew.csv', data.table = FALSE)
#for location EF
desplot::desplot(
  design%>% filter(loc == "EF"),
  block ~ range * row,
  out2 = block,
  out1 = set,
  cex = 0.7,
  ticks = T
)
#for location MW
desplot::desplot(
  design %>% filter(loc == "MW"),
  block ~ range * row,
  out2 = block,
  out1 = set,
  cex = 0.7,
  ticks = T, text= name2,
)

#visualizing different accuracy of each model of joint location for narea
corr_N_wholewave_joint <-read.csv("./output/wholewave/corr_N_wholewave_joint.csv") %>%
  rename(wholewave = result_N_wholewave_joint.ac)
corr_N_highh2_joint <-read.csv("./output/highh2/corr_N_highh2_joint.csv")%>%
  rename(high_heritablewave = result_N_highh2_joint.ac)
corr_N_nirs_joint <- read.csv("./output/NIRS/corr_N_nirs_joint.csv")%>% 
  rename(nirswave = result_N_nirs_joint.ac)
corr_N_joint <- read.csv("./output/GBLUP/corr_N_joint.csv")%>% 
  rename( GBLUP= result_N.ac)

narea_accuracy_joint <- bind_cols(corr_N_joint,
                                  corr_N_nirs_joint,
                                  corr_N_highh2_joint,
                                  corr_N_wholewave_joint) %>%
                              pivot_longer(cols = 1:4,
                              names_to = "models",
                              values_to = "accuracy")
# Create a boxplot
joint_narea_plot<- ggplot(narea_accuracy_joint, aes(x = models, y = accuracy, 
                            fill = models)) +geom_boxplot() +
  labs(title = "Boxplot of prediction Accuracy for joint loc for narea",
       x = "Models", y = "Accuracy")+scale_y_continuous("prediction accuracy",
                                limits=c(-0.2,0.55), breaks=c(0, 0.25, 0.5))


#visualizing different accuracy of each model of joint location for narea
corr_S_wholewave_joint<-read.csv("./output/wholewave/corr_S_wholewave_joint.csv") %>%
  rename(wholewave= result_S_wholewave_joint.ac)
corr_S_highh2_joint<- read.csv("./output/highh2/corr_S_highh2_joint.csv")%>%
  rename(high_heritablewave = result_S_highh2_joint.ac)
corr_S_nirs_joint<- read.csv("./output/NIRS/corr_S_nirs_joint.csv")%>% 
  rename(nirswave = result_S_nirs_joint.ac)
corr_S_joint<-read.csv("./output/GBLUP/corr_S_joint.csv")%>% 
  rename( GBLUP= result_S.ac)

narea_accuracy_joint<- bind_cols(corr_S_joint,corr_S_nirs_joint,
                               corr_S_highh2_joint,corr_S_wholewave_joint)%>% 
  pivot_longer(cols = 1:4,names_to = "models", values_to = "accuracy")
# Create a boxplot
joint_narea_plot<- ggplot(narea_accuracy_joint, aes(x = models, y = accuracy, 
                                                fill = models)) +
  geom_boxplot() +
  labs(title = "Boxplot of prediction Accuracy for joint loc for sla",
       x = "Models",
       y = "Accuracy")+
  scale_y_continuous("prediction accuracy",limits=c(-0.2,0.55), breaks=c(0, 0.25, 0.5))

#visualizing model accuracy training on ef and validating on mw location for narea
corr_N_wholewave_efmw<- read.csv("./output/wholewave/corr_N_wholewave_efmw.csv")%>%
  rename(wholewave= result_narea_wholewave_efmw.ac)
corr_N_highh2_efmw<- read.csv("./output/highh2/corr_N_highh2_efmw.csv")%>% 
  rename(high_heritablewave= result_narea_highh2_efmw.ac)
corr_N_nirs_efmw<- read.csv("./output/NIRS/ corr_N_nirs_efmw.csv")%>% 
  rename(nirswave= result_narea_nirs_efmw.ac)
corr_N_efmw<- read.csv("./output/GBLUP/corr_N_efmw.csv")%>% 
  rename(GBLUP = result_narea_efmw.ac)

narea_accuracy_efmw<- bind_cols(corr_N_efmw,corr_N_nirs_efmw,corr_N_highh2_efmw,
                                corr_N_wholewave_efmw )%>% 
  pivot_longer(cols = 1:4,names_to = "models", values_to = "accuracy")

# Create a boxplot
efmw_narea_plot<- ggplot(narea_accuracy_efmw, aes(x = models, 
                                                  y = accuracy, fill = models)) +
  geom_boxplot() +
  labs(title = "Boxplot of prediction Accuracy using ef training and mw val for narea",
       x = "Models",
       y = "Accuracy")+
  scale_y_continuous("prediction accuracy",limits=c(-0.2,0.55), breaks=c(0, 0.25, 0.5))

#visualizing model accuracy training on ef and validating on mw location for sla
corr_S_wholewave_efmw<- read.csv("./output/wholewave/corr_S_wholewave_efmw.csv")%>% 
  rename(wholewave= result_sla_wholewave_efmw.ac)
corr_S_highh2_efmw<- read.csv("./output/highh2/corr_S_highh2_efmw.csv")%>% 
  rename(high_heritablewave= result_sla_highh2_efmw.ac)
corr_S_nirs_efmw<- read.csv("./output/NIRS/ corr_sla_nirs_efmw.csv")%>% 
  rename(nirswave= result_sla_nirs_efmw.ac)
corr_S_efmw<- read.csv("./output/GBLUP/corr_S_efmw.csv")%>% 
  rename(GBLUP = result_sla_efmw.ac)

sla_accuracy_efmw<- bind_cols(corr_S_efmw,corr_S_nirs_efmw,corr_S_highh2_efmw,
                              corr_S_wholewave_efmw )%>% 
  pivot_longer(cols = 1:4,names_to = "models", values_to = "accuracy")

# Create a boxplot
efmw_sla_plot<- ggplot(sla_accuracy_efmw, aes(x = models, y = accuracy, 
                                              fill = models)) +
  geom_boxplot() +
  labs(title = "Boxplot of prediction Accuracy using ef training and mw val for sla",
       x = "Models",
       y = "Accuracy")+
  scale_y_continuous("prediction accuracy",limits=c(-0.2,0.55), breaks=c(0, 0.25, 0.5))

#visualizing model accuracy training on mw and validating on ef location for sla
corr_S_wholewave_mwef<- read.csv("./output/wholewave/corr_S_wholewave_mwef.csv")%>%
  rename(wholewave= result_sla_wholewave_mwef.ac)
corr_S_highh2_mwef<- read.csv("./output/highh2/corr_S_highh2_mwef.csv")%>% 
  rename(highly_heritablewave= result_sla_highh2_mwef.ac)
corr_S_nirs_mwef<- read.csv("./output/NIRS/ corr_sla_nirs_mwef.csv")%>% 
  rename(nirswave= result_sla_nirs_mwef.ac)
corr_S_mwef<- read.csv("./output/GBLUP/corr_S_mwef.csv")%>% 
  rename(GBLUP= result_sla_mwef.ac)
sla_accuracy_mwef<- bind_cols(corr_S_mwef,corr_S_nirs_mwef,corr_S_highh2_mwef,
                              corr_S_wholewave_mwef )%>% 
  pivot_longer(cols = 1:4, names_to = "models", values_to = "accuracy")
# Creating a boxplot
mwef_sla_plot<- ggplot(sla_accuracy_mwef, aes(x = models, y = accuracy, 
                                              fill = models)) +
  geom_boxplot() +
  labs(title = "Boxplot of prediction Accuracy using mw training and ef val for sla",
       x = "Models",
       y = "Accuracy")+
  scale_y_continuous("prediction accuracy",limits=c(-0.2,0.55), breaks=c(0, 0.25, 0.5))

#visualizing model accuracy training on mw and validating on ef location for sla
corr_N_wholewave_mwef<- read.csv("./output/wholewave/corr_N_wholewave_mwef.csv")%>%
  rename(wholewave= result_narea_wholewave_mwef.ac)
corr_N_highh2_mwef<- read.csv("./output/highh2/corr_N_highh2_mwef.csv")%>% 
  rename(high_heritablewave= result_narea_highh2_mwef.ac)
corr_N_nirs_mwef<- read.csv("./output/NIRS/ corr_N_nirs_mwef.csv")%>% 
  rename(nirswave= result_narea_nirs_mwef.ac)
corr_N_mwef<- read.csv("./output/GBLUP/ corr_N_mwef.csv")%>%
  rename(GBLUP = result_narea_mwef.ac)
narea_accuracy_mwef<- bind_cols(corr_N_mwef,corr_N_nirs_mwef,corr_N_highh2_mwef,
                                corr_N_wholewave_mwef )%>%
  pivot_longer(cols = 1:4,names_to = "models", values_to = "accuracy")
# Create a boxplot
mwef_narea_plot<- ggplot(narea_accuracy_mwef, aes(x = models, y = accuracy, 
                                                  fill = models)) +
  geom_boxplot() +
  labs(title = "Boxplot of prediction Accuracy using mw training and ef val for narea",
       x = "Models",
       y = "Accuracy")+
  scale_y_continuous("Prediction Accuracy",limits=c(-0.2,0.55), breaks=c(0, 0.25, 0.5))

#visulaizing different cv scenarios based on trait of interest
narea_plot<- plot_grid(joint_narea_plot,efmw_narea_plot,mwef_narea_plot, ncol=3)
sla_plot<- plot_grid(joint_sla_plot,efmw_sla_plot,mwef_sla_plot, ncol=3)

#creating boxplot comparing GBLUP, Gh2, GWW, Gnirs for 10% for narea joint loc
Gh2<- read.csv("./output/Gh2_joint_narea/corrN_Gh2_10j.csv")%>% rename(Gh2= Gh2_10)
GBLUP<- read.csv("./output/GBLUP_reduced/corrN_10j.csv") %>% rename(GBLUP= G_10)
GWW<- read.csv("./output/GWW_joint_narea/corrN_GWW_10j.csv") %>% rename(GWW= GWW_10)
Gnirs<- read.csv("./output/Gnirs_joint_narea/corrN_Gnirs_10j.csv") %>% rename(Gnirs= Gnirs_10)
narea_10_joint <- bind_cols(GBLUP,Gh2,GWW,Gnirs)%>%
  pivot_longer(cols = 1:4,names_to = "models", values_to = "accuracy") %>%  mutate(scheme= "10%")

#creating boxplot comparing GBLUP, Gh2, GWW, Gnirs for 25% for narea joint loc
Gh2<- read.csv("./output/Gh2_joint_narea/corrN_Gh2_25j.csv")%>% rename(Gh2= Gh2_25)
GBLUP<- read.csv("./output/GBLUP_reduced/corrN_25j.csv") %>% rename(GBLUP= G_25)
GWW<- read.csv("./output/GWW_joint_narea/corrN_GWW_25j.csv") %>% rename(GWW= GWW_25)
Gnirs<- read.csv("./output/Gnirs_joint_narea/corrN_Gnirs_25j.csv") %>% rename(Gnirs= Gnirs_25)
narea_25_joint <- bind_cols(GBLUP,Gh2,GWW,Gnirs)%>%
  pivot_longer(cols = 1:4,names_to = "models", values_to = "accuracy") %>%  mutate(scheme= "25%")

#creating boxplot comparing GBLUP, Gh2, GWW, Gnirs for 50% for narea joint loc
Gh2<- read.csv("./output/Gh2_joint_narea/corrN_Gh2_50j.csv")%>% rename(Gh2= Gh2_50)
GBLUP<- read.csv("./output/GBLUP_reduced/corrN_50j.csv") %>% rename(GBLUP= G_50)
GWW<- read.csv("./output/GWW_joint_narea/corrN_GWW_50j.csv") %>% rename(GWW= GWW_50)
Gnirs<- read.csv("./output/Gnirs_joint_narea/corrN_Gnirs_50j.csv") %>% rename(Gnirs= Gnirs_50)
narea_50_joint <- bind_cols(GBLUP,Gh2,GWW,Gnirs)%>%
  pivot_longer(cols = 1:4,names_to = "models", values_to = "accuracy") %>%  mutate(scheme= "50%")
narea_joint <- bind_rows(narea_10_joint,narea_25_joint,narea_50_joint)
#creating plot for three different scheme 
narea_joint <-ggplot(narea_joint, aes(x = models, y = accuracy, fill = models)) +
  geom_boxplot() +
  facet_wrap(~ scheme, scales = "free", nrow = 1) +
  labs(title = "Boxplot of Accuracy by Scheme",
       x = "Models",
       y = "Prediction Accuracy") +
  scale_y_continuous("Prediction Accuracy", limits = c(0, 0.55), breaks = seq(0, 0.5, 0.1)) +
  theme(aspect.ratio = 0.7)
ggsave("./figures/narea_joint.png", plot = narea_joint, width = 15, height = 6, units = "in")

#creating boxplot comparing GBLUP, Gh2, GWW, Gnirs for 10% for sla joint loc
Gh2<- read.csv("./output/Gh2_joint_sla/corrS_Gh2_10j.csv")%>% rename(Gh2= Gh2_10)
GBLUP<- read.csv("./output/GBLUP_reduced/corrS_10j.csv") %>% rename(GBLUP= G_10)
GWW<- read.csv("./output/GWW_joint_sla/corrS_GWW_10j.csv") %>% rename(GWW= GWW_10)
Gnirs<- read.csv("./output/Gnirs_joint_sla/corrS_Gnirs_10j.csv") %>% rename(Gnirs= Gnirs_10)
sla_10_joint <- bind_cols(GBLUP,Gh2,GWW,Gnirs)%>%
  pivot_longer(cols = 1:4,names_to = "models", values_to = "accuracy") %>%  mutate(scheme= "10%")

#creating boxplot comparing GBLUP, Gh2, GWW, Gnirs for 25% for sla joint loc
Gh2<- read.csv("./output/Gh2_joint_sla/corrS_Gh2_25j.csv")%>% rename(Gh2= Gh2_25)
GBLUP<- read.csv("./output/GBLUP_reduced/corrS_25j.csv") %>% rename(GBLUP= G_25)
GWW<- read.csv("./output/GWW_joint_sla/corrS_GWW_25j.csv") %>% rename(GWW= GWW_25)
Gnirs<- read.csv("./output/Gnirs_joint_sla/corrS_Gnirs_25j.csv") %>% rename(Gnirs= Gnirs_25)
sla_25_joint <- bind_cols(GBLUP,Gh2,GWW,Gnirs)%>%
  pivot_longer(cols = 1:4,names_to = "models", values_to = "accuracy") %>%  mutate(scheme= "25%")


#creating boxplot comparing GBLUP, Gh2, GWW, Gnirs for 50% for sla joint loc
Gh2<- read.csv("./output/Gh2_joint_sla/corrS_Gh2_50j.csv")%>% rename(Gh2= Gh2_50)
GBLUP<- read.csv("./output/GBLUP_reduced/corrS_50j.csv") %>% rename(GBLUP= G_50)
GWW<- read.csv("./output/GWW_joint_sla/corrS_GWW_50j.csv") %>% rename(GWW= GWW_50)
Gnirs<- read.csv("./output/Gnirs_joint_sla/corrS_Gnirs_50j.csv") %>% rename(Gnirs= Gnirs_50)
sla_50_joint <- bind_cols(GBLUP,Gh2,GWW,Gnirs)%>%
  pivot_longer(cols = 1:4,names_to = "models", values_to = "accuracy") %>%  mutate(scheme= "50%")

sla_joint <- bind_rows(sla_10_joint,sla_25_joint, sla_50_joint)

#creating plot for three different scheme 
sla_joint<-ggplot(sla_joint, aes(x = models, y = accuracy, fill = models)) +
  geom_boxplot() +
  facet_wrap(~ scheme, scales = "free", nrow = 1) +
  labs(title = "Boxplot of Accuracy by Scheme",
       x = "Models",
       y = "Prediction Accuracy") +
  scale_y_continuous("Prediction Accuracy", limits = c(0, 0.55), breaks = seq(0, 0.5, 0.1)) +
  theme(aspect.ratio = 0.7)

ggsave("./figures/sla_joint.png", plot = sla_joint, width = 15, height = 6, units = "in")


#creating boxplot comparing GBLUP, Gh2, GWW, Gnirs for 10% for sla efmw loc
Gh2<- read.csv("./output/Gh2_efmw_sla/corrS_Gh2_10efmw.csv")%>% rename(Gh2= Gh2_10)
GBLUP<- read.csv("./output/GBLUP_reduced/corrS_10efmw.csv") %>% rename(GBLUP= G_10)
GWW<- read.csv("./output/GWW_efmw_sla/corrS_GWW_10efmw.csv") %>% rename(GWW= GWW_10)
Gnirs<- read.csv("./output/Gnirs_efmw_sla/corrS_Gnirs_10efmw.csv") %>% rename(Gnirs= Gnirs_10)
sla_10_efmw <- bind_cols(GBLUP,Gh2,GWW,Gnirs)%>%
  pivot_longer(cols = 1:4,names_to = "models", values_to = "accuracy") %>%  mutate(scheme= "10%")

#creating boxplot comparing GBLUP, Gh2, GWW, Gnirs for 25% for sla joint loc
Gh2<- read.csv("./output/Gh2_efmw_sla/corrS_Gh2_25efmw.csv")%>% rename(Gh2= Gh2_25)
GBLUP<- read.csv("./output/GBLUP_reduced/corrS_25efmw.csv") %>% rename(GBLUP= G_25)
GWW<- read.csv("./output/GWW_efmw_sla/corrS_GWW_25efmw.csv") %>% rename(GWW= GWW_25)
Gnirs<- read.csv("./output/Gnirs_efmw_sla/corrS_Gnirs_25efmw.csv") %>% rename(Gnirs= Gnirs_25)
sla_25_efmw <- bind_cols(GBLUP,Gh2,GWW,Gnirs)%>%
  pivot_longer(cols = 1:4,names_to = "models", values_to = "accuracy") %>%  mutate(scheme= "25%")

#creating boxplot comparing GBLUP, Gh2, GWW, Gnirs for 50% for sla joint loc
Gh2<- read.csv("./output/Gh2_efmw_sla/corrS_Gh2_50efmw.csv")%>% rename(Gh2= Gh2_50)
GBLUP<- read.csv("./output/GBLUP_reduced/corrS_50efmw.csv") %>% rename(GBLUP= G_50)
GWW<- read.csv("./output/GWW_efmw_sla/corrS_GWW_50efmw.csv") %>% rename(GWW= GWW_50)
Gnirs<- read.csv("./output/Gnirs_efmw_sla/corrS_Gnirs_50efmw.csv") %>% rename(Gnirs= Gnirs_50)
sla_50_efmw <- bind_cols(GBLUP,Gh2,GWW,Gnirs)%>%
  pivot_longer(cols = 1:4,names_to = "models", values_to = "accuracy") %>%  mutate(scheme= "50%")

sla_efmw <- bind_rows(sla_10_efmw,sla_25_efmw, sla_50_efmw)

#creating plot for three different scheme 
sla_efmw<-ggplot(sla_efmw, aes(x = models, y = accuracy, fill = models)) +
  geom_boxplot() +
  facet_wrap(~ scheme, scales = "free", nrow = 1) +
  labs(title = "Boxplot of Accuracy by Scheme",
       x = "Models",
       y = "Prediction Accuracy") +
  scale_y_continuous("Prediction Accuracy", limits = c(0, 0.55), breaks = seq(0, 0.5, 0.1)) +
  theme(aspect.ratio = 0.7)

ggsave("./figures/sla_efmw.png", plot = sla_efmw, width = 15, height = 6, units = "in")

#creating boxplot comparing GBLUP, Gh2, GWW, Gnirs for 10% for sla mwef loc
Gh2<- read.csv("./output/Gh2_mwef_sla/corrS_Gh2_10mwef.csv")%>% rename(Gh2= Gh2_10)
GBLUP<- read.csv("./output/GBLUP_reduced/corrS_10mwef.csv") %>% rename(GBLUP= G_10)
GWW<- read.csv("./output/GWW_mwef_sla/corrS_GWW_10mwef.csv") %>% rename(GWW= GWW_10)
Gnirs<- read.csv("./output/Gnirs_mwef_sla/corrS_Gnirs_10mwef.csv") %>% rename(Gnirs= Gnirs_10)
sla_10_mwef <- bind_cols(GBLUP,Gh2,GWW,Gnirs)%>%
  pivot_longer(cols = 1:4,names_to = "models", values_to = "accuracy") %>%  mutate(scheme= "10%")

#creating boxplot comparing GBLUP, Gh2, GWW, Gnirs for 25% for sla joint loc
Gh2<- read.csv("./output/Gh2_mwef_sla/corrS_Gh2_25mwef.csv")%>% rename(Gh2= Gh2_25)
GBLUP<- read.csv("./output/GBLUP_reduced/corrS_25mwef.csv") %>% rename(GBLUP= G_25)
GWW<- read.csv("./output/GWW_mwef_sla/corrS_GWW_25mwef.csv") %>% rename(GWW= GWW_25)
Gnirs<- read.csv("./output/Gnirs_mwef_sla/corrS_Gnirs_25mwef.csv") %>% rename(Gnirs= Gnirs_25)
sla_25_mwef <- bind_cols(GBLUP,Gh2,GWW,Gnirs)%>%
  pivot_longer(cols = 1:4,names_to = "models", values_to = "accuracy") %>%  mutate(scheme= "25%")

#creating boxplot comparing GBLUP, Gh2, GWW, Gnirs for 50% for sla joint loc
Gh2<- read.csv("./output/Gh2_mwef_sla/corrS_Gh2_50mwef.csv")%>% rename(Gh2= Gh2_50)
GBLUP<- read.csv("./output/GBLUP_reduced/corrS_50mwef.csv") %>% rename(GBLUP= G_50)
GWW<- read.csv("./output/GWW_mwef_sla/corrS_GWW_50mwef.csv") %>% rename(GWW= GWW_50)
Gnirs<- read.csv("./output/Gnirs_mwef_sla/corrS_Gnirs_50mwef.csv") %>% rename(Gnirs= Gnirs_50)
sla_50_mwef <- bind_cols(GBLUP,Gh2,GWW,Gnirs)%>%
  pivot_longer(cols = 1:4,names_to = "models", values_to = "accuracy") %>%  mutate(scheme= "50%")

sla_mwef <- bind_rows(sla_10_mwef,sla_25_mwef, sla_50_mwef)

#creating plot for three different scheme 
sla_mwef<-ggplot(sla_mwef, aes(x = models, y = accuracy, fill = models)) +
  geom_boxplot() +
  facet_wrap(~ scheme, scales = "free", nrow = 1) +
  labs(title = "Boxplot of Accuracy by Scheme",
       x = "Models",
       y = "Prediction Accuracy") +
  scale_y_continuous("Prediction Accuracy", limits = c(0, 0.55), breaks = seq(0, 0.5, 0.1)) +
  theme(aspect.ratio = 0.7)

ggsave("./figures/sla_mwef.png", plot = sla_mwef, width = 15, height = 6, units = "in")

#creating boxplot comparing GBLUP, Gh2, GWW, Gnirs for 10% for narea efmw loc
Gh2<- read.csv("./output/Gh2_efmw_narea/corrN_Gh2_10efmw.csv")%>% rename(Gh2= Gh2_10)
GBLUP<- read.csv("./output/GBLUP_reduced/corrN_10efmw.csv") %>% rename(GBLUP= G_10)
GWW<- read.csv("./output/GWW_efmw_narea/corrN_GWW_10efmw.csv") %>% rename(GWW= GWW_10)
Gnirs<- read.csv("./output/Gnirs_efmw_narea/corrN_Gnirs_10efmw.csv") %>% rename(Gnirs= Gnirs_10)
narea_10_efmw <- bind_cols(GBLUP,Gh2,GWW,Gnirs)%>%
  pivot_longer(cols = 1:4,names_to = "models", values_to = "accuracy") %>%  mutate(scheme= "10%")

#creating boxplot comparing GBLUP, Gh2, GWW, Gnirs for 25% for narea efmw loc
Gh2<- read.csv("./output/Gh2_efmw_narea/corrN_Gh2_25efmw.csv")%>% rename(Gh2= Gh2_25)
GBLUP<- read.csv("./output/GBLUP_reduced/corrN_25efmw.csv") %>% rename(GBLUP= G_25)
GWW<- read.csv("./output/GWW_efmw_narea/corrN_GWW_25efmw.csv") %>% rename(GWW= GWW_25)
Gnirs<- read.csv("./output/Gnirs_efmw_narea/corrN_Gnirs_25efmw.csv") %>% rename(Gnirs= Gnirs_25)
narea_25_efmw <- bind_cols(GBLUP,Gh2,GWW,Gnirs)%>%
  pivot_longer(cols = 1:4,names_to = "models", values_to = "accuracy") %>%  mutate(scheme= "25%")

#creating boxplot comparing GBLUP, Gh2, GWW, Gnirs for 50% for narea efmw loc
Gh2<- read.csv("./output/Gh2_efmw_narea/corrN_Gh2_50efmw.csv")%>% rename(Gh2= Gh2_50)
GBLUP<- read.csv("./output/GBLUP_reduced/corrN_50efmw.csv") %>% rename(GBLUP= G_50)
GWW<- read.csv("./output/GWW_efmw_narea/corrN_GWW_50efmw.csv") %>% rename(GWW= GWW_50)
Gnirs<- read.csv("./output/Gnirs_efmw_narea/corrN_Gnirs_50efmw.csv") %>% rename(Gnirs= Gnirs_50)
narea_50_efmw <- bind_cols(GBLUP,Gh2,GWW,Gnirs)%>%
  pivot_longer(cols = 1:4,names_to = "models", values_to = "accuracy") %>%  mutate(scheme= "50%")
narea_efmw <- bind_rows(narea_10_efmw,narea_25_efmw,narea_50_efmw)
#creating plot for three different scheme 
narea_efmw <-ggplot(narea_efmw, aes(x = models, y = accuracy, fill = models)) +
  geom_boxplot() +
  facet_wrap(~ scheme, scales = "free", nrow = 1) +
  labs(title = "Boxplot of Accuracy by Scheme",
       x = "Models",
       y = "Prediction Accuracy") +
  scale_y_continuous("Prediction Accuracy", limits = c(0, 0.55), breaks = seq(0, 0.5, 0.1)) +
  theme(aspect.ratio = 0.7)
ggsave("./figures/narea_efmw.png", plot = narea_efmw, width = 15, height = 6, units = "in")

#creating boxplot comparing GBLUP, Gh2, GWW, Gnirs for 10% for narea mwef loc
Gh2<- read.csv("./output/Gh2_mwef_narea/corrN_Gh2_10mwef.csv")%>% rename(Gh2= Gh2_10)
GBLUP<- read.csv("./output/GBLUP_reduced/corrN_10mwef.csv") %>% rename(GBLUP= G_10)
GWW<- read.csv("./output/GWW_mwef_narea/corrN_GWW_10mwef.csv") %>% rename(GWW= GWW_10)
Gnirs<- read.csv("./output/Gnirs_mwef_narea/corrN_Gnirs_10mwef.csv") %>% rename(Gnirs= Gnirs_10)
narea_10_mwef <- bind_cols(GBLUP,Gh2,GWW,Gnirs)%>%
  pivot_longer(cols = 1:4,names_to = "models", values_to = "accuracy") %>%  mutate(scheme= "10%")

#creating boxplot comparing GBLUP, Gh2, GWW, Gnirs for 25% for narea mwef loc
Gh2<- read.csv("./output/Gh2_mwef_narea/corrN_Gh2_25mwef.csv")%>% rename(Gh2= Gh2_25)
GBLUP<- read.csv("./output/GBLUP_reduced/corrN_25mwef.csv") %>% rename(GBLUP= G_25)
GWW<- read.csv("./output/GWW_mwef_narea/corrN_GWW_25mwef.csv") %>% rename(GWW= GWW_25)
Gnirs<- read.csv("./output/Gnirs_mwef_narea/corrN_Gnirs_25mwef.csv") %>% rename(Gnirs= Gnirs_25)
narea_25_mwef <- bind_cols(GBLUP,Gh2,GWW,Gnirs)%>%
  pivot_longer(cols = 1:4,names_to = "models", values_to = "accuracy") %>%  mutate(scheme= "25%")

#creating boxplot comparing GBLUP, Gh2, GWW, Gnirs for 50% for narea mwef loc
Gh2<- read.csv("./output/Gh2_mwef_narea/corrN_Gh2_50mwef.csv")%>% rename(Gh2= Gh2_50)
GBLUP<- read.csv("./output/GBLUP_reduced/corrN_50mwef.csv") %>% rename(GBLUP= G_50)
GWW<- read.csv("./output/GWW_mwef_narea/corrN_GWW_50mwef.csv") %>% rename(GWW= GWW_50)
Gnirs<- read.csv("./output/Gnirs_mwef_narea/corrN_Gnirs_50mwef.csv") %>% rename(Gnirs= Gnirs_50)
narea_50_mwef <- bind_cols(GBLUP,Gh2,GWW,Gnirs)%>%
  pivot_longer(cols = 1:4,names_to = "models", values_to = "accuracy") %>%  mutate(scheme= "50%")
narea_mwef <- bind_rows(narea_10_mwef,narea_25_mwef,narea_50_mwef)
#creating plot for three different scheme 
narea_mwef <-ggplot(narea_mwef, aes(x = models, y = accuracy, fill = models)) +
  geom_boxplot() +
  facet_wrap(~ scheme, scales = "free", nrow = 1) +
  labs(title = "Boxplot of Accuracy by Scheme",
       x = "Models",
       y = "Prediction Accuracy") +
  scale_y_continuous("Prediction Accuracy", limits = c(0, 0.55), breaks = seq(0, 0.5, 0.1)) +
  theme(aspect.ratio = 0.7)
ggsave("./figures/narea_mwef.png", plot = narea_mwef, width = 15, height = 6, units = "in")



#compare model nrmse GBLUP, Gh2, GWW, Gnirs for 10% joint for narea
GBLUP_reduced10<- read.csv("./output/GBLUP_reduced/narea_joint10.csv")
Gh2_10<- read.csv("./output/Gh2_joint_narea/narea_Gh2_10_joint.csv")
Gnirs_10<- read.csv("./output/Gnirs_joint_narea/narea_Gnirs_10_joint.csv")
GWW_10<- read.csv("./output/GWW_joint_narea/narea_GWW_10_joint.csv")
nrmse_GBLUP<- calculate_nrmse(observed= GBLUP_reduced10$narea, predicted=GBLUP_reduced10$predicted.value)
nrmse_Gh2<- calculate_nrmse(observed= Gh2_10$narea, predicted=Gh2_10$predicted.value)
nrmse_Gnirs<- calculate_nrmse(observed= Gnirs_10$narea, predicted=Gnirs_10$predicted.value)
nrmse_GWW<- calculate_nrmse(observed= GWW_10$narea, predicted=GWW_10$predicted.value)

#25% for narea joint
GBLUP_reduced25<- read.csv("./output/GBLUP_reduced/narea_joint25.csv")
Gh2_25<- read.csv("./output/Gh2_joint_narea/narea_Gh2_25_joint.csv")
Gnirs_25<- read.csv("./output/Gnirs_joint_narea/narea_Gnirs_25_joint.csv")
GWW_25<- read.csv("./output/GWW_joint_narea/narea_GWW_25_joint.csv")
nrmse_GBLUP<- calculate_nrmse(observed= GBLUP_reduced25$narea, predicted=GBLUP_reduced25$predicted.value)
nrmse_Gh2<- calculate_nrmse(observed= Gh2_25$narea, predicted=Gh2_25$predicted.value)
nrmse_Gnirs<- calculate_nrmse(observed= Gnirs_25$narea, predicted=Gnirs_25$predicted.value)
nrmse_GWW<- calculate_nrmse(observed= GWW_25$narea, predicted=GWW_25$predicted.value)

#50% for narea joint
GBLUP_reduced50<- read.csv("./output/GBLUP_reduced/narea_joint50.csv")
Gh2_50<- read.csv("./output/Gh2_joint_narea/narea_Gh2_50_joint.csv")
Gnirs_50<- read.csv("./output/Gnirs_joint_narea/narea_Gnirs_50_joint.csv")
GWW_50<- read.csv("./output/GWW_joint_narea/narea_GWW_50_joint.csv")
nrmse_GBLUP<- calculate_nrmse(observed= GBLUP_reduced50$narea, predicted=GBLUP_reduced50$predicted.value)
nrmse_Gh2<- calculate_nrmse(observed= Gh2_50$narea, predicted=Gh2_50$predicted.value)
nrmse_Gnirs<- calculate_nrmse(observed= Gnirs_50$narea, predicted=Gnirs_50$predicted.value)
nrmse_GWW<- calculate_nrmse(observed= GWW_50$narea, predicted=GWW_50$predicted.value)

#selecting topk ind from loc ef and mw by reading blues from first stage
#top 20% indv for slablues ef and their mean
temp_data <- slablues_ef_correct %>% arrange(desc(sla)) %>% na.omit()
num_indv<- 165
n1 <- nrow(temp_data)  # total number of genotypes in data 
topk1 <- min(n1, num_indv)  # Select minimum of 165 and the number of genotypes
taxa_blues <- temp_data[1:topk1, ]
mean_value <- mean(taxa_blues[["sla"]], na.rm = TRUE)#309.0880
top_taxa_names <- taxa_blues$taxa
write.csv(top_taxa_names, "./output/selected/top_slablues_ef.csv")

#top 20% indv for slablues mw and their mean
temp_data <- slablues_mw_correct %>% arrange(desc(sla)) %>% na.omit()
num_indv<- 165
n1 <- nrow(temp_data)  # total number of genotypes in data 
topk1 <- min(n1, num_indv)  # Select minimum of 165 and the number of genotypes
taxa_blues <- temp_data[1:topk1, ]
mean_value <- mean(taxa_blues[["sla"]], na.rm = TRUE)#312.0316
top_taxa_names <- taxa_blues$taxa
write.csv(top_taxa_names, "./output/selected/top_slablues_mw.csv")

#top 20% indv for nareablues ef and their mean
temp_data <- nareablues_ef_correct %>% arrange(desc(narea)) %>% na.omit()
num_indv<- 165
n1 <- nrow(temp_data)  # total number of genotypes in data 
topk1 <- min(n1, num_indv)  # Select minimum of 165 and the number of genotypes
taxa_blues <- temp_data[1:topk1, ]
mean_value <- mean(taxa_blues[["narea"]], na.rm = TRUE)#1.5611
top_taxa_names <- taxa_blues$taxa
write.csv(top_taxa_names, "./output/selected/top_nareablues_ef.csv")

#top 20% indv for nareablues mw and their mean
temp_data <- nareablues_mw_correct %>% arrange(desc(narea)) %>% na.omit()
num_indv<- 165
n1 <- nrow(temp_data)  # total number of genotypes in data 
topk1 <- min(n1, num_indv)  # Select minimum of 165 and the number of genotypes
taxa_blues <- temp_data[1:topk1, ]
mean_value <- mean(taxa_blues[["narea"]], na.rm = TRUE)#1.56841
top_taxa_names <- taxa_blues$taxa
write.csv(top_taxa_names, "./output/selected/top_nareablues_mw.csv")

#selecting top20%indv and their means for narea ef from GBLUP_efmw
GBLUP_efmw_narea<- read.csv("./output/GBLUP/narea_efmw.csv") %>% na.omit()
GBLUP_efmw_narea<- top_selected(data= GBLUP_efmw_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_efmw_narea$result_top_taxa), each = length(GBLUP_efmw_narea$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_efmw_narea$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GBLUP_efmw_narea$result_mean_value),
                      mean_value = GBLUP_efmw_narea$result_mean_value)
GBLUP_topn_efmw <- merge(selected_df, mean_df, by = "rep")
write.csv(GBLUP_topn_efmw,"./output/selected/GBLUP_topn_efmw.csv", row.names=F)

#selecting top20%indv and their means for narea mw from GBLUP_mwef
GBLUP_mwef_narea<- read.csv("./output/GBLUP/narea_mwef.csv")
GBLUP_mwef_narea<- top_selected(data= GBLUP_mwef_narea)
#extracting mean and selected individual in a data frame per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_mwef_narea$result_top_taxa), each = length(GBLUP_mwef_narea$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_mwef_narea$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GBLUP_mwef_narea$result_mean_value),
                      mean_value = GBLUP_mwef_narea$result_mean_value)
GBLUP_topn_mwef <- merge(selected_df, mean_df, by = "rep")
write.csv(GBLUP_topn_mwef,"./output/selected/GBLUP_topn_mwef.csv", row.names=F)

#selecting top20%indv and their means for sla ef from GBLUP_efmw
GBLUP_efmw_sla<- read.csv("./output/GBLUP/sla_efmw.csv")
GBLUP_efmw_sla<- top_selected(data= GBLUP_efmw_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_efmw_sla$result_top_taxa), each = length(GBLUP_efmw_sla$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_efmw_sla$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GBLUP_efmw_sla$result_mean_value),
                      mean_value = GBLUP_efmw_sla$result_mean_value)
GBLUP_tops_efmw <- merge(selected_df, mean_df, by = "rep")
write.csv(GBLUP_tops_efmw,"./output/selected/GBLUP_tops_efmw.csv", row.names=F)

#selecting top20%indv and their means for sla mw from GBLUP_mwef
GBLUP_mwef_sla<- read.csv("./output/GBLUP/sla_mwef.csv")
GBLUP_mwef_sla<- top_selected(data= GBLUP_mwef_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_mwef_sla$result_top_taxa), each = length(GBLUP_mwef_sla$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_mwef_sla$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GBLUP_mwef_sla$result_mean_value),
                      mean_value = GBLUP_mwef_sla$result_mean_value)
GBLUP_tops_mwef <- merge(selected_df, mean_df, by = "rep")
write.csv(GBLUP_tops_mwef,"./output/selected/GBLUP_tops_mwef.csv", row.names=F)

#selecting top20%indv and their means for sla ef from Gh2_efmw
Gh2_efmw10_sla<- read.csv("./output/Gh2_efmw_sla/sla_Gh2_10_efmw.csv")
Gh2_efmw10_sla<- top_selected(data= Gh2_efmw10_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gh2_efmw10_sla$result_top_taxa), each = length(Gh2_efmw10_sla$result_top_taxa[[1]])),
                          taxa = unlist(Gh2_efmw10_sla$result_top_taxa))
mean_df <- data.frame(rep = seq_along(Gh2_efmw10_sla$result_mean_value),
                      mean_value = Gh2_efmw10_sla$result_mean_value)
Gh2_tops10_efmw <- merge(selected_df, mean_df, by = "rep")
write.csv(Gh2_tops10_efmw,"./output/selected/Gh2_tops10_efmw.csv", row.names=F)

#selecting top20%indv and their means for sla mw from Gh2_mwef
Gh2_mwef10_sla<- read.csv("./output/Gh2_mwef_sla/sla_Gh2_10_mwef.csv")
Gh2_mwef10_sla<- top_selected(data= Gh2_mwef10_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gh2_mwef10_sla$result_top_taxa), each = length(Gh2_mwef10_sla$result_top_taxa[[1]])),
                          taxa = unlist(Gh2_mwef10_sla$result_top_taxa))
mean_df <- data.frame(rep = seq_along(Gh2_mwef10_sla$result_mean_value),
                      mean_value = Gh2_efmw10_sla$result_mean_value)
Gh2_tops10_mwef <- merge(selected_df, mean_df, by = "rep")
write.csv(Gh2_tops10_mwef,"./output/selected/Gh2_tops10_mwef.csv", row.names=F)

#selecting top20%indv and their means for sla ef from Gh2_efmw
Gh2_efmw25_sla<- read.csv("./output/Gh2_efmw_sla/sla_Gh2_25_efmw.csv")
Gh2_efmw25_sla<- top_selected(data= Gh2_efmw25_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gh2_efmw25_sla$result_top_taxa), each = length(Gh2_efmw25_sla$result_top_taxa[[1]])),
                          taxa = unlist(Gh2_efmw25_sla$result_top_taxa))
mean_df <- data.frame(rep = seq_along(Gh2_efmw25_sla$result_mean_value),
                      mean_value = Gh2_efmw25_sla$result_mean_value)
Gh2_tops25_efmw <- merge(selected_df, mean_df, by = "rep")
write.csv(Gh2_tops25_efmw,"./output/selected/Gh2_tops25_efmw.csv", row.names=F)

#selecting top20%indv and their means for sla mw from Gh2_mwef
Gh2_mwef25_sla<- read.csv("./output/Gh2_mwef_sla/sla_Gh2_25_mwef.csv")
Gh2_mwef25_sla<- top_selected(data= Gh2_mwef25_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gh2_mwef25_sla$result_top_taxa), each = length(Gh2_mwef25_sla$result_top_taxa[[1]])),
                          taxa = unlist(Gh2_mwef25_sla$result_top_taxa))
mean_df <- data.frame(rep = seq_along(Gh2_mwef25_sla$result_mean_value),
                      mean_value = Gh2_efmw25_sla$result_mean_value)
Gh2_tops25_mwef <- merge(selected_df, mean_df, by = "rep")
write.csv(Gh2_tops25_mwef,"./output/selected/Gh2_tops25_mwef.csv", row.names=F)

#selecting top20%indv and their means for sla ef from Gh2_efmw
Gh2_efmw50_sla<- read.csv("./output/Gh2_efmw_sla/sla_Gh2_50_efmw.csv")
Gh2_efmw50_sla<- top_selected(data= Gh2_efmw50_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gh2_efmw50_sla$result_top_taxa), each = length(Gh2_efmw50_sla$result_top_taxa[[1]])),
                          taxa = unlist(Gh2_efmw50_sla$result_top_taxa))
mean_df <- data.frame(rep = seq_along(Gh2_efmw50_sla$result_mean_value),
                      mean_value = Gh2_efmw50_sla$result_mean_value)
Gh2_tops50_efmw <- merge(selected_df, mean_df, by = "rep")
write.csv(Gh2_tops50_efmw,"./output/selected/Gh2_tops50_efmw.csv", row.names=F)

#selecting top20%indv and their means for sla mw from Gh2_mwef
Gh2_mwef50_sla<- read.csv("./output/Gh2_mwef_sla/sla_Gh2_50_mwef.csv")
Gh2_mwef50_sla<- top_selected(data= Gh2_mwef50_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gh2_mwef50_sla$result_top_taxa), each = length(Gh2_mwef50_sla$result_top_taxa[[1]])),
                          taxa = unlist(Gh2_mwef50_sla$result_top_taxa))
mean_df <- data.frame(rep = seq_along(Gh2_mwef50_sla$result_mean_value),
                      mean_value = Gh2_efmw50_sla$result_mean_value)
Gh2_tops50_mwef <- merge(selected_df, mean_df, by = "rep")
write.csv(Gh2_tops50_mwef,"./output/selected/Gh2_tops50_mwef.csv", row.names=F)

#selecting top20%indv and their means for narea ef from Gh2_efmw
Gh2_efmw10_narea<- read.csv("./output/Gh2_efmw_narea/narea_Gh2_10_efmw.csv")
Gh2_efmw10_narea<- top_selected(data= Gh2_efmw10_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gh2_efmw10_narea$result_top_taxa), each = length(Gh2_efmw10_narea$result_top_taxa[[1]])),
                          taxa = unlist(Gh2_efmw10_narea$result_top_taxa))
mean_df <- data.frame(rep = seq_along(Gh2_efmw10_narea$result_mean_value),
                      mean_value = Gh2_efmw10_narea$result_mean_value)
Gh2_topn10_efmw <- merge(selected_df, mean_df, by = "rep")
write.csv(Gh2_topn10_efmw,"./output/selected/Gh2_topn10_efmw.csv", row.names=F)

#selecting top20%indv and their means for narea mw from Gh2_mwef
Gh2_mwef10_narea<- read.csv("./output/Gh2_mwef_narea/narea_Gh2_10_mwef.csv")
Gh2_mwef10_narea<- top_selected(data= Gh2_mwef10_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gh2_mwef10_narea$result_top_taxa), each = length(Gh2_mwef10_narea$result_top_taxa[[1]])),
                          taxa = unlist(Gh2_mwef10_narea$result_top_taxa))
mean_df <- data.frame(rep = seq_along(Gh2_mwef10_narea$result_mean_value),
                      mean_value = Gh2_efmw10_narea$result_mean_value)
Gh2_topn10_mwef <- merge(selected_df, mean_df, by = "rep")
write.csv(Gh2_topn10_mwef,"./output/selected/Gh2_topn10_mwef.csv", row.names=F)



#selecting top20%indv and their means for narea ef from Gh2_efmw
Gh2_efmw25_narea<- read.csv("./output/Gh2_efmw_narea/narea_Gh2_25_efmw.csv")
Gh2_efmw25_narea<- top_selected(data= Gh2_efmw25_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gh2_efmw25_narea$result_top_taxa), each = length(Gh2_efmw25_narea$result_top_taxa[[1]])),
                          taxa = unlist(Gh2_efmw25_narea$result_top_taxa))
mean_df <- data.frame(rep = seq_along(Gh2_efmw25_narea$result_mean_value),
                      mean_value = Gh2_efmw25_narea$result_mean_value)
Gh2_topn25_efmw <- merge(selected_df, mean_df, by = "rep")
write.csv(Gh2_topn25_efmw,"./output/selected/Gh2_topn25_efmw.csv", row.names=F)

#selecting top20%indv and their means for narea mw from Gh2_mwef
Gh2_mwef25_narea<- read.csv("./output/Gh2_mwef_narea/narea_Gh2_25_mwef.csv")
Gh2_mwef25_narea<- top_selected(data= Gh2_mwef25_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gh2_mwef25_narea$result_top_taxa), each = length(Gh2_mwef25_narea$result_top_taxa[[1]])),
                          taxa = unlist(Gh2_mwef25_narea$result_top_taxa))
mean_df <- data.frame(rep = seq_along(Gh2_mwef25_narea$result_mean_value),
                      mean_value = Gh2_efmw25_narea$result_mean_value)
Gh2_topn25_mwef <- merge(selected_df, mean_df, by = "rep")
write.csv(Gh2_topn25_mwef,"./output/selected/Gh2_topn25_mwef.csv", row.names=F)
#selecting top20%indv and their means for narea ef from Gh2_efmw
Gh2_efmw50_narea<- read.csv("./output/Gh2_efmw_narea/narea_Gh2_50_efmw.csv")
Gh2_efmw50_narea<- top_selected(data= Gh2_efmw50_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gh2_efmw50_narea$result_top_taxa), each = length(Gh2_efmw50_narea$result_top_taxa[[1]])),
                          taxa = unlist(Gh2_efmw50_narea$result_top_taxa))
mean_df <- data.frame(rep = seq_along(Gh2_efmw50_narea$result_mean_value),
                      mean_value = Gh2_efmw50_narea$result_mean_value)
Gh2_topn50_efmw <- merge(selected_df, mean_df, by = "rep")
write.csv(Gh2_topn50_efmw,"./output/selected/Gh2_topn50_efmw.csv", row.names=F)

#selecting top20%indv and their means for narea mw from Gh2_mwef
Gh2_mwef50_narea<- read.csv("./output/Gh2_mwef_narea/narea_Gh2_50_mwef.csv")
Gh2_mwef50_narea<- top_selected(data= Gh2_mwef50_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gh2_mwef50_narea$result_top_taxa), each = length(Gh2_mwef50_narea$result_top_taxa[[1]])),
                          taxa = unlist(Gh2_mwef50_narea$result_top_taxa))
mean_df <- data.frame(rep = seq_along(Gh2_mwef50_narea$result_mean_value),
                      mean_value = Gh2_efmw50_narea$result_mean_value)
Gh2_topn50_mwef <- merge(selected_df, mean_df, by = "rep")
write.csv(Gh2_topn50_mwef,"./output/selected/Gh2_topn50_mwef.csv", row.names=F)




#GWW selected individual
#selecting top20%indv and their means for sla ef from GWW efmw
GWW_efmw10_sla<- read.csv("./output/GWW_efmw_sla/sla_GWW_10_efmw.csv")
GWW_efmw10_sla<- top_selected(data= GWW_efmw10_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GWW_efmw10_sla$result_top_taxa), each = length(GWW_efmw10_sla$result_top_taxa[[1]])),
                          taxa = unlist(GWW_efmw10_sla$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GWW_efmw10_sla$result_mean_value),
                      mean_value = GWW_efmw10_sla$result_mean_value)
GWW_tops10_efmw <- merge(selected_df, mean_df, by = "rep")
write.csv(GWW_tops10_efmw,"./output/selected/GWW_tops10_efmw.csv", row.names=F)

#selecting top20%indv and their means for sla mw from GWW_mwef
GWW_mwef10_sla<- read.csv("./output/GWW_mwef_sla/sla_GWW_10_mwef.csv")
GWW_mwef10_sla<- top_selected(data= GWW_mwef10_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GWW_mwef10_sla$result_top_taxa), each = length(GWW_mwef10_sla$result_top_taxa[[1]])),
                          taxa = unlist(GWW_mwef10_sla$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GWW_mwef10_sla$result_mean_value),
                      mean_value = GWW_efmw10_sla$result_mean_value)
GWW_tops10_mwef <- merge(selected_df, mean_df, by = "rep")
write.csv(GWW_tops10_mwef,"./output/selected/GWW_tops10_mwef.csv", row.names=F)

#selecting top20%indv and their means for sla ef from GWW_efmw
GWW_efmw25_sla<- read.csv("./output/GWW_efmw_sla/sla_GWW_25_efmw.csv")
GWW_efmw25_sla<- top_selected(data= GWW_efmw25_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GWW_efmw25_sla$result_top_taxa), each = length(GWW_efmw25_sla$result_top_taxa[[1]])),
                          taxa = unlist(GWW_efmw25_sla$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GWW_efmw25_sla$result_mean_value),
                      mean_value = GWW_efmw25_sla$result_mean_value)
GWW_tops25_efmw <- merge(selected_df, mean_df, by = "rep")
write.csv(GWW_tops25_efmw,"./output/selected/GWW_tops25_efmw.csv", row.names=F)

#selecting top20%indv and their means for sla mw from GWW_mwef
GWW_mwef25_sla<- read.csv("./output/GWW_mwef_sla/sla_GWW_25_mwef.csv")
GWW_mwef25_sla<- top_selected(data= GWW_mwef25_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GWW_mwef25_sla$result_top_taxa), each = length(GWW_mwef25_sla$result_top_taxa[[1]])),
                          taxa = unlist(GWW_mwef25_sla$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GWW_mwef25_sla$result_mean_value),
                      mean_value = GWW_efmw25_sla$result_mean_value)
GWW_tops25_mwef <- merge(selected_df, mean_df, by = "rep")
write.csv(GWW_tops25_mwef,"./output/selected/GWW_tops25_mwef.csv", row.names=F)

#selecting top20%indv and their means for sla ef from GWW_efmw
GWW_efmw50_sla<- read.csv("./output/GWW_efmw_sla/sla_GWW_50_efmw.csv")
GWW_efmw50_sla<- top_selected(data= GWW_efmw50_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GWW_efmw50_sla$result_top_taxa), each = length(GWW_efmw50_sla$result_top_taxa[[1]])),
                          taxa = unlist(GWW_efmw50_sla$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GWW_efmw50_sla$result_mean_value),
                      mean_value = GWW_efmw50_sla$result_mean_value)
GWW_tops50_efmw <- merge(selected_df, mean_df, by = "rep")
write.csv(GWW_tops50_efmw,"./output/selected/GWW_tops50_efmw.csv", row.names=F)

#selecting top20%indv and their means for sla mw from GWW_mwef
GWW_mwef50_sla<- read.csv("./output/GWW_mwef_sla/sla_GWW_50_mwef.csv")
GWW_mwef50_sla<- top_selected(data= GWW_mwef50_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GWW_mwef50_sla$result_top_taxa), each = length(GWW_mwef50_sla$result_top_taxa[[1]])),
                          taxa = unlist(GWW_mwef50_sla$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GWW_mwef50_sla$result_mean_value),
                      mean_value = GWW_efmw50_sla$result_mean_value)
GWW_tops50_mwef <- merge(selected_df, mean_df, by = "rep")
write.csv(GWW_tops50_mwef,"./output/selected/GWW_tops50_mwef.csv", row.names=F)

#selecting top20%indv and their means for narea ef from GWW_efmw
GWW_efmw10_narea<- read.csv("./output/GWW_efmw_narea/narea_GWW_10_efmw.csv")
GWW_efmw10_narea<- top_selected(data= GWW_efmw10_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GWW_efmw10_narea$result_top_taxa), each = length(GWW_efmw10_narea$result_top_taxa[[1]])),
                          taxa = unlist(GWW_efmw10_narea$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GWW_efmw10_narea$result_mean_value),
                      mean_value = GWW_efmw10_narea$result_mean_value)
GWW_topn10_efmw <- merge(selected_df, mean_df, by = "rep")
write.csv(GWW_topn10_efmw,"./output/selected/GWW_topn10_efmw.csv", row.names=F)

#selecting top20%indv and their means for narea mw from GWW_mwef
GWW_mwef10_narea<- read.csv("./output/GWW_mwef_narea/narea_GWW_10_mwef.csv")
GWW_mwef10_narea<- top_selected(data= GWW_mwef10_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GWW_mwef10_narea$result_top_taxa), each = length(GWW_mwef10_narea$result_top_taxa[[1]])),
                          taxa = unlist(GWW_mwef10_narea$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GWW_mwef10_narea$result_mean_value),
                      mean_value = GWW_efmw10_narea$result_mean_value)
GWW_topn10_mwef <- merge(selected_df, mean_df, by = "rep")
write.csv(GWW_topn10_mwef,"./output/selected/GWW_topn10_mwef.csv", row.names=F)

#selecting top20%indv and their means for narea mw from GWW_mwef
GWW_mwef10_narea<- read.csv("./output/GWW_mwef_narea/narea_GWW_10_mwef.csv")
GWW_mwef10_narea<- top_selected(data= GWW_mwef10_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GWW_mwef10_narea$result_top_taxa), each = length(GWW_mwef10_narea$result_top_taxa[[1]])),
                          taxa = unlist(GWW_mwef10_narea$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GWW_mwef10_narea$result_mean_value),
                      mean_value = GWW_efmw10_narea$result_mean_value)
GWW_topn10_mwef <- merge(selected_df, mean_df, by = "rep")
write.csv(GWW_topn10_mwef,"./output/selected/GWW_topn10_mwef.csv", row.names=F)

#selecting top20%indv and their means for narea ef from GWW_efmw
GWW_efmw25_narea<- read.csv("./output/GWW_efmw_narea/narea_GWW_25_efmw.csv")
GWW_efmw25_narea<- top_selected(data= GWW_efmw25_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GWW_efmw25_narea$result_top_taxa), each = length(GWW_efmw25_narea$result_top_taxa[[1]])),
                          taxa = unlist(GWW_efmw25_narea$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GWW_efmw25_narea$result_mean_value),
                      mean_value = GWW_efmw25_narea$result_mean_value)
GWW_topn25_efmw <- merge(selected_df, mean_df, by = "rep")
write.csv(GWW_topn25_efmw,"./output/selected/GWW_topn25_efmw.csv", row.names=F)

#selecting top20%indv and their means for narea mw from GWW_mwef
GWW_mwef25_narea<- read.csv("./output/GWW_mwef_narea/narea_GWW_25_mwef.csv")
GWW_mwef25_narea<- top_selected(data= GWW_mwef25_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GWW_mwef25_narea$result_top_taxa), each = length(GWW_mwef25_narea$result_top_taxa[[1]])),
                          taxa = unlist(GWW_mwef25_narea$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GWW_mwef25_narea$result_mean_value),
                      mean_value = GWW_efmw25_narea$result_mean_value)
GWW_topn25_mwef <- merge(selected_df, mean_df, by = "rep")
write.csv(GWW_topn25_mwef,"./output/selected/GWW_topn25_mwef.csv", row.names=F)

#selecting top20%indv and their means for narea ef from GWW_efmw
GWW_efmw50_narea<- read.csv("./output/GWW_efmw_narea/narea_GWW_50_efmw.csv")
GWW_efmw50_narea<- top_selected(data= GWW_efmw50_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GWW_efmw50_narea$result_top_taxa), each = length(GWW_efmw50_narea$result_top_taxa[[1]])),
                          taxa = unlist(GWW_efmw50_narea$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GWW_efmw50_narea$result_mean_value),
                      mean_value = GWW_efmw50_narea$result_mean_value)
GWW_topn50_efmw <- merge(selected_df, mean_df, by = "rep")
write.csv(GWW_topn50_efmw,"./output/selected/GWW_topn50_efmw.csv", row.names=F)

#selecting top20%indv and their means for narea mw from GWW_mwef
GWW_mwef50_narea<- read.csv("./output/GWW_mwef_narea/narea_GWW_50_mwef.csv")
GWW_mwef50_narea<- top_selected(data= GWW_mwef50_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GWW_mwef50_narea$result_top_taxa), each = length(GWW_mwef50_narea$result_top_taxa[[1]])),
                          taxa = unlist(GWW_mwef50_narea$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GWW_mwef50_narea$result_mean_value),
                      mean_value = GWW_efmw50_narea$result_mean_value)
GWW_topn50_mwef <- merge(selected_df, mean_df, by = "rep")
write.csv(GWW_topn50_mwef,"./output/selected/GWW_topn50_mwef.csv", row.names=F)

#Gnirs individual selection
#selecting top20%indv and their means for sla ef from Gnirs efmw
Gnirs_efmw10_sla<- read.csv("./output/Gnirs_efmw_sla/sla_Gnirs_10_efmw.csv")
Gnirs_efmw10_sla<- top_selected(data= Gnirs_efmw10_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gnirs_efmw10_sla$result_top_taxa), each = length(Gnirs_efmw10_sla$result_top_taxa[[1]])),
                          taxa = unlist(Gnirs_efmw10_sla$result_top_taxa))
mean_df <- data.frame(rep = seq_along(Gnirs_efmw10_sla$result_mean_value),
                      mean_value = Gnirs_efmw10_sla$result_mean_value)
Gnirs_tops10_efmw <- merge(selected_df, mean_df, by = "rep")
write.csv(Gnirs_tops10_efmw,"./output/selected/Gnirs_tops10_efmw.csv", row.names=F)

#selecting top20%indv and their means for sla mw from Gnirs_mwef
Gnirs_mwef10_sla<- read.csv("./output/Gnirs_mwef_sla/sla_Gnirs_10_mwef.csv")
Gnirs_mwef10_sla<- top_selected(data= Gnirs_mwef10_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gnirs_mwef10_sla$result_top_taxa), each = length(Gnirs_mwef10_sla$result_top_taxa[[1]])),
                          taxa = unlist(Gnirs_mwef10_sla$result_top_taxa))
mean_df <- data.frame(rep = seq_along(Gnirs_mwef10_sla$result_mean_value),
                      mean_value = Gnirs_efmw10_sla$result_mean_value)
Gnirs_tops10_mwef <- merge(selected_df, mean_df, by = "rep")
write.csv(Gnirs_tops10_mwef,"./output/selected/Gnirs_tops10_mwef.csv", row.names=F)

#selecting top20%indv and their means for sla ef from Gnirs_efmw
Gnirs_efmw25_sla<- read.csv("./output/Gnirs_efmw_sla/sla_Gnirs_25_efmw.csv")
Gnirs_efmw25_sla<- top_selected(data= Gnirs_efmw25_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gnirs_efmw25_sla$result_top_taxa), each = length(Gnirs_efmw25_sla$result_top_taxa[[1]])),
                          taxa = unlist(Gnirs_efmw25_sla$result_top_taxa))
mean_df <- data.frame(rep = seq_along(Gnirs_efmw25_sla$result_mean_value),
                      mean_value = Gnirs_efmw25_sla$result_mean_value)
Gnirs_tops25_efmw <- merge(selected_df, mean_df, by = "rep")
write.csv(Gnirs_tops25_efmw,"./output/selected/Gnirs_tops25_efmw.csv", row.names=F)

#selecting top20%indv and their means for sla mw from Gnirs_mwef
Gnirs_mwef25_sla<- read.csv("./output/Gnirs_mwef_sla/sla_Gnirs_25_mwef.csv")
Gnirs_mwef25_sla<- top_selected(data= Gnirs_mwef25_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gnirs_mwef25_sla$result_top_taxa), each = length(Gnirs_mwef25_sla$result_top_taxa[[1]])),
                          taxa = unlist(Gnirs_mwef25_sla$result_top_taxa))
mean_df <- data.frame(rep = seq_along(Gnirs_mwef25_sla$result_mean_value),
                      mean_value = Gnirs_efmw25_sla$result_mean_value)
Gnirs_tops25_mwef <- merge(selected_df, mean_df, by = "rep")
write.csv(Gnirs_tops25_mwef,"./output/selected/Gnirs_tops25_mwef.csv", row.names=F)

#selecting top20%indv and their means for sla ef from Gnirs_efmw
Gnirs_efmw50_sla<- read.csv("./output/Gnirs_efmw_sla/sla_Gnirs_50_efmw.csv")
Gnirs_efmw50_sla<- top_selected(data= Gnirs_efmw50_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gnirs_efmw50_sla$result_top_taxa), each = length(Gnirs_efmw50_sla$result_top_taxa[[1]])),
                          taxa = unlist(Gnirs_efmw50_sla$result_top_taxa))
mean_df <- data.frame(rep = seq_along(Gnirs_efmw50_sla$result_mean_value),
                      mean_value = Gnirs_efmw50_sla$result_mean_value)
Gnirs_tops50_efmw <- merge(selected_df, mean_df, by = "rep")
write.csv(Gnirs_tops50_efmw,"./output/selected/Gnirs_tops50_efmw.csv", row.names=F)

#selecting top20%indv and their means for sla mw from Gnirs_mwef
Gnirs_mwef50_sla<- read.csv("./output/Gnirs_mwef_sla/sla_Gnirs_50_mwef.csv")
Gnirs_mwef50_sla<- top_selected(data= Gnirs_mwef50_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gnirs_mwef50_sla$result_top_taxa), each = length(Gnirs_mwef50_sla$result_top_taxa[[1]])),
                          taxa = unlist(Gnirs_mwef50_sla$result_top_taxa))
mean_df <- data.frame(rep = seq_along(Gnirs_mwef50_sla$result_mean_value),
                      mean_value = Gnirs_efmw50_sla$result_mean_value)
Gnirs_tops50_mwef <- merge(selected_df, mean_df, by = "rep")
write.csv(Gnirs_tops50_mwef,"./output/selected/Gnirs_tops50_mwef.csv", row.names=F)

#selecting top20%indv and their means for narea ef from Gnirs_efmw
Gnirs_efmw10_narea<- read.csv("./output/Gnirs_efmw_narea/narea_Gnirs_10_efmw.csv")
Gnirs_efmw10_narea<- top_selected(data= Gnirs_efmw10_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gnirs_efmw10_narea$result_top_taxa), each = length(Gnirs_efmw10_narea$result_top_taxa[[1]])),
                          taxa = unlist(Gnirs_efmw10_narea$result_top_taxa))
mean_df <- data.frame(rep = seq_along(Gnirs_efmw10_narea$result_mean_value),
                      mean_value = Gnirs_efmw10_narea$result_mean_value)
Gnirs_topn10_efmw <- merge(selected_df, mean_df, by = "rep")
write.csv(Gnirs_topn10_efmw,"./output/selected/Gnirs_topn10_efmw.csv", row.names=F)

#selecting top20%indv and their means for narea mw from Gnirs_mwef
Gnirs_mwef10_narea<- read.csv("./output/Gnirs_mwef_narea/narea_Gnirs_10_mwef.csv")
Gnirs_mwef10_narea<- top_selected(data= Gnirs_mwef10_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gnirs_mwef10_narea$result_top_taxa), each = length(Gnirs_mwef10_narea$result_top_taxa[[1]])),
                          taxa = unlist(Gnirs_mwef10_narea$result_top_taxa))
mean_df <- data.frame(rep = seq_along(Gnirs_mwef10_narea$result_mean_value),
                      mean_value = Gnirs_efmw10_narea$result_mean_value)
Gnirs_topn10_mwef <- merge(selected_df, mean_df, by = "rep")
write.csv(Gnirs_topn10_mwef,"./output/selected/Gnirs_topn10_mwef.csv", row.names=F)

#selecting top20%indv and their means for narea mw from Gnirs_mwef
Gnirs_mwef10_narea<- read.csv("./output/Gnirs_mwef_narea/narea_Gnirs_10_mwef.csv")
Gnirs_mwef10_narea<- top_selected(data= Gnirs_mwef10_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gnirs_mwef10_narea$result_top_taxa), each = length(Gnirs_mwef10_narea$result_top_taxa[[1]])),
                          taxa = unlist(Gnirs_mwef10_narea$result_top_taxa))
mean_df <- data.frame(rep = seq_along(Gnirs_mwef10_narea$result_mean_value),
                      mean_value = Gnirs_efmw10_narea$result_mean_value)
Gnirs_topn10_mwef <- merge(selected_df, mean_df, by = "rep")
write.csv(Gnirs_topn10_mwef,"./output/selected/Gnirs_topn10_mwef.csv", row.names=F)

#selecting top20%indv and their means for narea ef from Gnirs_efmw
Gnirs_efmw25_narea<- read.csv("./output/Gnirs_efmw_narea/narea_Gnirs_25_efmw.csv")
Gnirs_efmw25_narea<- top_selected(data= Gnirs_efmw25_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gnirs_efmw25_narea$result_top_taxa), each = length(Gnirs_efmw25_narea$result_top_taxa[[1]])),
                          taxa = unlist(Gnirs_efmw25_narea$result_top_taxa))
mean_df <- data.frame(rep = seq_along(Gnirs_efmw25_narea$result_mean_value),
                      mean_value = Gnirs_efmw25_narea$result_mean_value)
Gnirs_topn25_efmw <- merge(selected_df, mean_df, by = "rep")
write.csv(Gnirs_topn25_efmw,"./output/selected/Gnirs_topn25_efmw.csv", row.names=F)

#selecting top 165 indv and their means for narea mw from Gnirs_mwef
Gnirs_mwef25_narea<- read.csv("./output/Gnirs_mwef_narea/narea_Gnirs_25_mwef.csv")
Gnirs_mwef25_narea<- top_selected(data= Gnirs_mwef25_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gnirs_mwef25_narea$result_top_taxa), each = length(Gnirs_mwef25_narea$result_top_taxa[[1]])),
                          taxa = unlist(Gnirs_mwef25_narea$result_top_taxa))
mean_df <- data.frame(rep = seq_along(Gnirs_mwef25_narea$result_mean_value),
                      mean_value = Gnirs_efmw25_narea$result_mean_value)
Gnirs_topn25_mwef <- merge(selected_df, mean_df, by = "rep")
write.csv(Gnirs_topn25_mwef,"./output/selected/Gnirs_topn25_mwef.csv", row.names=F)

#
#selecting top20%indv and their means for narea ef from Gnirs_efmw
Gnirs_efmw50_narea<- read.csv("./output/Gnirs_efmw_narea/narea_Gnirs_50_efmw.csv")
Gnirs_efmw50_narea<- top_selected(data= Gnirs_efmw50_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gnirs_efmw50_narea$result_top_taxa), each = length(Gnirs_efmw50_narea$result_top_taxa[[1]])),
                          taxa = unlist(Gnirs_efmw50_narea$result_top_taxa))
mean_df <- data.frame(rep = seq_along(Gnirs_efmw50_narea$result_mean_value),
                      mean_value = Gnirs_efmw50_narea$result_mean_value)
Gnirs_topn50_efmw <- merge(selected_df, mean_df, by = "rep")
write.csv(Gnirs_topn50_efmw,"./output/selected/Gnirs_topn50_efmw.csv", row.names=F)

#selecting top20%indv and their means for narea mw from Gnirs_mwef
Gnirs_mwef50_narea<- read.csv("./output/Gnirs_mwef_narea/narea_Gnirs_50_mwef.csv")
Gnirs_mwef50_narea<- top_selected(data= Gnirs_mwef50_narea)
#extracting mean and selected individual in a data frame per replication
selected_df <- data.frame(rep = rep(seq_along(Gnirs_mwef50_narea$result_top_taxa), each = length(Gnirs_mwef50_narea$result_top_taxa[[1]])),
                          taxa = unlist(Gnirs_mwef50_narea$result_top_taxa))
mean_df <- data.frame(rep = seq_along(Gnirs_mwef50_narea$result_mean_value),
                      mean_value = Gnirs_efmw50_narea$result_mean_value)
Gnirs_topn50_mwef <- merge(selected_df, mean_df, by = "rep")
write.csv(Gnirs_topn50_mwef,"./output/selected/Gnirs_topn50_mwef.csv", row.names=F)
#selecting top 165 indv from GBLUP reduced 
GBLUP_redefmw10_narea<- read.csv("./output/GBLUP_reduced/narea_efmw10.csv") %>% na.omit()
GBLUP_redefmw10_narea<- top_selected(data= GBLUP_redefmw10_narea)
#extracting mean and selected individual in a data frame per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_redefmw10_narea$result_top_taxa), each = length(GBLUP_redefmw10_narea$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_redefmw10_narea$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GBLUP_redefmw10_narea$result_mean_value),
                      mean_value = GBLUP_redefmw10_narea$result_mean_value)
GBLUPred_topn10_efmw <- merge(selected_df, mean_df, by = "rep")
write.csv(GBLUPred_topn10_efmw,"./output/selected/GBLUPred_topn10_efmw.csv", row.names=F)
#selecting top 165 indv from GBLUP reduced for 25 and 50%
GBLUP_redefmw25_narea<- read.csv("./output/GBLUP_reduced/narea_efmw25.csv") %>% na.omit()
GBLUP_redefmw25_narea<- top_selected(data= GBLUP_redefmw25_narea)
#extracting mean and selected individual in a data frame per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_redefmw25_narea$result_top_taxa), each = length(GBLUP_redefmw25_narea$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_redefmw25_narea$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GBLUP_redefmw25_narea$result_mean_value),
                      mean_value = GBLUP_redefmw25_narea$result_mean_value)
GBLUPred_topn25_efmw <- merge(selected_df, mean_df, by = "rep")
write.csv(GBLUPred_topn25_efmw,"./output/selected/GBLUPred_topn25_efmw.csv", row.names=F)
#selecting top 165 indv from GBLUP reduced for ef and 50%
GBLUP_redefmw50_narea<- read.csv("./output/GBLUP_reduced/narea_efmw50.csv") %>% na.omit()
GBLUP_redefmw50_narea<- top_selected(data= GBLUP_redefmw50_narea)
#extracting mean and selected individual in a data frame per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_redefmw50_narea$result_top_taxa), each = length(GBLUP_redefmw50_narea$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_redefmw50_narea$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GBLUP_redefmw50_narea$result_mean_value),
                      mean_value = GBLUP_redefmw50_narea$result_mean_value)
GBLUPred_topn50_efmw <- merge(selected_df, mean_df, by = "rep")


#selecting top 165 indv from GBLUP reduced 
GBLUP_redmwef10_narea<- read.csv("./output/GBLUP_reduced/narea_mwef10.csv") %>% na.omit()
GBLUP_redmwef10_narea<- top_selected(data= GBLUP_redmwef10_narea)
#extracting mean and selected individual in a data frame per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_redmwef10_narea$result_top_taxa), each = length(GBLUP_redmwef10_narea$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_redmwef10_narea$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GBLUP_redmwef10_narea$result_mean_value),
                      mean_value = GBLUP_redmwef10_narea$result_mean_value)
GBLUPred_topn10_mwef <- merge(selected_df, mean_df, by = "rep")
write.csv(GBLUPred_topn10_mwef,"./output/selected/GBLUPred_topn10_mwef.csv", row.names=F)
#selecting top 165 indv from GBLUP reduced for 25 and 50%
GBLUP_redmwef25_narea<- read.csv("./output/GBLUP_reduced/narea_mwef25.csv") %>% na.omit()
GBLUP_redmwef25_narea<- top_selected(data= GBLUP_redmwef25_narea)
#extracting mean and selected individual in a data frame per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_redmwef25_narea$result_top_taxa), each = length(GBLUP_redmwef25_narea$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_redmwef25_narea$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GBLUP_redmwef25_narea$result_mean_value),
                      mean_value = GBLUP_redmwef25_narea$result_mean_value)
GBLUPred_topn25_mwef <- merge(selected_df, mean_df, by = "rep")
write.csv(GBLUPred_topn25_mwef,"./output/selected/GBLUPred_topn25_mwef.csv", row.names=F)
#selecting top 165 indv from GBLUP reduced for ef and 50%
GBLUP_redmwef50_narea<- read.csv("./output/GBLUP_reduced/narea_mwef50.csv") %>% na.omit()
GBLUP_redmwef50_narea<- top_selected(data= GBLUP_redmwef50_narea)
#extracting mean and selected individual in a data frame per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_redmwef50_narea$result_top_taxa), each = length(GBLUP_redmwef50_narea$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_redmwef50_narea$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GBLUP_redmwef50_narea$result_mean_value),
                      mean_value = GBLUP_redmwef50_narea$result_mean_value)
GBLUPred_topn50_mwef <- merge(selected_df, mean_df, by = "rep")
write.csv(GBLUPred_topn50_mwef,"./output/selected/GBLUPred_topn50_mwef.csv", row.names=F)

#selecting top 165 indv from GBLUP reduced 
GBLUP_redmwef10_sla<- read.csv("./output/GBLUP_reduced/sla_mwef10.csv") %>% na.omit()
GBLUP_redmwef10_sla<- top_selected(data= GBLUP_redmwef10_sla)
#extracting mean and selected individual in a data frame per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_redmwef10_sla$result_top_taxa), each = length(GBLUP_redmwef10_sla$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_redmwef10_sla$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GBLUP_redmwef10_sla$result_mean_value),
                      mean_value = GBLUP_redmwef10_sla$result_mean_value)
GBLUPred_tops10_mwef <- merge(selected_df, mean_df, by = "rep")
write.csv(GBLUPred_tops10_mwef,"./output/selected/GBLUPred_tops10_mwef.csv", row.names=F)
#selecting top 165 indv from GBLUP reduced for 25 and 50%
GBLUP_redmwef25_sla<- read.csv("./output/GBLUP_reduced/sla_mwef25.csv") %>% na.omit()
GBLUP_redmwef25_sla<- top_selected(data= GBLUP_redmwef25_sla)
#extracting mean and selected individual in a data frame per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_redmwef25_sla$result_top_taxa), each = length(GBLUP_redmwef25_sla$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_redmwef25_sla$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GBLUP_redmwef25_sla$result_mean_value),
                      mean_value = GBLUP_redmwef25_sla$result_mean_value)
GBLUPred_tops25_mwef <- merge(selected_df, mean_df, by = "rep")
write.csv(GBLUPred_tops25_mwef,"./output/selected/GBLUPred_tops25_mwef.csv", row.names=F)
#selecting top 165 indv from GBLUP reduced for ef and 50%
GBLUP_redmwef50_sla<- read.csv("./output/GBLUP_reduced/sla_mwef50.csv") %>% na.omit()
GBLUP_redmwef50_sla<- top_selected(data= GBLUP_redmwef50_sla)
#extracting mean and selected individual in a data frame per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_redmwef50_sla$result_top_taxa), each = length(GBLUP_redmwef50_sla$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_redmwef50_sla$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GBLUP_redmwef50_sla$result_mean_value),
                      mean_value = GBLUP_redmwef50_sla$result_mean_value)
GBLUPred_tops50_mwef <- merge(selected_df, mean_df, by = "rep")
write.csv(GBLUPred_tops50_mwef,"./output/selected/GBLUPred_tops50_mwef.csv", row.names=F)


#selecting top 165 indv from GBLUP reduced 
GBLUP_redefmw10_sla<- read.csv("./output/GBLUP_reduced/sla_efmw10.csv") %>% na.omit()
GBLUP_redefmw10_sla<- top_selected(data= GBLUP_redefmw10_sla)
#extracting mean and selected individual in a data frame per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_redefmw10_sla$result_top_taxa), each = length(GBLUP_redefmw10_sla$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_redefmw10_sla$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GBLUP_redefmw10_sla$result_mean_value),
                      mean_value = GBLUP_redefmw10_sla$result_mean_value)
GBLUPred_tops10_efmw <- merge(selected_df, mean_df, by = "rep")
write.csv(GBLUPred_tops10_efmw,"./output/selected/GBLUPred_tops10_efmw.csv", row.names=F)
#selecting top 165 indv from GBLUP reduced for 25 and 50%
GBLUP_redefmw25_sla<- read.csv("./output/GBLUP_reduced/sla_efmw25.csv") %>% na.omit()
GBLUP_redefmw25_sla<- top_selected(data= GBLUP_redefmw25_sla)
#extracting mean and selected individual in a data frame per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_redefmw25_sla$result_top_taxa), each = length(GBLUP_redefmw25_sla$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_redefmw25_sla$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GBLUP_redefmw25_sla$result_mean_value),
                      mean_value = GBLUP_redefmw25_sla$result_mean_value)
GBLUPred_tops25_efmw <- merge(selected_df, mean_df, by = "rep")
write.csv(GBLUPred_tops25_efmw,"./output/selected/GBLUPred_tops25_efmw.csv", row.names=F)
#selecting top 165 indv from GBLUP reduced for ef and 50%
GBLUP_redefmw50_sla<- read.csv("./output/GBLUP_reduced/sla_efmw50.csv") %>% na.omit()
GBLUP_redefmw50_sla<- top_selected(data= GBLUP_redefmw50_sla)
#extracting mean and selected individual in a data frame per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_redefmw50_sla$result_top_taxa), each = length(GBLUP_redefmw50_sla$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_redefmw50_sla$result_top_taxa))
mean_df <- data.frame(rep = seq_along(GBLUP_redefmw50_sla$result_mean_value),
                      mean_value = GBLUP_redefmw50_sla$result_mean_value)
GBLUPred_tops50_efmw <- merge(selected_df, mean_df, by = "rep")
write.csv(GBLUPred_tops50_efmw,"./output/selected/GBLUPred_tops50_efmw.csv", row.names=F)


# reading data for narea mw for all model and blues
top_nareablues_mw <- read.csv("./output/selected/top_nareablues_mw.csv")
GBLUP_topn_mwef <- read.csv("./output/selected/GBLUP_topn_mwef.csv")
Gh2_topn10_mwef <- read.csv("./output/selected/Gh2_topn10_mwef.csv")
GWW_topn10_mwef <- read.csv("./output/selected/GWW_topn10_mwef.csv")
Gnirs_topn10_mwef <- read.csv("./output/selected/Gnirs_topn10_mwef.csv")
GBLUPred_topn10_mwef <- read.csv("./output/selected/GBLUPred_topn10_mwef.csv")
Gh2_topn25_mwef <- read.csv("./output/selected/Gh2_topn25_mwef.csv")
GWW_topn25_mwef <- read.csv("./output/selected/GWW_topn25_mwef.csv")
Gnirs_topn25_mwef <- read.csv("./output/selected/Gnirs_topn25_mwef.csv")
GBLUPred_topn25_mwef <- read.csv("./output/selected/GBLUPred_topn25_mwef.csv")
Gh2_topn50_mwef <- read.csv("./output/selected/Gh2_topn50_mwef.csv")
GWW_topn50_mwef <- read.csv("./output/selected/GWW_topn50_mwef.csv")
Gnirs_topn50_mwef <- read.csv("./output/selected/Gnirs_topn50_mwef.csv")
GBLUPred_topn50_mwef <- read.csv("./output/selected/GBLUPred_topn50_mwef.csv")

# Combining the dataset for different schemes
mean_mw_narea <- bind_rows(
  mutate(Gh2_topn10_mwef, Model = "Gh2", Scheme = "10%"),
  mutate(Gnirs_topn10_mwef, Model = "Gnirs", Scheme = "10%"),
  mutate(GWW_topn10_mwef, Model = "WW", Scheme = "10%"),
  mutate(GBLUP_topn_mwef, Model = "GBLUP", Scheme = "10%"),
  mutate(GBLUPred_topn10_mwef, Model = "GBLUP_red", Scheme = "10%"),
  mutate(Gh2_topn25_mwef, Model = "Gh2", Scheme = "25%"),
  mutate(Gnirs_topn25_mwef, Model = "Gnirs", Scheme = "25%"),
  mutate(GWW_topn25_mwef, Model = "WW", Scheme = "25%"),
  mutate(GBLUP_topn_mwef, Model = "GBLUP", Scheme = "25%"),
  mutate(GBLUPred_topn25_mwef, Model = "GBLUP_red", Scheme = "25%"),
  mutate(Gh2_topn50_mwef, Model = "Gh2", Scheme = "50%"),
  mutate(Gnirs_topn50_mwef, Model = "Gnirs", Scheme = "50%"),
  mutate(GWW_topn50_mwef, Model = "WW", Scheme = "50%"),
  mutate(GBLUP_topn_mwef, Model = "GBLUP", Scheme = "50%"),
  mutate(GBLUPred_topn50_mwef, Model = "GBLUP_red", Scheme = "50%")
)
blues_values <- c(1.56841, 1.56841,1.56841)
model_colors_cb <- brewer.pal(5, "Set1")

# Plotting the combined data with facets for different schemes and colored boxplots 
ggplot(mean_mw_narea, aes(x = Model, y = mean_value, fill = Model)) +
  geom_boxplot() +
  geom_hline(yintercept = blues_values, linetype = "dashed", color = "blue") +
  labs(x = "Model", y = "Mean Value") +
  ggtitle("Model Performance; BLUES vs model schemes for Narea in loc mw") +
  theme_minimal() +
  scale_fill_manual(values = model_colors_cb) +  
  facet_wrap(~Scheme, scales = "free")


# reading data for narea mw for all model and blues
top_nareablues_ef <- read.csv("./output/selected/top_nareablues_ef.csv")
GBLUP_topn_efmw <- read.csv("./output/selected/GBLUP_topn_efmw.csv")
Gh2_topn10_efmw <- read.csv("./output/selected/Gh2_topn10_efmw.csv")
GWW_topn10_efmw <- read.csv("./output/selected/GWW_topn10_efmw.csv")
Gnirs_topn10_efmw <- read.csv("./output/selected/Gnirs_topn10_efmw.csv")
GBLUPred_topn10_efmw <- read.csv("./output/selected/GBLUPred_topn10_efmw.csv")
Gh2_topn25_efmw <- read.csv("./output/selected/Gh2_topn25_efmw.csv")
GWW_topn25_efmw <- read.csv("./output/selected/GWW_topn25_efmw.csv")
Gnirs_topn25_efmw <- read.csv("./output/selected/Gnirs_topn25_efmw.csv")
GBLUPred_topn25_efmw <- read.csv("./output/selected/GBLUPred_topn25_efmw.csv")
Gh2_topn50_efmw <- read.csv("./output/selected/Gh2_topn50_efmw.csv")
GWW_topn50_efmw <- read.csv("./output/selected/GWW_topn50_efmw.csv")
Gnirs_topn50_efmw <- read.csv("./output/selected/Gnirs_topn50_efmw.csv")
GBLUPred_topn50_efmw <- read.csv("./output/selected/GBLUPred_topn50_efmw.csv")

# Combining the dataset for different schemes
mean_ef_narea <- bind_rows(
  mutate(Gh2_topn10_efmw, Model = "Gh2", Scheme = "10%"),
  mutate(Gnirs_topn10_efmw, Model = "Gnirs", Scheme = "10%"),
  mutate(GWW_topn10_efmw, Model = "WW", Scheme = "10%"),
  mutate(GBLUP_topn_efmw, Model = "GBLUP", Scheme = "10%"),
  mutate(GBLUPred_topn10_efmw, Model = "GBLUP_red", Scheme = "10%"),
  mutate(Gh2_topn25_efmw, Model = "Gh2", Scheme = "25%"),
  mutate(Gnirs_topn25_efmw, Model = "Gnirs", Scheme = "25%"),
  mutate(GWW_topn25_efmw, Model = "WW", Scheme = "25%"),
  mutate(GBLUP_topn_efmw, Model = "GBLUP", Scheme = "25%"),
  mutate(GBLUPred_topn25_efmw, Model = "GBLUP_red", Scheme = "25%"),
  mutate(Gh2_topn50_efmw, Model = "Gh2", Scheme = "50%"),
  mutate(Gnirs_topn50_efmw, Model = "Gnirs", Scheme = "50%"),
  mutate(GWW_topn50_efmw, Model = "WW", Scheme = "50%"),
  mutate(GBLUP_topn_efmw, Model = "GBLUP", Scheme = "50%"),
  mutate(GBLUPred_topn50_efmw, Model = "GBLUP_red", Scheme = "50%")
)
blues_values <- c(1.5611, 1.5611,1.5611)
model_colors_cb <- brewer.pal(5, "Set1")

# Plotting the combined data with facets for different schemes and colored boxplots 
ggplot(mean_ef_narea, aes(x = Model, y = mean_value, fill = Model)) +
  geom_boxplot() +
  geom_hline(yintercept = blues_values, linetype = "dashed", color = "blue") +
  labs(x = "Model", y = "Mean Value") +
  ggtitle("Model Performance; BLUES vs model schemes for Narea in loc ef") +
  theme_minimal() +
  scale_fill_manual(values = model_colors_cb) +  
  facet_wrap(~Scheme, scales = "free")


# reading the data for different schemes for sla ef
top_slablues_ef <- read.csv("./output/selected/top_slablues_ef.csv")
GBLUP_tops_efmw <- read.csv("./output/selected/GBLUP_tops_efmw.csv")
Gh2_tops10_efmw <- read.csv("./output/selected/Gh2_tops10_efmw.csv")
GWW_tops10_efmw <- read.csv("./output/selected/GWW_tops10_efmw.csv")
Gnirs_tops10_efmw <- read.csv("./output/selected/Gnirs_tops10_efmw.csv")
GBLUPred_tops10_efmw <- read.csv("./output/selected/GBLUPred_tops10_efmw.csv")
Gh2_tops25_efmw <- read.csv("./output/selected/Gh2_tops25_efmw.csv")
GWW_tops25_efmw <- read.csv("./output/selected/GWW_tops25_efmw.csv")
Gnirs_tops25_efmw <- read.csv("./output/selected/Gnirs_tops25_efmw.csv")
GBLUPred_tops25_efmw <- read.csv("./output/selected/GBLUPred_tops25_efmw.csv")
Gh2_tops50_efmw <- read.csv("./output/selected/Gh2_tops50_efmw.csv")
GWW_tops50_efmw <- read.csv("./output/selected/GWW_tops50_efmw.csv")
Gnirs_tops50_efmw <- read.csv("./output/selected/Gnirs_tops50_efmw.csv")
GBLUPred_tops50_efmw <- read.csv("./output/selected/GBLUPred_tops50_efmw.csv")

# Combining the dataset for different schemes
mean_ef_sla <- bind_rows(
  mutate(Gh2_tops10_efmw, Model = "Gh2", Scheme = "10%"),
  mutate(Gnirs_tops10_efmw, Model = "Gnirs", Scheme = "10%"),
  mutate(GWW_tops10_efmw, Model = "WW", Scheme = "10%"),
  mutate(GBLUP_tops_efmw, Model = "GBLUP", Scheme = "10%"),
  mutate(GBLUPred_tops10_efmw, Model = "GBLUP_red", Scheme = "10%"),
  mutate(Gh2_tops25_efmw, Model = "Gh2", Scheme = "25%"),
  mutate(Gnirs_tops25_efmw, Model = "Gnirs", Scheme = "25%"),
  mutate(GWW_tops25_efmw, Model = "WW", Scheme = "25%"),
  mutate(GBLUP_tops_efmw, Model = "GBLUP", Scheme = "25%"),
  mutate(GBLUPred_tops25_efmw, Model = "GBLUP_red", Scheme = "25%"),
  mutate(Gh2_tops50_efmw, Model = "Gh2", Scheme = "50%"),
  mutate(Gnirs_tops50_efmw, Model = "Gnirs", Scheme = "50%"),
  mutate(GWW_tops50_efmw, Model = "WW", Scheme = "50%"),
  mutate(GBLUP_tops_efmw, Model = "GBLUP", Scheme = "50%"),
  mutate(GBLUPred_tops50_efmw, Model = "GBLUP_red", Scheme = "50%")
)
blues_values <- c(309.008,309.008,309.008 )
model_colors_cb <- brewer.pal(5, "Set1")

# Plotting the combined data with facets for different schemes and colored boxplots 
ggplot(mean_ef_sla, aes(x = Model, y = mean_value, fill = Model)) +
  geom_boxplot() +
  geom_hline(yintercept = blues_values, linetype = "dashed", color = "blue") +
  labs(x = "Model", y = "Mean Value") +
  ggtitle("Model Performance; BLUES vs model schemes for sla in loc ef") +
  theme_minimal() +
  scale_fill_manual(values = model_colors_cb) +  
  facet_wrap(~Scheme, scales = "free")

# reading the data for different schemes for sla mw
top_slablues_mw <- read.csv("./output/selected/top_slablues_mw.csv")
GBLUP_tops_mwef <- read.csv("./output/selected/GBLUP_tops_mwef.csv")
Gh2_tops10_mwef <- read.csv("./output/selected/Gh2_tops10_mwef.csv")
GWW_tops10_mwef <- read.csv("./output/selected/GWW_tops10_mwef.csv")
Gnirs_tops10_mwef <- read.csv("./output/selected/Gnirs_tops10_mwef.csv")
GBLUPred_tops10_mwef <- read.csv("./output/selected/GBLUPred_tops10_mwef.csv")
Gh2_tops25_mwef <- read.csv("./output/selected/Gh2_tops25_mwef.csv")
GWW_tops25_mwef <- read.csv("./output/selected/GWW_tops25_mwef.csv")
Gnirs_tops25_mwef <- read.csv("./output/selected/Gnirs_tops25_mwef.csv")
GBLUPred_tops25_mwef <- read.csv("./output/selected/GBLUPred_tops25_mwef.csv")
Gh2_tops50_mwef <- read.csv("./output/selected/Gh2_tops50_mwef.csv")
GWW_tops50_mwef <- read.csv("./output/selected/GWW_tops50_mwef.csv")
Gnirs_tops50_mwef <- read.csv("./output/selected/Gnirs_tops50_mwef.csv")
GBLUPred_tops50_mwef <- read.csv("./output/selected/GBLUPred_tops50_mwef.csv")

# Combining dataset for different schemes for sla mw
mean_mw_sla <- bind_rows(
  mutate(Gh2_tops10_mwef, Model = "Gh2", Scheme = "10%"),
  mutate(Gnirs_tops10_mwef, Model = "Gnirs", Scheme = "10%"),
  mutate(GWW_tops10_mwef, Model = "WW", Scheme = "10%"),
  mutate(GBLUP_tops_mwef, Model = "GBLUP", Scheme = "10%"),
  mutate(GBLUPred_tops10_mwef, Model = "GBLUP_red", Scheme = "10%"),
  mutate(Gh2_tops25_mwef, Model = "Gh2", Scheme = "25%"),
  mutate(Gnirs_tops25_mwef, Model = "Gnirs", Scheme = "25%"),
  mutate(GWW_tops25_mwef, Model = "WW", Scheme = "25%"),
  mutate(GBLUP_tops_mwef, Model = "GBLUP", Scheme = "25%"),
  mutate(GBLUPred_tops25_mwef, Model = "GBLUP_red", Scheme = "25%"),
  mutate(Gh2_tops50_mwef, Model = "Gh2", Scheme = "50%"),
  mutate(Gnirs_tops50_mwef, Model = "Gnirs", Scheme = "50%"),
  mutate(GWW_tops50_mwef, Model = "WW", Scheme = "50%"),
  mutate(GBLUP_tops_mwef, Model = "GBLUP", Scheme = "50%"),
  mutate(GBLUPred_tops50_mwef, Model = "GBLUP_red", Scheme = "50%")
)
blues_values <- c(312.0316,312.0316,312.0316 )
model_colors_cb <- brewer.pal(5, "Set1")

# Plotting the combined data with facets for different schemes and colored boxplots 
ggplot(mean_mw_sla, aes(x = Model, y = mean_value, fill = Model)) +
  geom_boxplot() +
  geom_hline(yintercept = blues_values, linetype = "dashed", color = "blue") +
  labs(x = "Model", y = "Mean Value") +
  ggtitle("Model Performance; BLUES vs model schemes for sla in loc mw") +
  theme_minimal() +
  scale_fill_manual(values = model_colors_cb) +  
  facet_wrap(~Scheme, scales = "free")



