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


#visualizing different accuracy of each model of joint location for sla
corr_S_wholewave_joint<-read.csv("./output/wholewave/corr_S_wholewave_joint.csv") %>%
  rename(wholewave= result_S_wholewave_joint.ac)
corr_S_highh2_joint<- read.csv("./output/highh2/corr_S_highh2_joint.csv")%>%
  rename(high_heritablewave = result_S_highh2_joint.ac)
corr_S_nirs_joint<- read.csv("./output/NIRS/corr_S_nirs_joint.csv")%>% 
  rename(nirswave = result_S_nirs_joint.ac)
corr_S_joint<-read.csv("./output/GBLUP/corr_S_joint.csv")%>% 
  rename( GBLUP= result_S.ac)

sla_accuracy_joint<- bind_cols(corr_S_joint,corr_S_nirs_joint,
                               corr_S_highh2_joint,corr_S_wholewave_joint)%>% 
  pivot_longer(cols = 1:4,names_to = "models", values_to = "accuracy")
# Create a boxplot
joint_sla_plot<- ggplot(sla_accuracy_joint, aes(x = models, y = accuracy, 
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
