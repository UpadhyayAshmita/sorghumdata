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


#combined models visulaization
#visualizing the GWW model compared to genomic model and wholewave model for sla efmw

corrS_GWW_10efmw<- read.csv("./output/GWW_efmw_sla/corrS_GWW_10efmw.csv")
corrS_GWW_25efmw<- read.csv("./output/GWW_efmw_sla/corrS_GWW_25efmw.csv")
corrS_GWW_50efmw<- read.csv("./output/GWW_efmw_sla/corrS_GWW_50efmw.csv")
corr_S_wholewave_efmw<- read.csv("./output/wholewave/corr_S_wholewave_efmw.csv")%>% 
  rename(wholewave= result_sla_wholewave_efmw.ac)
corr_S_efmw<- read.csv("./output/GBLUP/corr_S_efmw.csv")%>% 
  rename(GBLUP = result_sla_efmw.ac)

sla_GWW_efmw<- bind_cols(corr_S_efmw,corr_S_wholewave_efmw,corrS_GWW_50efmw,
                              corrS_GWW_25efmw,corrS_GWW_10efmw )%>%
  pivot_longer(cols = 1:5,names_to = "models", values_to = "accuracy")

#creating GWW vs whole wave Vs GBLUP model accuracy for sla efmw
p1<- ggplot(sla_GWW_efmw, aes(x = models, y = accuracy, 
                                              fill = models)) +
  geom_boxplot() +
  labs(title = "Genomic+ wholewave model vs GBLUP and wholewave model using ef training and mw val for sla",
       x = "Models",
       y = "Accuracy")+
  scale_y_continuous("prediction accuracy",limits=c(-0.2,0.55), breaks=c(0,0.1,0.2,0.3,0.4,0.5))+ theme(aspect.ratio = 0.7)

#visualizing the Gnirs model compared to genomic model and wholewave model for sla efmw
corrS_Gnirs_10efmw<- read.csv("./output/Gnirs_efmw_sla/corrS_Gnirs_10efmw.csv")
corrS_Gnirs_25efmw<- read.csv("./output/Gnirs_efmw_sla/corrS_Gnirs_25efmw.csv")
corrS_Gnirs_50efmw<- read.csv("./output/Gnirs_efmw_sla/corrS_Gnirs_50efmw.csv")
corr_S_nirs_efmw<- read.csv("./output/NIRS/corr_S_nirs_efmw.csv")%>% 
  rename(nirs= result_sla_nirs_efmw.ac)
corr_S_efmw<- read.csv("./output/GBLUP/corr_S_efmw.csv")%>% 
  rename(GBLUP = result_sla_efmw.ac)

sla_Gnirs_efmw <- bind_cols(corr_S_efmw,corr_S_nirs_efmw,corrS_Gnirs_50efmw,
                         corrS_Gnirs_25efmw,corrS_Gnirs_10efmw )%>%
  pivot_longer(cols = 1:5,names_to = "models", values_to = "accuracy")
#creating Gnirs Vs GBLUP and nirs model accuracy for sla efmw
p2<- ggplot(sla_Gnirs_efmw, aes(x = models, y = accuracy, 
                              fill = models)) +
  geom_boxplot() +
  labs(title = "Genomic+ nirs model vs GBLUP and nirs model using ef training and mw val for sla",
       x = "Models",
       y = "Accuracy")+
  scale_y_continuous("prediction accuracy",limits=c(-0.2,0.55), breaks=c(0,0.1,0.2,0.3,0.4,0.5))+ theme(aspect.ratio = 0.7)

#visualizing the Gh2 model compared to genomic model and highh2 model for sla efmw
corrS_Gh2_10efmw<- read.csv("./output/Gh2_efmw_sla/corrS_Gh2_10efmw.csv")
corrS_Gh2_25efmw<- read.csv("./output/Gh2_efmw_sla/corrS_Gh2_25efmw.csv")
corrS_Gh2_50efmw<- read.csv("./output/Gh2_efmw_sla/corrS_Gh2_50efmw.csv")
corr_S_highh2_efmw<- read.csv("./output/highh2/corr_S_highh2_efmw.csv")%>% 
  rename(highh2= result_sla_highh2_efmw.ac)
corr_S_efmw<- read.csv("./output/GBLUP/corr_S_efmw.csv")%>% 
  rename(GBLUP = result_sla_efmw.ac)

sla_Gh2_efmw <- bind_cols(corr_S_efmw,corr_S_highh2_efmw,corrS_Gh2_50efmw,
                            corrS_Gh2_25efmw,corrS_Gh2_10efmw )%>%
  pivot_longer(cols = 1:5,names_to = "models", values_to = "accuracy")
#creating Gh2 vs whole wave Vs GBLUP model accuracy for sla efmw
p3<- ggplot(sla_Gh2_efmw, aes(x = models, y = accuracy, 
                                fill = models)) +
  geom_boxplot() +
  labs(title = "Genomic+ h2 model vs GBLUP and highly heritable wave model using ef training and mw val for sla",
       x = "Models",
       y = "Accuracy")+
  scale_y_continuous("prediction accuracy",limits=c(-0.2,0.55), breaks=c(0,0.1,0.2,0.3,0.4,0.5))+ theme(aspect.ratio = 0.7)



#visualizing the GWW model compared to genomic model and wholewave model for sla mwef
corrS_GWW_10mwef<- read.csv("./output/GWW_mwef_sla/corrS_GWW_10mwef.csv")
corrS_GWW_25mwef<- read.csv("./output/GWW_mwef_sla/corrS_GWW_25mwef.csv")
corrS_GWW_50mwef<- read.csv("./output/GWW_mwef_sla/corrS_GWW_50mwef.csv")
corr_S_wholewave_mwef<- read.csv("./output/wholewave/corr_S_wholewave_mwef.csv")%>% 
  rename(wholewave= result_sla_wholewave_mwef.ac)
corr_S_mwef<- read.csv("./output/GBLUP/corr_S_mwef.csv")%>% 
  rename(GBLUP = result_sla_mwef.ac)

sla_GWW_mwef<- bind_cols(corr_S_mwef,corr_S_wholewave_mwef,corrS_GWW_50mwef,
                         corrS_GWW_25mwef,corrS_GWW_10mwef )%>%
  pivot_longer(cols = 1:5,names_to = "models", values_to = "accuracy")
#creating GWW vs whole wave Vs GBLUP model accuracy for sla mwef
p4<- ggplot(sla_GWW_mwef, aes(x = models, y = accuracy, 
                              fill = models)) +
  geom_boxplot() +
  labs(title = "Genomic+ wholewave model vs GBLUP and wholewave model using mw training and ef val for sla",
       x = "Models",
       y = "Accuracy")+
  scale_y_continuous("prediction accuracy",limits=c(-0.2,0.55), breaks=c(0,0.1,0.2,0.3,0.4,0.5))+ theme(aspect.ratio = 0.7)

#visualizing the Gnirs model compared to genomic model and wholewave model for sla efmw
corrS_Gnirs_10mwef<- read.csv("./output/Gnirs_mwef_sla/corrS_Gnirs_10mwef.csv")
corrS_Gnirs_25mwef<- read.csv("./output/Gnirs_mwef_sla/corrS_Gnirs_25mwef.csv")
corrS_Gnirs_50mwef<- read.csv("./output/Gnirs_mwef_sla/corrS_Gnirs_50mwef.csv")
corr_S_nirs_mwef<- read.csv("./output/NIRS/corr_sla_nirs_mwef.csv")%>% 
  rename(nirs= result_sla_nirs_mwef.ac)
corr_S_mwef<- read.csv("./output/GBLUP/corr_S_mwef.csv")%>% 
  rename(GBLUP = result_sla_mwef.ac)

sla_Gnirs_mwef <- bind_cols(corr_S_mwef,corr_S_nirs_mwef,corrS_Gnirs_50mwef,
                            corrS_Gnirs_25mwef,corrS_Gnirs_10mwef )%>%
  pivot_longer(cols = 1:5,names_to = "models", values_to = "accuracy")
#creating Gnirs Vs GBLUP and nirs model accuracy for sla efmw
p5<- ggplot(sla_Gnirs_mwef, aes(x = models, y = accuracy, 
                                fill = models)) +
  geom_boxplot() +
  labs(title = "Genomic+ nirs model vs GBLUP and nirs model using mw training and ef val for sla",
       x = "Models",
       y = "Accuracy")+
  scale_y_continuous("prediction accuracy",limits=c(-0.2,0.55), breaks=c(0,0.1,0.2,0.3,0.4,0.5))+ theme(aspect.ratio = 0.7)

#visualizing the Gh2 model compared to genomic model and highh2 model for sla efmw
corrS_Gh2_10mwef<- read.csv("./output/Gh2_mwef_sla/corrS_Gh2_10mwef.csv")
corrS_Gh2_25mwef<- read.csv("./output/Gh2_mwef_sla/corrS_Gh2_25mwef.csv")
corrS_Gh2_50mwef<- read.csv("./output/Gh2_mwef_sla/corrS_Gh2_50mwef.csv")
corr_S_highh2_mwef<- read.csv("./output/highh2/corr_S_highh2_mwef.csv")%>% 
  rename(highh2= result_sla_highh2_mwef.ac)
corr_S_mwef<- read.csv("./output/GBLUP/corr_S_mwef.csv")%>% 
  rename(GBLUP = result_sla_mwef.ac)

sla_Gh2_mwef <- bind_cols(corr_S_mwef,corr_S_highh2_mwef,corrS_Gh2_50mwef,
                          corrS_Gh2_25mwef,corrS_Gh2_10mwef )%>%
  pivot_longer(cols = 1:5,names_to = "models", values_to = "accuracy")
#creating Gh2 vs whole wave Vs GBLUP model accuracy for sla efmw
p6 <- ggplot(sla_Gh2_mwef, aes(x = models, y = accuracy, 
                              fill = models)) +
  geom_boxplot() +
  labs(title = "Genomic+ h2 model vs GBLUP and highly heritable wave model using mw training and ef val for sla",
       x = "Models",
       y = "Accuracy")+
  scale_y_continuous("prediction accuracy",limits=c(-0.2,0.55), breaks=c(0,0.1,0.2,0.3,0.4,0.5))+ theme(aspect.ratio = 0.7)

#visualizing the GWW model compared to genomic model and wholewave model for narea mwef
corrN_GWW_10mwef<- read.csv("./output/GWW_mwef_narea/corrN_GWW_10mwef.csv")
corrN_GWW_25mwef<- read.csv("./output/GWW_mwef_narea/corrN_GWW_25mwef.csv")
corrN_GWW_50mwef<- read.csv("./output/GWW_mwef_narea/corrN_GWW_50mwef.csv")
corr_N_wholewave_mwef<- read.csv("./output/wholewave/corr_N_wholewave_mwef.csv")%>% 
  rename(wholewave= result_narea_wholewave_mwef.ac)

corr_N_mwef <- read.csv("./output/GBLUP/ corr_N_mwef.csv")%>% 
  rename(GBLUP = result_narea_mwef.ac)

narea_GWW_mwef<- bind_cols(corr_N_mwef,corr_N_wholewave_mwef,corrN_GWW_50mwef,
                         corrN_GWW_25mwef,corrN_GWW_10mwef )%>%
  pivot_longer(cols = 1:5,names_to = "models", values_to = "accuracy")
#creating GWW vs whole wave Vs GBLUP model accuracy for narea mwef
p7<- ggplot(narea_GWW_mwef, aes(x = models, y = accuracy, 
                              fill = models)) +
  geom_boxplot() +
  labs(title = "Genomic+ wholewave model vs GBLUP and wholewave model using mw training and ef val for narea",
       x = "Models",
       y = "Accuracy")+
  scale_y_continuous("prediction accuracy",limits=c(-0.2,0.55), breaks=c(0,0.1,0.2,0.3,0.4,0.5))+ theme(aspect.ratio = 0.7)

#visualizing the Gnirs model compared to genomic model and wholewave model for narea mwef
corrN_Gnirs_10mwef<- read.csv("./output/Gnirs_mwef_narea/corrN_Gnirs_10mwef.csv")
corrN_Gnirs_25mwef<- read.csv("./output/Gnirs_mwef_narea/corrN_Gnirs_25mwef.csv")
corrN_Gnirs_50mwef<- read.csv("./output/Gnirs_mwef_narea/corrN_Gnirs_50mwef.csv")
corr_N_nirs_mwef<- read.csv("./output/NIRS/corr_N_nirs_mwef.csv")%>% 
  rename(nirs= result_narea_nirs_mwef.ac)
corr_N_mwef<- read.csv("./output/GBLUP/ corr_N_mwef.csv")%>% 
  rename(GBLUP = result_narea_mwef.ac)

narea_Gnirs_mwef <- bind_cols(corr_N_mwef,corr_N_nirs_mwef,corrN_Gnirs_50mwef,
                            corrN_Gnirs_25mwef,corrN_Gnirs_10mwef )%>%
  pivot_longer(cols = 1:5,names_to = "models", values_to = "accuracy")
#creating Gnirs Vs GBLUP and nirs model accuracy for narea mwef
p8<- ggplot(narea_Gnirs_mwef, aes(x = models, y = accuracy, 
                                fill = models)) +
  geom_boxplot() +
  labs(title = "Genomic+ nirs model vs GBLUP and nirs model using mw training and ef val for narea",
       x = "Models",
       y = "Accuracy")+
  scale_y_continuous("prediction accuracy",limits=c(-0.2,0.55), breaks=c(0,0.1,0.2,0.3,0.4,0.5))+ theme(aspect.ratio = 0.7)

#visualizing the Gh2 model compared to genomic model and highh2 model for narea efmw
corrN_Gh2_10mwef<- read.csv("./output/Gh2_mwef_narea/corrN_Gh2_10mwef.csv")
corrN_Gh2_25mwef<- read.csv("./output/Gh2_mwef_narea/corrN_Gh2_25mwef.csv")
corrN_Gh2_50mwef<- read.csv("./output/Gh2_mwef_narea/corrN_Gh2_50mwef.csv")
corr_N_highh2_mwef<- read.csv("./output/highh2/corr_N_highh2_mwef.csv")%>% 
  rename(highh2= result_narea_highh2_mwef.ac)
corr_N_mwef<- read.csv("./output/GBLUP/ corr_N_mwef.csv")%>% 
  rename(GBLUP = result_narea_mwef.ac)

narea_Gh2_mwef <- bind_cols(corr_N_mwef,corr_N_highh2_mwef,corrN_Gh2_50mwef,
                          corrN_Gh2_25mwef,corrN_Gh2_10mwef )%>%
  pivot_longer(cols = 1:5,names_to = "models", values_to = "accuracy")
#creating Gh2 vs whole wave Vs GBLUP model accuracy for narea mwef
p9 <- ggplot(narea_Gh2_mwef, aes(x = models, y = accuracy, 
                               fill = models)) +
  geom_boxplot() +
  labs(title = "Genomic+ h2 model vs GBLUP and highly heritable wave model using mw training and ef val for narea",
       x = "Models",
       y = "Accuracy")+
  scale_y_continuous("prediction accuracy",limits=c(-0.2,0.55), breaks=c(0,0.1,0.2,0.3,0.4,0.5))+ theme(aspect.ratio = 0.7)

#visualizing the GWW model compared to genomic model and wholewave model for narea efmw
corrN_GWW_10efmw<- read.csv("./output/GWW_efmw_narea/corrN_GWW_10efmw.csv")
corrN_GWW_25efmw<- read.csv("./output/GWW_efmw_narea/corrN_GWW_25efmw.csv")
corrN_GWW_50efmw<- read.csv("./output/GWW_efmw_narea/corrN_GWW_50efmw.csv")
corr_N_wholewave_efmw<- read.csv("./output/wholewave/corr_N_wholewave_efmw.csv")%>% 
  rename(wholewave= result_narea_wholewave_efmw.ac)

corr_N_efmw <- read.csv("./output/GBLUP/corr_N_efmw.csv")%>% 
  rename(GBLUP = result_narea_efmw.ac)

narea_GWW_efmw<- bind_cols(corr_N_efmw,corr_N_wholewave_efmw,corrN_GWW_50efmw,
                           corrN_GWW_25efmw,corrN_GWW_10efmw )%>%
  pivot_longer(cols = 1:5,names_to = "models", values_to = "accuracy")
#creating GWW vs whole wave Vs GBLUP model accuracy for narea efmw
p10<- ggplot(narea_GWW_efmw, aes(x = models, y = accuracy, 
                                fill = models)) +
  geom_boxplot() +
  labs(title = "Genomic+ wholewave model vs GBLUP and wholewave model using ef training and mw val for narea",
       x = "Models",
       y = "Accuracy")+
  scale_y_continuous("prediction accuracy",limits=c(-0.2,0.55), breaks=c(0,0.1,0.2,0.3,0.4,0.5))+ theme(aspect.ratio = 0.7)

#visualizing the Gnirs model compared to genomic model and wholewave model for narea efmw
corrN_Gnirs_10efmw<- read.csv("./output/Gnirs_efmw_narea/corrN_Gnirs_10efmw.csv")
corrN_Gnirs_25efmw<- read.csv("./output/Gnirs_efmw_narea/corrN_Gnirs_25efmw.csv")
corrN_Gnirs_50efmw<- read.csv("./output/Gnirs_efmw_narea/corrN_Gnirs_50efmw.csv")
corr_N_nirs_efmw<- read.csv("./output/NIRS/corr_N_nirs_efmw.csv")%>% 
  rename(nirs= result_narea_nirs_efmw.ac)
corr_N_efmw<- read.csv("./output/GBLUP/corr_N_efmw.csv")%>% 
  rename(GBLUP = result_narea_efmw.ac)
narea_Gnirs_efmw <- bind_cols(corr_N_efmw,corr_N_nirs_efmw,corrN_Gnirs_50efmw,
                              corrN_Gnirs_25efmw,corrN_Gnirs_10efmw )%>%
  pivot_longer(cols = 1:5,names_to = "models", values_to = "accuracy")
#creating Gnirs Vs GBLUP and nirs model accuracy for narea efmw
p11<- ggplot(narea_Gnirs_efmw, aes(x = models, y = accuracy, 
                                  fill = models)) +
  geom_boxplot() +
  labs(title = "Genomic+ nirs model vs GBLUP and nirs model using ef training and mw val for narea",
       x = "Models",
       y = "Accuracy")+
  scale_y_continuous("prediction accuracy",limits=c(-0.2,0.55), breaks=c(0,0.1,0.2,0.3,0.4,0.5))+ theme(aspect.ratio = 0.7)

#visualizing the Gh2 model compared to genomic model and highh2 model for narea efmw
corrN_Gh2_10efmw<- read.csv("./output/Gh2_efmw_narea/corrN_Gh2_10efmw.csv")
corrN_Gh2_25efmw<- read.csv("./output/Gh2_efmw_narea/corrN_Gh2_25efmw.csv")
corrN_Gh2_50efmw<- read.csv("./output/Gh2_efmw_narea/corrN_Gh2_50efmw.csv")
corr_N_highh2_efmw<- read.csv("./output/highh2/corr_N_highh2_efmw.csv")%>% 
  rename(highh2= result_narea_highh2_efmw.ac)
corr_N_efmw<- read.csv("./output/GBLUP/corr_N_efmw.csv")%>% 
  rename(GBLUP = result_narea_efmw.ac)
narea_Gh2_efmw <- bind_cols(corr_N_efmw,corr_N_highh2_efmw,corrN_Gh2_50efmw,
                            corrN_Gh2_25efmw,corrN_Gh2_10efmw )%>%
  pivot_longer(cols = 1:5,names_to = "models", values_to = "accuracy")
#creating Gh2 vs whole wave Vs GBLUP model accuracy for narea efmw
p12 <- ggplot(narea_Gh2_efmw, aes(x = models, y = accuracy, 
                                 fill = models)) +
  geom_boxplot() +
  labs(title = "Genomic+ h2 model vs GBLUP and highly heritable wave model using ef training and mw val for narea",
       x = "Models",
       y = "Accuracy")+
  scale_y_continuous("prediction accuracy",limits=c(-0.2,0.55), breaks=c(0,0.1,0.2,0.3,0.4,0.5))+ theme(aspect.ratio = 0.7)

#visualizing the GWW model compared to genomic model and wholewave model for narea joint
corrN_GWW_10j<- read.csv("./output/GWW_joint_narea/corrN_GWW_10j.csv")
corrN_GWW_25j<- read.csv("./output/GWW_joint_narea/corrN_GWW_25j.csv")
corrN_GWW_50j<- read.csv("./output/GWW_joint_narea/corrN_GWW_50j.csv")
corr_N_wholewave_joint<- read.csv("./output/wholewave/corr_N_wholewave_joint.csv")%>% 
  rename(wholewave= result_N_wholewave_joint.ac)
corr_N_joint <- read.csv("./output/GBLUP/corr_N_joint.csv")%>% 
  rename(GBLUP = result_N.ac)
narea_GWW_joint<- bind_cols(corr_N_joint,corr_N_wholewave_joint,corrN_GWW_50j,
                           corrN_GWW_25j,corrN_GWW_10j )%>%
  pivot_longer(cols = 1:5,names_to = "models", values_to = "accuracy")
#creating GWW vs whole wave Vs GBLUP model accuracy for narea joint
p13<- ggplot(narea_GWW_joint, aes(x = models, y = accuracy, 
                                 fill = models)) +
  geom_boxplot() +
  labs(title = "Genomic+ wholewave model vs GBLUP and wholewave model for narea joint loc",
       x = "Models",
       y = "Accuracy")+
  scale_y_continuous("prediction accuracy",limits=c(-0.2,0.55), breaks=c(0,0.1,0.2,0.3,0.4,0.5))+ theme(aspect.ratio = 0.7)

#visualizing the Gnirs model compared to genomic model and wholewave model for narea joint
corrN_Gnirs_10j<- read.csv("./output/Gnirs_joint_narea/corrN_Gnirs_10j.csv")
corrN_Gnirs_25j<- read.csv("./output/Gnirs_joint_narea/corrN_Gnirs_25j.csv")
corrN_Gnirs_50j<- read.csv("./output/Gnirs_joint_narea/corrN_Gnirs_50j.csv")
corr_N_nirs_joint<- read.csv("./output/NIRS/corr_N_nirs_joint.csv")%>% 
  rename(nirs= result_N_nirs_joint.ac)
corr_N_joint<- read.csv("./output/GBLUP/corr_N_joint.csv")%>% 
  rename(GBLUP = result_N.ac)
narea_Gnirs_joint <- bind_cols(corr_N_joint,corr_N_nirs_joint,corrN_Gnirs_50j,
                              corrN_Gnirs_25j,corrN_Gnirs_10j )%>%
  pivot_longer(cols = 1:5,names_to = "models", values_to = "accuracy")
#creating Gnirs Vs GBLUP and nirs model accuracy for narea joint
p14<- ggplot(narea_Gnirs_joint, aes(x = models, y = accuracy, 
                                   fill = models)) +
  geom_boxplot() +
  labs(title = "Genomic+ nirs model vs GBLUP and nirs model for narea joint",
       x = "Models",
       y = "Accuracy")+
  scale_y_continuous("prediction accuracy",
                     limits=c(-0.2,0.55), 
                     breaks=c(0,0.1,0.2,0.3,0.4,0.5))+ theme(aspect.ratio = 0.7)

#visualizing the Gh2 model compared to genomic model and highh2 model for narea joint
corrN_Gh2_10j<- read.csv("./output/Gh2_joint_narea/corrN_Gh2_10j.csv")
corrN_Gh2_25j<- read.csv("./output/Gh2_joint_narea/corrN_Gh2_25j.csv")
corrN_Gh2_50j<- read.csv("./output/Gh2_joint_narea/corrN_Gh2_50j.csv")
corr_N_highh2_joint<- read.csv("./output/highh2/corr_N_highh2_joint.csv")%>% 
  rename(highh2= result_N_highh2_joint.ac)
corr_N_joint<- read.csv("./output/GBLUP/corr_N_joint.csv")%>% 
  rename(GBLUP = result_N.ac)
narea_Gh2_joint <- bind_cols(corr_N_joint,corr_N_highh2_joint,corrN_Gh2_50j,
                            corrN_Gh2_25j,corrN_Gh2_10j)%>%
  pivot_longer(cols = 1:5,names_to = "models", values_to = "accuracy")
#creating Gh2 vs whole wave Vs GBLUP model accuracy for narea joint
p15 <- ggplot(narea_Gh2_joint, aes(x = models, y = accuracy, 
                                  fill = models)) +
  geom_boxplot() +
  labs(title = "Genomic+ h2 model vs GBLUP and highly heritable wave model for narea joint loc",
       x = "Models",
       y = "Accuracy")+
  scale_y_continuous("prediction accuracy",limits=c(-0.2,0.55), breaks=c(0,0.1,0.2,0.3,0.4,0.5))+ theme(aspect.ratio = 0.7)

#visualizing the GWW model compared to genomic model and wholewave model for sla joint
corrS_GWW_10j<- read.csv("./output/GWW_joint_sla/corrS_GWW_10j.csv")
corrS_GWW_25j<- read.csv("./output/GWW_joint_sla/corrS_GWW_25j.csv")
corrS_GWW_50j<- read.csv("./output/GWW_joint_sla/corrS_GWW_50j.csv")
corr_S_wholewave_joint<- read.csv("./output/wholewave/corr_S_wholewave_joint.csv")%>% 
  rename(wholewave= result_S_wholewave_joint.ac)
corr_S_joint <- read.csv("./output/GBLUP/corr_S_joint.csv")%>% 
  rename(GBLUP = result_S.ac)
sla_GWW_joint<- bind_cols(corr_S_joint,corr_S_wholewave_joint,corrS_GWW_50j,
                            corrS_GWW_25j,corrS_GWW_10j )%>%
  pivot_longer(cols = 1:5,names_to = "models", values_to = "accuracy")
#creating GWW vs whole wave Vs GBLUP model accuracy for sla joint
p16<- ggplot(sla_GWW_joint, aes(x = models, y = accuracy, 
                                  fill = models)) +
  geom_boxplot() +
  labs(title = "Genomic+ wholewave model vs GBLUP and wholewave model for sla joint loc",
       x = "Models",
       y = "Accuracy")+
  scale_y_continuous("prediction accuracy",limits=c(-0.2,0.55), breaks=c(0,0.1,0.2,0.3,0.4,0.5))+ theme(aspect.ratio = 0.7)

#visualizing the Gnirs model compared to genomic model and wholewave model for sla joint
corrS_Gnirs_10j<- read.csv("./output/Gnirs_joint_sla/corrS_Gnirs_10j.csv")
corrS_Gnirs_25j<- read.csv("./output/Gnirs_joint_sla/corrS_Gnirs_25j.csv")
corrS_Gnirs_50j<- read.csv("./output/Gnirs_joint_sla/corrS_Gnirs_50j.csv")
corr_S_nirs_joint<- read.csv("./output/NIRS/corr_S_nirs_joint.csv")%>% 
  rename(nirs= result_S_nirs_joint.ac)
corr_S_joint<- read.csv("./output/GBLUP/corr_S_joint.csv")%>% 
  rename(GBLUP = result_S.ac)
sla_Gnirs_joint <- bind_cols(corr_S_joint,corr_S_nirs_joint,corrS_Gnirs_50j,
                               corrS_Gnirs_25j,corrS_Gnirs_10j )%>%
  pivot_longer(cols = 1:5,names_to = "models", values_to = "accuracy")
#creating Gnirs Vs GBLUP and nirs model accuracy for sla joint
p17<- ggplot(sla_Gnirs_joint, aes(x = models, y = accuracy, 
                                    fill = models)) +
  geom_boxplot() +
  labs(title = "Genomic+ nirs model vs GBLUP and nirs model for sla joint",
       x = "Models",
       y = "Accuracy")+
  scale_y_continuous("prediction accuracy",limits=c(-0.2,0.55), breaks=c(0,0.1,0.2,0.3,0.4,0.5))+ theme(aspect.ratio = 0.7)

#visualizing the Gh2 model compared to genomic model and highh2 model for narea joint
corrS_Gh2_10j<- read.csv("./output/Gh2_joint_sla/corrS_Gh2_10j.csv")
corrS_Gh2_25j<- read.csv("./output/Gh2_joint_sla/corrS_Gh2_25j.csv")
corrS_Gh2_50j<- read.csv("./output/Gh2_joint_sla/corrS_Gh2_50j.csv")
corr_S_highh2_joint<- read.csv("./output/highh2/corr_S_highh2_joint.csv")%>% 
  rename(highh2= result_S_highh2_joint.ac)
corr_S_joint<- read.csv("./output/GBLUP/corr_S_joint.csv")%>% 
  rename(GBLUP = result_S.ac)
sla_Gh2_joint <- bind_cols(corr_S_joint,corr_S_highh2_joint,corrS_Gh2_50j,
                             corrS_Gh2_25j,corrS_Gh2_10j)%>%
  pivot_longer(cols = 1:5,names_to = "models", values_to = "accuracy")
#creating Gh2 vs whole wave Vs GBLUP model accuracy for sla joint
p18 <- ggplot(sla_Gh2_joint, aes(x = models, y = accuracy, 
                                   fill = models)) +
  geom_boxplot() +
  labs(title = "Genomic+ h2 model vs GBLUP and highly heritable wave model for sla joint loc",
       x = "Models",
       y = "Accuracy")+
  scale_y_continuous("prediction accuracy",limits=c(-0.2,0.55), breaks=c(0,0.1,0.2,0.3,0.4,0.5))+ theme(aspect.ratio = 0.7)





