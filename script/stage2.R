# ---------------------loading packages--------------------
library(data.table)
library(cvTools)
library(ASRgenomics)
library(asreml)
library(dplyr)
library(ggplot2)
library(tidyverse)
source("./script/aux_functions.R")
#reading blues from first stage analysis for narea 
selected_lines<- read.csv("./data/selected_lines.csv")
bluesjoint_narea_correct<- fread("./output/traitsoutput/bluesjoint_narea_correct.csv",
                                 data.table = F) %>% 
  filter(bluesjoint_narea_correct$taxa %in% selected_lines$x) %>% 
  mutate(taxa = factor(taxa))%>% rename("narea"= predicted.value)
bluesjoint_sla_correct<- fread("./output/traitsoutput/bluesjoint_sla_correct.csv",
                               data.table = F)%>%
  filter(bluesjoint_sla_correct$taxa %in% selected_lines$x) %>% 
  mutate(taxa = factor(taxa))%>% rename("sla"= predicted.value)
slablues_mw_correct<- fread("./output/traitsoutput/slablues_mw_correct.csv", data.table = F)%>% 
  filter(slablues_mw_correct$taxa %in% selected_lines$x) %>% 
  mutate(taxa = factor(taxa))%>% rename("sla"= predicted.value)
nareablues_mw_correct<- fread("./output/traitsoutput/nareablues_mw_correct.csv", data.table = F)%>% 
  filter(nareablues_mw_correct$taxa %in% selected_lines$x) %>% 
  mutate(taxa = factor(taxa))%>% rename("narea"= predicted.value)
slablues_ef_correct<- fread("./output/traitsoutput/slablues_ef_correct.csv", data.table = F)%>%
  filter(slablues_ef_correct$taxa %in% selected_lines$x) %>% 
  mutate(taxa = factor(taxa))%>% rename("sla"= predicted.value)
nareablues_ef_correct<- fread("./output/traitsoutput/nareablues_ef_correct.csv", data.table = F)%>% 
  filter(nareablues_ef_correct$taxa %in% selected_lines$x) %>% 
  mutate(taxa = factor(taxa)) %>% rename("narea"= predicted.value)

# ---------------------fitting second stage model for narea trait and training
# the model with joint loc blues i.e nitroblues and validating on same data ---------------------
# ---------------------creating list with 5 fold and 20 reps----------------------
GINV <- readRDS(file = "./data/relmatrices/GINV.rds")
sort<- create_folds( individuals= bluesjoint_narea_correct$taxa, nfolds= 5, 
                     reps = 20, seed = 123)
# ---------------------narea-running cross validation for joint loc---------------------
result_N <-
  crossv(
    sort = sort,
    train = na.omit(bluesjoint_narea_correct),
    validation =  bluesjoint_narea_correct,
    mytrait = "narea",
    GINV = GINV,
    data= "GBLUP",
    scheme= "joint"
  )
corr_N_joint<- data.frame(result_N$ac)
fwrite(corr_N_joint, "./output/GBLUP/corr_N_joint.csv")
# ---------------------SLA-running cross validation for joint loc---------------------
result_S <- crossv(sort = sort,
                  train = na.omit(bluesjoint_sla_correct),
                  validation = bluesjoint_sla_correct, 
                  mytrait = "sla", 
                  GINV = GINV,
                  data = "GBLUP",
                  scheme= "joint")
corr_S_joint<- data.frame(result_S$ac)
fwrite(corr_S_joint, "./output/GBLUP/corr_S_joint.csv")
# ---------------------SLA-running cross validation for ...---------------------
#SLA- CV with ef as a train and mw validation
result_sla_efmw<- crossv(sort = sort,
                  train = na.omit(slablues_ef_correct),
                  validation = slablues_mw_correct, 
                  mytrait = "sla", 
                  GINV = GINV,
                  data = "GBLUP",
                  scheme = "efmw")
corr_S_efmw<- data.frame(result_sla_efmw$ac)
fwrite(corr_S_efmw, "./output/GBLUP/corr_S_efmw.csv")
#narea training the model with ef and validating on mw--------------------
result_narea_efmw<- crossv(sort = sort,
                         train = na.omit(nareablues_ef_correct),
                         validation = nareablues_mw_correct, 
                         mytrait = "narea", 
                         GINV = GINV,
                         data = "GBLUP",
                         scheme = "efmw"
                         )
corr_N_efmw<- data.frame(result_narea_efmw$ac)
fwrite(corr_N_efmw, "./output/GBLUP/corr_N_efmw.csv")

#narea training the model with mw and validating on ef--------------------
result_narea_mwef<- crossv(sort = sort,
                           train = na.omit(nareablues_mw_correct),
                           validation = nareablues_ef_correct, 
                           mytrait = "narea", 
                           GINV = GINV,
                           data = "GBLUP",
                           scheme = "mwef")
corr_N_mwef<- data.frame(result_narea_mwef$ac)
fwrite(corr_N_mwef, "./output/GBLUP/ corr_N_mwef.csv")

#SLA training the model with mw and validating on ef--------------------
result_sla_mwef<- crossv(sort = sort,
                           train = na.omit(slablues_mw_correct),
                           validation = slablues_ef_correct, 
                           mytrait = "sla", 
                           GINV= GINV,
                         data = "GBLUP",
                         scheme = "mwef")

corr_S_mwef<- data.frame(result_sla_mwef$ac)
fwrite(corr_S_mwef, "./output/GBLUP/corr_S_mwef.csv")

#NIRS model
nirs_inv<- readRDS(file = "./data/relmatrices/nirs_inv.rds")
# fitting model for narea using train nareajoint and val on nareajoint using NIRS matrix--------------------
result_N_nirs_joint <-
  crossv(
    sort = sort,
    train = na.omit(bluesjoint_narea_correct),
    validation =  bluesjoint_narea_correct,
    mytrait = "narea",
    GINV = nirs_inv,
    data= "NIRS",
    scheme= "nirs_joint"
  )
corr_N_nirs_joint<- data.frame(result_N_nirs_joint$ac)
fwrite(corr_N_nirs_joint, "./output/NIRS/corr_N_nirs_joint.csv")
#joint sla nirs
result_S_nirs_joint <-
  crossv(
    sort = sort,
    train = na.omit(bluesjoint_sla_correct),
    validation =  bluesjoint_sla_correct,
    mytrait = "sla",
    GINV = nirs_inv,
    data= "NIRS",
    scheme= "nirs_joint"
  )
corr_S_nirs_joint<- data.frame(result_S_nirs_joint$ac)
fwrite(corr_S_nirs_joint, "./output/NIRS/corr_S_nirs_joint.csv")

#narea training ef validating mw using NIRS ef
nirs_inv_ef<- readRDS(file = "./data/relmatrices/nirs_inv_ef.rds")
result_narea_nirs_efmw<- crossv(sort = sort,
                           train = na.omit(nareablues_ef_correct),
                           validation = nareablues_mw_correct, 
                           mytrait = "narea", 
                           GINV = nirs_inv_ef,
                           data = "NIRS",
                           scheme = "nirs_efmw"
)
corr_N_nirs_efmw<- data.frame(result_narea_nirs_efmw$ac)
fwrite(corr_N_nirs_efmw, "./output/NIRS/corr_N_nirs_efmw.csv")
#narea training mw validating ef using NIRSmw
nirs_inv_mw <- readRDS(file = "./data/relmatrices/nirs_inv_mw.rds")
result_narea_nirs_mwef<- crossv(sort = sort,
                                train = na.omit(nareablues_mw_correct),
                                validation = nareablues_ef_correct, 
                                mytrait = "narea", 
                                GINV = nirs_inv_mw,
                                data = "NIRS",
                                scheme = "nirs_mwef"
)
corr_N_nirs_mwef<- data.frame(result_narea_nirs_mwef$ac)
fwrite(corr_N_nirs_mwef, "./output/NIRS/corr_N_nirs_mwef.csv")
#sla training ef validating mw using nirs ef matrix
result_sla_nirs_efmw<- crossv(sort = sort,
                                train = na.omit(slablues_ef_correct),
                                validation = slablues_mw_correct, 
                                mytrait = "sla", 
                                GINV = nirs_inv_ef,
                                data = "NIRS",
                                scheme = "nirs_efmw")
corr_S_nirs_efmw<- data.frame(result_sla_nirs_efmw$ac)
fwrite(corr_S_nirs_efmw, "./output/NIRS/corr_S_nirs_efmw.csv")
#sla training mw validating ef using nirs mw matrix
result_sla_nirs_mwef<- crossv(sort = sort,
                              train = na.omit(slablues_mw_correct),
                              validation = slablues_ef_correct, 
                              mytrait = "sla", 
                              GINV = nirs_inv_mw,
                              data = "NIRS",
                              scheme = "nirs_mwef")
corr_sla_nirs_mwef<- data.frame(result_sla_nirs_mwef$ac)
fwrite(corr_sla_nirs_mwef, "./output/NIRS/corr_sla_nirs_mwef.csv")

#whole wave model
re_w_joint_inv <- readRDS(file = "./data/relmatrices/re_w_joint_inv.rds")
# -----narea-running cross validation for joint loc using whole wave information---------------------
result_N_wholewave_joint<-
  crossv(
    sort = sort,
    train = na.omit(bluesjoint_narea_correct),
    validation =  bluesjoint_narea_correct,
    mytrait = "narea",
    GINV = re_w_joint_inv,
    data= "wholewave",
    scheme= "wholewave_joint"
  )
corr_N_wholewave_joint<- data.frame(result_N_wholewave_joint$ac)
fwrite(corr_N_wholewave_joint, "./output/wholewave/corr_N_wholewave_joint.csv")

# -----SLA-running cross validation for joint loc using whole wave info---------------------
result_S_wholewave_joint<- crossv(sort = sort,
                  train = na.omit(bluesjoint_sla_correct),
                  validation = bluesjoint_sla_correct, 
                  mytrait = "sla", 
                  GINV = re_w_joint_inv,
                  data = "wholewave",
                  scheme= "wholewave_joint")
corr_S_wholewave_joint<- data.frame(result_S_wholewave_joint$ac)
fwrite(corr_S_wholewave_joint, "./output/wholewave/corr_S_wholewave_joint.csv")

#------narea training ef validating mw using whole wave ef
re_w_ef_inv <- readRDS(file = "./data/relmatrices/re_w_ef_inv.rds")
result_narea_wholewave_efmw<- crossv(sort = sort,
                                train = na.omit(nareablues_ef_correct),
                                validation = nareablues_mw_correct, 
                                mytrait = "narea", 
                                GINV = re_w_ef_inv,
                                data = "wholewave",
                                scheme = "wholewave_efmw"
)
corr_N_wholewave_efmw<- data.frame(result_narea_wholewave_efmw$ac)
fwrite(corr_N_wholewave_efmw, "./output/wholewave/corr_N_wholewave_efmw.csv")
#narea training mw validating ef using wholewave mw
re_w_mw_inv <- readRDS(file = "./data/relmatrices/re_w_mw_inv.rds")
result_narea_wholewave_mwef<- crossv(sort = sort,
                                     train = na.omit(nareablues_mw_correct),
                                     validation = nareablues_ef_correct, 
                                     mytrait = "narea", 
                                     GINV = re_w_mw_inv,
                                     data = "wholewave",
                                     scheme = "wholewave_mwef"
)
corr_N_wholewave_mwef<- data.frame(result_narea_wholewave_mwef$ac)
fwrite(corr_N_wholewave_mwef, "./output/wholewave/corr_N_wholewave_mwef.csv")
#sla training ef and validating mw using wholewave ef
re_w_ef_inv <- readRDS(file = "./data/relmatrices/re_w_ef_inv.rds")
result_sla_wholewave_efmw<- crossv(sort = sort,
                                     train = na.omit(slablues_ef_correct),
                                     validation = slablues_mw_correct, 
                                     mytrait = "sla", 
                                     GINV = re_w_ef_inv,
                                     data = "wholewave",
                                     scheme = "wholewave_efmw"
)
corr_S_wholewave_efmw<- data.frame(result_sla_wholewave_efmw$ac)
fwrite(corr_S_wholewave_efmw, "./output/wholewave/corr_S_wholewave_efmw.csv")
#sla training mw validating ef using wholewave mw
re_w_mw_inv <- readRDS(file = "./data/relmatrices/re_w_mw_inv.rds")
result_sla_wholewave_mwef<- crossv(sort = sort,
                                     train = na.omit(slablues_mw_correct),
                                     validation = slablues_ef_correct, 
                                     mytrait = "sla", 
                                     GINV = re_w_mw_inv,
                                     data = "wholewave",
                                     scheme = "wholewave_mwef"
)
corr_S_wholewave_mwef<- data.frame(result_sla_wholewave_mwef$ac)
fwrite(corr_S_wholewave_mwef, "./output/wholewave/corr_S_wholewave_mwef.csv")

#--highly heritable models using high heritable relationship matrix for jointloc
re_h2_joint_inv <- readRDS(file = "./data/relmatrices/re_h2_joint_inv.rds")
# -----narea-running cross validation for joint loc using highh2 matrix--------------------
result_N_highh2_joint<-
  crossv(
    sort = sort,
    train = na.omit(bluesjoint_narea_correct),
    validation =  bluesjoint_narea_correct,
    mytrait = "narea",
    GINV = re_h2_joint_inv,
    data= "highh2",
    scheme= "highh2_joint"
  )
corr_N_highh2_joint<- data.frame(result_N_highh2_joint$ac)
fwrite(corr_N_highh2_joint, "./output/highh2/corr_N_highh2_joint.csv")
# -----SLA-running cv for joint loc using highh2 matrix--------------------
result_S_highh2_joint<- crossv(sort = sort,
                            train = na.omit(bluesjoint_sla_correct),
                            validation = bluesjoint_sla_correct, 
                            mytrait = "sla", 
                            GINV = re_h2_joint_inv,
                            data = "highh2",
                            scheme= "highh2_joint")
corr_S_highh2_joint<- data.frame(result_S_highh2_joint$ac)
fwrite(corr_S_highh2_joint, "./output/highh2/corr_S_highh2_joint.csv")
#------narea training ef validating mw using high h2 ef loc matrix
re_h2_ef_inv <- readRDS(file = "./data/relmatrices/re_h2_ef_inv.rds")
result_narea_highh2_efmw<- crossv(sort = sort,
                                     train = na.omit(nareablues_ef_correct),
                                     validation = nareablues_mw_correct, 
                                     mytrait = "narea", 
                                     GINV = re_h2_ef_inv,
                                     data = "highh2",
                                     scheme = "highh2_efmw"
)
corr_N_highh2_efmw<- data.frame(result_narea_highh2_efmw$ac)
fwrite(corr_N_highh2_efmw, "./output/highh2/corr_N_highh2_efmw.csv")
#narea training mw validating ef using highh2 mw matrix
re_h2_mw_inv <- readRDS(file = "./data/relmatrices/re_h2_mw_inv.rds")
result_narea_highh2_mwef<- crossv(sort = sort,
                                     train = na.omit(nareablues_mw_correct),
                                     validation = nareablues_ef_correct, 
                                     mytrait = "narea", 
                                     GINV = re_h2_mw_inv,
                                     data = "highh2",
                                     scheme = "highh2_mwef"
)
corr_N_highh2_mwef<- data.frame(result_narea_highh2_mwef$ac)
fwrite(corr_N_highh2_mwef, "./output/highh2/corr_N_highh2_mwef.csv")
#sla training ef and validating mw using highh2 ef matrix
re_h2_ef_inv <- readRDS(file = "./data/relmatrices/re_h2_ef_inv.rds")
result_sla_highh2_efmw<- crossv(sort = sort,
                                   train = na.omit(slablues_ef_correct),
                                   validation = slablues_mw_correct, 
                                   mytrait = "sla", 
                                   GINV = re_h2_ef_inv,
                                   data = "highh2",
                                   scheme = "highh2_efmw"
)
corr_S_highh2_efmw<- data.frame(result_sla_highh2_efmw$ac)
fwrite(corr_S_highh2_efmw, "./output/highh2/corr_S_highh2_efmw.csv")
#sla training mw validating ef using highh2 mw
re_h2_mw_inv <- readRDS(file = "./data/relmatrices/re_h2_mw_inv.rds")
result_sla_highh2_mwef<- crossv(sort = sort,
                                   train = na.omit(slablues_mw_correct),
                                   validation = slablues_ef_correct, 
                                   mytrait = "sla", 
                                   GINV = re_h2_mw_inv,
                                   data = "highh2",
                                   scheme = "highh2_mwef"
)
corr_S_highh2_mwef<- data.frame(result_sla_highh2_mwef$ac)
fwrite(corr_S_highh2_mwef, "./output/highh2/corr_S_highh2_mwef.csv")




#GBLUP_reduced model after removing 10,25 and 50% indvidual
G10inv_joint <- readRDS(file = "./data/relmatrices/GBLUP/joint/G10inv_joint.rds")
result_N_G10_j<- crossv(sort = sort,
                        train = na.omit(bluesjoint_narea_correct),
                        validation = bluesjoint_narea_correct, 
                        mytrait = "narea", 
                        GINV = G10inv_joint,
                        data = "GBLUP_reduced",
                        scheme = "joint10"
)
corrN_G_10j<- data.frame(G_10 = rep(NA,20))
corrN_G_10j$G_10 <- result_N_G10_j$ac
fwrite(corrN_G_10j, "./output/GBLUP_reduced/corrN_10j.csv")

G25inv_joint <- readRDS(file = "./data/relmatrices/GBLUP/joint/G25inv_joint.rds")
result_N_G25_j<- crossv(sort = sort,
                        train = na.omit(bluesjoint_narea_correct),
                        validation = bluesjoint_narea_correct, 
                        mytrait = "narea", 
                        GINV = G25inv_joint,
                        data = "GBLUP_reduced",
                        scheme = "joint25"
)
corrN_G_25j<- data.frame(G_25 = rep(NA,20))
corrN_G_25j$G_25 <- result_N_G25_j$ac
fwrite(corrN_G_25j, "./output/GBLUP_reduced/corrN_25j.csv")

G50inv_joint <- readRDS(file = "./data/relmatrices/GBLUP/joint/G50inv_joint.rds")
result_N_G50_j<- crossv(sort = sort,
                        train = na.omit(bluesjoint_narea_correct),
                        validation = bluesjoint_narea_correct, 
                        mytrait = "narea", 
                        GINV = G50inv_joint,
                        data = "GBLUP_reduced",
                        scheme = "joint50"
)
corrN_G_50j<- data.frame(G_50 = rep(NA,20))
corrN_G_50j$G_50 <- result_N_G50_j$ac
fwrite(corrN_G_50j, "./output/GBLUP_reduced/corrN_50j.csv")

G10inv_joint <- readRDS(file = "./data/relmatrices/GBLUP/joint/G10inv_joint.rds")
result_S_G10_j<- crossv(sort = sort,
                        train = na.omit(bluesjoint_sla_correct),
                        validation = bluesjoint_sla_correct, 
                        mytrait = "sla", 
                        GINV = G10inv_joint,
                        data = "GBLUP_reduced",
                        scheme = "joint10"
)
corrS_G_10j<- data.frame(G_10 = rep(NA,20))
corrS_G_10j$G_10 <- result_S_G10_j$ac
fwrite(corrS_G_10j, "./output/GBLUP_reduced/corrS_10j.csv")
G25inv_joint <- readRDS(file = "./data/relmatrices/GBLUP/joint/G25inv_joint.rds")
result_S_G25_j<- crossv(sort = sort,
                        train = na.omit(bluesjoint_sla_correct),
                        validation = bluesjoint_sla_correct, 
                        mytrait = "sla", 
                        GINV = G25inv_joint,
                        data = "GBLUP_reduced",
                        scheme = "joint25"
)
corrS_G_25j<- data.frame(G_25 = rep(NA,20))
corrS_G_25j$G_25 <- result_S_G25_j$ac
fwrite(corrS_G_25j, "./output/GBLUP_reduced/corrS_25j.csv")

G50inv_joint <- readRDS(file = "./data/relmatrices/GBLUP/joint/G50inv_joint.rds")
result_S_G50_j<- crossv(sort = sort,
                        train = na.omit(bluesjoint_sla_correct),
                        validation = bluesjoint_sla_correct, 
                        mytrait = "sla", 
                        GINV = G50inv_joint,
                        data = "GBLUP_reduced",
                        scheme = "joint50"
)
corrS_G_50j<- data.frame(G_50 = rep(NA,20))
corrS_G_50j$G_50 <- result_S_G50_j$ac
fwrite(corrS_G_50j, "./output/GBLUP_reduced/corrS_50j.csv")
G10inv_ef <- readRDS(file = "./data/relmatrices/GBLUP/ef/G10inv_ef.rds")
result_N_G10_efmw<- crossv(sort = sort,
                           train = na.omit(nareablues_ef_correct),
                           validation = nareablues_mw_correct, 
                           mytrait = "narea", 
                           GINV = G10inv_ef,
                           data = "GBLUP_reduced",
                           scheme = "efmw10"
)
corrN_G_10efmw<- data.frame(G_10 = rep(NA,20))
corrN_G_10efmw$G_10 <- result_N_G10_efmw$ac
fwrite(corrN_G_10efmw, "./output/GBLUP_reduced/corrN_10efmw.csv")


G25inv_ef <- readRDS(file = "./data/relmatrices/GBLUP/ef/G25inv_ef.rds")
result_N_G25_efmw<- crossv(sort = sort,
                           train = na.omit(nareablues_ef_correct),
                           validation = nareablues_mw_correct, 
                           mytrait = "narea", 
                           GINV = G25inv_ef,
                           data = "GBLUP_reduced",
                           scheme = "efmw25"
)
corrN_G_25efmw<- data.frame(G_25 = rep(NA,20))
corrN_G_25efmw$G_25 <- result_N_G25_efmw$ac
fwrite(corrN_G_25efmw, "./output/GBLUP_reduced/corrN_25efmw.csv")
G50inv_ef <- readRDS(file = "./data/relmatrices/GBLUP/ef/G50inv_ef.rds")
result_N_G50_efmw<- crossv(sort = sort,
                           train = na.omit(nareablues_ef_correct),
                           validation = nareablues_mw_correct, 
                           mytrait = "narea", 
                           GINV = G50inv_ef,
                           data = "GBLUP_reduced",
                           scheme = "efmw50"
)
corrN_G_50efmw<- data.frame(G_50 = rep(NA,20))
corrN_G_50efmw$G_50 <- result_N_G50_efmw$ac
fwrite(corrN_G_50efmw, "./output/GBLUP_reduced/corrN_50efmw.csv")
G10inv_mw <- readRDS(file = "./data/relmatrices/GBLUP/mw/G10inv_mw.rds")
result_N_G10_mwef<- crossv(sort = sort,
                           train = na.omit(nareablues_mw_correct),
                           validation = nareablues_ef_correct, 
                           mytrait = "narea", 
                           GINV = G10inv_mw,
                           data = "GBLUP_reduced",
                           scheme = "mwef10"
)
corrN_G_10mwef<- data.frame(G_10 = rep(NA,20))
corrN_G_10mwef$G_10 <- result_N_G10_mwef$ac
fwrite(corrN_G_10mwef, "./output/GBLUP_reduced/corrN_10mwef.csv")

G25inv_mw <- readRDS(file = "./data/relmatrices/GBLUP/mw/G25inv_mw.rds")
result_N_G25_mwef<- crossv(sort = sort,
                           train = na.omit(nareablues_mw_correct),
                           validation = nareablues_ef_correct, 
                           mytrait = "narea", 
                           GINV = G25inv_mw,
                           data = "GBLUP_reduced",
                           scheme = "mwef25"
)
corrN_G_25mwef<- data.frame(G_25 = rep(NA,20))
corrN_G_25mwef$G_25 <- result_N_G25_mwef$ac
fwrite(corrN_G_25mwef, "./output/GBLUP_reduced/corrN_25mwef.csv")
G50inv_mw <- readRDS(file = "./data/relmatrices/GBLUP/mw/G50inv_mw.rds")
result_N_G50_mwef<- crossv(sort = sort,
                           train = na.omit(nareablues_mw_correct),
                           validation = nareablues_ef_correct, 
                           mytrait = "narea", 
                           GINV = G50inv_mw,
                           data = "GBLUP_reduced",
                           scheme = "mwef50"
)
corrN_50mwef<- data.frame(G_50 = rep(NA,20))
corrN_50mwef$G_50 <- result_N_G50_mwef$ac
fwrite(corrN_50mwef, "./output/GBLUP_reduced/corrN_50mwef.csv")

G10inv_mw <- readRDS(file = "./data/relmatrices/GBLUP/mw/G10inv_mw.rds")
result_S_G10_mwef<- crossv(sort = sort,
                           train = na.omit(slablues_mw_correct),
                           validation = slablues_ef_correct, 
                           mytrait = "sla", 
                           GINV = G10inv_mw,
                           data = "GBLUP_reduced",
                           scheme = "mwef10"
)
corrS_G_10mwef<- data.frame(G_10 = rep(NA,20))
corrS_G_10mwef$G_10 <- result_S_G10_mwef$ac
fwrite(corrS_G_10mwef, "./output/GBLUP_reduced/corrS_10mwef.csv")

G25inv_mw <- readRDS(file = "./data/relmatrices/GBLUP/mw/G25inv_mw.rds")
result_S_G25_mwef<- crossv(sort = sort,
                           train = na.omit(slablues_mw_correct),
                           validation = slablues_ef_correct, 
                           mytrait = "sla", 
                           GINV = G25inv_mw,
                           data = "GBLUP_reduced",
                           scheme = "mwef25"
)
corrS_G_25mwef<- data.frame(G_25 = rep(NA,20))
corrS_G_25mwef$G_25 <- result_S_G25_mwef$ac
fwrite(corrS_G_25mwef, "./output/GBLUP_reduced/corrS_25mwef.csv")
G50inv_mw <- readRDS(file = "./data/relmatrices/GBLUP/mw/G50inv_mw.rds")
result_S_G50_mwef<- crossv(sort = sort,
                           train = na.omit(slablues_mw_correct),
                           validation = slablues_ef_correct, 
                           mytrait = "sla", 
                           GINV = G50inv_mw,
                           data = "GBLUP_reduced",
                           scheme = "mwef50"
)
corrS_G_50mwef<- data.frame(G_50 = rep(NA,20))
corrS_G_50mwef$G_50 <- result_N_G50_mwef$ac
fwrite(corrS_G_50mwef, "./output/GBLUP_reduced/corrS_50mwef.csv")

G10inv_ef <- readRDS(file = "./data/relmatrices/GBLUP/ef/G10inv_ef.rds")
result_S_G10_efmw<- crossv(sort = sort,
                           train = na.omit(slablues_ef_correct),
                           validation = slablues_mw_correct, 
                           mytrait = "sla", 
                           GINV = G10inv_ef,
                           data = "GBLUP_reduced",
                           scheme = "efmw10"
)
corrS_G_10efmw<- data.frame(G_10 = rep(NA,20))
corrS_G_10efmw$G_10 <- result_S_G10_efmw$ac
fwrite(corrS_G_10efmw, "./output/GBLUP_reduced/corrS_10efmw.csv")

G25inv_ef <- readRDS(file = "./data/relmatrices/GBLUP/ef/G25inv_ef.rds")
result_S_G25_efmw<- crossv(sort = sort,
                           train = na.omit(slablues_ef_correct),
                           validation = slablues_mw_correct, 
                           mytrait = "sla", 
                           GINV = G25inv_ef,
                           data = "GBLUP_reduced",
                           scheme = "efmw25"
)
corrS_G_25efmw<- data.frame(G_25 = rep(NA,20))
corrS_G_25efmw$G_25 <- result_S_G25_efmw$ac
fwrite(corrS_G_25efmw, "./output/GBLUP_reduced/corrS_25efmw.csv")
G50inv_ef <- readRDS(file = "./data/relmatrices/GBLUP/ef/G50inv_ef.rds")
result_S_G50_efmw<- crossv(sort = sort,
                           train = na.omit(slablues_ef_correct),
                           validation = slablues_mw_correct, 
                           mytrait = "sla", 
                           GINV = G50inv_ef,
                           data = "GBLUP_reduced",
                           scheme = "efmw50"
)
corrS_G_50efmw<- data.frame(G_50 = rep(NA,20))
corrS_G_50efmw$G_50 <- result_S_G50_efmw$ac
fwrite(corrS_G_50efmw, "./output/GBLUP_reduced/corrS_50efmw.csv")


#Gh2 model
# --cv for narea joint location using Gh2 matrix with 10%,25%, 50% scheme ---------------------
Gh2inv_joint <- readRDS(file = "./data/relmatrices/Gh2/joint/10/Gh2inv_joint.rds")
result_N_Gh2_j<- crossv(sort = sort,
                                train = na.omit(bluesjoint_narea_correct),
                                validation = bluesjoint_narea_correct, 
                                mytrait = "narea", 
                                GINV = Gh2inv_joint,
                                data = "Gh2_joint_narea",
                                scheme = "Gh2_10_joint"
)
corrN_Gh2_10j<- data.frame(Gh2_10 = rep(NA,20))
corrN_Gh2_10j$Gh2_10 <- result_N_Gh2_j$ac
fwrite(corrN_Gh2_10j, "./output/Gh2_joint_narea/corrN_Gh2_10j.csv")
# ---------------------25%---------------------
Gh2inv_joint <- readRDS(file = "./data/relmatrices/Gh2/joint/25/Gh2inv_joint.rds")
result_N_Gh2_j<- crossv(sort = sort,
                        train = na.omit(bluesjoint_narea_correct),
                        validation = bluesjoint_narea_correct, 
                        mytrait = "narea", 
                        GINV = Gh2inv_joint,
                        data = "Gh2_joint_narea",
                        scheme = "Gh2_25_joint"
)
corrN_Gh2_25j<- data.frame(Gh2_25 = rep(NA,20))
corrN_Gh2_25j$Gh2_25 <- result_N_Gh2_j$ac
fwrite(corrN_Gh2_25j, "./output/Gh2_joint_narea/corrN_Gh2_25j.csv")
# ---------------------50%---------------------
Gh2inv_joint <- readRDS(file = "./data/relmatrices/Gh2/joint/50/Gh2inv_joint.rds")
result_N_Gh2_j<- crossv(sort = sort,
                        train = na.omit(bluesjoint_narea_correct),
                        validation = bluesjoint_narea_correct, 
                        mytrait = "narea", 
                        GINV = Gh2inv_joint,
                        data = "Gh2_joint_narea",
                        scheme = "Gh2_50_joint"
)
corrN_Gh2_50j<- data.frame(Gh2_50 = rep(NA,20))
corrN_Gh2_50j$Gh2_50 <- result_N_Gh2_j$ac
fwrite(corrN_Gh2_50j, "./output/Gh2_joint_narea/corrN_Gh2_50j.csv")
# -------cv for sla joint location using Gh2 matrix with 10%,25%, 50% scheme ---------------------
Gh2inv_joint <- readRDS(file = "./data/relmatrices/Gh2/joint/10/Gh2inv_joint.rds")
result_S_Gh2_j<- crossv(sort = sort,
                        train = na.omit(bluesjoint_sla_correct),
                        validation = bluesjoint_sla_correct, 
                        mytrait = "sla", 
                        GINV = Gh2inv_joint,
                        data = "Gh2_joint_sla",
                        scheme = "Gh2_10_joint"
)
corrS_Gh2_10j<- data.frame(Gh2_10 = rep(NA,20))
corrS_Gh2_10j$Gh2_10 <- result_S_Gh2_j$ac
fwrite(corrS_Gh2_10j, "./output/Gh2_joint_sla/corrS_Gh2_10j.csv")
# ---------------------25%---------------------
Gh2inv_joint <- readRDS(file = "./data/relmatrices/Gh2/joint/25/Gh2inv_joint.rds")
result_S_Gh2_j<- crossv(sort = sort,
                        train = na.omit(bluesjoint_sla_correct),
                        validation = bluesjoint_sla_correct, 
                        mytrait = "sla", 
                        GINV = Gh2inv_joint,
                        data = "Gh2_joint_sla",
                        scheme = "Gh2_25_joint"
)
corrS_Gh2_25j<- data.frame(Gh2_25 = rep(NA,20))
corrS_Gh2_25j$Gh2_25 <- result_S_Gh2_j$ac
fwrite(corrS_Gh2_25j, "./output/Gh2_joint_sla/corrS_Gh2_25j.csv")
# ---------------------50%---------------------
Gh2inv_joint <- readRDS(file = "./data/relmatrices/Gh2/joint/50/Gh2inv_joint.rds")
result_S_Gh2_j<- crossv(sort = sort,
                        train = na.omit(bluesjoint_sla_correct),
                        validation = bluesjoint_sla_correct, 
                        mytrait = "sla", 
                        GINV = Gh2inv_joint,
                        data = "Gh2_joint_sla",
                        scheme = "Gh2_50_joint"
)
corrS_Gh2_50j<- data.frame(Gh2_50 = rep(NA,20))
corrS_Gh2_50j$Gh2_50 <- result_S_Gh2_j$ac
fwrite(corrS_Gh2_50j, "./output/Gh2_joint_sla/corrS_Gh2_50j.csv")
# ----cv for narea train ef val on mw using Gh2 matrix with 10%,25%, 50% scheme ---------------------
Gh2inv_ef <- readRDS(file = "./data/relmatrices/Gh2/ef/10/Gh2inv_ef.rds")
result_N_Gh2_efmw<- crossv(sort = sort,
                        train = na.omit(nareablues_ef_correct),
                        validation = nareablues_mw_correct, 
                        mytrait = "narea", 
                        GINV = Gh2inv_ef,
                        data = "Gh2_efmw_narea",
                        scheme = "Gh2_10_efmw"
)
corrN_Gh2_10efmw<- data.frame(Gh2_10 = rep(NA,20))
corrN_Gh2_10efmw$Gh2_10 <- result_N_Gh2_efmw$ac
fwrite(corrN_Gh2_10efmw, "./output/Gh2_efmw_narea/corrN_Gh2_10efmw.csv")
# ---------------------25%---------------------
Gh2inv_ef <- readRDS(file = "./data/relmatrices/Gh2/ef/25/Gh2inv_ef.rds")
result_N_Gh2_efmw<- crossv(sort = sort,
                           train = na.omit(nareablues_ef_correct),
                           validation = nareablues_mw_correct, 
                           mytrait = "narea", 
                           GINV = Gh2inv_ef,
                           data = "Gh2_efmw_narea",
                           scheme = "Gh2_25_efmw"
)
corrN_Gh2_25efmw<- data.frame(Gh2_25 = rep(NA,20))
corrN_Gh2_25efmw$Gh2_25 <- result_N_Gh2_efmw$ac
fwrite(corrN_Gh2_25efmw, "./output/Gh2_efmw_narea/corrN_Gh2_25efmw.csv")
# ---------------------50%---------------------
Gh2inv_ef <- readRDS(file = "./data/relmatrices/Gh2/ef/50/Gh2inv_ef.rds")
result_N_Gh2_efmw<- crossv(sort = sort,
                           train = na.omit(nareablues_ef_correct),
                           validation = nareablues_mw_correct, 
                           mytrait = "narea", 
                           GINV = Gh2inv_ef,
                           data = "Gh2_efmw_narea",
                           scheme = "Gh2_50_efmw"
)
corrN_Gh2_50efmw<- data.frame(Gh2_50 = rep(NA,20))
corrN_Gh2_50efmw$Gh2_50 <- result_N_Gh2_efmw$ac
fwrite(corrN_Gh2_50efmw, "./output/Gh2_efmw_narea/corrN_Gh2_50efmw.csv")


# ----cv for narea train mw val on ef using Gh2 matrix with 10%,25%, 50% scheme ---------------------
Gh2inv_mw <- readRDS(file = "./data/relmatrices/Gh2/mw/10/Gh2inv_mw.rds")
result_N_Gh2_mwef<- crossv(sort = sort,
                           train = na.omit(nareablues_mw_correct),
                           validation = nareablues_ef_correct, 
                           mytrait = "narea", 
                           GINV = Gh2inv_mw,
                           data = "Gh2_mwef_narea",
                           scheme = "Gh2_10_mwef"
)
corrN_Gh2_10mwef<- data.frame(Gh2_10 = rep(NA,20))
corrN_Gh2_10mwef$Gh2_10 <- result_N_Gh2_mwef$ac
fwrite(corrN_Gh2_10mwef, "./output/Gh2_mwef_narea/corrN_Gh2_10mwef.csv")
# ---------------------25%---------------------
Gh2inv_mw <- readRDS(file = "./data/relmatrices/Gh2/mw/25/Gh2inv_mw.rds")
result_N_Gh2_mwef<- crossv(sort = sort,
                           train = na.omit(nareablues_mw_correct),
                           validation = nareablues_ef_correct, 
                           mytrait = "narea", 
                           GINV = Gh2inv_mw,
                           data = "Gh2_mwef_narea",
                           scheme = "Gh2_25_mwef"
)
corrN_Gh2_25mwef<- data.frame(Gh2_25 = rep(NA,20))
corrN_Gh2_25mwef$Gh2_25 <- result_N_Gh2_mwef$ac
fwrite(corrN_Gh2_25mwef, "./output/Gh2_mwef_narea/corrN_Gh2_25mwef.csv")
# ---------------------50%---------------------
Gh2inv_mw <- readRDS(file = "./data/relmatrices/Gh2/mw/50/Gh2inv_mw.rds")
result_N_Gh2_mwef<- crossv(sort = sort,
                           train = na.omit(nareablues_mw_correct),
                           validation = nareablues_ef_correct, 
                           mytrait = "narea", 
                           GINV = Gh2inv_mw,
                           data = "Gh2_mwef_narea",
                           scheme = "Gh2_50_mwef"
)
corrN_Gh2_50mwef<- data.frame(Gh2_50 = rep(NA,20))
corrN_Gh2_50mwef$Gh2_50 <- result_N_Gh2_mwef$ac
fwrite(corrN_Gh2_50mwef, "./output/Gh2_mwef_narea/corrN_Gh2_50mwef.csv")       
# ---------------------cv for sla train mw val on ef using Gh2 matrix with 10%,25%, 50% scheme ---------------------
Gh2inv_mw <- readRDS(file = "./data/relmatrices/Gh2/mw/10/Gh2inv_mw.rds")
result_S_Gh2_mwef<- crossv(sort = sort,
                           train = na.omit(slablues_mw_correct),
                           validation = slablues_ef_correct, 
                           mytrait = "sla", 
                           GINV = Gh2inv_mw,
                           data = "Gh2_mwef_sla",
                           scheme = "Gh2_10_mwef"
)
corrS_Gh2_10mwef<- data.frame(Gh2_10 = rep(NA,20))
corrS_Gh2_10mwef$Gh2_10 <- result_S_Gh2_mwef$ac
fwrite(corrS_Gh2_10mwef, "./output/Gh2_mwef_sla/corrS_Gh2_10mwef.csv")
# ---------------------25%---------------------
Gh2inv_mw <- readRDS(file = "./data/relmatrices/Gh2/mw/25/Gh2inv_mw.rds")
result_S_Gh2_mwef<- crossv(sort = sort,
                           train = na.omit(slablues_mw_correct),
                           validation = slablues_ef_correct, 
                           mytrait = "sla", 
                           GINV = Gh2inv_mw,
                           data = "Gh2_mwef_sla",
                           scheme = "Gh2_25_mwef"
)
corrS_Gh2_25mwef<- data.frame(Gh2_25 = rep(NA,20))
corrS_Gh2_25mwef$Gh2_25 <- result_S_Gh2_mwef$ac
fwrite(corrS_Gh2_25mwef, "./output/Gh2_mwef_sla/corrS_Gh2_25mwef.csv")
# ---------------------50%---------------------
Gh2inv_mw <- readRDS(file = "./data/relmatrices/Gh2/mw/50/Gh2inv_mw.rds")
result_S_Gh2_mwef<- crossv(sort = sort,
                           train = na.omit(slablues_mw_correct),
                           validation = slablues_ef_correct, 
                           mytrait = "sla", 
                           GINV = Gh2inv_mw,
                           data = "Gh2_mwef_sla",
                           scheme = "Gh2_50_mwef"
)
corrS_Gh2_50mwef<- data.frame(Gh2_50 = rep(NA,20))
corrS_Gh2_50mwef$Gh2_50 <- result_S_Gh2_mwef$ac
fwrite(corrS_Gh2_50mwef, "./output/Gh2_mwef_sla/corrS_Gh2_50mwef.csv")       
# -----cv for narea train ef val on mw using Gh2 matrix with 10%,25%, 50% scheme ---------------------
Gh2inv_ef <- readRDS(file = "./data/relmatrices/Gh2/ef/10/Gh2inv_ef.rds")
result_S_Gh2_efmw<- crossv(sort = sort,
                           train = na.omit(slablues_ef_correct),
                           validation = slablues_mw_correct, 
                           mytrait = "sla", 
                           GINV = Gh2inv_ef,
                           data = "Gh2_efmw_sla",
                           scheme = "Gh2_10_efmw"
)
corrS_Gh2_10efmw<- data.frame(Gh2_10 = rep(NA,20))
corrS_Gh2_10efmw$Gh2_10 <- result_S_Gh2_efmw$ac
fwrite(corrS_Gh2_10efmw, "./output/Gh2_efmw_sla/corrS_Gh2_10efmw.csv")
# ---------------------25%---------------------
Gh2inv_ef <- readRDS(file = "./data/relmatrices/Gh2/ef/25/Gh2inv_ef.rds")
result_S_Gh2_efmw<- crossv(sort = sort,
                           train = na.omit(slablues_ef_correct),
                           validation = slablues_mw_correct, 
                           mytrait = "sla", 
                           GINV = Gh2inv_ef,
                           data = "Gh2_efmw_sla",
                           scheme = "Gh2_25_efmw"
)
corrS_Gh2_25efmw<- data.frame(Gh2_25 = rep(NA,20))
corrS_Gh2_25efmw$Gh2_25 <- result_S_Gh2_efmw$ac
fwrite(corrS_Gh2_25efmw, "./output/Gh2_efmw_sla/corrS_Gh2_25efmw.csv")

# ---------------------50%---------------------
Gh2inv_ef <- readRDS(file = "./data/relmatrices/Gh2/ef/50/Gh2inv_ef.rds")
result_S_Gh2_efmw<- crossv(sort = sort,
                           train = na.omit(slablues_ef_correct),
                           validation = slablues_mw_correct, 
                           mytrait = "sla", 
                           GINV = Gh2inv_ef,
                           data = "Gh2_efmw_sla",
                           scheme = "Gh2_50_efmw"
)
corrS_Gh2_50efmw<- data.frame(Gh2_50 = rep(NA,20))
corrS_Gh2_50efmw$Gh2_50 <- result_S_Gh2_efmw$ac
fwrite(corrS_Gh2_50efmw, "./output/Gh2_efmw_sla/corrS_Gh2_50efmw.csv")     
#Gnirs model
# ---------------------cv for narea joint location using Gnirs matrix with 10%,25%, 50% scheme ---------------------
Gnirsinv_joint <- readRDS(file = "./data/relmatrices/Gnirs/joint/10/Gnirsinv_joint.rds")
result_N_Gnirs_j<- crossv(sort = sort,
                        train = na.omit(bluesjoint_narea_correct),
                        validation = bluesjoint_narea_correct, 
                        mytrait = "narea", 
                        GINV = Gnirsinv_joint,
                        data = "Gnirs_joint_narea",
                        scheme = "Gnirs_10_joint"
)
corrN_Gnirs_10j<- data.frame(Gnirs_10 = rep(NA,20))
corrN_Gnirs_10j$Gnirs_10 <- result_N_Gnirs_j$ac
fwrite(corrN_Gnirs_10j, "./output/Gnirs_joint_narea/corrN_Gnirs_10j.csv")
# ---------------------25%---------------------
Gnirsinv_joint <- readRDS(file = "./data/relmatrices/Gnirs/joint/25/Gnirsinv_joint.rds")
result_N_Gnirs_j<- crossv(sort = sort,
                        train = na.omit(bluesjoint_narea_correct),
                        validation = bluesjoint_narea_correct, 
                        mytrait = "narea", 
                        GINV = Gnirsinv_joint,
                        data = "Gnirs_joint_narea",
                        scheme = "Gnirs_25_joint"
)
corrN_Gnirs_25j<- data.frame(Gnirs_25 = rep(NA,20))
corrN_Gnirs_25j$Gnirs_25 <- result_N_Gnirs_j$ac
fwrite(corrN_Gnirs_25j, "./output/Gnirs_joint_narea/corrN_Gnirs_25j.csv")
# ---------------------50%---------------------
Gnirsinv_joint <- readRDS(file = "./data/relmatrices/Gnirs/joint/50/Gnirsinv_joint.rds")
result_N_Gnirs_j<- crossv(sort = sort,
                        train = na.omit(bluesjoint_narea_correct),
                        validation = bluesjoint_narea_correct, 
                        mytrait = "narea", 
                        GINV = Gnirsinv_joint,
                        data = "Gnirs_joint_narea",
                        scheme = "Gnirs_50_joint"
)
corrN_Gnirs_50j<- data.frame(Gnirs_50 = rep(NA,20))
corrN_Gnirs_50j$Gnirs_50 <- result_N_Gnirs_j$ac
fwrite(corrN_Gnirs_50j, "./output/Gnirs_joint_narea/corrN_Gnirs_50j.csv")
# ---------------------cv for sla joint location using Gnirs matrix with 10%,25%, 50% scheme ---------------------
Gnirsinv_joint <- readRDS(file = "./data/relmatrices/Gnirs/joint/10/Gnirsinv_joint.rds")
result_S_Gnirs_j<- crossv(sort = sort,
                        train = na.omit(bluesjoint_sla_correct),
                        validation = bluesjoint_sla_correct, 
                        mytrait = "sla", 
                        GINV = Gnirsinv_joint,
                        data = "Gnirs_joint_sla",
                        scheme = "Gnirs_10_joint"
)
corrS_Gnirs_10j<- data.frame(Gnirs_10 = rep(NA,20))
corrS_Gnirs_10j$Gnirs_10 <- result_S_Gnirs_j$ac
fwrite(corrS_Gnirs_10j, "./output/Gnirs_joint_sla/corrS_Gnirs_10j.csv")
# ---------------------25%---------------------
Gnirsinv_joint <- readRDS(file = "./data/relmatrices/Gnirs/joint/25/Gnirsinv_joint.rds")
result_S_Gnirs_j<- crossv(sort = sort,
                        train = na.omit(bluesjoint_sla_correct),
                        validation = bluesjoint_sla_correct, 
                        mytrait = "sla", 
                        GINV = Gnirsinv_joint,
                        data = "Gnirs_joint_sla",
                        scheme = "Gnirs_25_joint"
)
corrS_Gnirs_25j<- data.frame(Gnirs_25 = rep(NA,20))
corrS_Gnirs_25j$Gnirs_25 <- result_S_Gnirs_j$ac
fwrite(corrS_Gnirs_25j, "./output/Gnirs_joint_sla/corrS_Gnirs_25j.csv")


# ---------------------50%---------------------
Gnirsinv_joint <- readRDS(file = "./data/relmatrices/Gnirs/joint/50/Gnirsinv_joint.rds")
result_S_Gnirs_j<- crossv(sort = sort,
                        train = na.omit(bluesjoint_sla_correct),
                        validation = bluesjoint_sla_correct, 
                        mytrait = "sla", 
                        GINV = Gnirsinv_joint,
                        data = "Gnirs_joint_sla",
                        scheme = "Gnirs_50_joint"
)
corrS_Gnirs_50j<- data.frame(Gnirs_50 = rep(NA,20))
corrS_Gnirs_50j$Gnirs_50 <- result_S_Gnirs_j$ac
fwrite(corrS_Gnirs_50j, "./output/Gnirs_joint_sla/corrS_Gnirs_50j.csv")
# ---------------------cv for narea train ef val on mw using Gnirs matrix with 10%,25%, 50% scheme ---------------------
Gnirsinv_ef <- readRDS(file = "./data/relmatrices/Gnirs/ef/10/Gnirsinv_ef.rds")
result_N_Gnirs_efmw<- crossv(sort = sort,
                           train = na.omit(nareablues_ef_correct),
                           validation = nareablues_mw_correct, 
                           mytrait = "narea", 
                           GINV = Gnirsinv_ef,
                           data = "Gnirs_efmw_narea",
                           scheme = "Gnirs_10_efmw"
)
corrN_Gnirs_10efmw<- data.frame(Gnirs_10 = rep(NA,20))
corrN_Gnirs_10efmw$Gnirs_10 <- result_N_Gnirs_efmw$ac
fwrite(corrN_Gnirs_10efmw, "./output/Gnirs_efmw_narea/corrN_Gnirs_10efmw.csv")
# ---------------------25%---------------------
Gnirsinv_ef <- readRDS(file = "./data/relmatrices/Gnirs/ef/25/Gnirsinv_ef.rds")
result_N_Gnirs_efmw<- crossv(sort = sort,
                           train = na.omit(nareablues_ef_correct),
                           validation = nareablues_mw_correct, 
                           mytrait = "narea", 
                           GINV = Gnirsinv_ef,
                           data = "Gnirs_efmw_narea",
                           scheme = "Gnirs_25_efmw"
)
corrN_Gnirs_25efmw<- data.frame(Gnirs_25 = rep(NA,20))
corrN_Gnirs_25efmw$Gnirs_25 <- result_N_Gnirs_efmw$ac
fwrite(corrN_Gnirs_25efmw, "./output/Gnirs_efmw_narea/corrN_Gnirs_25efmw.csv")
# ---------------------50%---------------------
Gnirsinv_ef <- readRDS(file = "./data/relmatrices/Gnirs/ef/50/Gnirsinv_ef.rds")
result_N_Gnirs_efmw<- crossv(sort = sort,
                           train = na.omit(nareablues_ef_correct),
                           validation = nareablues_mw_correct, 
                           mytrait = "narea", 
                           GINV = Gnirsinv_ef,
                           data = "Gnirs_efmw_narea",
                           scheme = "Gnirs_50_efmw"
)
corrN_Gnirs_50efmw<- data.frame(Gnirs_50 = rep(NA,20))
corrN_Gnirs_50efmw$Gnirs_50 <- result_N_Gnirs_efmw$ac
fwrite(corrN_Gnirs_50efmw, "./output/Gnirs_efmw_narea/corrN_Gnirs_50efmw.csv")     

# ---------------------cv for narea train mw val on ef using Gnirs matrix with 10%,25%, 50% scheme ---------------------
Gnirsinv_mw <- readRDS(file = "./data/relmatrices/Gnirs/mw/10/Gnirsinv_mw.rds")
result_N_Gnirs_mwef<- crossv(sort = sort,
                           train = na.omit(nareablues_mw_correct),
                           validation = nareablues_ef_correct, 
                           mytrait = "narea", 
                           GINV = Gnirsinv_mw,
                           data = "Gnirs_mwef_narea",
                           scheme = "Gnirs_10_mwef"
)
corrN_Gnirs_10mwef<- data.frame(Gnirs_10 = rep(NA,20))
corrN_Gnirs_10mwef$Gnirs_10 <- result_N_Gnirs_mwef$ac
fwrite(corrN_Gnirs_10mwef, "./output/Gnirs_mwef_narea/corrN_Gnirs_10mwef.csv")
# ---------------------25%---------------------
Gnirsinv_mw <- readRDS(file = "./data/relmatrices/Gnirs/mw/25/Gnirsinv_mw.rds")
result_N_Gnirs_mwef<- crossv(sort = sort,
                           train = na.omit(nareablues_mw_correct),
                           validation = nareablues_ef_correct, 
                           mytrait = "narea", 
                           GINV = Gnirsinv_mw,
                           data = "Gnirs_mwef_narea",
                           scheme = "Gnirs_25_mwef"
)
corrN_Gnirs_25mwef<- data.frame(Gnirs_25 = rep(NA,20))
corrN_Gnirs_25mwef$Gnirs_25 <- result_N_Gnirs_mwef$ac
fwrite(corrN_Gnirs_25mwef, "./output/Gnirs_mwef_narea/corrN_Gnirs_25mwef.csv")

# ---------------------50%---------------------
Gnirsinv_mw <- readRDS(file = "./data/relmatrices/Gnirs/mw/50/Gnirsinv_mw.rds")
result_N_Gnirs_mwef<- crossv(sort = sort,
                           train = na.omit(nareablues_mw_correct),
                           validation = nareablues_ef_correct, 
                           mytrait = "narea", 
                           GINV = Gnirsinv_mw,
                           data = "Gnirs_mwef_narea",
                           scheme = "Gnirs_50_mwef"
)
corrN_Gnirs_50mwef<- data.frame(Gnirs_50 = rep(NA,20))
corrN_Gnirs_50mwef$Gnirs_50 <- result_N_Gnirs_mwef$ac
fwrite(corrN_Gnirs_50mwef, "./output/Gnirs_mwef_narea/corrN_Gnirs_50mwef.csv") 

# ---------------------cv for sla train mw val on ef using Gnirs matrix with 10%,25%, 50% scheme ---------------------
Gnirsinv_mw <- readRDS(file = "./data/relmatrices/Gh2/mw/10/Gh2inv_mw.rds")
result_S_Gnirs_mwef<- crossv(sort = sort,
                           train = na.omit(slablues_mw_correct),
                           validation = slablues_ef_correct, 
                           mytrait = "sla", 
                           GINV = Gnirsinv_mw,
                           data = "Gnirs_mwef_sla",
                           scheme = "Gnirs_10_mwef"
)
corrS_Gnirs_10mwef<- data.frame(Gnirs_10 = rep(NA,20))
corrS_Gnirs_10mwef$Gnirs_10 <- result_S_Gnirs_mwef$ac
fwrite(corrS_Gnirs_10mwef, "./output/Gnirs_mwef_sla/corrS_Gnirs_10mwef.csv")

# ---------------------25%---------------------
Gnirsinv_mw <- readRDS(file = "./data/relmatrices/Gnirs/mw/25/Gnirsinv_mw.rds")
result_S_Gnirs_mwef<- crossv(sort = sort,
                           train = na.omit(slablues_mw_correct),
                           validation = slablues_ef_correct, 
                           mytrait = "sla", 
                           GINV = Gnirsinv_mw,
                           data = "Gnirs_mwef_sla",
                           scheme = "Gnirs_25_mwef"
)
corrS_Gnirs_25mwef<- data.frame(Gnirs_25 = rep(NA,20))
corrS_Gnirs_25mwef$Gnirs_25 <- result_S_Gnirs_mwef$ac
fwrite(corrS_Gnirs_25mwef, "./output/Gnirs_mwef_sla/corrS_Gnirs_25mwef.csv")
# ---------------------50%---------------------
Gnirsinv_mw <- readRDS(file = "./data/relmatrices/Gnirs/mw/50/Gnirsinv_mw.rds")
result_S_Gnirs_mwef<- crossv(sort = sort,
                           train = na.omit(slablues_mw_correct),
                           validation = slablues_ef_correct, 
                           mytrait = "sla", 
                           GINV = Gnirsinv_mw,
                           data = "Gnirs_mwef_sla",
                           scheme = "Gnirs_50_mwef"
)
corrS_Gnirs_50mwef<- data.frame(Gnirs_50 = rep(NA,20))
corrS_Gnirs_50mwef$Gnirs_50 <- result_S_Gnirs_mwef$ac
fwrite(corrS_Gnirs_50mwef, "./output/Gnirs_mwef_sla/corrS_Gnirs_50mwef.csv")    

# ---------------------cv for sla train ef val on mw using Gnirs matrix with 10%,25%, 50% scheme ---------------------
Gnirsinv_ef <- readRDS(file = "./data/relmatrices/Gnirs/ef/10/Gnirsinv_ef.rds")
result_S_Gnirs_efmw<- crossv(sort = sort,
                           train = na.omit(slablues_ef_correct),
                           validation = slablues_mw_correct, 
                           mytrait = "sla", 
                           GINV = Gnirsinv_ef,
                           data = "Gnirs_efmw_sla",
                           scheme = "Gnirs_10_efmw"
)
corrS_Gnirs_10efmw<- data.frame(Gnirs_10 = rep(NA,20))
corrS_Gnirs_10efmw$Gnirs_10 <- result_S_Gnirs_efmw$ac
fwrite(corrS_Gnirs_10efmw, "./output/Gnirs_efmw_sla/corrS_Gnirs_10efmw.csv")
# ---------------------25%---------------------
Gnirsinv_ef <- readRDS(file = "./data/relmatrices/Gnirs/ef/25/Gnirsinv_ef.rds")
result_S_Gnirs_efmw<- crossv(sort = sort,
                           train = na.omit(slablues_ef_correct),
                           validation = slablues_mw_correct, 
                           mytrait = "sla", 
                           GINV = Gnirsinv_ef,
                           data = "Gnirs_efmw_sla",
                           scheme = "Gnirs_25_efmw"
)
corrS_Gnirs_25efmw<- data.frame(Gnirs_25 = rep(NA,20))
corrS_Gnirs_25efmw$Gnirs_25 <- result_S_Gnirs_efmw$ac
fwrite(corrS_Gnirs_25efmw, "./output/Gnirs_efmw_sla/corrS_Gnirs_25efmw.csv")
# ---------------------50%---------------------
Gnirsinv_ef <- readRDS(file = "./data/relmatrices/Gnirs/ef/50/Gnirsinv_ef.rds")
result_S_Gnirs_efmw<- crossv(sort = sort,
                           train = na.omit(slablues_ef_correct),
                           validation = slablues_mw_correct, 
                           mytrait = "sla", 
                           GINV = Gnirsinv_ef,
                           data = "Gnirs_efmw_sla",
                           scheme = "Gnirs_50_efmw"
)
corrS_Gnirs_50efmw<- data.frame(Gnirs_50 = rep(NA,20))
corrS_Gnirs_50efmw$Gnirs_50 <- result_S_Gnirs_efmw$ac
fwrite(corrS_Gnirs_50efmw, "./output/Gnirs_efmw_sla/corrS_Gnirs_50efmw.csv")     

#GWW(Genomic Wholewave) model
# ---------------------CV for narea joint location using Gnirs matrix with 10%,25%, 50% scheme ---------------------
GWWinv_joint <- readRDS(file = "./data/relmatrices/GWW/joint/10/GWWinv_joint.rds")
result_N_GWW_j<- crossv(sort = sort,
                          train = na.omit(bluesjoint_narea_correct),
                          validation = bluesjoint_narea_correct, 
                          mytrait = "narea", 
                          GINV = GWWinv_joint,
                          data = "GWW_joint_narea",
                          scheme = "GWW_10_joint"
)
corrN_GWW_10j<- data.frame(GWW_10 = rep(NA,20))
corrN_GWW_10j$GWW_10 <- result_N_GWW_j$ac
fwrite(corrN_GWW_10j, "./output/GWW_joint_narea/corrN_GWW_10j.csv")

# ---------------------25%---------------------
GWWinv_joint <- readRDS(file = "./data/relmatrices/GWW/joint/25/GWWinv_joint.rds")
result_N_GWW_j<- crossv(sort = sort,
                          train = na.omit(bluesjoint_narea_correct),
                          validation = bluesjoint_narea_correct, 
                          mytrait = "narea", 
                          GINV = GWWinv_joint,
                          data = "GWW_joint_narea",
                          scheme = "GWW_25_joint"
)
corrN_GWW_25j<- data.frame(GWW_25 = rep(NA,20))
corrN_GWW_25j$GWW_25 <- result_N_GWW_j$ac
fwrite(corrN_GWW_25j, "./output/GWW_joint_narea/corrN_GWW_25j.csv")
# ---------------------50%---------------------
GWWinv_joint <- readRDS(file = "./data/relmatrices/GWW/joint/50/GWWinv_joint.rds")
result_N_GWW_j<- crossv(sort = sort,
                          train = na.omit(bluesjoint_narea_correct),
                          validation = bluesjoint_narea_correct, 
                          mytrait = "narea", 
                          GINV = GWWinv_joint,
                          data = "GWW_joint_narea",
                          scheme = "GWW_50_joint"
)
corrN_GWW_50j<- data.frame(GWW_50 = rep(NA,20))
corrN_GWW_50j$GWW_50 <- result_N_GWW_j$ac
fwrite(corrN_GWW_50j, "./output/GWW_joint_narea/corrN_GWW_50j.csv")
# ---------------------cv for sla joint location using GWW matrix with 10%,25%, 50% scheme ---------------------
GWWinv_joint <- readRDS(file = "./data/relmatrices/GWW/joint/10/GWWinv_joint.rds")
result_S_GWW_j<- crossv(sort = sort,
                          train = na.omit(bluesjoint_sla_correct),
                          validation = bluesjoint_sla_correct, 
                          mytrait = "sla", 
                          GINV = GWWinv_joint,
                          data = "GWW_joint_sla",
                          scheme = "GWW_10_joint"
)
corrS_GWW_10j<- data.frame(GWW_10 = rep(NA,20))
corrS_GWW_10j$GWW_10 <- result_S_GWW_j$ac
fwrite(corrS_GWW_10j, "./output/GWW_joint_sla/corrS_GWW_10j.csv")
# ---------------------25%---------------------
GWWinv_joint <- readRDS(file = "./data/relmatrices/GWW/joint/25/GWWinv_joint.rds")
result_S_GWW_j<- crossv(sort = sort,
                          train = na.omit(bluesjoint_sla_correct),
                          validation = bluesjoint_sla_correct, 
                          mytrait = "sla", 
                          GINV = GWWinv_joint,
                          data = "GWW_joint_sla",
                          scheme = "GWW_25_joint"
)
corrS_GWW_25j<- data.frame(GWW_25 = rep(NA,20))
corrS_GWW_25j$GWW_25 <- result_S_GWW_j$ac
fwrite(corrS_GWW_25j, "./output/GWW_joint_sla/corrS_GWW_25j.csv")

# ---------------------50%---------------------
GWWinv_joint <- readRDS(file = "./data/relmatrices/Gnirs/joint/50/Gnirsinv_joint.rds")
result_S_GWW_j<- crossv(sort = sort,
                          train = na.omit(bluesjoint_sla_correct),
                          validation = bluesjoint_sla_correct, 
                          mytrait = "sla", 
                          GINV = GWWinv_joint,
                          data = "GWW_joint_sla",
                          scheme = "GWW_50_joint"
)
corrS_GWW_50j<- data.frame(GWW_50 = rep(NA,20))
corrS_GWW_50j$GWW_50 <- result_S_GWW_j$ac
fwrite(corrS_GWW_50j, "./output/GWW_joint_sla/corrS_GWW_50j.csv")

# ----------cv for narea train ef val on mw using GWW matrix with 10%,25%, 50% scheme ---------------------
GWWinv_ef <- readRDS(file = "./data/relmatrices/GWW/ef/10/GWWinv_ef.rds")
result_N_GWW_efmw<- crossv(sort = sort,
                             train = na.omit(nareablues_ef_correct),
                             validation = nareablues_mw_correct, 
                             mytrait = "narea", 
                             GINV = GWWinv_ef,
                             data = "GWW_efmw_narea",
                             scheme = "GWW_10_efmw"
)
corrN_GWW_10efmw<- data.frame(GWW_10 = rep(NA,20))
corrN_GWW_10efmw$GWW_10 <- result_N_GWW_efmw$ac
fwrite(corrN_GWW_10efmw, "./output/GWW_efmw_narea/corrN_GWW_10efmw.csv")
# ---------------------25%---------------------
GWWinv_ef <- readRDS(file = "./data/relmatrices/GWW/ef/25/GWWinv_ef.rds")
result_N_GWW_efmw<- crossv(sort = sort,
                             train = na.omit(nareablues_ef_correct),
                             validation = nareablues_mw_correct, 
                             mytrait = "narea", 
                             GINV = GWWinv_ef,
                             data = "GWW_efmw_narea",
                             scheme = "GWW_25_efmw"
)
corrN_GWW_25efmw<- data.frame(GWW_25 = rep(NA,20))
corrN_GWW_25efmw$GWW_25 <- result_N_GWW_efmw$ac
fwrite(corrN_GWW_25efmw, "./output/GWW_efmw_narea/corrN_GWW_25efmw.csv")
# ---------------------50%---------------------
GWWinv_ef <- readRDS(file = "./data/relmatrices/GWW/ef/50/GWWinv_ef.rds")
result_N_GWW_efmw<- crossv(sort = sort,
                             train = na.omit(nareablues_ef_correct),
                             validation = nareablues_mw_correct, 
                             mytrait = "narea", 
                             GINV = GWWinv_ef,
                             data = "GWW_efmw_narea",
                             scheme = "GWW_50_efmw"
)
corrN_GWW_50efmw<- data.frame(GWW_50 = rep(NA,20))
corrN_GWW_50efmw$GWW_50 <- result_N_GWW_efmw$ac
fwrite(corrN_GWW_50efmw, "./output/GWW_efmw_narea/corrN_GWW_50efmw.csv")  


# ---------------------cv for narea train mw val on ef using GWW matrix with 10%,25%, 50% scheme ---------------------
GWWinv_mw <- readRDS(file = "./data/relmatrices/GWW/mw/10/GWWinv_mw.rds")
result_N_GWW_mwef<- crossv(sort = sort,
                             train = na.omit(nareablues_mw_correct),
                             validation = nareablues_ef_correct, 
                             mytrait = "narea", 
                             GINV = GWWinv_mw,
                             data = "GWW_mwef_narea",
                             scheme = "GWW_10_mwef"
)
corrN_GWW_10mwef<- data.frame(GWW_10 = rep(NA,20))
corrN_GWW_10mwef$GWW_10 <- result_N_GWW_mwef$ac
fwrite(corrN_GWW_10mwef, "./output/GWW_mwef_narea/corrN_GWW_10mwef.csv")

# ---------------------25%---------------------
GWWinv_mw <- readRDS(file = "./data/relmatrices/GWW/mw/25/GWWinv_mw.rds")
result_N_GWW_mwef<- crossv(sort = sort,
                             train = na.omit(nareablues_mw_correct),
                             validation = nareablues_ef_correct, 
                             mytrait = "narea", 
                             GINV = GWWinv_mw,
                             data = "GWW_mwef_narea",
                             scheme = "GWW_25_mwef"
)
corrN_GWW_25mwef<- data.frame(GWW_25 = rep(NA,20))
corrN_GWW_25mwef$GWW_25 <- result_N_GWW_mwef$ac
fwrite(corrN_GWW_25mwef, "./output/GWW_mwef_narea/corrN_GWW_25mwef.csv")

# ---------------------50%---------------------
GWWinv_mw <- readRDS(file = "./data/relmatrices/GWW/mw/50/GWWinv_mw.rds")
result_N_GWW_mwef<- crossv(sort = sort,
                             train = na.omit(nareablues_mw_correct),
                             validation = nareablues_ef_correct, 
                             mytrait = "narea", 
                             GINV = GWWinv_mw,
                             data = "GWW_mwef_narea",
                             scheme = "GWW_50_mwef"
)
corrN_GWW_50mwef<- data.frame(GWW_50 = rep(NA,20))
corrN_GWW_50mwef$GWW_50 <- result_N_GWW_mwef$ac
fwrite(corrN_GWW_50mwef, "./output/GWW_mwef_narea/corrN_GWW_50mwef.csv") 

# -------cv for sla train mw val on ef using GWW matrix with 10%,25%, 50% scheme ---------------------
GWWinv_mw <- readRDS(file = "./data/relmatrices/GWW/mw/10/GWWinv_mw.rds")
result_S_GWW_mwef<- crossv(sort = sort,
                             train = na.omit(slablues_mw_correct),
                             validation = slablues_ef_correct, 
                             mytrait = "sla", 
                             GINV = GWWinv_mw,
                             data = "GWW_mwef_sla",
                             scheme = "GWW_10_mwef"
)
corrS_GWW_10mwef<- data.frame(GWW_10 = rep(NA,20))
corrS_GWW_10mwef$GWW_10 <- result_S_GWW_mwef$ac
fwrite(corrS_GWW_10mwef, "./output/GWW_mwef_sla/corrS_GWW_10mwef.csv")
# ---------------------25%---------------------
GWWinv_mw <- readRDS(file = "./data/relmatrices/GWW/mw/25/GWWinv_mw.rds")
result_S_GWW_mwef<- crossv(sort = sort,
                             train = na.omit(slablues_mw_correct),
                             validation = slablues_ef_correct, 
                             mytrait = "sla", 
                             GINV = GWWinv_mw,
                             data = "GWW_mwef_sla",
                             scheme = "GWW_25_mwef"
)
corrS_GWW_25mwef<- data.frame(GWW_25 = rep(NA,20))
corrS_GWW_25mwef$GWW_25 <- result_S_GWW_mwef$ac
fwrite(corrS_GWW_25mwef, "./output/GWW_mwef_sla/corrS_GWW_25mwef.csv")

# ---------------------50%---------------------
GWWinv_mw <- readRDS(file = "./data/relmatrices/GWW/mw/50/GWWinv_mw.rds")
result_S_GWW_mwef<- crossv(sort = sort,
                             train = na.omit(slablues_mw_correct),
                             validation = slablues_ef_correct, 
                             mytrait = "sla", 
                             GINV = GWWinv_mw,
                             data = "GWW_mwef_sla",
                             scheme = "GWW_50_mwef"
)
corrS_GWW_50mwef<- data.frame(GWW_50 = rep(NA,20))
corrS_GWW_50mwef$GWW_50 <- result_S_GWW_mwef$ac
fwrite(corrS_GWW_50mwef, "./output/GWW_mwef_sla/corrS_GWW_50mwef.csv")       

# ---------------------cv for sla train ef val on mw using GWW matrix with 10%,25%, 50% scheme ---------------------
GWWinv_ef <- readRDS(file = "./data/relmatrices/GWW/ef/10/GWWinv_ef.rds")
result_S_GWW_efmw<- crossv(sort = sort,
                             train = na.omit(slablues_ef_correct),
                             validation = slablues_mw_correct, 
                             mytrait = "sla", 
                             GINV = GWWinv_ef,
                             data = "GWW_efmw_sla",
                             scheme = "GWW_10_efmw"
)
corrS_GWW_10efmw<- data.frame(GWW_10 = rep(NA,20))
corrS_GWW_10efmw$GWW_10 <- result_S_GWW_efmw$ac
fwrite(corrS_GWW_10efmw, "./output/GWW_efmw_sla/corrS_GWW_10efmw.csv")

# ---------------------25%---------------------
GWWinv_ef <- readRDS(file = "./data/relmatrices/GWW/ef/25/GWWinv_ef.rds")
result_S_GWW_efmw<- crossv(sort = sort,
                             train = na.omit(slablues_ef_correct),
                             validation = slablues_mw_correct, 
                             mytrait = "sla", 
                             GINV = GWWinv_ef,
                             data = "GWW_efmw_sla",
                             scheme = "GWW_25_efmw"
)
corrS_GWW_25efmw<- data.frame(GWW_25 = rep(NA,20))
corrS_GWW_25efmw$GWW_25 <- result_S_GWW_efmw$ac
fwrite(corrS_GWW_25efmw, "./output/GWW_efmw_sla/corrS_GWW_25efmw.csv")
# ---------------------50%---------------------
GWWinv_ef <- readRDS(file = "./data/relmatrices/GWW/ef/50/GWWinv_ef.rds")
result_S_GWW_efmw<- crossv(sort = sort,
                             train = na.omit(slablues_ef_correct),
                             validation = slablues_mw_correct, 
                             mytrait = "sla", 
                             GINV = GWWinv_ef,
                             data = "GWW_efmw_sla",
                             scheme = "GWW_50_efmw"
)
corrS_GWW_50efmw<- data.frame(GWW_50 = rep(NA,20))
corrS_GWW_50efmw$GWW_50 <- result_S_GWW_efmw$ac
fwrite(corrS_GWW_50efmw, "./output/GWW_efmw_sla/corrS_GWW_50efmw.csv")     
