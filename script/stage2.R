library(data.table)
library(cvTools)
library(ASRgenomics)
library(asreml)
library(dplyr)
library(ggplot2)
source("./script/aux_functions.R")
#reading blues from first stage analysis for narea 
selected_lines<- read.csv("./data/selected_lines.csv")
bluesjoint_narea_correct<- fread("./output/traitsoutput/bluesjoint_narea_correct.csv",data.table = F) %>% 
  filter(bluesjoint_narea_correct$taxa %in% selected_lines$x) %>% 
  mutate(taxa = factor(taxa))%>% rename("narea"= predicted.value)
bluesjoint_sla_correct<- fread("./output/traitsoutput/bluesjoint_sla_correct.csv" , data.table = F)%>%
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
sort<- create_folds( individuals= bluesjoint_narea_correct$taxa, nfolds= 5, reps = 20, seed = 1)
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
# ---------------------SLA-running cross validation for joint loc---------------------
result_S<- crossv(sort = sort,
                  train = bluesjoint_sla_correct,
                  validation = bluesjoint_sla_correct, 
                  mytrait = "sla", 
                  GINV = GINV,
                  data = "GBLUP",
                  scheme= "joint")
corr_S_joint<- data.frame(result_S$ac)

# ---------------------SLA-running cross validation for ...---------------------
#SLA- CV with ef as a train and mw validation
result_sla_efmw<- crossv(sort = sort,
                  train = slablues_ef_correct,
                  validation = slablues_mw_correct, 
                  mytrait = "sla", 
                  GINV = GINV,
                  data = "GBLUP",
                  scheme = "efmw")
corr_S_efmw<- data.frame(result_sla_efmw$ac)
# ---------------------SLA-running cross validation for joint loc---------------------
#narea training the model withef and validating on nareabluesmw--------------------
result_narea_efmw<- crossv(sort = sort,
                         train = nareablues_ef_correct,
                         validation = nareablues_mw_correct, 
                         mytrait = "narea", 
                         GINV = GINV,
                         data = "GBLUP",
                         scheme = "efmw"
                         )
corr_N_efmw<- data.frame(result_narea_efmw$ac)
boxplot(corr_N_efmw$result_narea_efmw.ac)
#narea training the model with mw and validating on ef--------------------
result_narea_mwef<- crossv(sort = sort,
                           train = nareablues_mw_correct,
                           validation = nareablues_ef_correct, 
                           mytrait = "narea", 
                           GINV = GINV,
                           data = "GBLUP",
                           scheme = "mwef")
corr_N_mwef<- data.frame(result_narea_mwef$ac)
#SLA training the model with mw and validating on ef--------------------
result_sla_mwef<- crossv(sort = sort,
                           train = slablues_mw_correct,
                           validation = slablues_ef_correct, 
                           mytrait = "sla", 
                           kin = kin,
                         data = "GBLUP",
                         scheme = "mwef")

corr_S_mwef<- data.frame(result_sla_mwef$ac)

#NIRS model
sort<- create_folds( individuals= bluesjoint_narea_correct$taxa, nfolds= 5, reps = 20, seed = 135)
nirs_inv<- readRDS(file = "./data/relmatrices/nirs_inv.rds")
# fitting model for narea using train nareajoint and val on nareajoint using NIRS matrix--------------------
result_N_nirs <-
  crossv(
    sort = sort,
    train = bluesjoint_narea_correct,
    validation =  bluesjoint_narea_correct,
    mytrait = "narea",
    GINV = nirs_inv,
    data= "NIRS",
    scheme= "nirs_joint"
  )
corr_N_joint<- data.frame(result_N_nirs$ac)
#joint sla nirs
result_S_nirs <-
  crossv(
    sort = sort,
    train = bluesjoint_sla_correct,
    validation =  bluesjoint_sla_correct,
    mytrait = "sla",
    GINV = nirs_inv,
    data= "NIRS",
    scheme= "nirs_joint"
  )
corr_S_nirs_joint<- data.frame(result_S_nirs$ac)
sla_NIRS_joint<- fread("./output/NIRS/sla_joint.csv")
#narea training ef validating mw using NIRS ef
re_nirs_ef<- read.csv('./data/relmatrices/re_nirs_ef.csv')
rownames(re_nirs_ef) <- colnames(re_nirs_ef)
Gb <- G.tuneup(G = as.matrix(re_nirs_ef), bend = TRUE, eig.tol = 1e-06)$Gb
nirs_inv_ef<- G.inverse(G = Gb , sparseform = T)
result_narea_nirs_efmw<- crossv(sort = sort,
                           train = nareablues_ef_correct,
                           validation = nareablues_mw_correct, 
                           mytrait = "narea", 
                           GINV = nirs_inv_ef,
                           data = "NIRS",
                           scheme = "nirs_efmw"
)

#narea training mw validating ef using NIRSmw
nirs_inv_mw <- readRDS(file = "./data/relmatrices/nirs_inv_mw.rds")
result_narea_nirs_mwef<- crossv(sort = sort,
                                train = nareablues_mw_correct,
                                validation = nareablues_ef_correct, 
                                mytrait = "narea", 
                                GINV = nirs_inv_mw,
                                data = "NIRS",
                                scheme = "nirs_mwef"
)
corr_N_nirs_mwef<- data.frame(result_narea_nirs_mwef$ac)
fwrite(corr_N_nirs_mwef, "./output/NIRS/ corr_N_nirs_mwef.csv")
#sla training ef validating mw using nirs ef matrix
result_sla_nirs_efmw<- crossv(sort = sort,
                                train = slablues_ef_correct,
                                validation = slablues_mw_correct, 
                                mytrait = "sla", 
                                GINV = nirs_inv_ef,
                                data = "NIRS",
                                scheme = "nirs_efmw")
#sla training mw validating ef using nirs mw matrix
result_sla_nirs_mwef<- crossv(sort = sort,
                              train = slablues_mw_correct,
                              validation = slablues_ef_correct, 
                              mytrait = "sla", 
                              GINV = nirs_inv_mw,
                              data = "NIRS",
                              scheme = "nirs_mwef")
corr_sla_nirs_mwef<- data.frame(result_sla_nirs_mwef$ac)
fwrite(corr_sla_nirs_mwef, "./output/NIRS/ corr_sla_nirs_mwef.csv")

#whole wave model

# re_w_joint<- read.csv("./data/relmatrices/re_w_joint.csv")
# rownames(re_w_joint) <- colnames(re_w_joint)
# Gb <- G.tuneup(G = as.matrix(re_w_joint), bend = TRUE, eig.tol = 1e-06)$Gb
# re_w_joint_inv<- G.inverse(G = Gb , sparseform = T)
#saveRDS(re_w_joint_inv, file = "./data/relmatrices/re_w_joint_inv.rds")
re_w_joint_inv <- readRDS(file = "./data/relmatrices/re_w_joint_inv.rds")
sort<- create_folds( individuals= bluesjoint_narea_correct$taxa, nfolds= 5, reps = 20, seed = 129)
# -----narea-running cross validation for joint loc using whole wave information---------------------
result_N_wholewave<-
  crossv(
    sort = sort,
    train = na.omit(bluesjoint_narea_correct),
    validation =  bluesjoint_narea_correct,
    mytrait = "narea",
    GINV = re_w_joint_inv,
    data= "wholewave",
    scheme= "wholewave_joint"
  )
corr_N_wholewave_joint<- data.frame(result_N_wholewave$ac)
fwrite(corr_N_wholewave_joint, "./output/wholewave/corr_N_wholewave_joint.csv")

# -----SLA-running cross validation for joint loc using whole wave info---------------------
result_S_wholewave<- crossv(sort = sort,
                  train = bluesjoint_sla_correct,
                  validation = bluesjoint_sla_correct, 
                  mytrait = "sla", 
                  GINV = re_w_joint_inv,
                  data = "wholewave",
                  scheme= "wholewave_joint")
corr_S_wholewave_joint<- data.frame(result_S_wholewave$ac)
fwrite(corr_S_wholewave_joint, "./output/wholewave/corr_S_wholewave_joint.csv")

#------narea training ef validating mw using whole wave ef
# re_w_ef<- read.csv('./data/relmatrices/re_w_ef.csv')
# rownames(re_w_ef) <- colnames(re_w_ef)
# Gb <- G.tuneup(G = as.matrix(re_w_ef), bend = TRUE, eig.tol = 1e-06)$Gb
# re_w_ef_inv<- G.inverse(G = Gb , sparseform = T)
# saveRDS(re_w_ef_inv, file = "./data/relmatrices/re_w_ef_inv.rds")
re_w_ef_inv <- readRDS(file = "./data/relmatrices/re_w_ef_inv.rds")
result_narea_wholewave_efmw<- crossv(sort = sort,
                                train = nareablues_ef_correct,
                                validation = nareablues_mw_correct, 
                                mytrait = "narea", 
                                GINV = re_w_ef_inv,
                                data = "wholewave",
                                scheme = "wholewave_efmw"
)
corr_N_wholewave_efmw<- data.frame(result_narea_wholewave_efmw$ac)
fwrite(corr_N_wholewave_efmw, "./output/wholewave/corr_N_wholewave_efmw.csv")
#narea training mw validating ef using wholewave mw
# re_w_mw<- read.csv('./data/relmatrices/re_w_MW.csv')
# rownames(re_w_mw) <- colnames(re_w_mw)
# Gb <- G.tuneup(G = as.matrix(re_w_mw), bend = TRUE, eig.tol = 1e-06)$Gb
# re_w_mw_inv<- G.inverse(G = Gb , sparseform = T)
# saveRDS(re_w_mw_inv, file = "./data/relmatrices/re_w_mw_inv.rds")
re_w_mw_inv <- readRDS(file = "./data/relmatrices/re_w_mw_inv.rds")
result_narea_wholewave_mwef<- crossv(sort = sort,
                                     train = nareablues_mw_correct,
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
                                     train = slablues_ef_correct,
                                     validation = slablues_mw_correct, 
                                     mytrait = "sla", 
                                     GINV = re_w_ef_inv,
                                     data = "wholewave",
                                     scheme = "wholewave_efmw"
)
corr_S_wholewave_efmw<- data.frame(result_sla_wholewave_efmw$ac)
#fwrite(corr_S_wholewave_efmw, "./output/wholewave/corr_S_wholewave_efmw.csv")


#sla training mw validating ef using wholewave mw
re_w_mw_inv <- readRDS(file = "./data/relmatrices/re_w_mw_inv.rds")
result_sla_wholewave_mwef<- crossv(sort = sort,
                                     train = slablues_mw_correct,
                                     validation = slablues_ef_correct, 
                                     mytrait = "sla", 
                                     GINV = re_w_mw_inv,
                                     data = "wholewave",
                                     scheme = "wholewave_mwef"
)
corr_S_wholewave_mwef<- data.frame(result_sla_wholewave_mwef$ac)
fwrite(corr_S_wholewave_mwef, "./output/wholewave/corr_S_wholewave_mwef.csv")

#highly heritable models using high heritable relationship matrix for joint loc
# # re_h2_joint<- read.csv("./data/relmatrices/re_h2_joint.csv")
# # rownames(re_h2_joint) <- colnames(re_h2_joint)
# # Gb <- G.tuneup(G = as.matrix(re_h2_joint), bend = TRUE, eig.tol = 1e-06)$Gb
# re_h2_joint_inv<- G.inverse(G = Gb , sparseform = T)
# saveRDS(re_h2_joint_inv, file = "./data/relmatrices/re_h2_joint_inv.rds")
re_h2_joint_inv <- readRDS(file = "./data/relmatrices/re_h2_joint_inv.rds")
sort<- create_folds( individuals= bluesjoint_narea_correct$taxa, nfolds= 5, reps = 20, seed = 147)
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
#fwrite(corr_N_highh2_joint, "./output/highh2/corr_N_highh2_joint.csv")

# -----SLA-running cv for joint loc using highh2 matrix--------------------
result_S_highh2_joint<- crossv(sort = sort,
                            train = bluesjoint_sla_correct,
                            validation = bluesjoint_sla_correct, 
                            mytrait = "sla", 
                            GINV = re_h2_joint_inv,
                            data = "highh2",
                            scheme= "highh2_joint")
corr_S_highh2_joint<- data.frame(result_S_highh2_joint$ac)
fwrite(corr_S_highh2_joint, "./output/highh2/corr_S_highh2_joint.csv")

#------narea training ef validating mw using high h2 ef loc matrix
# re_h2_ef<-read.csv("./data/relmatrices/re_h2_ef.csv")
# rownames(re_h2_ef) <- colnames(re_h2_ef)
# Gb <- G.tuneup(G = as.matrix(re_h2_ef), bend = TRUE, eig.tol = 1e-06)$Gb
# re_h2_ef_inv<- G.inverse(G = Gb , sparseform = T)
# saveRDS(re_h2_ef_inv, file = "./data/relmatrices/re_h2_ef_inv.rds")
re_h2_ef_inv <- readRDS(file = "./data/relmatrices/re_h2_ef_inv.rds")
result_narea_highh2_efmw<- crossv(sort = sort,
                                     train = nareablues_ef_correct,
                                     validation = nareablues_mw_correct, 
                                     mytrait = "narea", 
                                     GINV = re_h2_ef_inv,
                                     data = "highh2",
                                     scheme = "highh2_efmw"
)
corr_N_highh2_efmw<- data.frame(result_narea_highh2_efmw$ac)
fwrite(corr_N_highh2_efmw, "./output/highh2/corr_N_highh2_efmw.csv")
#narea training mw validating ef using highh2 mw matrix
# re_h2_mw<-read.csv("./data/relmatrices/re_h2_mw.csv")
# rownames(re_h2_mw) <- colnames(re_h2_mw)
# Gb <- G.tuneup(G = as.matrix(re_h2_mw), bend = TRUE, eig.tol = 1e-06)$Gb
# re_h2_mw_inv<- G.inverse(G = Gb , sparseform = T)
# saveRDS(re_h2_mw_inv, file = "./data/relmatrices/re_h2_mw_inv.rds")
re_h2_mw_inv <- readRDS(file = "./data/relmatrices/re_h2_mw_inv.rds")
result_narea_h2_mwef<- crossv(sort = sort,
                                     train = nareablues_mw_correct,
                                     validation = nareablues_ef_correct, 
                                     mytrait = "narea", 
                                     GINV = re_h2_mw_inv,
                                     data = "highh2",
                                     scheme = "highh2_mwef"
)
corr_N_highh2_mwef<- data.frame(result_narea_h2_mwef$ac)
#fwrite(corr_N_highh2_mwef, "./output/highh2/corr_N_highh2_mwef.csv")

#sla training ef and validating mw using highh2 ef matrix
re_h2_ef_inv <- readRDS(file = "./data/relmatrices/re_h2_ef_inv.rds")
result_sla_highh2_efmw<- crossv(sort = sort,
                                   train = slablues_ef_correct,
                                   validation = slablues_mw_correct, 
                                   mytrait = "sla", 
                                   GINV = re_h2_ef_inv,
                                   data = "highh2",
                                   scheme = "highh2_efmw"
)
corr_S_highh2_efmw<- data.frame(result_sla_highh2_efmw$ac)
#fwrite(corr_S_highh2_efmw, "./output/highh2/corr_S_highh2_efmw.csv")


#sla training mw validating ef using highh2 mw
re_h2_mw_inv <- readRDS(file = "./data/relmatrices/re_h2_mw_inv.rds")
result_sla_highh2_mwef<- crossv(sort = sort,
                                   train = slablues_mw_correct,
                                   validation = slablues_ef_correct, 
                                   mytrait = "sla", 
                                   GINV = re_h2_mw_inv,
                                   data = "highh2",
                                   scheme = "highh2_mwef"
)
corr_S_highh2_mwef<- data.frame(result_sla_highh2_mwef$ac)
#fwrite(corr_S_highh2_mwef, "./output/highh2/corr_S_highh2_mwef.csv")



