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
sort<- create_folds( individuals= bluesjoint_narea_correct$taxa, nfolds= 5, reps =20)

# ---------------------narea-running cross validation for joint loc---------------------
result_N <-
  crossv(
    sort = sort,
    train = bluesjoint_narea_correct,
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
colnames(gebv_S_efmw) <- c("taxa", "predicted.value", "sla", "repetition")
gebv_S_efmw<- read.csv("./GBLUP_output/gebv_sla/trainef_valmw/gebv_S_efmw.csv")
corr_S_efmw<- read.csv("./GBLUP_output/gebv_sla/trainef_valmw/corr_S_efmw.csv")
boxplot(corr_S_efmw$X)
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
boxplot(corr_N_efmw$X)
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
                           GINV = GINV,
                         data = "GBLUP",
                         scheme = "mwef")

corr_S_mwef<- data.frame(result_sla_mwef$ac)

# ---------------------fitting model for trait narea using train set from wavebluesef and validation on waveblues mw---------------------
wavebluesEF_wider_corrected<- fread("./output/wavebluesEF_wider_corrected.csv")%>% left_join(nareabluesef,by=c("taxa"))
wavebluesEF_wider_corrected<- transform(wavebluesEF_wider_corrected,
                                        taxa= factor(taxa))
#for validation set from mw location
wavebluesMW_wider_corrected<- fread("./output/wavebluesMW_wider_corrected.csv")%>% left_join(nareabluesmw,by=c("taxa"))
wavebluesMW_wider_corrected<- transform(wavebluesMW_wider_corrected,
                                        taxa= factor(taxa))

# ---------------------cross validation set---------------------
sort<- create_folds(individuals = wavebluesEF_wider_corrected$taxa, nfolds= 5, reps= 20, seed= 148)
result<- crossv(sort = sort,train = wavebluesEF_wider_corrected, validation = wavebluesMW_wider_corrected, mytrait = "narea", GINV = GINV)

# ---------------------fitting model and cross validation using wavebluesmw and wavebluesef for narea---------------------
wavebluesEF_wider_corrected<- fread("./output/wavebluesEF_wider_corrected.csv")%>% left_join(nareabluesef,by=c("taxa"))
wavebluesEF_wider_corrected<- transform(wavebluesEF_wider_corrected,
                                        taxa= factor(taxa))
#for validation set from mw location
wavebluesMW_wider_corrected<- fread("./output/wavebluesMW_wider_corrected.csv")%>% left_join(nareabluesmw,by=c("taxa"))
wavebluesMW_wider_corrected<- transform(wavebluesMW_wider_corrected,
                                        taxa= factor(taxa))
# ---------------------cross validation set---------------------
sort<- create_folds(individuals=wavebluesMW_wider_corrected$taxa, nfolds= 5, reps= 20, seed= 157)
result<- crossv(sort = sort,train = slabluesef, validation = slabluesmw, mytrait = "sla", GINV = GINV)
# --------------------fitting model---------------------
acN <- c()
gebv_mwvsef <- list()
for(j in 1:length(sort)){
  rN<-list()
  for(i in 1:5){
    test<-wavebluesmw_wider_corrected
    test[test$taxa %in% sort[[j]][[i]],"narea"] <- NA
    cat(i, j,"\n")
    # ---------------------fitting model---------------------
    
    modelN <- asreml(
      fixed = narea ~ 1,
      random =  ~  vm(taxa, source = GINV$Ginv.sparse),
      data = wavebluesMW_wider_corrected, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "taxa"))
    cat("\n")
    
    # Update if necessary
    if(!modelN$converge){ modelN <- update.asreml(modelN) }
    if(!modelN$converge){ modelN <- update.asreml(modelN) }
    
    rN[[i]] <- modelN$predictions$pvals[modelN$predictions$pvals$taxa %in% sort[[j]][[i]],1:2]
  }
  
  raN<- Reduce(rbind, rN)
  gebv_efvsmw[[j]] <- raN
  raN <- raN %>% left_join(wavebluesef_wider_corrected[,c("taxa", "narea")])
  acN[j]<-cor(raN[,2], raN[,3], use = "complete.obs")
}
fwrite( raN, "./output/gebv_mwvsef.csv", sep = "\t", col.names = T)
gebv_mwvsef<- fread("./output/gebv_mwvsef.csv")
cor<- cor.test(gebv_mwvsef$predicted.value, gebv_mwvsef$narea)

#fitting second stage for sla using wavebluesEF as a train set and wavebluesMW as a validation set 
# ---------------------fitting model and cross validation using wavebluesef and waveblues mw for sla trait of ineterst---------------------
wavebluesEF_wider_corrected<- fread("./output/wavebluesEF_wider_corrected.csv")%>% left_join(slabluesef,by=c("taxa"))
wavebluesEF_wider_corrected<- transform(wavebluesEF_wider_corrected,
                                        taxa= factor(taxa))
