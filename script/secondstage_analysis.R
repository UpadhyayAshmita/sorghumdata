library(data.table)
library(cvTools)
library(bestNormalize)
library(fs)
library(ASRgenomics)
library(asreml)
library(dplyr)
#-----reading the blues data file-----
blues<- fread("./output/blues.txt")

#reading blues from first stage analysis for Narea 
# blues for Narea trait
Nitroblues <- subset(blues, select = c("TAXA", "Narea"))

# blues for SLA trait
SLAblues <- subset(blues, select = c("TAXA", "SLA"))
SLAblues<- transform(SLAblues,
                       TAXA= factor(TAXA))

nitroblues<- fread("./output/nitroblues.csv")
nitroblues<- transform(nitroblues,
                       TAXA= factor(TAXA))

#second stage analysis 

#### creating folders 5-fold cross-validation ####
set.seed(134)
sort<-list()
nl <- length(unique(nitroblues$TAXA))
for(a in 1:20){
  folds <- cvFolds(nl,type ="random", K=5, R = 1)
  Sample<-cbind(folds$which,folds$subsets)
  cv<-split(levels(nitroblues$TAXA)[Sample[,2]], f=Sample[,1])#a list of subsets (folds) of TAXA
  sort[[a]]<-cv
}
rm(a, folds, cv, Sample)

#----- Second stage for trait Narea --------
acN <- c()
gebvN <- list()

for(j in 1:length(sort)){
  rN<-list()
  for(i in 1:5){
    test<-nitroblues
    test[test$TAXA %in% sort[[j]][[i]],"Narea"] <- NA
    cat(i, j,"\n")
    modelN <- asreml(
      fixed = Narea ~ 1,
      random =  ~  vm(TAXA, source = GINV$Ginv.sparse),
      data = test, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "TAXA"))
    cat("\n")
    
    # Update if necessary
    if(!modelN$converge){ modelN <- update.asreml(modelN) }
    if(!modelN$converge){ modelN <- update.asreml(modelN) }
    
    rN[[i]] <- modelN$predictions$pvals[modelN$predictions$pvals$TAXA %in% sort[[j]][[i]],1:2]
  }
  
  raN<- Reduce(rbind, rN)
  gebvN[[j]] <- raN
  fwrite(raN, "gebvN.txt", sep = "\t", append = T, col.names = F)
  raN <- raN %>% left_join(nitroblues[,c("TAXA", "Narea")])
  acN[j]<-cor(raN[,2], raN[,3], use = "complete.obs")
}
fwrite(data.frame("TAXA","GEBV"), "gebvN.txt", sep = "\t", col.names = F)
gebv<- fread("./output/gebvN.txt")

#cross validation for SLA
set.seed(135)
sort<-list()
nl <- length(unique(SLAblues$TAXA))
for(a in 1:20){
  folds <- cvFolds(nl,type ="random", K=5, R = 1)
  Sample<-cbind(folds$which,folds$subsets)
  cv<-split(levels(SLAblues$TAXA)[Sample[,2]], f=Sample[,1])#a list of subsets (folds) of TAXA
  sort[[a]]<-cv
}
rm(a, folds, cv, Sample)

#----- second stage for trait SLA --------
acS <- c()
gebvSLA <- list()

for(j in 1:length(sort)){
  rS<-list()
  for(i in 1:5){
    testsla<-SLAblues               
    testsla[testsla$TAXA %in% sort[[j]][[i]],"SLA"] <- NA
    cat(i, j,"\n")
    modelSLA <- asreml(
      fixed = SLA ~ 1,
      random =  ~  vm(TAXA, source = GINV$Ginv.sparse),
      data = testsla, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "TAXA"))
    cat("\n")
    
    # Update if necessary
    if(!modelSLA$converge){ modelSLA <- update.asreml(modelSLA) }
    if(!modelSLA$converge){ modelSLA <- update.asreml(modelSLA) }
    
    rS[[i]] <- modelSLA$predictions$pvals[modelSLA$predictions$pvals$TAXA %in% sort[[j]][[i]],1:2]
  }
  
  raS<- Reduce(rbind, rS)
  gebvSLA[[j]] <- raS
  fwrite(raS, "gebvSLA.txt", sep = "\t", append = T, col.names = F)
  raS <- raS %>% left_join(SLAblues[,c("TAXA", "SLA")])
  acS[j]<-cor(raS[,2], raS[,3], use = "complete.obs")
}
       
