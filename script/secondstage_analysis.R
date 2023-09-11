library(data.table)
library(cvTools)
library(bestNormalize)
library(fs)
library(ASRgenomics)
library(asreml)
library(dplyr)
library(janitor)
#-----reading the blues data file-----
blues<- fread("./output/blues.csv")

#reading blues from first stage analysis for narea 
# blues for narea trait
Nitroblues <- subset(blues, select = c("taxa", "narea"))

fwrite(Nitroblues, "./output/nitroblues.csv",row.names = FALSE)

# blues for sla trait
slablues <- subset(blues, select = c("taxa", "sla"))
slablues<- transform(slablues,
                       taxa= factor(taxa))

nitroblues<- fread("./output/nitroblues.csv")
nitroblues<- transform(nitroblues,
                       taxa= factor(taxa))

#second stage analysis 

#### creating folders 5-fold cross-validation ####
set.seed(134)
sort<-list()
nl <- length(unique(nitroblues$taxa))
for(a in 1:20){
  folds <- cvFolds(nl,type ="random", K=5, R = 1)
  Sample<-cbind(folds$which,folds$subsets)
  cv<-split(levels(nitroblues$taxa)[Sample[,2]], f=Sample[,1])#a list of subsets (folds) of taxa
  sort[[a]]<-cv
}
rm(a, folds, cv, Sample)

#----- Second stage for trait narea --------
acN <- c()
gebvN <- list()

for(j in 1:length(sort)){
  rN<-list()
  for(i in 1:5){
    test<-nitroblues
    test[test$taxa %in% sort[[j]][[i]],"narea"] <- NA
    cat(i, j,"\n")
    # ---------------------fitting model---------------------
    
    modelN <- asreml(
      fixed = narea ~ 1,
      random =  ~  vm(taxa, source = GINV$Ginv.sparse),
      data = test, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "taxa"))
    cat("\n")
    
    # Update if necessary
    if(!modelN$converge){ modelN <- update.asreml(modelN) }
    if(!modelN$converge){ modelN <- update.asreml(modelN) }
    
    rN[[i]] <- modelN$predictions$pvals[modelN$predictions$pvals$taxa %in% sort[[j]][[i]],1:2]
  }
  
  raN<- Reduce(rbind, rN)
  gebvN[[j]] <- raN
  raN <- raN %>% left_join(nitroblues[,c("taxa", "narea")])
  acN[j]<-cor(raN[,2], raN[,3], use = "complete.obs")
}
fwrite( raN, "./output/gebvN.csv", sep = "\t", col.names = T)
gebvN<- fread("./output/gebvN.csv")

#cross validation for sla
set.seed(135)
sort<-list()
nl <- length(unique(slablues$taxa))
for(a in 1:20){
  folds <- cvFolds(nl,type ="random", K=5, R = 1)
  Sample<-cbind(folds$which,folds$subsets)
  cv<-split(levels(slablues$taxa)[Sample[,2]], f=Sample[,1])#a list of subsets (folds) of taxa
  sort[[a]]<-cv
}
rm(a, folds, cv, Sample)

#----- second stage for trait sla --------
acS <- c()
gebvsla <- list()

for(j in 1:length(sort)){
  rS<-list()
  for(i in 1:5){
    testsla<-slablues               
    testsla[testsla$taxa %in% sort[[j]][[i]],"sla"] <- NA
    cat(i, j,"\n")
    
    modelsla <- asreml(
      fixed = sla ~ 1,
      random =  ~  vm(taxa, source = GINV$Ginv.sparse),
      data = testsla, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "taxa"))
    cat("\n")
    
    # Update if necessary
    if(!modelsla$converge){ modelsla <- update.asreml(modelsla) }
    if(!modelsla$converge){ modelsla <- update.asreml(modelsla) }
    
    rS[[i]] <- modelsla$predictions$pvals[modelsla$predictions$pvals$taxa %in% sort[[j]][[i]],1:2]
  }
  
  raS<- Reduce(rbind, rS)
  gebvsla[[j]] <- raS
  raS <- raS %>% left_join(slablues[,c("taxa", "sla")])
  fwrite(raS, "./output/gebvsla.csv", sep = "\t", append = T, col.names = F)
  acS[j]<-cor(raS[,2], raS[,3], use = "complete.obs")
}
       
gebvsla<- fread("./output/gebvsla.csv")
