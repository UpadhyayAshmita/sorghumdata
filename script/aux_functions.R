# ---------------------cross validation function---------------------
crossv <- function(sort,
                   train,
                   validation,
                   GINV,
                   mytrait,
                   data,
                   scheme
                   ) {
  library(asreml)
  library(dplyr)
  library(data.table)
  ac <- c()
  gebv <- list()
  
  for(j in 1:length(sort)){
    r<-list()
    for(i in 1:5){
      test<-train
      test[test$taxa %in% sort[[j]][[i]], mytrait ] <- NA
      cat(i, j,"\n")
      # ---------------------fitting model---------------------
      
      model <- asreml(
        fixed = get(mytrait) ~ 1, 
        random =  ~ vm(taxa, source = GINV$Ginv.sparse),
        data = test, na.action = na.method(x = "include"),
        predict = predict.asreml(classify = "taxa"))
      
      # Update if necessary
      if(!model$converge){ model <- update.asreml(model) }
      if(!model$converge){ model <- update.asreml(model) }
      
      r[[i]] <- model$predictions$pvals[model$predictions$pvals$taxa %in% sort[[j]][[i]],1:2]
    }
    
    gebv[[j]] <- Reduce(rbind, r)
    gebv[[j]] <- gebv[[j]] %>%
      left_join(validation[,c("taxa", mytrait)]) %>%
      mutate(rep = j)
    ac[j]<-cor(gebv[[j]][,2], gebv[[j]][,3], use = "complete.obs")
  }
  gebv<- bind_rows(gebv, .id = "rep")
 write.csv(gebv, paste0("./output/", data, "/",  mytrait, "_", scheme, ".csv"),row.names=F)
 result<- list(gebv= gebv,ac= ac)
 return(result)
}

#relationship matrix function
calculate_relationship <-
  function(data) {
    data <- data %>% 
      tibble::column_to_rownames("taxa") %>% 
      as.matrix()
    W <- data %*% t(data) / (ncol(data) - 1)
    return(W)
  }

#creating sort
# ---------------------creating list with 5 fold and 20 reps----------------------
create_folds<- function(individuals, nfolds, reps, seed = 123){
  library(cvTools)
  library(dplyr)
  set.seed(seed)
  sort<-list()
  individuals <- as.factor(individuals)
  nl <- length(unique(individuals))
  for(a in 1:reps){
  folds <- cvFolds(nl,type ="random", K= nfolds)
  Sample<-cbind(folds$which,folds$subsets)
  cv<-split(levels(individuals)[Sample[,2]], f=Sample[,1])#a list of subsets (folds) of taxa
  sort[[a]]<-cv
 }
 return(sort)
}
