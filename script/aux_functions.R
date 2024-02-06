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
      test<-na.omit(train)
      test[test$taxa %in% sort[[j]][[i]], mytrait ] <- NA
      cat(i, j,"\n")
      # ---------------------fitting model---------------------
      
      model <- asreml(
        fixed = get(mytrait) ~ 1, 
        random =  ~ vm(taxa,source = GINV$Ginv.sparse),#source= kin, singG= "NSD")
        data = test, na.action = na.method(x = "include"),
        family = asr_gaussian (dispersion = 1), #those are the important things
        weights = W, #those are the important things
        maxit = 50,
        trace = FALSE,
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

#function for combine relationship matrix
library(data.table)
# Function to perform operations on kinship matrix and wholewave matrix
process_kinship <- function(kin_matrix = kin,
                            wave_matrix = re_nirs_joint,
                            remove_percent,
                            seed=123,
                            rep,
                            type= NULL) {
  set.seed(seed)
  N <- nrow(kin_matrix)
  remove <- sample(1:N, size = round(N * remove_percent))
  greduced <- kin_matrix[-remove, -remove]
  kin_matrix[remove, remove] <- wave_matrix[remove, remove]
  if (is.null(type)) {
    scheme <- deparse(substitute(wave_matrix))
    scheme <- unlist(strsplit(scheme, "_"))
    loc <- scheme[3]
    type <- scheme[2]
  }
  write.csv(kin_matrix, paste0("./data/relmatrices/","G",type, "/",loc , "/",remove_percent * 100,"/", "G", type,"_",rep, ".csv"))
  rownames(kin_matrix) <- colnames(kin_matrix)
  Gb <- G.tuneup(G = as.matrix(kin_matrix), bend = TRUE, eig.tol = 1e-06)$Gb
  comb_matrixinv <- G.inverse(G = Gb, sparseform = TRUE)
   saveRDS(comb_matrixinv,paste0("./data/relmatrices/","G",type, "/",loc , "/",remove_percent * 100,"/", "G", type,"_",rep,  ".rds"))
    return(process_kinship)
   }

#function to calculate NRMSE
calculate_nrmse <- function(observed, predicted) {
  data <- data.frame(observed = observed, predicted = predicted)
  #data$observed[is.na(data$observed)] <- data$predicted[is.na(data$observed)]
  data<- na.omit(data)
  rmse <- sqrt(mean((data$predicted - data$observed)^2))
  range_response <- max(data$observed) - min(data$observed)
  nrmse <- rmse / range_response
  return(nrmse)
}
