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


kin<- fread('./data/relmatrices/kinship_additive.txt', data.table= F)
re_nirs_joint<- read.csv("./data/relmatrices/re_nirs_joint.csv")
# Function to perform operations on kinship matrix and wholewave matrix
process_kinship <- function(kin_matrix, 
                            wave_matrix, 
                            remove_percent,
                            comb_matrix, 
                            seed=123,
                            scheme,
                            loc,
                            percent) {
  set.seed(seed)
  N <- nrow(kin_matrix)
  remove <- sample(1:N, size = round(N * remove_percent))
  wave_matrix <- kin_matrix
  (paste0(comb_matrix))[remove, ] <- wave_matrix[remove, ]
  comb_matrix[, remove] <- wave_matrix[, remove]
  #write comb_matrix matrix to file
  write.csv(comb_matrix, paste0("./data/relmatrices/", scheme, "/",loc, "/",percent,remove_percent * 100, ".csv"))
  # Compute inverse of comb matrix
  rownames(comb_matrix) <- colnames(comb_matrix)
  Gb <- G.tuneup(G = as.matrix(comb_matrix), bend = TRUE, eig.tol = 1e-06)$Gb
  comb_matrixinv <- G.inverse(G = Gb, sparseform = TRUE)
  # Save inverse of comb matrix to file
   saveRDS(comb_matrixinv,paste0("./data/relmatrices/", scheme, "/",loc, "/",percent,remove_percent * 100, ".rds"))
}
