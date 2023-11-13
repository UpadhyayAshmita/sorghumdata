# ---------------------loading package--------------------

library(tidyverse)
library(asreml)
library(data.table)
library(ASRgenomics)
library(dplyr)
library(janitor)

# ---------------------loading data---------------------
design<- fread('./data/design.csv',data.table= FALSE)

# ---------------------processing of data---------------------

design<-
  design |> clean_names() |> mutate(
    name2 = factor(name2),
    taxa = factor(taxa),
    loc = factor(loc),
    set = factor(set),
    block = factor(block),
    range = factor(range),
    row = factor(row),
    uni = c(1:960, 1:960)
  ) |>   arrange(loc, range, row)


#calculating blues for each wavelength 
variables <- colnames(design)[10:2160]
models <- vector("list",length(variables))
wavebluesEF <-data.frame()

# ---------------------fitting model---------------------

for (i in 1:length(variables)) {
  cat(variables[i], '\n')
  tryCatch({
    model <- asreml(
      fixed = get(variables[i]) ~set + name2,
      random = ~block,
      residual =  ~ar1(range):ar1(row),
      data = design,subset= loc== "EF", na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "name2", sed = TRUE)
    )
    
    if (!model$converge) {
      
      model <- update.asreml(model)
      
    }
    
    if (model$converge) {
      
      models[[i]] <- model
      
      #---------------------storing prediction---------------------
      temp <- models[[i]]$predictions$pvals[, 1:2]
      temp$wave <- variables[i]
      wavebluesEF <- bind_rows(wavebluesEF, temp)
      cat('\n')
      
    }  else {
      wavebluesEF <- bind_rows(wavebluesEF, data.frame(name2 = NA,
                                                       predicted.value = NA,
                                                       wave = variables[i]))
      cat('\n')
    }
    
    
  },
  
  error = function(err) {
    message("An error occured")
    print(err)
  })
  
}


fwrite( wavebluesEF, "./output/wavebluesEF.csv", row.names = FALSE)




