# ---------------------loading package--------------------

library(tidyverse)
library(asreml)
library(data.table)
library(ASRgenomics)
library(dplyr)
library(janitor)

asreml.options(
  workspace = '8gb',
  pworkspace = '8gb'
)

# ---------------------loading data---------------------
design<- fread('./data/design.csv', data.table = F)

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
variables <- colnames(design)[10:13]
models <- vector("list",length(variables))
waveblues <-data.frame()

# ---------------------fitting model---------------------

for (i in 1:length(variables)) {
  cat(variables[i], '\n')
  tryCatch({
    model <- asreml(
      fixed = get(variables[i]) ~loc + set + loc:set + name2 + name2:loc,
      random = ~loc:block,
      residual =  ~id(loc):ar1(range):ar1(row),
      data = design, na.action = na.method(x = "include"),
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
      waveblues <- bind_rows(waveblues, temp)
      cat('\n')
      
    }  else {
      waveblues <- bind_rows(waveblues, data.frame(name2 = NA,
                                                       predicted.value = NA,
                                                       wave = variables[i]))
      
    
      cat('\n')
    }
  },
  
  error = function(err) {
    message("An error occured")
    print(err)
  })

  rm(model)
  gc()
}


fwrite( waveblues, "./output/waveblues.csv", row.names = FALSE)
#waveblues<- fread("./output/waveblues.csv")
