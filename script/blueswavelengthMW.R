# ---------------------loading package--------------------

library(tidyverse)
library(asreml)
library(data.table)
library(ASRgenomics)
library(dplyr)
library(janitor)

#calculating blues for MW location for each wavelength 
asreml.options(
  workspace = '8gb',
  pworkspace = '8gb'
)


phenotypic_data_joint<- fread('./data/phenotypic_data_joint.csv',data.table= FALSE)
# ---------------------processing of data---------------------
phenotypic_data_joint<-
  phenotypic_data_joint |> clean_names() |> mutate(
    name2 = factor(name2),
    taxa = factor(taxa),
    loc = factor(loc),
    set = factor(set),
    block = factor(block),
    range = factor(range),
    row = factor(row),
  ) |>   arrange(loc, range, row)

#calculating blues for each wavelength 
variables <- colnames(phenotypic_data_joint)[12:2162]
models <- vector("list",length(variables))
wavebluesMW <-data.frame()
# ---------------------fitting model---------------------
for (i in 1:length(variables)) {
  cat(variables[i], '\n')
  tryCatch({
    model <- asreml(
      fixed = get(variables[i]) ~set + name2,
      random = ~block,
      residual =  ~ar1(range):ar1(row),
      data = phenotypic_data_joint,subset= loc== "MW", na.action = na.method(x = "include"),
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
      wavebluesMW <- bind_rows(wavebluesMW, temp)
      cat('\n')
      
    }  else {
      wavebluesMW <- bind_rows(wavebluesMW, data.frame(name2 = NA,
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
fwrite( wavebluesMW, "./output/wavebluesMW.csv", row.names = FALSE)
