# ---------------------library---------------------
library(tidyverse)
library(asreml)
library(data.table)
library(ASRgenomics)
library(janitor)
source('./script/outlier.R')


asreml.options(
  workspace = '8gb',
  pworkspace = '8gb'
)


# reading arguments from command line
args <- commandArgs(trailingOnly = TRUE)
loc <- args[1]


# ---------------------load data---------------------
designnew<- fread('./data/designnew.csv',data.table= FALSE)


# filter data if loc is not NA
if (!is.na(loc)) {
  
  # subset for the model
  subs <- designnew$LOC == loc
  
  # formulas for one-location model
  fixed_terms <- c("name2 + set")
  random_formula <- as.formula("~ block")
  residual_formula <- as.formula("~ ar1(range):ar1(row)")
  
  # output file
  output_file <- paste0("./output/waveblues", loc, ".csv")
  
  # formulas for all locations model
} else {
  
  # subset for the model (all locations)
  subs <- designnew$LOC %in% c('EF', 'MW')
  
  fixed_terms <- c("loc + set + loc:set + name2 + name2:loc")
  random_formula <- as.formula("~ loc:block")
  residual_formula <- as.formula("~ id(loc):ar1(range):ar1(row)")
  
  # output file
  output_file <- "./output/waveblues.csv"
  
}


# ---------------------processing data---------------------
designnew<-
  designnew %>% clean_names() %>% mutate(
    name2 = factor(name2),
    taxa = factor(taxa),
    loc = factor(loc),
    set = factor(set),
    block = factor(block),
    range = factor(range),
    row = factor(row),
    uni = c(1:960, 1:960)
  ) %>%   arrange(loc, range, row)


#calculating blues for each wavelength 

models <- list()
variables <- colnames(designnew)[11:2161]
waveblues <- data.frame()

for(i in 1:length(variables)) {
  
  # ---------------------fitting model---------------------
  cat(variables[i], '\n')
  
  fixed_formula <- as.formula(paste(c(variables[i], "~", fixed_terms), collapse = " "))
  
  models[[i]] <- asreml(
    fixed = fixed_formula,
    random = random_formula,
    residual = residual_formula,
    data = designnew, 
    subset = subs,
    na.action = na.method(x = "include"),
    predict = predict.asreml(classify = "name2", sed = TRUE)
  )

  # ---------------------check outlier---------------------
  out <- outlier(models[[i]]$residuals)
  
  # ---------------------fit the model without outlier---------------------
  if (length(out) > 0) {
    designnew[out, variables[i]] <- NA
    models[[i]] <- asreml(
      fixed = fixed_formula,
      random = random_formula,
      residual = residual_formula,
      data = designnew, 
      subset = subs,
      na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "name2", sed = TRUE)
    )
  }

  # ---------------------storing prediction---------------------
  temp <- models[[i]]$predictions$pvals[, 1:2]
  temp$wave <- variables[i]
  waveblues <- rbind(waveblues, temp)
  cat('\n')
  
}

cat('Writing file...', '\n')
fwrite(waveblues, output_file, row.names= FALSE)
cat('Done!')
