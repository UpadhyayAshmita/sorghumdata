# ---------------------loading library---------------------
library(tidyverse)
library(asreml)
library(data.table)
library(ASRgenomics)
library(janitor)
source('./script/outlier.R')

# ---------------------loading data---------------------

designnew<- fread('./data/designnew.csv',data.table= FALSE)

# ---------------------processing of data---------------------

designnew<-
  designnew %>% clean_names()%>% mutate(
    name2 = factor(name2),
    taxa = factor(taxa),
    loc = factor(loc),
    set = factor(set),
    block = factor(block),
    range = factor(range),
    row = factor(row),
    uni = c(1:960, 1:960)
  ) %>%   arrange(loc, range, row)

# ---------------------calculating heritability---------------------

models <- list()
variables <- colnames(designnew)[11:2161]
h2 <- data.frame(matrix(NA, length(variables), 2))  #n rows and 2 columns
h2[,1] <- variables

colnames(h2) <- c("wave", "h2")


#for(i in variables)
  for (i in 1:length(variables)) {
    wavelength<- variables[i]
    
  
  # ---------------------fitting model---------------------

  cat( wavelength, '\n')
  
  models[[wavelength]] <- asreml(
  fixed = get(wavelength) ~ set + loc + loc: set,
  random = ~name2 + loc:block + loc:name2,
  residual =  ~id(loc):ar1(range):ar1(row),
  data = designnew, na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)

  # ---------------------checking outliers---------------------
  
out <- outlier(models[[wavelength]]$residuals)

if(length(out) > 0) {
  designnew[,wavelength][designnew$LOC][out] <- NA
    models[[wavelength]] <- asreml(
      fixed = get(wavelength) ~ set + loc + loc:set,
      random =  ~ name2 + loc:block + loc:name2,
      residual =  ~ id(loc):ar1(range):ar1(row),
      data = designnew, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "name2", sed = TRUE)
      
  )
}

# ---------------------calculating heritability from formula---------------------

h2[i,2]  <- (1 - ((models[[wavelength]]$predictions$avsed["mean"] ^ 2) /(2 * summary(models[[wavelength]])$varcomp["name2", "component"])))
}

cat('\n')

fwrite(h2, "./output/h2.csv",row.names = F)
h2<- fread("./output/h2.csv", data.table = FALSE)

cat('done!')







