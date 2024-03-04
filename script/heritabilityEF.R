# ---------------------loading library---------------------
library(tidyverse)
library(asreml)
library(data.table)
library(ASRgenomics)
library(janitor)
source('./script/outlier.R')

# ---------------------load data---------------------
phenotypic_data<- fread('./data/phenotypic_data.csv',data.table= FALSE)
# ---------------------processing of data---------------------
phenotypic_data<-
  phenotypic_data %>% mutate(
    name2 = factor(name2),
    taxa = factor(taxa),
    loc = factor(loc),
    set = factor(set),
    block = factor(block),
    range = factor(range),
    row = factor(row)
  ) %>%   arrange(loc, range, row)
# ---------------------fitting model--------------------
models <- list()
variables <- colnames(phenotypic_data)[12:2162]
h2 <- data.frame(matrix(NA, length(variables), 2))  #n rows and 2 columns
h2[,1] <- variables
colnames(h2) <- c("wave", "h2")
for (i in 1:length(variables)) {
  wavelength<- variables[i]
  # ---------------------fitting model--------------------
  cat( wavelength, '\n')
  temp <- asreml(
    fixed = get(wavelength) ~ set,
    random = ~name2 + block,
    residual =  ~ar1(range):ar1(row),
    data = phenotypic_data, subset= loc== "EF",na.action = na.method(x = "include"),
    predict = predict.asreml(classify = "name2", sed = TRUE)
  )
  if (!temp$converge) {
    temp <- update.asreml(temp)
  }  
  if (!temp$converge) {
    models[[wavelength]] <- NA
    h2[i,2]  <- NA   
  } else {
    # ---------------------checking outliers---------------------
    out <- outlier(temp$residuals)
    if(length(out) > 0) {
      phenotypic_data[out, wavelength] <- NA
      
      temp <- asreml(
        fixed = get(wavelength) ~ set,
        random = ~name2 + block,
        residual =  ~ar1(range):ar1(row),
        data = phenotypic_data, subset= loc== "EF",na.action = na.method(x = "include"),
        predict = predict.asreml(classify = "name2", sed = TRUE)
      )
    }
    # ---------------------calculating heritability from Cullis formula---------------------
    models[[wavelength]] <- temp
    h2[i,2]  <- (1 - ((models[[wavelength]]$predictions$avsed["mean"] ^ 2) /(2 * summary(models[[wavelength]])$varcomp["name2", "component"])))
  }
}
cat('\n')
# ---------------------save data with no outlier---------------------
fwrite(phenotypic_data, './data/phenotypic_dataef_no_outlier.csv', sep = '\t', row.names = F)
fwrite(h2, "./output/h2EF.csv",row.names = F)

cat('done!')