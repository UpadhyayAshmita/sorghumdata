# ---------------------loading library---------------------
library(tidyverse)
library(asreml)
library(data.table)
library(ASRgenomics)
library(janitor)
source('./script/outlier.R')

# ---------------------loading data---------------------

design_ef<- fread('./data/designEF_no_outlier.csv',data.table= FALSE) |>
  filter(loc == "EF")
design_mw<- fread('./data/designMW_no_outlier.csv',data.table= FALSE)|>
  filter(loc == "MW")

design <- bind_rows(design_ef, design_mw)
fwrite(design, "./data/design.csv", row.names = F)
fread("./data/design.csv")

# ---------------------processing of data---------------------

design<-
  design %>% clean_names()%>% mutate(
    name2 = factor(name2),
    taxa = factor(taxa),
    loc = factor(loc),
    set = factor(set),
    block = factor(block),
    range = factor(range),
    row = factor(row),
    uni = c(1:960, 1:960)
  ) %>%   arrange(loc, range, row)
design<- design%>% select(-c(v1))

# ---------------------calculating heritability---------------------


models <- list()
variables <- colnames(design)[10:2160]
h2 <- data.frame(matrix(NA, length(variables), 2))  #n rows and 2 columns
h2[,1] <- variables

colnames(h2) <- c("wave", "h2")


#for(i in variables)
for (i in 1:length(variables)) {
  wavelength<- variables[i]
  
  
  # ---------------------fitting model---------------------
  
  cat( wavelength, '\n')
  
  temp <- asreml(
    fixed = get(wavelength) ~ set + loc + loc: set,
    random = ~name2 + loc:block + loc:name2,
    residual =  ~id(loc):ar1(range):ar1(row),
    data = design, na.action = na.method(x = "include"),
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
      design[out, wavelength] <- NA
      
      temp <- asreml(
        fixed = get(wavelength) ~ set + loc + loc: set,
        random = ~name2 + loc:block + loc:name2,
        residual =  ~id(loc):ar1(range):ar1(row),
        data = design, na.action = na.method(x = "include"),
        predict = predict.asreml(classify = "name2", sed = TRUE)
      )
      
      
    }
    
    
    # ---------------------calculating heritability from formula---------------------
    models[[wavelength]] <- temp
    
    h2[i,2]  <- (1 - ((models[[wavelength]]$predictions$avsed["mean"] ^ 2) /(2 * summary(models[[wavelength]])$varcomp["name2", "component"])))
    
  }
}

cat('\n')

# ---------------------save data with no outlier---------------------

fwrite(design, './data/design_no_outlier.csv', sep = '\t', row.names = F)
fwrite(h2, "./output/h2.csv",row.names = F)
#h2<- fread("./output/h2.csv", data.table = FALSE)

cat('done!')







