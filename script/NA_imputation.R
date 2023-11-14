# ---------------------loading package--------------------
library(tidyverse)
library(asreml)
library(data.table)
library(ASRgenomics)
library(dplyr)
library(janitor)
# ---------------------loading data---------------------
phenotypic_data<- fread('./data/phenotypic_data.csv',data.table= FALSE)
wavebluesEF<- fread("./output/wavebluesEF.csv")
na_wave <- subset(wavebluesEF, is.na(wavebluesEF$predicted.value))
phenotypic_data_filtered_ef<- phenotypic_data %>% select(c(1:11),"wave_1140", "wave_1148", "wave_1149","wave_1150", "wave_1154", "wave_1157","wave_1159", 
  "wave_1162", "wave_1163", "wave_1164", "wave_1166", "wave_1167"
  ,"wave_1168","wave_1170","wave_1172", "wave_1173","wave_1176", "wave_1178","wave_1179","wave_1182","wave_1183")
fwrite(phenotypic_data_filtered_ef, "./intermediate/phenotypic_data_filtered_ef.csv",row.names = F)
phenotypic_data_filtered_ef<- fread("./intermediate/phenotypic_data_filtered_ef.csv")
# ---------------------processing of data---------------------
phenotypic_data_filtered_ef<-
  phenotypic_data_filtered_ef |> mutate(
    name2 = factor(name2),
    taxa = factor(taxa),
    loc = factor(loc),
    set = factor(set),
    block = factor(block),
    range = factor(range),
    row = factor(row),
  ) |>   arrange(loc, range, row)


#calculating blues for each wavelength 
variables <- colnames(phenotypic_data_filtered_ef)[12:32]
models <- vector("list",length(variables))
wavebluesEF_na <-data.frame()
# ---------------------fitting model---------------------
for (i in 1:length(variables)) {
  cat(variables[i], '\n')
  tryCatch({
    model <- asreml(
      fixed = get(variables[i]) ~set + name2,
      random = ~block,
      #residual =  ~ar1(range):ar1(row),
      data = phenotypic_data_filtered_ef,subset= loc== "EF", na.action = na.method(x = "include"),
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
      wavebluesEF <- bind_rows(wavebluesEF_na, temp)
      cat('\n')
      
    }  else {
      wavebluesEF_na <- bind_rows(wavebluesEF_na, data.frame(name2 = NA,
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
fwrite( wavebluesEF_na, "./intermediate/nawavebluesef.csv", row.names = FALSE)
nawavebluesef<- fread("./intermediate/nawavebluesef.csv")
wavebluesEF_binded<- rbind(wavebluesEF,nawavebluesef)
wavebluesEF_binded<- wavebluesEF_binded %>% mutate(temp = as.integer(str_sub(wave , start = 6)))%>% 
  arrange(temp, name2) %>% 
  select(-temp)
wavebluesef<- wavebluesEF_binded %>% distinct(wave, predicted.value, name2)
fwrite(wavebluesef, "./output/wavebluesef.csv", row.names = F)

# ---------------------calculating blues for missing waves that has na in orginal blues data for location "mw"---------------------
wavebluesMW<- fread("./output/wavebluesMW.csv")
na_wavemw <- subset(wavebluesMW, is.na(wavebluesMW$predicted.value))
phenotypic_data_filtered_mw<- design %>% select(c(1:11),"wave_390", "wave_391")

phenotypic_data_filtered_mw<-
  phenotypic_data_filtered_mw |> clean_names() |> mutate(
    name2 = factor(name2),
    taxa = factor(taxa),
    loc = factor(loc),
    set = factor(set),
    block = factor(block),
    range = factor(range),
    row = factor(row),
  ) |>   arrange(loc, range, row)
# ---------------------fitting model---------------------
variables <- colnames(phenotypic_data_filtered_mw)[10:11]
models <- vector("list",length(variables))
wavebluesMW <-data.frame()
# ---------------------fitting model---------------------
for (i in 1:length(variables)) {
  cat(variables[i], '\n')
  tryCatch({
    model <- asreml(
      fixed = get(variables[i]) ~set + name2,
      random = ~block,
      #residual =  ~ar1(range):ar1(row),
      data = phenotypic_data_filtered_mw,subset= loc== "MW", na.action = na.method(x = "include"),
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
fwrite(wavebluesMW, "./intermediate/nawavebluesmw.csv", row.names= F)
nawavebluesmw<- fread("./intermediate/nawavebluesmw.csv", data.table= F)
wavebluesMW_binded<- rbind(wavebluesMW,nawavebluesmw)
wavebluesMW_binded<- wavebluesMW_binded %>% mutate(temp = as.integer(str_sub(wave , start = 6)))%>% 
  arrange(temp, name2) %>% 
  select(-temp)
wavebluesmw<- wavebluesmw_binded %>% distinct(wave, predicted.value, name2)
fwrite(wavebluesmw, "./output/wavebluesmw.csv", row.names = F)
# ---------------------dealing with 33 na observation in h2 mw data and calculating h2 through different model---------------------
# ---------------------calculating heritability---------------------
# ---------------------loading library---------------------
library(tidyverse)
library(asreml)
library(data.table)
library(ASRgenomics)
source('./script/outlier.R')

# ---------------------loading data---------------------

phenotypic_data<- fread('./data/designnew.csv',data.table= FALSE)

designnew_long<- pivot_longer(designnew, cols = 11:2161,names_to = "wave",values_to = "blues")

wave_h2na<-na_h2mw$wave


designnew_na<- designnew_long %>% filter(wave %in% wave_h2na)

designfiltered<- pivot_wider(designnew_na, names_from = "wave", values_from = "blues")

# ---------------------processing of data---------------------

designfiltered<-
  designfiltered %>% clean_names()%>% mutate(
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
variables <- colnames(designfiltered)[12:ncol(designfiltered)]
h2 <- data.frame(matrix(NA, length(variables), 2))  #n rows and 2 columns
h2[,1] <- variables

colnames(h2) <- c("wave", "h2")


#for(i in variables)
for (i in 1:length(variables)) {
  wavelength<- variables[i]
  
  
  # ---------------------fitting model---------------------
  
  cat( wavelength, '\n')
  
  temp <- asreml(
    fixed = get(wavelength) ~ set,
    random = ~name2 + block,
    #residual =  ~ar1(range):ar1(row),
    data = designfiltered, subset= loc== "MW",na.action = na.method(x = "include"),
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
      designfiltered[out, wavelength] <- NA
      
      temp <- asreml(
        fixed = get(wavelength) ~ set,
        random = ~name2 + block,
        #residual =  ~ar1(range):ar1(row),
        data = designfiltered, subset= loc== "MW",na.action = na.method(x = "include"),
        predict = predict.asreml(classify = "name2", sed = TRUE)
      )
      
      
    }
    
    
    # ---------------------calculating heritiability from formula---------------------
    models[[wavelength]] <- temp
    
    h2[i,2]  <- (1 - ((models[[wavelength]]$predictions$avsed["mean"] ^ 2) /(2 * summary(models[[wavelength]])$varcomp["name2", "component"])))
    
  }
}

cat('\n')

# ---------------------save data with no outlier---------------------

#fwrite(designnew, './data/designMW_no_outlier.csv', sep = '\t', row.names = F)

fwrite(h2, "./intermediate/nah2MW.csv",row.names = F)

nah2mw<- fread("./intermediate/nah2MW.csv", data.table = F)
h2MW<- bind_rows(h2mw,nah2mw)

fwrite(h2MW,"./output/h2MW.csv",row.names = F)
h2MW<- fread("./output/h2MW.csv")
