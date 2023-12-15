# ---------------------loading package--------------------
library(tidyverse)
library(asreml)
library(data.table)
library(ASRgenomics)
library(dplyr)
source("./script/outlier.R")
# ---------------------loading data---------------------
phenotypic_data<- fread('./data/phenotypic_data.csv',data.table= FALSE)
wavebluesEF<- fread("./output/wavebluesEF.csv")
missingingEF <- wavebluesEF %>%
  filter(is.na(predicted.value))%>% 
  select(wave)
phenotypic_data_filtered_ef<- phenotypic_data %>% 
  select(c(1:11),missingingEF$wave)
#fwrite(phenotypic_data_filtered_ef, "./output/intermediate/phenotypic_data_filtered_ef.csv",row.names = F)
phenotypic_data_filtered_ef<- fread("./output/intermediate/phenotypic_data_filtered_ef.csv", data.table= FALSE)
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
variables <- colnames(phenotypic_data_filtered_ef)[grepl("wave", colnames(phenotypic_data_filtered_ef))]
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
      wavebluesEF_na <- bind_rows(wavebluesEF_na, temp)
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
fwrite( wavebluesEF_na, "./output/intermediate/nawavebluesef.csv", row.names = FALSE)
nawavebluesef<- fread("./output/intermediate/nawavebluesef.csv" ,data.table= FALSE)
wavebluesEF_binded <-
  wavebluesEF %>%
  filter(!is.na(predicted.value)) %>%
  bind_rows(nawavebluesef)
wavebluesEF_binded<- wavebluesEF_binded %>% mutate(temp = as.integer(str_sub(wave , start = 6)))%>% 
  arrange(temp, name2) %>% 
  select(-temp)
wavebluesEF_binded<- wavebluesEF_binded %>% distinct(wave, predicted.value, name2)
#fwrite(wavebluesEF_binded, "./output/wavebluesEF_binded.csv")
wavebluesEF_wider<- wavebluesEF_binded %>% pivot_wider(names_from = wave, values_from = predicted.value,names_repair = "check_unique")
#fwrite(wavebluesEF_wider, "./output/wavebluesEF_wider.csv")
# ---------------------calculating blues for missing waves that has na in orginal blues data for location "mw"---------------------
wavebluesMW<- fread("./output/wavebluesMW.csv",data.table= FALSE)
na_wavemw <- subset(wavebluesMW, is.na(wavebluesMW$predicted.value))
phenotypic_data_filtered_mw<- phenotypic_data %>% select(c(1:11),"wave_390")
phenotypic_data_filtered_mw<-
  phenotypic_data_filtered_mw |> mutate(
    name2 = factor(name2),
    taxa = factor(taxa),
    loc = factor(loc),
    set = factor(set),
    block = factor(block),
    range = factor(range),
    row = factor(row),
  ) |>   arrange(loc, range, row)
# ---------------------fitting model---------------------
variables <- colnames(phenotypic_data_filtered_mw)[grepl("wave", colnames(phenotypic_data_filtered_mw))]
models <- vector("list",length(variables))
wavebluesMW_na <- data.frame()
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
      wavebluesMW_na <- bind_rows(wavebluesMW_na, temp)
      cat('\n')
      
    }  else {
      wavebluesMW_na <- bind_rows(wavebluesMW_na, data.frame(name2 = NA,
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

fwrite(wavebluesMW_na, "./output/intermediate/nawavebluesmw.csv", row.names= F)
nawavebluesmw<- fread("./output/intermediate/nawavebluesmw.csv", data.table= F)
wavebluesMW_binded <-
  wavebluesMW %>%
  filter(!is.na(predicted.value)) %>%
  bind_rows(nawavebluesmw)
wavebluesMW_binded<- wavebluesMW_binded %>% mutate(temp = as.integer(str_sub(wave , start = 6)))%>% #esma feri run garera na vako wavebluesMW le analyse garna prxa tara ya maile already na impute gareko suruko bluesMW mai 
  arrange(temp, name2) %>% 
  select(-temp)
wavebluesMW_binded<- wavebluesMW_binded %>% distinct(wave, predicted.value, name2)
#fwrite(wavebluesMW_binded, "./output/wavebluesMW_binded.csv", row.names = F)
wavebluesMW_wider<- wavebluesMW_binded %>% pivot_wider(names_from= wave, values_from= predicted.value, names_repair= "check_unique")
#fwrite(wavebluesMW_wider, "./output/wavebluesMW_wider.csv")


#calculating blues for joint location for each wavelength that has missing values from stage 1 
asreml.options(
  workspace = '8gb',
  pworkspace = '8gb'
)
# ---------------------loading data---------------------
phenotypic_data_joint <- fread('./data/phenotypic_data_joint.csv', data.table = F)
wavejoint<- fread("./output/waveblues_sbf.csv", data.table = F)
wavejoint_long<- pivot_longer(wavejoint, cols = 2:2152,names_to = "wave",values_to = "predicted.value")%>% filter(!is.na(name2))
missingingjoint <- wavejoint_long %>%
  filter(is.na(predicted.value))%>% 
  select(wave)
phenotypic_data_filtered_joint<- phenotypic_data_joint %>% 
  select(c(1:11),missingingjoint$wave)
phenotypic_data_filtered_joint<-
  phenotypic_data_filtered_joint |> mutate(
    name2 = factor(name2),
    taxa = factor(taxa),
    loc = factor(loc),
    set = factor(set),
    block = factor(block),
    range = factor(range),
    row = factor(row),
  ) |>   arrange(loc, range, row)
# ---------------------fitting model---------------------
models <- list()
variables <- colnames(phenotypic_data_filtered_joint)[grepl("wave", colnames(phenotypic_data_filtered_joint))]
models <- vector("list",length(variables))
waveblues_na <- data.frame()
for (i in 1:length(variables)) {
  cat(variables[i], '\n')
  tryCatch({
    model <- asreml(
      fixed = get(variables[i]) ~loc + set + loc:set + name2 + name2:loc,
      random = ~loc:block,
      #residual =  ~id(loc):ar1(range):ar1(row),
      data = phenotypic_data_filtered_joint, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "name2", sed = TRUE)
    )
    if (!model$converge) {
      model <- update.asreml(model)
    }
    if (model$converge) {
      models[[i]] <- model
      #---------------------storing prediction---------------------
      temp <- data.frame(name2 = models[[i]]$predictions$pvals$name2, predicted.value = models[[i]]$predictions$pvals$predicted.value)
      temp$wave <- variables[i]
      waveblues_na <- bind_rows(waveblues_na, temp)
      cat('\n')
      
    }  else {
      waveblues_na <- bind_rows(waveblues_na, data.frame(name2 = NA,
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
fwrite( waveblues_na, "./output/intermediate/nawaveblues.csv", row.names = FALSE)
nawaveblues<- fread("./output/intermediate/nawaveblues.csv")
wavejoint_binded <-
  wavejoint_long %>%
  filter(!is.na(predicted.value)) %>%
  bind_rows(nawaveblues)
wavejoint_binded<- wavejoint_binded %>% mutate(temp = as.integer(str_sub(wave , start = 6)))%>% #esma feri run garera na vako wavebluesMW le analyse garna prxa tara ya maile already na impute gareko suruko bluesMW mai 
  arrange(temp, name2) %>% 
  select(-temp)
wavejoint_binded<- wavejoint_binded %>% distinct(wave, predicted.value, name2)
fwrite(wavejoint_binded, "./output/wavejoint_binded.csv", row.names = F)
# ---------------------imputing na in heritability data for mw location---------------------
phenotypic_data<- fread('./data/phenotypic_data.csv',data.table= FALSE)
h2MW<- fread("./output/h2MW.csv")
nah2_mw <- h2MW %>%
  filter(is.na(h2))%>% 
  select(wave)
phenotypic_data_filtered_h2mw<- phenotypic_data %>% 
  select(c(1:11),nah2_mw$wave)
# ---------------------processing of data---------------------

phenotypic_data_filtered_h2mw<-
  phenotypic_data_filtered_h2mw%>% mutate(
    name2 = factor(name2),
    taxa = factor(taxa),
    loc = factor(loc),
    set = factor(set),
    block = factor(block),
    range = factor(range),
    row = factor(row),
  ) %>%   arrange(loc, range, row)

# ---------------------calculating heritability---------------------

models <- list()
variables <- colnames(phenotypic_data_filtered_h2mw)[12:ncol(phenotypic_data_filtered_h2mw)]
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
    data = phenotypic_data_filtered_h2mw, subset= loc== "MW",na.action = na.method(x = "include"),
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
      phenotypic_data_filtered_h2mw[out, wavelength] <- NA
      
      temp <- asreml(
        fixed = get(wavelength) ~ set,
        random = ~name2 + block,
        #residual =  ~ar1(range):ar1(row),
        data = phenotypic_data_filtered_h2mw, subset= loc== "MW",na.action = na.method(x = "include"),
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
#fwrite(h2, "./output/intermediate/nah2MW.csv",row.names = F)
nah2MW<- fread("./output/intermediate/nah2MW.csv", data.table = F)
h2MW_binded<- bind_rows(h2MW,nah2MW) %>% drop_na() %>% 
  mutate(tempo = as.integer(str_sub(wave , start = 6)))%>% 
  arrange(tempo, h2) %>% 
  select(-tempo)
fwrite(h2MW_binded,"./output/h2MW_binded.csv",row.names = F)
h2MW_binded<- fread("./output/h2MW_binded.csv")
