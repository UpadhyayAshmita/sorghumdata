# ---------------------loading package--------------------
library(tidyverse)
library(asreml)
library(data.table)
library(ASRgenomics)
source('./script/outlier.R')
# ---------------------loading data---------------------
phenotypic_data<- fread('./data/phenotypic_data.csv',data.table= FALSE)
wavebluesEF<- fread("./output/wavebluesEF.csv")
na_wave <- subset(wavebluesEF, is.na(wavebluesEF$predicted.value))
phenotypic_data_filtered_ef<- phenotypic_data %>% select(c(1:11),"wave_1140", "wave_1148", "wave_1149","wave_1150", "wave_1154", "wave_1157","wave_1159", 
  "wave_1162", "wave_1163", "wave_1164", "wave_1166", "wave_1167"
  ,"wave_1168","wave_1170","wave_1172", "wave_1173","wave_1176", "wave_1178","wave_1179","wave_1182","wave_1183")
#fwrite(phenotypic_data_filtered_ef, "./intermediate/phenotypic_data_filtered_ef.csv",row.names = F)
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
wavebluesEF_binded<- wavebluesEF_binded %>% distinct(wave, predicted.value, name2)
#fwrite(wavebluesEF_binded, "./output/wavebluesEF_binded.csv")
wavebluesEF_wider<- wavebluesEF_binded %>% pivot_wider(names_from = wave, values_from = predicted.value,names_repair = "check_unique")
#fwrite(wavebluesEF_wider, "./output/wavebluesEF_wider.csv")
# ---------------------calculating blues for missing waves that has na in orginal blues data for location "mw"---------------------
wavebluesMW<- fread("./output/wavebluesMW.csv")
na_wavemw <- subset(wavebluesMW, is.na(wavebluesMW$predicted.value))
phenotypic_data_filtered_mw<- phenotypic_data %>% select(c(1:11),"wave_390", "wave_391")

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
variables <- colnames(phenotypic_data_filtered_mw)[12:13]
models <- vector("list",length(variables))
wavebluesM <-data.frame()
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
      wavebluesM <- bind_rows(wavebluesM, temp)
      cat('\n')
      
    }  else {
      wavebluesM <- bind_rows(wavebluesM, data.frame(name2 = NA,
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
fwrite(wavebluesM, "./intermediate/nawavebluesmw.csv", row.names= F)
nawavebluesmw<- fread("./intermediate/nawavebluesmw.csv", data.table= F)
wavebluesMW_binded<- rbind(wavebluesMW,nawavebluesmw)
wavebluesMW_binded<- wavebluesMW %>% mutate(temp = as.integer(str_sub(wave , start = 6)))%>% #esma feri run garera na vako wavebluesMW le analyse garna prxa tara ya maile already na impute gareko suruko bluesMW mai 
  arrange(temp, name2) %>% 
  select(-temp)
wavebluesMW_binded<- wavebluesMW_binded %>% distinct(wave, predicted.value, name2)
fwrite(wavebluesMW_binded, "./output/wavebluesMW_binded.csv", row.names = F)
wavebluesMW_wider<- wavebluesMW_binded %>% pivot_wider(names_from= wave, values_from= predicted.value, names_repair= "check_unique")
#fwrite(wavebluesMW_wider, "./output/wavebluesMW_wider.csv")
# ---------------------dealing with 33 na observation in h2 mw data and calculating h2 through different model---------------------
# ---------------------loading data---------------------
phenotypic_data<- fread('./data/phenotypic_data.csv',data.table= FALSE)
phenotypic_data_longer<- pivot_longer(phenotypic_data, cols = 12:2162,names_to = "wave",values_to = "blues")
h2MW<- fread("./output/h2MW.csv")
na_h2mw <- subset(h2MW, is.na(h2MW$h2))
wave_h2na<-na_h2mw$wave
phenotypicdata_na<-phenotypic_data_longer %>% filter(wave %in% wave_h2na)
phenotypic_data_filter_h2mw<- pivot_wider(phenotypicdata_na, names_from = "wave", values_from = "blues")
# ---------------------processing of data---------------------
phenotypic_data_filter_h2mw<-
  phenotypic_data_filter_h2mw%>% mutate(
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
variables <- colnames(phenotypic_data_filter_h2mw)[12:ncol(phenotypic_data_filter_h2mw)]
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
    data = phenotypic_data_filter_h2mw, subset= loc== "MW",na.action = na.method(x = "include"),
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
      phenotypic_data_filter_h2mw[out, wavelength] <- NA
      temp <- asreml(
        fixed = get(wavelength) ~ set,
        random = ~name2 + block,
        #residual =  ~ar1(range):ar1(row),
        data = phenotypic_data_filter_h2mw, subset= loc== "MW",na.action = na.method(x = "include"),
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
h2MW_binded<- bind_rows(h2MW,nah2mw)
h2MW_binded<- drop_na(h2MW_binded)
#fwrite(h2MW_binded,"./output/h2MW_binded.csv",row.names = F)
h2MW_binded<- fread("./output/h2MW_binded.csv")
