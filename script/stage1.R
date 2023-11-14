library(tidyverse)
library(asreml)
source('./script/outlier.R')
library(data.table)
library(ASRgenomics)
# ---------------------allocating workspace--------------------
asreml.options(
  workspace = '8gb',
  pworkspace = '8gb'
)
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
# ---------------------calculating heritability of each wavelength for location MW---------------------
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
    data = phenotypic_data, subset= loc== "MW",na.action = na.method(x = "include"),
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
        data = phenotypic_data, subset= loc== "MW",na.action = na.method(x = "include"),
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
fwrite(phenotypic_data, './data/phenotypic_datamw_no_outlier.csv', sep = '\t', row.names = F)
fwrite(h2, "./output/h2MW.csv",row.names = F)
#h2<- fread("./output/h2MW.csv", data.table = FALSE)
cat('done!')
# ---------------------calculating heritability of wavelength for loc EF---------------------
# ---------------------loading data---------------------
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
fwrite(h2, "./output/h2ef.csv",row.names = F)
#h2<- fread("./output/h2MW.csv", data.table = FALSE)
cat('done!')

# ---------------------heritability for joint location---------------------
# ---------------------loading data---------------------
phenotypic_datamw<- fread('./data/phenotypic_datamw_no_outlier.csv',data.table= FALSE) |>
  filter(loc == "MW")
phenotypic_dataef<- fread('./data/phenotypic_dataef_no_outlier.csv',data.table= FALSE)|>
  filter(loc == "EF")
phenotypic_data_joint<- bind_rows(phenotypic_datamw, phenotypic_dataef)
fwrite(phenotypic_data_joint, "./data/phenotypic_data_joint.csv", row.names = F)

# ---------------------processing of data---------------------
phenotypic_data_joint<-
  phenotypic_data_joint %>% clean_names()%>% mutate(
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
variables <- colnames(design)[12:2162]
h2 <- data.frame(matrix(NA, length(variables), 2))  #n rows and 2 columns
h2[,1] <- variables
colnames(h2) <- c("wave", "h2")
for (i in 1:length(variables)) {
  wavelength<- variables[i]
  # ---------------------fitting model---------------------
  cat( wavelength, '\n')
  temp <- asreml(
    fixed = get(wavelength) ~ set + loc + loc: set,
    random = ~name2 + loc:block + loc:name2,
    residual =  ~id(loc):ar1(range):ar1(row),
    data = phenotypic_data_joint, na.action = na.method(x = "include"),
    predict = predict.asreml(classify = "name2", sed = TRUE)
  )
  if (!temp$converge) {
    temp <- update.asreml(temp)
  }  
  if (!temp$converge) {
    models[[wavelength]] <- NA
    h2[i,2]  <- NA   
  } else {
    # ---------------------checking outliers--------------------
    out <- outlier(temp$residuals)
    if(length(out) > 0) {
      design[out, wavelength] <- NA
      temp <- asreml(
        fixed = get(wavelength) ~ set + loc + loc: set,
        random = ~name2 + loc:block + loc:name2,
        residual =  ~id(loc):ar1(range):ar1(row),
        data = phenotypic_data_joint, na.action = na.method(x = "include"),
        predict = predict.asreml(classify = "name2", sed = TRUE)
      )
    }
    # ---------------------calculating heritability from cullis formula---------------------
    models[[wavelength]] <- temp
    h2[i,2]  <- (1 - ((models[[wavelength]]$predictions$avsed["mean"] ^ 2) /(2 * summary(models[[wavelength]])$varcomp["name2", "component"])))
  }
}
cat('\n')
# ---------------------save data with no outlier---------------------
fwrite(phenotypic_data_joint, './data/phenotypic_data_joint_no_outlier.csv', sep = '\t', row.names = F)
fwrite(h2, "./output/h2.csv",row.names = F)
#h2<- fread("./output/h2.csv", data.table = FALSE)
cat('done!')

#calculating blues for joint location for each wavelength 
asreml.options(
  workspace = '8gb',
  pworkspace = '8gb'
)
# ---------------------loading data---------------------
phenotypic_data_joint <- fread('./data/phenotypic_data_joint.csv', data.table = F)
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
# ---------------------fitting model---------------------
models <- list()
variables <- colnames(phenotypic_data_joint)[12:2162]
waveblues <- data.frame()
for (i in 1:length(variables)) {
  cat(variables[i], '\n')
  tryCatch({
    model <- asreml(
      fixed = get(variables[i]) ~loc + set + loc:set + name2 + name2:loc,
      random = ~loc:block,
      residual =  ~id(loc):ar1(range):ar1(row),
      data = phenotypic_data_joint, na.action = na.method(x = "include"),
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
fwrite( waveblues, "waveblues.csv", row.names = FALSE)
#waveblues<- fread("./output/waveblues.csv")
cat('Done with Joint blues')

# ---------------------calculating blues for loc MW---------------------
# ---------------------loading data---------------------
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

# ---------------------calculating blues for EF location---------------------
# ---------------------loading data---------------------
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
variables <- colnames(phenotypic_data_joint)[12:2162]
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
      data = phenotypic_data_joint,subset= loc== "EF", na.action = na.method(x = "include"),
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
# ---------------------Obtaining BLUEs for N and SLA---------------------
# -------------------calculating heritability for narea for each location---------------------
# ---------------------fitting model---------------------
nareaMW <- asreml(
  fixed = narea ~ set,
  random = ~name2 + block,
  residual =  ~ar1(range):ar1(row),
  data = phenotypic_data, subset= loc== "MW",na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)
nareaMW<- update.asreml(nareaMW)
# ---------------------checking outliers---------------------
out <- outlier(nareaMW$residuals)
if(length(out) > 0) {
  phenotypic_data[out, "narea"] <- NA
} else{
  nareaMW <- asreml(
    fixed = narea ~ set,
    random = ~name2 + block,
    residual =  ~ar1(range):ar1(row),
    data = phenotypic_data, subset= loc== "MW",na.action = na.method(x = "include"),
    predict = predict.asreml(classify = "name2", sed = TRUE)
  )
  nareaMW<- update.asreml(nareaMW)
}
# ---------------------calculating and storing the heritability---------------------
temp <- (1 - ((nareaMW$predictions$avsed["mean"] ^ 2) /(2 * summary(nareaMW)$varcomp["name2", "component"])))

# ---------------------h2 was found 0.587 for narea for location mw---------------------
# ---------------------storing outlier filtered data for location mw---------------------
fwrite(phenotypic_data, "./data/phenotypedataMW.csv", row.names= FALSE)

# ---------------------reading data---------------------
phenotypedataMW<-  fread("./data/phenotypedataMW.csv")
# ---------------------checking residual plot for the model---------------------
#plot(h2narea)

# ---------------------calculating heritability for narea for ef location---------------------
# ---------------------fitting model---------------------
nareaEF <- asreml(
  fixed = narea ~ set,
  random = ~name2 + block,
  residual =  ~ar1(range):ar1(row),
  data = phenotypic_data, subset= loc== "EF",na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)

nareaEF<- update.asreml(nareaEF)

# ---------------------checking outliers---------------------
out <- outlier(nareaEF$residuals)

if(length(out) > 0) {
  phenotypic_data[out, "narea"] <- NA
}else {
  
  nareaEF <- asreml(
    fixed = narea ~ set,
    random = ~name2 + block,
    residual =  ~ar1(range):ar1(row),
    data = phenotypic_data, subset= loc== "EF",na.action = na.method(x = "include"),
    predict = predict.asreml(classify = "name2", sed = TRUE)
  )
  nareaEF<- update.asreml(nareaEF)
}
# ---------------------calculating and storing heritability---------------------
temp2<- (1 - ((nareaEF$predictions$avsed["mean"] ^ 2) /(2 * summary(nareaEF)$varcomp["name2", "component"])))

# ---------------------saving outliers filtered data for location ef for narea---------------------
fwrite(phenotypic_data, "./data/phenotypedataEF.csv", row.names= FALSE)

# ---------------------h2 for narea for loc ef was 0.355---------------------
# ---------------------checking residual plot---------------------
#plot(nareaEF)
# ---------------------calculating heritability for sla for MW location---------------------
# ---------------------fitting model---------------------
slaMW <- asreml(
  fixed = sla ~ set,
  random = ~name2 + block,
  residual =  ~ar1(range):ar1(row),
  data = phenotypic_data, subset= loc== "MW",na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)

slaMW<- update.asreml(slaMW)

# ---------------------checking outliers---------------------
out <- outlier(slaMW$residuals)

if(length(out) > 0) {
  phenotypic_data[out, "sla"] <- NA
}else{
  
  slaMW <- asreml(
    fixed = sla ~ set,
    random = ~name2 + block,
    residual =  ~ar1(range):ar1(row),
    data = phenotypic_data, subset= loc== "MW",na.action = na.method(x = "include"),
    predict = predict.asreml(classify = "name2", sed = TRUE)
  )
  slaMW<- update.asreml(slaMW)
}

# ---------------------calculating and storing heritability---------------------
temp3 <- (1 - ((slaMW$predictions$avsed["mean"] ^ 2) /(2 * summary(slaMW)$varcomp["name2", "component"])))

# ---------------------storing outliers filtered data for sla for loc mw---------------------
fwrite(phenotypic_data, "./data/phenotypedataslaMW.csv", row.names= FALSE)

# ---------------------checking residual plot---------------------
#plot(slaMW)

# ---------------------calculating heritability for sla for EF location ---------------------
# ---------------------fitting model---------------------
slaEF <- asreml(
  fixed = sla ~ set,
  random = ~name2 + block,
  residual =  ~ar1(range):ar1(row),
  data = phenotypic_data, subset= loc== "EF",na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)

slaEF<- update.asreml(slaEF)

# ---------------------checking outliers---------------------
out <- outlier(slaMW$residuals)

if(length(out) > 0) {
  phenotypic_data[out, "sla"] <- NA
} else{
  slaEF <- asreml(
    fixed = sla ~ set,
    random = ~name2 + block,
    residual =  ~ar1(range):ar1(row),
    data = phenotypic_data, subset= loc== "EF",na.action = na.method(x = "include"),
    predict = predict.asreml(classify = "name2", sed = TRUE)
  )
  slaEF<- update.asreml(slaEF)
}

# ---------------------calculating and storing h2 for sla for loc ef---------------------
temp4 <- (1 - ((slaEF$predictions$avsed["mean"] ^ 2) /(2 * summary(slaEF)$varcomp["name2", "component"])))

# ---------------------storing outliers filtered data for sla for loc ef---------------------
fwrite(phenotypic_data, "./data/phenotypedataslaEF.csv", row.names= FALSE)

# ---------------------residual plot---------------------
#plot(slaEF)

# ---------------------joint heritability model for narea ---------------------
# ---------------------binding both location narea outliers filtered data to pass into joint model---------------------

narea_ef<- fread('./data/phenotypedataEF.csv',data.table= FALSE) |>
  filter(loc == "EF")
narea_mw<- fread('./data/phenotypedataMW.csv',data.table= FALSE)|>
  filter(loc == "MW")
narea<- bind_rows(narea_ef, narea_mw)

# ---------------------storing narea outliers filtered both location data as narea---------------------
fwrite(narea, "./data/narea.csv",row.names = F)

# ---------------------reading the data ---------------------
#narea<- fread("./data/narea.csv", data.table= F)
# ---------------------fitting the model ---------------------
nareajointh2 <- asreml(
  fixed = narea ~ set + loc + loc: set,
  random = ~name2 + loc:block + loc:name2,
  residual =  ~id(loc):ar1(range):ar1(row),
  data = narea, na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)

nareajointh2<- update.asreml(nareajointh2)

# ---------------------calulating h2 for narea joint loc---------------------
temp5 <- (1 - ((nareajointh2$predictions$avsed["mean"] ^ 2) /(2 * summary(nareajointh2)$varcomp["name2", "component"])))

# -----------------------binding both location outliers filtered sla data into one-------------------

sla_ef<- fread('./data/phenotypedataslaEF.csv',data.table= FALSE) |>
  filter(loc == "EF")
sla_mw<- fread('./data/phenotypedataslaMW.csv',data.table= FALSE)|>
  filter(loc == "MW")
sla<- bind_rows(sla_ef, sla_mw)
# ---------------------storing binded both location outliers filtered data into one file sla---------------------
fwrite(sla, "./data/sla.csv",row.names = F)
# ---------------------reading the sla data for combined location---------------------
#sla<- fread("./data/sla.csv", data.table = F)
# ---------------------preprocessing of data---------------------
# ---------------------model fitting---------------------

slajointh2 <- asreml(
  fixed = sla ~ set + loc + loc: set,
  random = ~name2 + loc:block + loc:name2,
  residual =  ~id(loc):ar1(range):ar1(row),
  data = sla, na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)

slajointh2<- update.asreml(slajointh2)

#calculating h2 for joint location for sla

temp6 <- (1 - ((slajointh2$predictions$avsed["mean"] ^ 2) /(2 * summary(slajointh2)$varcomp["name2", "component"])))

# ---------------------calculating blues for narea and sla with outliers filtered data narea and sla---------------------
# ---------------------binding both trait both location filtered data---------------------
narea<- subset(narea, select= -sla)
sla<- subset(sla, select= -narea)
multitrait<- full_join(narea, sla)
fwrite(multitrait,"./data/multitrait.csv", row.names = F)#multitrait is outlier filtered data for narea and sla combined
#calculating BLUES of narea and sla for  EF location 
# ---------------------loading data---------------------
multitrait<- fread("./data/multitrait.csv", data.table = F)
models <- list()
bluesEF <- list()
a <- 1
for(i in c("narea","sla")){
  # ---------------------fitting model---------------------
  
  models[[i]] <- asreml(
    fixed = get(i) ~ name2 +set,
    random =  ~  block,
    #residual =  ~ id(range):id(row),
    data = multitrait, subset = loc == "EF", na.action = na.method(x = c("include")),
    predict = predict.asreml(classify = "name2",sed = TRUE))
  # ---------------------storing prediction---------------------
  
  bluesEF[[i]] <-   models[[i]]$predictions$pvals[,1:2]
  colnames(bluesEF[[i]]) <- c("name2", i)
  a <- a + 1
}
#calculating BLUES of narea and sla for MW location
models <- list()
bluesMW <- list()
a <- 1
for(i in c("narea","sla")){
  
  # ---------------------fitting model---------------------
  models[[i]] <- asreml(
    fixed = get(i) ~ set+ name2,
    random =  ~  block,
    #residual =  ~ id(range):id(row),
    data = multitrait, subset = loc == "MW", na.action = na.method(x = c("include")),
    predict = predict.asreml(classify = "name2",sed = TRUE))
  
  # ---------------------storing prediction---------------------
  
  bluesMW[[i]] <-   models[[i]]$predictions$pvals[,1:2]  
  colnames(bluesMW[[i]]) <- c("name2", i)
  a <- a + 1}
blues <- 
  bind_rows(
    "EF" = bluesEF %>% reduce(left_join, by = "name2"),
    "MW" = bluesMW %>% reduce(left_join, by = "name2"),
    .id = "loc") %>% left_join(multitrait %>% filter(!duplicated(name2)) %>% select(name2))
# ---------------------filtering name2 basede on corrected name from name_west file---------------------
Names_WEST<- fread("./data/Names_WEST.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
blues <- blues |> left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>% mutate(taxa = ifelse((name2 != Corrected_names) & grepl("PI", Corrected_names), NA, Corrected_names))
blues<- subset(blues, !is.na(taxa))
# length(unique(blues$taxa))
# blues<- drop_na(blues)
# ---------------------reading the blues data file---------------------
blues<- blues %>% select(-name2,-Corrected_names)
# ---------------------storing blues value for both trait and both location---------------------
fwrite(blues, "./traitsoutput/blues.csv")

# ---------------------loading blues for sla and narea for two location---------------------
blues<- fread("./traitsoutput/blues.csv")
nareabluesef<- subset(blues, loc== "EF") %>% select(-sla)
nareabluesmw<- subset(blues, loc== "MW") %>% select(-sla)
slabluesef<- subset(blues, loc== "EF") %>% select(-narea)
slabluesmw<- subset(blues, loc== "MW") %>% select(-narea)
fwrite(nareabluesef,"./traitsoutput/nareabluesef.csv", row.names = F)
fwrite(nareabluesmw,"./traitsoutput/nareabluesmw.csv", row.names = F)

fwrite(slabluesef,"./traitsoutput/slabluesef.csv", row.names = F)
fwrite(slabluesmw,"./traitsoutput/slabluesmw.csv", row.names = F)

