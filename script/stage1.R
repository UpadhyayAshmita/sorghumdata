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
  phenotypic_data_joint %>% mutate(
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
waveblues<- fread("./output/waveblues.csv")
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
fread("./output/wavebluesMW.csv")
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
multitrait<-
  multitrait|> mutate(
    name2 = factor(name2),
    taxa = factor(taxa),
    loc = factor(loc),
    set = factor(set),
    block = factor(block),
    range = factor(range),
    row = factor(row),
  ) |>   arrange(loc, range, row)


# ---------------------NareaEF blues fitting model---------------------
nareablues_ef_simpler <- asreml(
  fixed = narea ~ name2 +set,
  random =  ~  block,
  #residual =  ~ ar1(range):ar1(row),
  data = multitrait, subset = loc == "EF", 
  na.action = na.method(x = c("include")))
  nareablues_ef_simpler <- update.asreml(nareablues_ef_simpler)
nareablues_ef <- asreml(
  fixed = narea ~ name2 +set,
  random =  ~  block,
  residual =  ~ ar1(range):ar1(row),
  data = multitrait, subset = loc == "EF", 
  na.action = na.method(x = c("include")),
  G.param = nareablues_ef_simpler$G.param,
  R.param = nareablues_ef_simpler$R.param,
  predict = predict.asreml(classify = "name2",vcov = TRUE)
)
nareablues_ef<- update.asreml(nareablues_ef)
nareablues_ef<- update.asreml(nareablues_ef)
  # ---------------------storing prediction---------------------
w <- as.matrix(nareablues_ef$predictions$vcov)
not_missing <- !is.na(nareablues_ef$predictions$pvals$predicted.value)
w <- w[not_missing , not_missing]
w <- round(diag(solve(w)), 3) #Smith et al 2001
nareablues_ef = data.frame(name2 = nareablues_ef$predictions$pvals$name2,
                        predicted.value = round(nareablues_ef$predictions$pvals$predicted.value, 3),
                        W = ifelse(not_missing, w, NA))
fwrite(nareablues_ef,"./output/traitsoutput/nareablues_ef.csv", row.names = F)
# ---------------------SlaEF blues fitting model---------------------
slablues_ef_simpler <- asreml(
  fixed = sla ~ name2 +set,
  random =  ~  block,
  #residual =  ~ ar1(range):ar1(row),
  data = multitrait, subset = loc == "EF", 
  na.action = na.method(x = c("include")))
  slablues_ef_simpler<- update.asreml(slablues_ef_simpler)
slablues_ef <- asreml(
  fixed = sla ~ name2 +set,
  random =  ~  block,
  residual =  ~ ar1(range):ar1(row),
  data = multitrait, subset = loc == "EF", 
  na.action = na.method(x = c("include")),
  G.param = slablues_ef_simpler$G.param,
  R.param = slablues_ef_simpler$R.param,
  predict = predict.asreml(classify = "name2",vcov = TRUE))
slablues_ef<- update.asreml(slablues_ef)

# ---------------------storing prediction---------------------
w <- as.matrix(slablues_ef$predictions$vcov)
not_missing <- !is.na(slablues_ef$predictions$pvals$predicted.value)
w <- w[not_missing , not_missing]
w <- round(diag(solve(w)), 3) #Smith et al 2001
slablues_ef = data.frame(name2 = slablues_ef$predictions$pvals$name2,
                           predicted.value = round(slablues_ef$predictions$pvals$predicted.value, 3),
                           W = ifelse(not_missing, w, NA))
fwrite(slablues_ef,"./output/traitsoutput/slablues_ef.csv", row.names = F)
#BLUES of narea and sla for MW location
# ---------------------fitting model---------------------
  nareablues_mw_simpler <- asreml(
    fixed = narea ~ set+ name2,
    random =  ~  block,
    #residual =  ~ ar1(range):ar1(row),
    data = multitrait, subset = loc == "MW", 
    na.action = na.method(x = c("include")))

nareablues_mw_simpler<- update.asreml(nareablues_mw_simpler)
nareablues_mw <- asreml(
  fixed = narea ~ set+ name2,
  random =  ~  block,
  residual =  ~ ar1(range):ar1(row),
  data = multitrait, subset = loc == "MW", 
  na.action = na.method(x = c("include")),
  G.param = nareablues_mw_simpler$G.param,
  R.param = nareablues_mw_simpler$R.param,
  predict = predict.asreml(classify = "name2",vcov = TRUE))
nareablues_mw<- update.asreml(nareablues_mw)

  # ---------------------storing prediction---------------------
w <- as.matrix(nareablues_mw$predictions$vcov)
not_missing <- !is.na(nareablues_mw$predictions$pvals$predicted.value)
w <- w[not_missing , not_missing]
w <- round(diag(solve(w)), 3) #Smith et al 2001
nareablues_mw = data.frame(name2 = nareablues_mw$predictions$pvals$name2,
                           predicted.value = round(nareablues_mw$predictions$pvals$predicted.value, 3),
                           W = ifelse(not_missing, w, NA)) 
fwrite(nareablues_mw,"./output/traitsoutput/nareablues_mw.csv", row.names = F)
# ---------------------SLA-EF model fitting for blues---------------------
slablues_mw <- asreml(
  fixed = sla ~ set+ name2,
  random =  ~  block,
  residual =  ~ ar1(range):ar1(row),
  data = multitrait, subset = loc == "MW", 
  na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2",vcov = TRUE))
slablues_mw<- update.asreml(slablues_mw)
# ---------------------storing prediction---------------------
w <- as.matrix(slablues_mw$predictions$vcov)
not_missing <- !is.na(slablues_mw$predictions$pvals$predicted.value)
w <- w[not_missing , not_missing]
w <- round(diag(solve(w)), 3) #Smith et al 2001
slablues_mw = data.frame(name2 = slablues_mw$predictions$pvals$name2,
                         predicted.value = round(slablues_mw$predictions$pvals$predicted.value, 3),
                         W = ifelse(not_missing, w, NA))
fwrite(slablues_mw,"./output/traitsoutput/slablues_mw.csv", row.names = F)
# ---------------------calculating blues for joint location for narea and sla
# ---------------------SLAjoint - fitting model---------------------
bluesjoint_sla <- asreml(
  fixed = sla ~ loc + set + loc:set + name2 + name2:loc,
  random =  ~ loc:block,
  residual =  ~ corgh(loc):ar1(range):ar1(row),
  data = multitrait, na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2",vcov = TRUE))
  # ---------------------storing prediction---------------------
w <- as.matrix(bluesjoint_sla$predictions$vcov)
not_missing <- !is.na(bluesjoint_sla$predictions$pvals$predicted.value)
w <- w[not_missing , not_missing]
w <- round(diag(solve(w)), 3) #Smith et al 2001
bluesjoint_sla= data.frame(name2 = bluesjoint_sla$predictions$pvals$name2,
                         predicted.value = round(bluesjoint_sla$predictions$pvals$predicted.value, 3),
                         W = ifelse(not_missing, w, NA))
fwrite(bluesjoint_sla,"./output/traitsoutput/bluesjoint_sla.csv", row.names = F)

# ---------------------N-joint model fitting for blues---------------------
bluesjoint_narea <- asreml(
  fixed = narea ~ loc + set + loc:set + name2 + name2:loc,
  random =  ~ loc:block,
  residual =  ~ corgh(loc):ar1(range):ar1(row),
  data = multitrait, na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "name2",vcov = TRUE))
# ---------------------storing prediction---------------------
w <- as.matrix(bluesjoint_narea$predictions$vcov)
not_missing <- !is.na(bluesjoint_narea$predictions$pvals$predicted.value)
w <- w[not_missing , not_missing]
w <- round(diag(solve(w)), 3) #Smith et al 2001
bluesjoint_narea= data.frame(name2 = bluesjoint_narea$predictions$pvals$name2,
                           predicted.value = round(bluesjoint_narea$predictions$pvals$predicted.value, 3),
                           W = ifelse(not_missing, w, NA))
write.csv(bluesjoint_narea, "./output/traitsoutput/bluesjoint_narea.csv", row.names = F)
# -------------filtering name2 based on corrected name from name_west file-----
Names_WEST<- fread("./data/Names_WEST_SF.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
bluesjoint_narea_correct <-fread("./output/traitsoutput/bluesjoint_narea.csv")|>
    left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>% 
                      mutate(taxa = ifelse((name2 != Corrected_names) & 
                      grepl("PI", Corrected_names), NA, Corrected_names))
bluesjoint_narea_correct<- subset(bluesjoint_narea_correct, !is.na(taxa))%>% 
  select(-name2,-Corrected_names)
#bluesjoint_narea_correct$predicted.value[is.na(bluesjoint_narea_correct$predicted.value)]<- median(bluesjoint_narea_correct$predicted.value, na.rm = T )

bluesjoint_sla_correct <- fread("./output/traitsoutput/bluesjoint_sla.csv")|>
  left_join(Names_WEST %>% 
  dplyr::select(name2, Corrected_names)) %>% 
  mutate(taxa = ifelse((name2 != Corrected_names) & 
                         grepl("PI", Corrected_names), NA, Corrected_names))
bluesjoint_sla_correct<- subset(bluesjoint_sla_correct, !is.na(taxa))%>% 
  select(-name2,-Corrected_names)
#bluesjoint_sla_correct$predicted.value[is.na(bluesjoint_sla_correct$predicted.value)]<- median(bluesjoint_sla_correct$predicted.value, na.rm = T )

nareablues_mw_correct<- fread("./output/traitsoutput/nareablues_mw.csv") %>% 
  left_join(Names_WEST %>%dplyr::select(name2, Corrected_names)) %>% 
  mutate(taxa = ifelse((name2 != Corrected_names) & 
 grepl("PI", Corrected_names), NA, Corrected_names))
nareablues_mw_correct<- subset(nareablues_mw_correct, !is.na(taxa))%>% select(-name2,-Corrected_names)
# nareablues_mw_correct$predicted.value[is.na(nareablues_mw_correct$predicted.value)]<- median(nareablues_mw_correct$predicted.value, na.rm = T )

nareablues_ef_correct<- fread("./output/traitsoutput/nareablues_ef.csv") %>% 
  left_join(Names_WEST %>%dplyr::select(name2, Corrected_names)) %>% 
  mutate(taxa = ifelse((name2 != Corrected_names) & 
                         grepl("PI", Corrected_names), NA, Corrected_names))
nareablues_ef_correct<- subset(nareablues_ef_correct, !is.na(taxa))%>% select(-name2,-Corrected_names)
# nareablues_ef_correct$predicted.value[is.na(nareablues_ef_correct$predicted.value)]<- median(nareablues_ef_correct$predicted.value, na.rm = T )

slablues_ef_correct<- fread("./output/traitsoutput/slablues_ef.csv") %>% 
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>% 
  mutate(taxa = ifelse((name2 != Corrected_names) & 
                         grepl("PI", Corrected_names), NA, Corrected_names))
 slablues_ef_correct<- subset(slablues_ef_correct, !is.na(taxa))%>% select(-name2,-Corrected_names)
# slablues_ef_correct$predicted.value[is.na(slablues_ef_correct$predicted.value)]<- median(slablues_ef_correct$predicted.value, na.rm = T )

slablues_mw_correct<- fread("./output/traitsoutput/slablues_mw.csv")%>% 
  left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>% 
  mutate(taxa = ifelse((name2 != Corrected_names) & 
                         grepl("PI", Corrected_names), NA, Corrected_names))
 slablues_mw_correct<- subset(slablues_mw_correct, !is.na(taxa))%>% select(-name2,-Corrected_names)
# slablues_mw_correct$predicted.value[is.na(slablues_mw_correct$predicted.value)]<- median(slablues_mw_correct$predicted.value, na.rm = T )



