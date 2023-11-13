library(tidyverse)
library(asreml)
source('./script/outlier.R')
library(data.table)
library(ASRgenomics)
#library(flextable)
#library(SmartEDA)
#library(explore)
#library(patchwork)
#library(janitor)

asreml.options(
  workspace = '8gb',
  pworkspace = '8gb'
)
# ---------------------load data---------------------
phenotypic_data<- fread('./data/phenotypic_data.csv',data.table= FALSE)

# ---------------------processing data---------------------
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
# ---------------------calculating heritability of each wavelength--------------------

models <- list()
variables <- colnames(phenotypic_data)[12:2162]
h2 <- data.frame(matrix(NA, length(variables), 2))  #n rows and 2 columns
h2[,1] <- variables

colnames(h2) <- c("wave", "h2")

for (i in 1:length(variables)) {
  wavelength<- variables[i]
  # ---------------------fitting model---------------------
  cat( wavelength, '\n')
  
  models[[wavelength]] <- asreml(
    fixed = get(wavelength) ~ set + loc + loc: set,
    random = ~name2 + loc:block + loc:name2,
    residual =  ~id(loc):ar1(range):ar1(row),
    data = phenotypic_data, na.action = na.method(x = "include"),
    predict = predict.asreml(classify = "name2", sed = TRUE)
  )
  
  # ---------------------checking outliers---------------------
  out <- outlier(models[[wavelength]]$residuals)
  
  if(length(out) > 0) {
    phenotypic_data[,wavelength][phenotypic_data$loc][out] <- NA
    models[[wavelength]] <- asreml(
      fixed = get(wavelength) ~ set + loc + loc:set,
      random =  ~ name2 + loc:block + loc:name2,
      residual =  ~ id(loc):ar1(range):ar1(row),
      data = phenotypic_data, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "name2", sed = TRUE)
    )
  }
  
  # ---------------------calculating heritability from Cullis---------------------
  h2[i,2]  <- (1 - ((models[[wavelength]]$predictions$avsed["mean"] ^ 2) /(2 * summary(models[[wavelength]])$varcomp["name2", "component"])))
}

fwrite(h2, "./output/h2.csv",row.names = F)
#h2<- fread("./output/h2.csv", data.table = FALSE)

#calculating blues for joint location for each wavelength 
models <- list()
variables <- colnames(phenotypic_data)[12:2162]
waveblues <- data.frame()

# ---------------------fitting model---------------------
for (i in 1:length(variables)) {
  cat(variables[i], '\n')
  
  models[[i]] <- tryCatch({
    asreml(
      fixed = get(variables[i]) ~ loc + set + loc:set + name2 + name2:loc,
      random = ~loc:block,
      residual =  ~ id(loc):ar1(range):ar1(row),
      data = phenotypic_data, na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "name2", sed = TRUE)
    )
  }, error = function(err) {
    message("An error occured")
    print(err)
  })
  
  # ---------------------check outlier---------------------
  
  out <- outlier(models[[i]]$residuals)
  # ---------------------fit the model without outlier---------------------
  if (length(out) > 0) {
    phenotypic_data[out, variables[i]] <- NA
    models[[i]] <- asreml(
      fixed = get(variables[i]) ~ loc + set + loc:set + name2 + name2:loc,
      random = ~ loc:block,
      residual = ~ id(loc):ar1(range):ar1(row),
      data = phenotypic_data, 
      na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "name2", sed = TRUE)
    )
  }
  
  # ---------------------storing prediction---------------------
  temp <- models[[i]]$predictions$pvals[, 1:2]
  temp$wave <- variables[i]
  waveblues <- rbind(waveblues, temp)
  cat(i '\n')
}

fwrite(waveblues, "./output/waveblues.csv", row.names= FALSE)
cat('Done with Joint blues')

#blues for each location MW
models <- list()
variables <- colnames(phenotypic_data)[12:2162]
wavebluesEF <- data.frame()

# ---------------------fitting model---------------------
for(i in 1:length(variables)) {
  models[[i]] <- asreml(
  fixed = get(variables[i]) ~ name2 + set,
  random = ~ block,
  residual =  ~ ar1(range):ar1(row),
  data = phenotypic_data, subset = loc== "MW", na.action= na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE))

# ---------------------check outlier---------------------
out <- outlier(models[[i]]$residuals)

# ---------------------fit the model without outlier---------------------
if (length(out) > 0) {
  phenotypic_data[out, variables[i]] <- NA
  models[[i]] <- asreml(
    fixed = get(variables[i]) ~ name2 + set
    random = ~block,
    residual = ~ar1(range):ar1(row),
    data = phenotypic_data,
    na.action = na.method(x = "include"),
    predict = predict.asreml(classify = "name2", sed = TRUE)
  )
}

# ---------------------storing the prediction in an object---------------------

temp1 <- models[[i]]$predictions$pvals[, 1:2]
temp1$wave <- variables[i]
wavebluesMW <- rbind(wavebluesMW, temp1)
}

fwrite(wavebluesMW, "./output/wavebluesMW.csv", row.names=FALSE)

#blues for EF location 
models <- list()
variables <- colnames(phenotypic_data)[12:2162]
wavebluesEF <- data.frame()

# ---------------------fitting model---------------------
for(i in 1:length(variables)) {
  models[[i]] <- asreml(
  fixed = get(variables[i]) ~ name2 + set,
  random = ~ block,
  residual =  ~ ar1(range):ar1(row),
  data = phenotypic_data, subset =loc== "EF", na.action= na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE))

# ---------------------check outlier---------------------
out <- outlier(models[[i]]$residuals)

# ---------------------fit the model without outlier---------------------
if (length(out) > 0) {
  phenotypic_data[out, variables[i]] <- NA
  models[[i]] <- asreml(
    fixed = get(variables[i]) ~ name2 + set
    random = ~block,
    residual = ~ar1(range):ar1(row),
    data = phenotypic_data,
    na.action = na.method(x = "include"),
    predict = predict.asreml(classify = "name2", sed = TRUE)
  )
}
# ---------------------storing the prediction in an object---------------------

temp2 <- models[[i]]$predictions$pvals[, 1:2]
temp2$wave <- variables[i]
wavebluesEF <- rbind(wavebluesEF, temp2)
}

fwrite(wavebluesEF,"./output/wavebluesEF.csv")

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
fwrite(multitrait,"./data/multitrait.csv", row.names = F)

#multitrait<- fread("./data/multitrait.csv", data.table = F)
Names_WEST<- fread("./data/Names_WEST_SF.csv", data.table = F)
#calculating BLUES of narea and sla for  EF location 
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

length(unique(blues$name2))

Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
blues <- blues |> left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>% mutate(taxa = ifelse((name2 != Corrected_names) & grepl("PI", Corrected_names), NA, Corrected_names))

blues<- subset(blues, !is.na(taxa))

#blues_taxa <- droplevels(blues[blues$taxa %in% rownames(kin), ]) 
length(unique(blues$taxa))
blues<- drop_na(blues)

# ---------------------storing blues value for both trait and both location---------------------
fwrite(blues, "./traitsoutput/blues.csv")

# ---------------------reading the blues data file---------------------
blues<- fread("./traitsoutput/blues.csv")

nareabluesef<- subset(blues, loc== "EF") %>% select(-4)
nareabluesmw<- subset(blues, loc== "MW") %>% select(-4)
slabluesef<- subset(blues, loc== "EF") %>% select(-3)
slabluesmw<- subset(blues, loc== "MW") %>% select(-3)


fwrite(nareabluesef,"./traitsoutput/nareabluesef.csv", row.names = F)
fwrite(nareabluesmw,"./traitsoutput/nareabluesmw.csv", row.names = F)

fwrite(slabluesef,"./traitsoutput/slabluesef.csv", row.names = F)
fwrite(slabluesmw,"./traitsoutput/slabluesmw.csv", row.names = F)
