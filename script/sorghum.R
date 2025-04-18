
library(tidyverse)
library(asreml)
source('./script/outlier.R')
library(data.table)
library(ASRgenomics)
library(flextable)
library(dplyr)
library(SmartEDA)
library(explore)
library(patchwork)
library(janitor)

#reading field design and wavelength data and filtering year 16 ("wavelength data available for year 16)

design<- read.csv('./data/design.csv') %>% filter(year!= 16)
averagedspectra <- read.csv('./data/Averaged_Spectra.csv')

averagedspectra <- averagedspectra %>%
  separate_wider_delim(
    Spectra,
    delim = "_",
    names = c("Spectra", "plotnum"),
    cols_remove = T
  ) %>%
  select(-plotnum) %>%
  group_by(Spectra) %>%
  summarise_all(mean) %>%
  ungroup() %>%
  rename(PlotID = Spectra)

designnew<- left_join(design,averagedspectra, by = "PlotID") # Apply left_join dplyr function

#finding 7 different check varieties that are replicated in all block inside each set in location "EF"

is_check_ef <- table(design[design$loc == 'EF', ]$name2) > 1
checks_ef <- names(is_check_ef[is_check_ef == TRUE])


##finding 7 different check varieties that are replicated in all block inside each set in location "MW"

is_check_mw<- table(design[design$loc == 'MW', ]$name2) > 1
checks_mw <- names(is_check_mw[is_check_mw == TRUE])


is_check_efsingle <- table(design[design$loc == 'EF',]$name2) >16
check_efsingle<- names(is_check_efsingle[is_check_efsingle == TRUE]) #so we can see check var "spx" is replicated 17 times while other 6 check varieties are replicated 16th times one in each block inside four sets in each location "EF" and " MW".

#write.csv(designnew, './data/designnew.csv')

designnew<- fread('./data/designnew.csv', data.table = FALSE)


#for location EF

desplot::desplot(
  design%>% filter(loc == "EF"),
  block ~ range * row,
  out2 = block,
  out1 = set,
  cex = 0.7,
  ticks = T
)


#for location MW
desplot::desplot(
  design %>% filter(loc == "MW"),
  block ~ range * row,
  out2 = block,
  out1 = set,
  cex = 0.7,
  ticks = T, text= name2,
)


asreml.options(
  workspace = '8gb',
  pworkspace = '8gb'
)



# ---------------------load data---------------------
designnew<- fread('./data/designnew.csv',data.table= FALSE)



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





# ---------------------calculating heritability of each wavelength--------------------

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
    designnew[,wavelength][designnew$loc][out] <- NA
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

#fwrite(h2, "./output/h2.csv",row.names = F)

h2<- fread("./output/h2.csv", data.table = FALSE)



#calculating blues for joint location for each wavelength 
#calculating blues for each wavelength 
models <- list()
variables <- colnames(designnew)[11:2161]
waveblues <- data.frame()

# ---------------------fitting model---------------------

for (i in 1:length(variables)) {
  cat(variables[i], '\n')
  
  models[[i]] <- tryCatch({
    asreml(
      fixed = get(variables[i]) ~ loc + set + loc:set + name2 + name2:loc,
      random = ~loc:block,
      residual =  ~ id(loc):ar1(range):ar1(row),
      data = designnew, na.action = na.method(x = "include"),
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
    designnew[out, variables[i]] <- NA
    models[[i]] <- asreml(
      fixed = get(variables[i]) ~ loc + set + loc:set + name2 + name2:loc,
      random = ~ loc:block,
      residual = ~ id(loc):ar1(range):ar1(row),
      data = designnew, 
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

fwrite(waveblues, "./output/waveblues.csv", row.names= FALSE)
cat('Done!')


#blues for each location MW

models <- list()
variables <- colnames(designnew)[11:2161]
wavebluesEF <- data.frame()

# ---------------------fitting model---------------------

  for(i in 1:length(variables)) {models[[i]] <- asreml(
    fixed = get(variables[i]) ~ name2 + set,
    random = ~ block,
    residual =  ~ ar1(range):ar1(row),
    data = designnew, subset = loc== "MW", na.action= na.method(x = "include"),
    predict = predict.asreml(classify = "name2", sed = TRUE))
  
  
  # ---------------------check outlier---------------------
  out <- outlier(models[[i]]$residuals)
  
  # ---------------------fit the model without outlier---------------------
  if (length(out) > 0) {
    designnew[out, variables[i]] <- NA
    models[[i]] <- asreml(
      fixed = get(variables[i]) ~ name2 + set
      random = ~block,
      residual = ~ar1(range):ar1(row),
      data = designnew,
      na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "name2", sed = TRUE)
    )
  }
  
  # ---------------------storing the prediction in an object---------------------
  
  temp1 <- models[[i]]$predictions$pvals[, 1:2]
  temp1$wave <- variables[i]
  wavebluesMW <- rbind(wavebluesMW, temp1)
  }

#fwrite( wavebluesMW"./output/wavebluesMW.csv", row.names=FALSE)

fread("./output/wavebluesMW.csv")


#blues for EF location 
models <- list()
variables <- colnames(designnew)[11:2161]
wavebluesEF <- data.frame()

# ---------------------fitting model---------------------

  for(i in 1:length(variables)) {models[[i]] <- asreml(
    fixed = get(variables[i]) ~ name2 + set,
    random = ~ block,
    residual =  ~ ar1(range):ar1(row),
    data = designnew, subset =loc== "EF", na.action= na.method(x = "include"),
    predict = predict.asreml(classify = "name2", sed = TRUE))
  
  
  # ---------------------check outlier---------------------
  out <- outlier(models[[i]]$residuals)

  # ---------------------fit the model without outlier---------------------
  if (length(out) > 0) {
    designnew[out, variables[i]] <- NA
    models[[i]] <- asreml(
      fixed = get(variables[i]) ~ name2 + set
      random = ~block,
      residual = ~ar1(range):ar1(row),
      data = designnew,
      na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "name2", sed = TRUE)
    )
  }
  # ---------------------storing the prediction in an object---------------------
  
  temp2 <- models[[i]]$predictions$pvals[, 1:2]
  temp2$wave <- variables[i]
  wavebluesEF <- rbind(wavebluesEF, temp2)
  }

#fwrite(wavebluesEF,"./output/wavebluesEF.csv")

#creating relationship matrix for each location and joint location from blues obtained for wavelength

#generating relationship matrix for whole spectra
  matrix.wholewave <- as.matrix(blueswave)
  rownames(matrix.wholewave) <- matrix.wholewave[, 1]
  matrix.wholewave <- matrix.wholewave[, -1]
  class(matrix.wholewave) <- 'numeric'

#creating relationship matrix by multiplying matrix with its transpose

  relationshipmatrix <- matrix.wholewave %*% t(matrix.wholewave) / (ncol(blueswave) - 1)
  dim(relationshipmatrix)
  hist(as.vector(relationshipmatrix))
  
  write.csv(relationshipmatrix,'wholewaverel.csv')
  whole.relation.matrix.csv<- read.csv('./output/wholewaverel.csv')
  
  #creating relationship matrix for nirs data
  
  dfnirs<- blueswave %>% select(name2,"Wave_392":"Wave_850")
  matrix_nirs<- as.matrix(dfnirs)
  rownames(matrix_nirs)<- matrix_nirs[, 1]
  matrix_nirs<- matrix_nirs[,-1]
  class(matrix_nirs)<- 'numeric'


#nirsmatrix* nirsmatrix(transpose/(total 392-850)wavelength))

nirswaverel<- matrix_nirs %*% t(matrix_nirs) / (ncol(dfnirs)-1)
dim(nirswaverel)
hist(as.vector(nirswaverel))

write.csv(nirswaverel,'nirswaverel.csv')
nirswaverel<- read.csv('./output/nirswaverel.csv')

#for selecting genotype with higher heritability
h2<-read.csv('./output/h2.csv') #reading h2 data



#visualizing trend in plot with wavelegth vs heritability 

h2<- h2 %>% mutate(wavenumber= as.numeric(str_replace(wave,"wave_", "")))
ggplot(data= h2, aes(x = wavenumber, y = h2)) + geom_point(color= "purple", size= 0.1) + labs(title= "heritability per wavelength", x= "wavenumber",y= "heritability")
ggplot(data = h2[h2$h2 > 0.4, ], aes(x = wavenumber, y = h2)) + geom_point(color= "purple", size= 0.1) + labs(title= "heritability per wavelength", x= "wavenumber",y= "heritability")
highh2wave<- h2[h2$h2 > 0.4712, c("wavenumber", "h2")]

fwrite(highh2wave,"./output/highh2wave.csv")
highh2wave<- fread("./output/highh2wave.csv")


#creating relationship matrix for wavelength with high heritability 
highh2matrix<- as.matrix(highh2wave)
rownames(highh2matrix) <- highh2matrix[, 1]
highh2matrix <- highh2matrix[, -1]
class(highh2matrix) <- 'numeric'

#creating relationship matrix by multiplying matrix with its transpose

highh2_relationmat <- highh2matrix %*% t(highh2matrix) / (ncol(highh2wave) - 1)
dim(highh2_relationmat)
hist(as.vector(highh2_relationmat))


write.csv(highh2_relationmat,'highh2_relationmat.csv')
highh2_relationmat<- read.csv("./output/highh2_relationmat.csv", row.names= "T")

#reading phenotypic data
phenotypedata<- read.csv('./data/phenotypes.csv')
phenotypedata<-
  phenotypedata %>% clean_names()%>% mutate(
    name2 = factor(name2),
    taxa = factor(taxa),
    loc = factor(loc),
    set = factor(set),
    block = factor(block),
    range = factor(range),
    row = factor(row),
    uni = c(1:960, 1:960)
  ) %>%   arrange(loc, range, row)


#some EDA

head(phenotypedata)

summary(phenotypedata)
table(phenotypedata$taxa,phenotypedata$set) #checks are replicated 8 times in each four set while other genotypes are replicated only 2 times in each set
round(cor(phenotypedata[,10:11],use="pairwise.complete.obs"),5)# narea and sla is negatively correlated i.e -0.6298

#desplot for location "EF"


desplot::desplot(
  phenotypedata%>% filter(loc == "EF"),
  block ~ range * row,
  out2 = block,
  out1 = set,
  cex = 0.7,
  ticks = T
)


#desplot for location "MW"
desplot::desplot(phenotypedata%>% filter(loc == "MW"),
  block ~ range * row,
  out2 = block,
  out1 = set,
  cex = 0.7,
  ticks = T
)

head(phenotypedata)
str(phenotypedata)

#exploring more through visual plot 

ExpNumStat(waveblues, by="A", round= 2) %>%  flextable()

designnew%>% explore()



# ---------------------reading phenotypic data for trait of interest---------------------



phenotypedata<- read.csv('./data/phenotypes.csv')

# ---------------------preprocessing of data---------------------

phenotypedata<-
  phenotypedata %>% clean_names()%>% mutate(
    name2 = factor(name2),
    taxa = factor(taxa),
    loc = factor(loc),
    set = factor(set),
    block = factor(block),
    range = factor(range),
    row = factor(row),
    uni = c(1:960, 1:960)
  ) %>%   arrange(loc, range, row)

# -------------------calculating heritability for narea for each location---------------------

# ---------------------fitting model---------------------

nareaMW <- asreml(
  fixed = narea ~ set,
  random = ~name2 + block,
  residual =  ~ar1(range):ar1(row),
  data = phenotypedata, subset= loc== "MW",na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "name2", sed = TRUE)
)

 nareaMW<- update.asreml(nareaMW)

  # ---------------------checking outliers---------------------
  
  
  out <- outlier(nareaMW$residuals)
  
  if(length(out) > 0) {
    phenotypedata[out, "narea"] <- NA
    
    } else{
    
    nareaMW <- asreml(
      fixed = narea ~ set,
      random = ~name2 + block,
      residual =  ~ar1(range):ar1(row),
      data = phenotypedata, subset= loc== "MW",na.action = na.method(x = "include"),
      predict = predict.asreml(classify = "name2", sed = TRUE)
    )
    nareaMW<- update.asreml(nareaMW)
    }
 
 
 # ---------------------calculating and storing the heritability---------------------
 
  temp <- (1 - ((nareaMW$predictions$avsed["mean"] ^ 2) /(2 * summary(nareaMW)$varcomp["name2", "component"])))
 
 # ---------------------h2 was found 0.587 for narea for location mw---------------------
 
# ---------------------storing outlier filtered data for location mw---------------------
 
#fwrite(phenotypedata, "./data/phenotypedataMW.csv", row.names= FALSE)

 
 
# ---------------------reading data---------------------
 
phenotypedataMW<-  fread("./data/phenotypedataMW.csv")

 # ---------------------checking residual plot for the model---------------------
  
plot(h2narea)
 
 
 
 
# ---------------------calculating heritability for narea for ef location---------------------

# ---------------------fitting model---------------------

nareaEF <- asreml(
   fixed = narea ~ set,
   random = ~name2 + block,
   residual =  ~ar1(range):ar1(row),
   data = phenotypedata, subset= loc== "EF",na.action = na.method(x = "include"),
   predict = predict.asreml(classify = "name2", sed = TRUE)
 )
 
nareaEF<- update.asreml(nareaEF)
 
 # ---------------------checking outliers---------------------
 
 
 out <- outlier(nareaEF$residuals)
 
 if(length(out) > 0) {
   phenotypedata[out, "narea"] <- NA
   }else {
   
  nareaEF <- asreml(
     fixed = narea ~ set,
     random = ~name2 + block,
     residual =  ~ar1(range):ar1(row),
     data = phenotypedata, subset= loc== "EF",na.action = na.method(x = "include"),
     predict = predict.asreml(classify = "name2", sed = TRUE)
   )
   nareaEF<- update.asreml(nareaEF)
   }

# ---------------------calculating and storing heritability---------------------

 temp2<- (1 - ((nareaEF$predictions$avsed["mean"] ^ 2) /(2 * summary(nareaEF)$varcomp["name2", "component"])))

# ---------------------saving outliers filtered data for location ef for narea---------------------

 fwrite(phenotypedata, "./data/phenotypedataEF.csv", row.names= FALSE)


# ---------------------h2 for narea for loc ef was 0.355---------------------


# ---------------------checking residual plot---------------------

  plot(nareaEF)





# ---------------------calculating heritability for sla for MW location---------------------
 

# ---------------------fitting model---------------------

 slaMW <- asreml(
   fixed = sla ~ set,
   random = ~name2 + block,
   residual =  ~ar1(range):ar1(row),
   data = phenotypedata, subset= loc== "MW",na.action = na.method(x = "include"),
   predict = predict.asreml(classify = "name2", sed = TRUE)
 )
 
 slaMW<- update.asreml(slaMW)
 
 # ---------------------checking outliers---------------------
 
 
 out <- outlier(slaMW$residuals)
 
 if(length(out) > 0) {
   phenotypedata[out, "sla"] <- NA
   }else{
   
   slaMW <- asreml(
     fixed = sla ~ set,
     random = ~name2 + block,
     residual =  ~ar1(range):ar1(row),
     data = phenotypedata, subset= loc== "MW",na.action = na.method(x = "include"),
     predict = predict.asreml(classify = "name2", sed = TRUE)
   )
   slaMW<- update.asreml(slaMW)
   }
 
 # ---------------------calculating and storing heritability---------------------
 
 temp3 <- (1 - ((slaMW$predictions$avsed["mean"] ^ 2) /(2 * summary(slaMW)$varcomp["name2", "component"])))
 
 # ---------------------storing outliers filtered data for sla for loc mw---------------------
 
 fwrite(phenotypedata, "./data/phenotypedataslaMW.csv", row.names= FALSE)
 
 
 
 # ---------------------checking residual plot---------------------
 
 plot(slaMW)
 
# ---------------------calculating heritability for sla for EF location ---------------------

 
# ---------------------fitting model---------------------

 slaEF <- asreml(
   fixed = sla ~ set,
   random = ~name2 + block,
   residual =  ~ar1(range):ar1(row),
   data = phenotypedata, subset= loc== "EF",na.action = na.method(x = "include"),
   predict = predict.asreml(classify = "name2", sed = TRUE)
 )
 
 slaEF<- update.asreml(slaEF)
 
 # ---------------------checking outliers---------------------
 
 
 out <- outlier(slaMW$residuals)
 
 if(length(out) > 0) {
   phenotypedata[out, "sla"] <- NA
 } else{
   slaEF <- asreml(
     fixed = sla ~ set,
     random = ~name2 + block,
     residual =  ~ar1(range):ar1(row),
     data = phenotypedata, subset= loc== "EF",na.action = na.method(x = "include"),
     predict = predict.asreml(classify = "name2", sed = TRUE)
   )
   slaEF<- update.asreml(slaEF)
 }
 
 
 # ---------------------calculating and storing h2 for sla for loc ef---------------------
 
 temp4 <- (1 - ((slaEF$predictions$avsed["mean"] ^ 2) /(2 * summary(slaEF)$varcomp["name2", "component"])))
 
 
 
 # ---------------------storing outliers filtered data for sla for loc ef---------------------
 
 fwrite(phenotypedata, "./data/phenotypedataslaEF.csv", row.names= FALSE)
 
 
 # ---------------------residual plot---------------------
 
 plot(slaEF)
 
 
 
 
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
 narea<- fread("./data/narea.csv", data.table= F)
 
 # ---------------------preprocessing of data---------------------
 
 narea<-
   narea %>% clean_names()%>% mutate(
     name2 = factor(name2),
     taxa = factor(taxa),
     loc = factor(loc),
     set = factor(set),
     block = factor(block),
     range = factor(range),
     row = factor(row),
     uni = c(1:960, 1:960)
   ) %>%   arrange(loc, range, row)
 
 
 
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
 sla<- fread("./data/sla.csv", data.table = F)
 
 
 # ---------------------preprocessing of data---------------------
 
 sla<-
   sla %>% clean_names()%>% mutate(
     name2 = factor(name2),
     taxa = factor(taxa),
     loc = factor(loc),
     set = factor(set),
     block = factor(block),
     range = factor(range),
     row = factor(row),
     uni = c(1:960, 1:960)
   ) %>%   arrange(loc, range, row)
 
 
 
 
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

 
# ---------------------finding best model---------------------

 
narea<- subset(narea, select= -sla)
sla<- subset(sla, select= -narea)


 
multitrait<- multitrait %>% clean_names()%>% mutate(
  name2 = factor(name2),
  taxa = factor(taxa),
  loc = factor(loc),
  set = factor(set),
  block = factor(block),
  range = factor(range),
  row = factor(row),
  uni = c(1:960, 1:960)
) %>%   arrange(loc, range, row)
 
 

#blues for location EF for trait narea 
 
 
modar1<-  asreml(fixed = narea ~ name2 +set,
                           random =  ~ block,
                           residual = ~ar1(range):ar1(row),
                           data = narea, subset= loc== "EF",na.action = na.method(x = c("include")),
                           predict = predict.asreml(classify = "name2",sed = TRUE))
modar1<- update.asreml(modar1)
 
summary(modar1)$aic
plot(modar1)
 
 
 modeliden <- asreml(
   fixed = narea ~ name2+ set,
   random =  ~block,
   residual =  ~ id(range):id(row),
   data = narea,subset= loc== "EF", na.action = na.method(x = c("include")),
   predict = predict.asreml(classify = "name2",sed = TRUE))
 
 modeliden<- update.asreml(modeliden)
 
 summary(modeliden)$aic
 plot(modeliden)

 
 
 # ---------------------binding both trait both location filtered data---------------------
 
 narea<- subset(narea, select= -sla)
 sla<- subset(sla, select= -narea)
 
 multitrait<- full_join(narea, sla)
 #fwrite(multitrait,"./data/multitrait.csv", row.names = F)
 
 
 # ---------------------reading data---------------------
 library(janitor)
 library(dplyr)
 multitrait<- fread("./data/multitrait.csv", data.table = F)
 Names_WEST<- fread("./data/Names_WEST_SF.csv", data.table = F)
 
 length(unique(multitrait$taxa))
  length(unique(Names_WEST$Name2))
 
 
 
 # ---------------------preprocessing data---------------------
 
 
 multitrait<- multitrait %>% clean_names()%>% mutate(
   name2 = factor(name2),
   taxa = factor(taxa),
   loc = factor(loc),
   set = factor(set),
   block = factor(block),
   range = factor(range),
   row = factor(row),
   uni = c(1:960, 1:960)
 ) %>%   arrange(loc, range, row)
 
 
 
 
 
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
  
  
  
#---------combining"bluesEF" and "bluesMW" based on a common "name2" column and creating new "loc" column.Then left joining the blues data with "phenotypedata" adding info of "name2",
#filtering to include only rows with "taxa" present in the "kin" matrix, grouped the data by "taxa",obtained the mean for each column (excluding "name2" and "loc"), 
#filtering out rows with Na's values in the "narea" column, obtained final "blues" data frame with the mean "blues" values for each unique "taxa".


blues <- 
  bind_rows(
  "EF" = bluesEF %>% reduce(left_join, by = "name2"),
  "MW" = bluesMW %>% reduce(left_join, by = "name2"),
  .id = "loc") %>% left_join(multitrait %>% filter(!duplicated(name2)) %>% select(name2))
    
length(unique(blues$name2))

Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
blues <- blues |> left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>% mutate(taxa = ifelse((name2 != Corrected_names) & grepl("PI", Corrected_names), NA, Corrected_names))



blues<- droplevels(blues[blues$taxa %in% rownames(kin), ])
 



length(unique(blues_taxa$taxa))




# ---------------------storing blues value for both trait and both location---------------------


  fwrite(blues, "./traitsoutput/blues.csv")


# ---------------------reading the blues data file---------------------


blues<- fread("./traitsoutput/blues.csv")

nareabluesef<- subset(blues, loc== "EF") %>% select(-3)
nareabluesmw<- subset(blues, loc== "MW") %>% select(-3)
slabluesef<- subset(blues, loc== "EF") %>% select(-2)
slabluesmw<- subset(blues, loc== "MW") %>% select(-2)


fwrite(nareabluesef,"./traitsoutput/nareabluesef.csv", row.names = F)
fwrite(nareabluesmw,"./traitsoutput/nareabluesmw.csv", row.names = F)

fwrite(slabluesef,"./traitsoutput/slabluesef.csv", row.names = F)
fwrite(slabluesmw,"./traitsoutput/slabluesmw.csv", row.names = F)

# alternative model for first stage 
#Aic and Bic of identity residual model was the lowest than ar1 and diag so using the blues obtained from previous for loop for second stage analysis
  


  