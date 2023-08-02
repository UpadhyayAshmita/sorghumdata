library(here)
library(tidyverse)
library(asreml)
source('./script/outlier.R')
library(data.table)
library(purrr)
library(tidyverse)
library(ASRgenomics)
library(asreml)
library(here)
library(flextable)
library(dplyr)
library(SmartEDA)
library(explore)
library(patchwork)

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

is_check_ef <- table(design[design$LOC == 'EF', ]$Name2) > 1
checks_ef <- names(is_check_ef[is_check_ef == TRUE])


##finding 7 different check varieties that are replicated in all block inside each set in location "EF"

is_check_mw<- table(design[design$LOC == 'MW', ]$Name2) > 1
checks_mw <- names(is_check_mw[is_check_mw == TRUE])


is_check_efsingle <- table(design[design$LOC == 'EF',]$Name2) >16
check_efsingle<- names(is_check_efsingle[is_check_efsingle == TRUE]) #so we can see check var "spx" is replicated 17 times while other 6 check varieties are replicated 16th times one in each block inside four sets in each location "EF" and " MW".







#for location EF

desplot::desplot(
  design %>% filter(LOC == "EF"),
  block ~ range * row,
  out2 = block,
  out1 = Set,
  cex = 0.7,
  ticks = T
)


#for location MW
desplot::desplot(
  design %>% filter(LOC == "MW"),
  block ~ range * row,
  out2 = block,
  out1 = Set,
  cex = 0.7,
  ticks = T, text= Name2,
)




table(design[(design$LOC== "EF") & (design$Set == "Q1")& (design$block== 1), ])





#fitting the model 

# y= mu + genotype+ location+ set+ block:set + e



designnew<- transform(designnew,
                       Name2= factor(Name2),
                       LOC= factor(LOC),
                       Set= factor(Set),
                       block= factor(block),
                       range= factor(range),
                       row= factor(row))

write.csv(designnew, './data/designnew.csv')

designnew<- read.csv('./data/designnew.csv')

models <- list()
variables <- colnames(designnew)[11:2161]
h2 <- data.frame(matrix(NA, length(variables), 2))  #n rows and 2 columns
h2[,1] <- variables
colnames(h2) <- c("wave", "h2")
rownames(h2) <- variables

for(i in variables) {
  models[[i]] <- asreml(
    fixed = get(i) ~ 1 + Set+ LOC,
    random = ~Name2 + LOC:block + LOC:Name2,
    #residual =  ~ corh(rep):id(range):(row),
    data = designnew, na.action = na.method(x = "include"),
    predict = predict.asreml(classify = "Name2", sed = TRUE)
  )
 
  out <- outlier(models[[i]]$residuals)
  
   if(length(out) > 0) {
      designnew[,i][designnew$LOC == "EF"][out] <- NA
     models[[i]] <- asreml(
       fixed = get(i) ~ 1 + Set + LOC,
       random =  ~ Name2 + LOC:block + LOC:Name2,
       #residual =  ~ id(rep):id(range):(row),
       data = designnew, na.action = na.method(x = "include"),
       predict = predict.asreml(classify = "Name2", sed = TRUE)
   )}
  
  h2[i,2]  <- (1 - ((models[[i]]$predictions$avsed["mean"] ^ 2) /(2 * summary(models[[i]])$varcomp["Name2", "component"])))
}

h2<- fwrite(h2, "h2.csv",row.names = F)
h2<- fread("./output/h2.csv")
summary(h2)



#calculating blues for each wavelength 

models <- list()
variables <- colnames(designnew)[11:2161]
waveblues <- data.frame()  #n rows and 2 columns
# waveblues[,1] <- variables
# colnames(waveblues) <- c("wave", "blues")
# rownames(waveblues) <- variables

for(i in 1:length(variables)) {
  models[[i]] <- asreml(
    fixed = get(variables[i]) ~ Name2,
    random = ~ LOC:block + Set,
    #residual =  ~ corh(rep):id(range):(row),
    data = designnew, na.action = na.method(x = "include"),
    predict = predict.asreml(classify = "Name2", sed = TRUE))
  temp <- models[[i]]$predictions$pvals[, 1:2]
  temp$wave <- variables[i]
  waveblues <- rbind(waveblues, temp)
}
 

frwite(waveblues, "waveblues.csv")

waveblues<- fread("./output/waveblues.csv")


blueswave<- pivot_wider(waveblues, 
                          id_cols = Name2, 
                          names_from = wave, 
                          values_from = predicted.value)


fwrite(blueswave, "blueswave.csv")
blueswave<- fread("./output/blueswave.csv")



#creating relationship matrix for whole spectra
matrix.wholewave <- as.matrix(df_prediction)
rownames(matrix.wholewave) <- matrix.wholewave[, 1]
matrix.wholewave <- matrix.wholewave[, -1]
class(matrix.wholewave) <- 'numeric'
 
#creating relationship matrix by multiplying matrix with its transpose
relationshipmatrix <- matrix.wholewave %*% t(matrix.wholewave) / (ncol(df_prediction) - 1)
dim(relationshipmatrix)
hist(as.vector(relationshipmatrix))
write.csv(relationshipmatrix,'whole.relation.matrix.csv')
whole.relation.matrix.csv<- read.csv('./output/whole.relation.matrix.csv')

#creating relationship matrix for nirs data

library(dplyr)
df_nirs<- df_prediction %>% select(Name2,"Wave_392":"Wave_850")

matrix_nirs<- as.matrix(df_nirs)
rownames(matrix_nirs)<- matrix_nirs[,1]
matrix_nirs<- matrix_nirs[,-1]
class(matrix_nirs)<- 'numeric'

#nirsmatrix* nirsmatrix(transpose/(total 392-850)wavelength))

relation.mat.nirs<- matrix_nirs %*% t(matrix_nirs) / (ncol(df_nirs)-1)
dim(relation.mat.nirs)
hist(as.vector(relation.mat.nirs))
write.csv(relation.mat.nirs,'NIRS.R.matrix.csv')
NIRS.R.matrix.csv<- read.csv('./output/NIRS.R.matrix.csv')

#for selecting genotype with higher heritability
h2<-read.csv('./output/h2.csv')


library(tidyverse)

h2<- h2 %>% mutate(wavenumber= as.numeric(str_replace(X1,"Wave_", "")))

library(ggplot2)
ggplot(data= h2, aes(x = wavenumber, y = X2)) + geom_point(color= "purple", size= 0.1) + labs(title= "heritability per wavelength", x= "wavenumber",y= "heritability")



ggplot(data = h2[h2$X2 > 0.4, ], aes(x = wave, y = h2)) + 
  geom_point(color= "purple", size= 0.1) + 
  labs(title= "heritability per wavelength", x= "wavenumber",y= "heritability")

h2_maxwavelength<- h2[h2$X2 > 0.4712, c("X1", "X2")]





#reading phenotypic data
phenotypedata<- read.csv('./data/phenotypes.csv')

phenotypedata<- transform(phenotypedata,
                          TAXA= factor(TAXA),
                          LOC= factor(LOC),
                          Name2= factor(Name2),
                          Set= factor(Set),
                          block= factor(block),
                          range= factor(range),
                          row= factor(row))


#some EDA
head(phenotypedata)
summary(phenotypedata)
table(phenotypedata$TAXA,phenotypedata$Set) #checks are replicated 8 times in each four set while other genotypes are replicated only 2 times in each set

round(cor(phenotypedata[,10:11],use="pairwise.complete.obs"),5)# Narea and SLA is negatively correlated i.e -0.6298


#desplot for location "EF"
desplot::desplot(
  phenotypedata%>% filter(LOC == "EF"),
  block ~ range * row,
  out2 = block,
  out1 = Set,
  cex = 0.7,
  ticks = T
)

#desplot for location "MW"
desplot::desplot(
  phenotypedata%>% filter(LOC == "MW"),
  block ~ range * row,
  out2 = block,
  out1 = Set,
  cex = 0.7,
  ticks = T
)

head(phenotypedata)
str(phenotypedata)


#exploring more through visual plot 
ExpNumStat(phenotypedata, by="A", round= 2) %>%  
  flextable()

phenotypedata%>% explore()

#calculating BLUES for  EF location 
models <- list()
bluesEF <- list()
a <- 1
for(i in c("Narea","SLA")){
  models[[i]] <- asreml(
    fixed = get(i) ~ Name2,
    random =  ~  Set+ block,
    #residual =  ~ id(range):id(row),
    data = phenotypedata, subset = LOC == "EF", na.action = na.method(x = c("include")),
    predict = predict.asreml(classify = "Name2",sed = TRUE))
  bluesEF[[i]] <-   models[[i]]$predictions$pvals[,1:2]
  colnames(bluesEF[[i]]) <- c("Name2", i)
  a <- a + 1
}


#calculating BLUES for MW location
models <- list()
bluesMW <- list()
a <- 1
for(i in c("Narea","SLA")){
  models[[i]] <- asreml(
    fixed = get(i) ~ Name2,
    random =  ~  Set + block,
    #residual =  ~ id(range):id(row),
    data = phenotypedata, subset = LOC == "MW", na.action = na.method(x = c("include")),
    predict = predict.asreml(classify = "Name2",sed = TRUE))
  bluesMW[[i]] <-   models[[i]]$predictions$pvals[,1:2]0  
  colnames(bluesMW[[i]]) <- c("Name2", i)
  a <- a + 1
}

#---------combining"bluesEF" and "bluesMW" based on a common "Name2" column and creating new "LOC" column.Then left joining the blues data with "phenotypedata" adding info of "Name2",
#filtering to include only rows with "TAXA" present in the "kin" matrix, grouped the data by "TAXA",obtained the mean for each column (excluding "Name2" and "LOC"), 
#filtering out rows with Na's values in the "Narea" column, obtained final "blues" data frame with the mean "blues" values for each unique "TAXA".


blues <-
  bind_rows(
    "EF" = bluesEF %>% reduce(left_join, by = "Name2"),
    "MW" = bluesMW %>% reduce(left_join, by = "Name2"),
    .id = "LOC"
  ) %>% left_join(phenotypedata %>% filter(!duplicated(Name2)) %>% select(Name2, TAXA))
blues <- droplevels(blues[blues$TAXA %in%rownames(kin), ])

blues <- blues %>%
  group_by(TAXA) %>% select(-Name2, -LOC) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE))) %>%
  filter(!is.nan(Narea))

fwrite(blues, "blues.txt", quote = F, sep = "\t")


#-----reading the blues data file-----

blues<- fread("./output/blues.txt")


# blues for Narea trait
Nitroblues <- subset(blues, select = c("TAXA", "Narea"))

# blues for SLA trait
SLAblues <- subset(blues, select = c("TAXA", "SLA"))



# alternative model for first stage 
# separating phenotypedata for two different data 
# 
#  Assuming your original data frame is named "df"
# 
# Create a data frame for location EF
# phenoEF <- subset(phenotypedata, LOC == "EF")
# fwrite(phenoEF, "phenoEF.csv")
# 
# Create a data frame for location MW
# phenoMW<- subset(phenotypedata, LOC == "MW")
# fwrite(phenoMW, "phenoMW.csv")
# 
# reading two location phenotypedata
# phenoEF<- fread("./data/phenoEF.csv")
# phenoMW<- fread("./data/phenoMW.csv")
# 
# 
# ---------------- Creating a vector of ascending and descending sequences
# ascending_values <- 1:40
# 
#  -------------------Calculating the number of rows in each ascending sequence
# rows_per_sequence <- 40
# 
#  -------------------Determine the row indices to apply ascending order
# ascending_rows_indices <- c(41:80, 121:160, 201:240, 281:320, 361:400, 441:480, 521:560, 601:640,681:720, 761:800, 841:880, 921:960)
# 
# ----------Creating a function to repeat ascending order values and keeping the rest unchanged
# apply_ascending_order <- function(ascend) {
#   ascending_values_rep <- rep(ascending_values, length.out = length(ascend))
#   ascend[ascending_rows_indices] <- ascending_values_rep[ascending_rows_indices]
#   return(ascend)
# }
# 
#  ------------------Applying ascending order on specified rows in the "row" column
# phenoEF$row<- apply_ascending_order(phenoEF$row)
# phenoMW$row<- apply_ascending_order(phenoMW$row)
# 
# fwrite(phenoEF, "./data/phenotypeEF.csv")
# fwrite(phenoMW, "./data/phenotypeMW.csv")



#reading each location phenotype data
phenotypeEF<- fread("./data/phenotypeEF.csv")
phenotypeMW<- fread("./data/phenotypeMW.csv")

#changing numerical variable into factor for the model
phenotypeEF<- transform(phenoEF,
                    TAXA= factor(TAXA),
                    LOC= factor(LOC),
                    Name2= factor(Name2),
                    Set= factor(Set),
                    block= factor(block),
                    range= factor(range),
                    row= factor(row))



phenotypeMW<- transform(phenoMW,
                    TAXA= factor(TAXA),
                    LOC= factor(LOC),
                    Name2= factor(Name2),
                    Set= factor(Set),
                    block= factor(block),
                    range= factor(range),
                    row= factor(row))
hist(phenotypeEF$SLA)

#---- Alternately comparing the identity,ar1 and diagonal model for the blues of Narea and SLA to find the best model than the general identity residual model

#----for Narea in location EF
modeliden <- asreml(
  fixed = Narea ~ Name2,
  random =  ~  Set + block,
  residual =  ~ id(range):id(row),
  data = phenotypeEF, na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "Name2",sed = TRUE))

summary(modeliden)$aic
plot(modeliden)
modeliden$loglik


blueiden<-summary(modeliden,coef=TRUE)$coef.fixed
head(blueiden)


modelar<- asreml(fixed = Narea ~ Name2,
                 random =  ~ Set + block,
                 residual = ~ar1(range) :ar1(row),
                 data = phenotypeEF, na.action = na.method(x = c("include")),
                 predict = predict.asreml(classify = "Name2",sed = TRUE))
summary(modelar)$aic
plot(modelar)
bluesar<-summary(modelar,coef=TRUE)$coef.fixed
head(bluesar)
modelar$loglik




#--for SLA in LOC EF

modeliden1 <- asreml(
  fixed = SLA ~ Name2,
  random =  ~  Set + block,
  residual =  ~ id(range):id(row),
  data = phenotypeEF, na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "Name2",sed = TRUE))

summary(modeliden1)$aic
plot(modeliden1)
modeliden1$loglik


blueiden1<-summary(modeliden1,coef=TRUE)$coef.fixed
head(blueiden1)


modelar1<- asreml(fixed = SLA ~ Name2,
                 random =  ~ Set + block,
                 residual = ~ar1(range) :ar1(row),
                 data = phenotypeEF, na.action = na.method(x = c("include")),
                 predict = predict.asreml(classify = "Name2",sed = TRUE))
summary(modelar1)$aic
plot(modelar1)
bluesar1<-summary(modelar1,coef=TRUE)$coef.fixed
head(bluesar1)
modelar1$loglik




#------------for location MW
#--for Narea


modeliden2 <- asreml(
  fixed = Narea ~ Name2,
  random =  ~  Set + block,
  residual =  ~ id(range):id(row),
  data = phenotypeMW, na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "Name2",sed = TRUE))

summary(modeliden2)$aic
plot(modeliden2)


modelar2<- asreml(fixed = Narea ~ Name2,
                  random =  ~ Set + block,
                  residual = ~ar1(range) :ar1(row),
                  data = phenotypeMW, na.action = na.method(x = c("include")),
                  predict = predict.asreml(classify = "Name2",sed = TRUE))
summary(modelar2)$aic
plot(modelar2)


#-----for SLA
modeliden3 <- asreml(
  fixed = SLA ~ Name2,
  random =  ~  Set + block,
  residual =  ~ id(range):id(row),
  data = phenotypeMW, na.action = na.method(x = c("include")),
  predict = predict.asreml(classify = "Name2",sed = TRUE))

summary(modeliden3)$aic
plot(modeliden3)


modelar3<- asreml(fixed = SLA ~ Name2,
                  random =  ~ Set + block,
                  residual = ~ar1(range) :ar1(row),
                  data = phenotypeMW, na.action = na.method(x = c("include")),
                  predict = predict.asreml(classify = "Name2",sed = TRUE))
summary(modelar3)$aic
plot(modelar3)

#Aic and Bic of identity residual model was the lowest than ar1 and diag so using the blues obtained from previous for loop for second stage analysis
















