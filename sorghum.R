library(here)
#install ASREML
install.packages("../asreml-4.1.0.176-CentOS-7-R4.2.tar.gz", repos = NULL)
# activate ASREML license
library(asreml)
asreml.license.activate()


source('outlier.R')


library(tidyverse)
design<- read.csv('design.csv') %>% filter(year!= 16)

averagedspectra <- read.csv('Averaged_Spectra.csv')
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

write.csv(designnew, 'designnew.csv')

designnew<- read.csv('designnew.csv')

models <- list()
variables <- colnames(designnew)[11:2161]
h2 <- data.frame(matrix(NA, length(variables), 2))  #n rows and 2 columns
h2[,1] <- variables
colnames(h2) <- c("wave", "h2")
rownames(h2) <- variables
for(i in variables) {
  models[[i]] <- asreml(
    fixed = get(i) ~ 1 + Set+ LOC,
    random = ~Name2 + Set:block + LOC:Name2,
    #residual =  ~ corh(rep):id(range):(row),
    data = designnew, na.action = na.method(x = "include"),
    predict = predict.asreml(classify = "Name2", sed = TRUE)
  )
 
  out <- outlier(models[[i]]$residuals)
  
   if(length(out) > 0) {
      designnew[,i][designnew$LOC == "EF"][out] <- NA
     models[[i]] <- asreml(
       fixed = get(i) ~ 1 + Set + LOC,
       random =  ~ Name2 + Set:block + LOC:Name2,
       #residual =  ~ id(rep):id(range):(row),
       data = designnew, na.action = na.method(x = "include"),
       predict = predict.asreml(classify = "Name2", sed = TRUE)
     ) }
  
  h2[i,2]  <- (1 - ((models[[i]]$predictions$avsed["mean"] ^ 2) /(2 * summary(models[[i]])$varcomp["Name2", "component"])))
}




create_df_predictions<- function(var, model) {
  df <- model$predictions$pvals[, 1:2]
  df$wavelength <- var
  return(df)
}
df_predictions <- purrr::map2(variables, models, create_df_predictions) %>% list_rbind()

df_predictions<- df_predictions %>% pivot_wider(names_from = wavelength, values_from = predicted.value)
write.csv(df_predictions,"prediction.csv")

df_prediction<- read.csv('prediction.csv')


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
whole.relation.matrix.csv<- read.csv('whole.relation.matrix.csv')

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
NIRS.R.matrix.csv<- read.csv('NIRS.R.matrix.csv')

#for selecting genotype with higher heritability
h2<-read.csv(paste0(here(), '/sorghumdata/h2.csv'))


library(tidyverse)

h2<- h2 %>% mutate(wavenumber= as.numeric(str_replace(X1,"Wave_", "")))

library(ggplot2)
ggplot(data= h2, aes(x = wavenumber, y = X2)) + geom_point(color= "purple", size= 0.1) + labs(title= "heritability per wavelength", x= "wavenumber",y= "heritability")



ggplot(data = h2[h2$X2 > 0.4, ], aes(x = wave, y = h2)) + 
  geom_point(color= "purple", size= 0.1) + 
  labs(title= "heritability per wavelength", x= "wavenumber",y= "heritability")

h2_maxwavelength<- h2[h2$X2 > 0.4712, c("X1", "X2")]




