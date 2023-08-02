# install ASREML
# install.packages("../asreml-4.1.0.176-CentOS-7-R4.2.tar.gz", repos = NULL)

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



models <- list()
variables <- colnames(designnew)[10:2160]
h2 <- data.frame(matrix(NA, length(variables), 2))  #n rows and 2 columns
h2[,1] <- variables
a <- 1
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
  
  h2[a,2]  <- (1 - ((models[[i]]$predictions$avsed["mean"] ^ 2) /(2 * summary(models[[i]])$varcomp["Name2", "component"])))
  h2[a,1]  <- i
  a <- a + 1
}


predicted.value<- models[[i]]$predictions$pvals

       

