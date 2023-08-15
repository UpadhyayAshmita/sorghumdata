# ---------------------library---------------------
library(tidyverse)
library(asreml)
library(data.table)
library(purrr)
library(tidyverse)
library(ASRgenomics)
library(janitor)

# ---------------------load data---------------------
designnew<- fread('./data/designnew.csv',data.table= FALSE)


# ---------------------processing data---------------------

designnew<-
  designnew %>% mutate(
    Name2 = factor(Name2),
    TAXA = factor(TAXA),
    LOC = factor(LOC),
    Set = factor(Set),
    block = factor(block),
    range = factor(range),
    row = factor(row),
    uni = c(1:960, 1:960)
  ) %>%   arrange(LOC, range, row)



designnew$LOC <- factor(designnew$LOC, levels = c('EF', 'MW'))
designnew <- designnew[order(designnew$LOC), ]


#calculating blues for each wavelength 

models <- list()
variables <- colnames(designnew)[11:2161]
waveblues <- data.frame()  #n rows and 2 columns
waveblues[,1] <- variables
colnames(waveblues) <- c("wave", "blues")
rownames(waveblues) <- variables

for(i in 1:length(variables)) {models[[i]] <- asreml(
  fixed = get(variables[i]) ~ Name2 + LOC+ Name2:LOC + Set + LOC:Set,
  random = ~ LOC:block + Set,
  residual =  ~ id(LOC):ar1(range):ar1(row),
  data = designnew, na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "Name2", sed = TRUE))
temp <- models[[i]]$predictions$pvals[, 1:2]
temp$wave <- variables[i]
waveblues <- rbind(waveblues, temp)}


frwite(waveblues, "./output/waveblues.csv", row.names= FALSE)





fwrite(dt[1:10,1:10], "./output/test.csv")

# 
# # ---------------------h2 calculation---------------------
# 
# modelh2 <-
#   asreml(
#     fixed =  wave_1 ~ loc + set + loc:set,
#     random =  ~  name2 + loc:name2 + loc:block,
#     residual =  ~ id(loc):ar1(range):ar1(row),
#     data = phenotypes,
#     na.action = na.method(x = "include")
#   )
# #calculating blues for each wavelength 
# 
# 
# # ---------------------blues joint---------------------
# modelh2 <-
#   asreml(
#     fixed =  wave_1 ~ loc + name2 + loc:name2 + set + loc:set,
#     random =  ~   loc:block,
#     residual =  ~ id(loc):ar1(range):ar1(row),
#     data = phenotypes,
#     na.action = na.method(x = "include")
#   )
# #calculating blues for each wavelength 
# 
# 
# # ---------------------blues MW---------------------
# 
# modelmw <-
#   asreml(
#     fixed =  wave_1 ~ name2 + set,
#     random =  ~   block,
#     residual =  ~ ar1(range):ar1(row),
#     data = phenotypes, subset = loc == "MW",
#     na.action = na.method(x = "include")
#   )
# # ---------------------blues EF---------------------
# 
# modelmw <-
#   asreml(
#     fixed =  wave_1 ~ name2 + set,
#     random =  ~   block,
#     residual =  ~ ar1(range):ar1(row),
#     data = phenotypes, subset = loc == "EF",
#     na.action = na.method(x = "include")
#   )
# 
# models <- list()
# variables <- colnames(designnew)[11:2161]
# waveblues <- data.frame()  #n rows and 2 columns
# waveblues[,1] <- variables
# colnames(waveblues) <- c("wave", "blues")
# rownames(waveblues) <- variables
# 
# for(i in 1:length(variables)) {models[[i]] <- asreml(
#   fixed = get(variables[i]) ~ Name2,
#   random = ~ LOC:block + Set,
#   residual =  ~ id(loc):ar1(range):ar1(row),
#   data = designnew, na.action = na.method(x = "include"),
#   predict = predict.asreml(classify = "Name2", sed = TRUE))
# temp <- models[[i]]$predictions$pvals[, 1:2]
# temp$wave <- variables[i]
# waveblues <- rbind(waveblues, temp)}
# 
# 
# fwrite(waveblues, "./wavebluejoints.csv")
# 
# # ---------------------for loop task 2---------------------
# 
# 
# 
# # ---------------------for loop task 3---------------------
