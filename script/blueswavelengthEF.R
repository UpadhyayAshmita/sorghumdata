# ---------------------loading package--------------------

library(tidyverse)
library(asreml)
library(data.table)
library(purrr)
library(tidyverse)
library(ASRgenomics)


# ---------------------loading data---------------------
designnew<- fread('./data/designnew.csv',data.table= FALSE)

# ---------------------processing of data---------------------

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

# ---------------------calculating wavelength blues for locationEF---------------------

models <- list()
variables <- colnames(designnew)[11:2161]
wavebluesEF <- data.frame()  #n rows and 2 columns
wavebluesEF[,1] <- variables
colnames(wavebluesEF) <- c("wave", "blues")
rownames(wavebluesEF) <- variables
for(i in 1:length(variables)) {models[[i]] <- asreml(
  fixed = get(variables[i]) ~ Name2 + Set,
  random = ~ block,
  residual =  ~ ar1(range):ar1(row),
  data = designnew, subset = LOC== "EF", na.action= na.method(x = "include"),
  predict = predict.asreml(classify = "Name2", sed = TRUE))
temp2 <- models[[i]]$predictions$pvals[, 1:2]
temp2$wave <- variables[i]
wavebluesEF <- rbind(wavebluesMW, temp2)}

fwrite("./output/wavebluesMW.csv", row.names = FALSE)