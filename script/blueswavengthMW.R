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

# ---------------------calculating wavelength blues for locationMW---------------------

models <- list()
variables <- colnames(designnew)[11:2161]
wavebluesMW <- data.frame()  #n rows and 2 columns
wavebluesMW[,1] <- variables
colnames(wavebluesMW) <- c("wave", "blues")
rownames(wavebluesMW) <- variables
for(i in 1:length(variables)) {models[[i]] <- asreml(
  fixed = get(variables[i]) ~ Name2 + Set,
  random = ~ block,
  residual =  ~ ar1(range):ar1(row),
  data = designnew, subset = LOC== "MW", na.action= na.method(x = "include"),
  predict = predict.asreml(classify = "Name2", sed = TRUE))
temp1 <- models[[i]]$predictions$pvals[, 1:2]
temp1$wave <- variables[i]
wavebluesMW <- rbind(wavebluesMW, temp1)}

fwrite("./output/wavebluesMW.csv", row.names = FALSE)





