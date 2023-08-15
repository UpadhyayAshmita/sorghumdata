# ---------------------loading library---------------------
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

# ---------------------calculating heritability---------------------

models <- list()
variables <- colnames(designnew)[11:2161]
h2 <- data.frame(matrix(NA, length(variables), 2))  #n rows and 2 columns
h2[,1] <- variables
colnames(h2) <- c("wave", "h2")
rownames(h2) <- variables

for(i in variables) { models[[i]] <- asreml(
  fixed = get(i) ~ 1 + Set + LOC + LOC: Set,
  random = ~Name2 + LOC:block + LOC:Name2,
  residual =  ~id(LOC):ar1(range):ar1(row),
  data = designnew, na.action = na.method(x = "include"),
  predict = predict.asreml(classify = "Name2", sed = TRUE)
)

out <- outlier(models[[i]]$residuals)

for(i in variables) {
  designnew[,i][designnew$LOC == "EF"][out] <- NA
  models[[i]] <- asreml(
    fixed = get(i) ~ 1 + Set + LOC + LOC:Set,
    random =  ~ Name2 + LOC:block + LOC:Name2,
    residual =  ~ id(LOC):ar1(range):ar1(row),
    data = designnew, na.action = na.method(x = "include"),
    predict = predict.asreml(classify = "Name2", sed = TRUE)
    
  )
}

h2[i,2]  <- (1 - ((models[[i]]$predictions$avsed["mean"] ^ 2) /(2 * summary(models[[i]])$varcomp["Name2", "component"])))}


fwrite(h2, "h2.csv",row.names = F)







