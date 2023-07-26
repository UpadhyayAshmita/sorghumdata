library(simplePHENOTYPES)
library(data.table)
library(AGHmatrix)
library(tidyverse)
library(ASRgenomics)
library(asreml)
library(here)


# converting from vcf to numeric...some files can be removed afterwards (only look at the *_numeric.txt file)
create_phenotypes(
  geno_file = "GBS002.pruned.vcf",
  add_QTN_num = 1,
  add_effect = 0.2,
  big_add_QTN_effect = 0.9,
  rep = 10,
  h2 = 0.7,
  model = "A",
  home_dir = getwd(),
  out_geno = "numeric",
)

# loading the numeric file
# SNPs are -1, 0, or 1
dt_num <- fread("./data/vcffiles/GBS002.pruned_numeric.txt", data.table = F)
dt_num[1:10, 1:10]


# the package AGHmatrix requires SNPs to be 0, 1, or 2
# transposing to have individuals in rows and SNPs in columns
# adding 1 to have it as 0, 1, 2
dt <- t(dt_num[, -1:-5]) + 1
rownames(dt) <- substr(rownames(dt), 1, nchar(rownames(dt)) / 2) # fix row names
dt[1:10,1:10]

# create an Additive relationship matrix
kin_A <- Gmatrix(dt)
kin_A_dt <- as.data.table(kin_A)
fwrite(kin_A_dt, "kinship_additive.txt", sep = "\t", quote = FALSE)

#reading data file and bending matrix and obtaining inverse of kinship matrix 

kin<- fread('./output/kinship_additive.txt', data.table= F)
rownames(kin) <- colnames(kin)
Gb <- G.tuneup(G = as.matrix(kin), bend = TRUE, eig.tol = 1e-06)$Gb
GINV <- G.inverse(G = Gb , sparseform = T)



#reading phenotypic data
phenotypedata<- read.csv('./data/phenotypes.csv')
column_name = 'row'
phenotypedata[[column_name]] <- 1:nrow(phenotypedata)
phenotypedata<- transform(phenotypedata,
                          TAXA= factor(TAXA),
                          LOC= factor(LOC),
                          Name2= factor(Name2),
                          Set= factor(Set),
                          block= factor(block),
                          range= factor(range),
                          row= factor(row))


#calculating BLUES for  EF location 
models <- list()
bluesEF <- list()
a <- 1
for(i in c("Narea","SLA")){
  models[[i]] <- asreml(
    fixed = get(i) ~ Name2,
    random =  ~  Set + Set:block,
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
    random =  ~  Set + Set:block,
    #residual =  ~ id(range):id(row),
    data = phenotypedata, subset = LOC == "MW", na.action = na.method(x = c("include")),
    predict = predict.asreml(classify = "Name2",sed = TRUE))
  bluesMW[[i]] <-   models[[i]]$predictions$pvals[,1:2]
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


#-------fitting GBLUP model with blues and kinship matrix obtained 


            
