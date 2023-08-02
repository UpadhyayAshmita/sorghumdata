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

<<<<<<< HEAD:script/genomic_kinship .R
=======
names <- phenotypedata %>% select(TAXA, Name2)%>% filter(!duplicated(Name2))

df_prediction <- df_prediction[,2:3] %>% left_join(names)

# pd <- df_prediction[df_prediction$TAXA %in% rownames(dt),] 




>>>>>>> de10e065c9c127d41c2c456c8b850cde80db3108:genomic_kinship .R
# create an Additive relationship matrix
kin_A <- Gmatrix(dt)
kin_A_dt <- as.data.table(kin_A)
fwrite(kin_A_dt, "kinship_additive.txt", sep = "\t", quote = FALSE)
<<<<<<< HEAD:script/genomic_kinship .R

#reading data file and bending matrix and obtaining inverse of kinship matrix 

kin<- fread('./output/kinship_additive.txt', data.table= F)
rownames(kin) <- colnames(kin)
Gb <- G.tuneup(G = as.matrix(kin), bend = TRUE, eig.tol = 1e-06)$Gb
GINV <- G.inverse(G = Gb , sparseform = T)





=======





models[[i]] <- asreml( fixed = get(i) ~ Name2, random = ~ Set + Set:block, 
                       residual = ~ id(range):id(row), data = wave, 
                       subset = env == "17EF", na.action = na.method(x = "include"), 
                       predict = predict.asreml(classify = "Name2") )


# Get G matrix.
#G <- G.matrix(M = as.matrix(your marker), method = "VanRaden")$G

#reading phenotypic data
phenotypedata<- read.csv('phenotypes.csv')
>>>>>>> de10e065c9c127d41c2c456c8b850cde80db3108:genomic_kinship .R




<<<<<<< HEAD:script/genomic_kinship .R
            
=======

Gb <- G.tuneup(G = G, bend = TRUE, eig.tol = 1e-03)$Gb


GINV <- G.inverse(G = Gb, sparseform = T)


#predicted value for nitrogen leaf area with vector having two column and 800 individual rows
nitrogendata<- filter(phenotypedata, LOC== "EF" &  year == "17")





phenotypedata<- transform(phenotypedata,
                      TAXA= factor(TAXA),
                      LOC= factor(LOC),
                      Set= factor(Set),
                      block= factor(block),
                      range= factor(range),
                      row= factor(row))





Nblups <- asreml( 
  fixed = Narea ~ TAXA, random = ~ Set + Set:block, 
  # residual = ~ id(range):id(row), 
  data = phenotypedata, 
  subset = (year == 17) & (LOC == "EF"),
  na.action = na.method(x = "include",y= "include"), 
  predict = predict.asreml(classify = "TAXA")
)




library(tidyverse)
>>>>>>> de10e065c9c127d41c2c456c8b850cde80db3108:genomic_kinship .R
