
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


