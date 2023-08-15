library(simplePHENOTYPES)
library(data.table)
library(AGHmatrix)

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
dt_num <- fread("GBS002.pruned_numeric.txt", data.table = F)
dt_num[1:10, 1:10]


# the package AGHmatrix requires SNPs to be 0, 1, or 2
# transposing to have individuals in rows and SNPs in columns
# adding 1 to have it as 0, 1, 2
dt <- t(dt_num[, -1:-5]) + 1
rownames(dt) <- substr(rownames(dt), 1, nchar(rownames(dt)) / 2) # fix row names
dt[1:10,1:10]

names <- phenotypedata %>% select(TAXA, Name2)%>% filter(!duplicated(Name2))

df_prediction <- df_prediction[,2:3] %>% left_join(names)

# pd <- df_prediction[df_prediction$TAXA %in% rownames(dt),] 




# create an Additive relationship matrix
kin_A <- Gmatrix(dt)
library(data.table)
kin_A_dt <- as.data.table(kin_A)
fwrite(kin_A_dt, "kinship_additive.txt", sep = "\t", quote = FALSE)



# Load the data.table packages
library(data.table)

# Read the VCF file using fread
vcf_data <- fread("/scrfs/storage/ashmitau/home/Documents/repos/sorghumdata/GBS002.pruned.vcf")




markers <- data.frame(hapmap$GD - 1)
colnames(markers) <- hapmap$GI$SNP
rownames(markers) <- hapmap$GT
rm(hapmap)

# Get G matrix.
G <- G.matrix(M = as.matrix(your marker), method = "VanRaden")$G

#reading phenotypic data
phenotypedata<- read.csv('phenotypes.csv')
Gb <- G.tuneup(G = G, bend = TRUE, eig.tol = 1e-03)$Gb


GINV <- G.inverse(G = Gb, sparseform = T)

model <- asreml(fixed = cbind(SLA, N) ~ trait,
                random = ~ corv(trait):vm(genotype, source = GINV$Ginv.sparse),
                residual = ~ id(units):id(trait),
                #start.values=T,
                data = means)





