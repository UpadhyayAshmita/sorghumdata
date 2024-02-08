# ---------------------loading packages---------------------
library(data.table)
library(tidyverse)
library(asreml)
library(ASRgenomics)
library(simplePHENOTYPES)
library(AGHmatrix)
source("./script/aux_functions.R")

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
colnames(dt_num) <- gsub("_.*","", colnames(dt_num))
nitroblues<- fread("./output/traitsoutput/bluesjoint_narea_correct.csv", data.table = F)
selected_lines <- colnames(dt_num)[-1:-5][colnames(dt_num)[-1:-5] %in% nitroblues$taxa]
dt_num <- dt_num[,c(colnames(dt_num)[1:5], selected_lines)]

# the package AGHmatrix requires SNPs to be 0, 1, or 2
# transposing to have individuals in rows and SNPs in columns
# adding 1 to have it as 0, 1, 2
dt <- t(dt_num[, -1:-5]) + 1
#rownames(dt) <- substr(rownames(dt), 1, nchar(rownames(dt)) / 2) # fix row names
dt[1:10,1:10]

# create an Additive relationship matrix
kin_A <- Gmatrix(dt)
kin_A_dt <- as.data.table(kin_A)
fwrite(kin_A_dt, "./data/relmatrices/kinship_additive.txt", sep = "\t", quote = FALSE)

#reading data file and bending matrix and obtaining inverse of kinship matrix 
kin<- fread('./data/relmatrices/kinship_additive.txt', data.table= F)
rownames(kin) <- colnames(kin)
kin<- as.matrix(kin)
# Gb <- G.tuneup(G = as.matrix(kin), bend = TRUE, eig.tol = 1e-06)$Gb
# GINV <- G.inverse(G = Gb , sparseform = T)
# saveRDS(GINV, file = "./data/relmatrices/GINV.rds")
# GINV <- readRDS(file = "./data/relmatrices/GINV.rds")
wavebluesEF_wider<- fread("./output/wavebluesEF_wider.csv")
# left joining corrected name with name2 from wavelength blues to filter individual that are also present in genotypic data ---------------------
Names_WEST<- fread("./data/Names_WEST_SF.csv", data.table = F)
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[2]<- "name2"
wavebluesEF_wider <- 
  wavebluesEF_wider |>
  left_join(Names_WEST) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) & grepl("PI", Corrected_names), NA, Corrected_names)) %>% 
  filter(!is.na(taxa)) %>% select(-name2, -Corrected_names ) %>% select(taxa, everything()) 
wavebluesEF_wider<- wavebluesEF_wider %>% filter(wavebluesEF_wider$taxa %in% selected_lines)
fwrite(wavebluesEF_wider, "./output/wavebluesEF_wider_corrected.csv")
#reading wavebluesef data
wavebluesEF_wider_corrected<- fread("./output/wavebluesEF_wider_corrected.csv")
# ---------------------matrix for whole spectra---------------------
re_w_ef<- calculate_relationship(wavebluesEF_wider_corrected)
# ---------------------saving relationship matrix for whole location---------------------
write.csv(re_w_ef,'./data/relmatrices/re_w_ef.csv',row.names= F)
re_w_ef<- read.csv('./data/relmatrices/re_w_ef.csv')
rownames(re_w_ef) <- colnames(re_w_ef)
Gb <- G.tuneup(G = as.matrix(re_w_ef), bend = TRUE, eig.tol = 1e-06)$Gb
re_w_ef_inv<- G.inverse(G = Gb , sparseform = T)
saveRDS(re_w_ef_inv, file = "./data/relmatrices/re_w_ef_inv.rds")
# ---------------------creating relationship matrix for nirs data---------------------
wavenirs<-wavebluesEF_wider_corrected %>%select(taxa,"wave_860":"wave_1660")
re_nirs_ef<- calculate_relationship(wavenirs)
write.csv(re_nirs_ef,'./data/relmatrices/re_nirs_ef.csv', row.names= F)
re_nirs_ef<- read.csv('./data/relmatrices/re_nirs_ef.csv')
rownames(re_nirs_ef) <- colnames(re_nirs_ef)
Gb <- G.tuneup(G = as.matrix(re_nirs_ef), bend = TRUE, eig.tol = 1e-06)$Gb
nirs_inv_ef<- G.inverse(G = Gb , sparseform = T)
saveRDS(nirs_inv_ef, file = "./data/relmatrices/nirs_inv_ef.rds")
nirs_inv_ef <- readRDS(file = "./data/relmatrices/nirs_inv_ef.rds")

# ---------------------for selecting genotype with higher heritability---------------------
# ---------------------reading h2 for ef location---------------------
h2ef<-fread('./output/h2EF.csv',data.table = F) #reading h2 data
#h2ef<- h2ef[-2152,]
#visualizing trend in plot with wavelegth vs heritability 
h2ef<- h2ef %>% mutate(wavenumber= as.numeric(str_replace(wave,"wave_", "")))
ggplot(data= h2ef, aes(x = wavenumber, y = h2)) + geom_point(color= "purple", size= 0.1) + labs(title= "heritability per wavelength", x= "wavenumber",y= "heritability")
ggplot(data = h2ef[h2ef$h2 > 0.4, ], aes(x = wavenumber, y = h2)) + geom_point(color= "purple", size= 0.1) + labs(title= "heritability per wavelength", x= "wavenumber",y= "heritability")
# ---------------------filtering wavelength that are highly heritable according to the trend and the threshold---------------------
highh2wave<- h2ef[h2ef$h2 > 0.4, c( "h2","wave")]
# ---------------------filtering bluesef data for wavelength that are highly heritable---------------------
wavebluesEF_wider_corrected<- fread("./output/wavebluesEF_wider_corrected.csv")
wavebluesEF_wider_corrected<- wavebluesEF_wider_corrected %>% pivot_longer(cols= 12:2152,names_to = "wave", values_to = "predicted.value")
highwaveblues<- wavebluesEF_wider_corrected %>%  filter(wave %in% highh2wave$wave)
highwaveblues<- highwaveblues %>% pivot_wider(names_from = "wave",values_from = "predicted.value")
# ---------------------saving the high heritable waveblues ---------------------
fwrite(highwaveblues, "./output/intermediate/highwaveblues.csv", row.names = F)
# ---------------------high heritable relationship matrix---------------------
highwaveblues<- fread("./output/intermediate/highwaveblues.csv")
highwaveblues_filtered<- highwaveblues%>%
  filter(taxa %in% wavebluesEF_wider_corrected$taxa)
re_h2_ef<- calculate_relationship(highwaveblues_filtered)
write.csv(re_h2_ef, "./data/relmatrices/re_h2_ef.csv", row.names = F)
re_h2_ef<-read.csv("./data/relmatrices/re_h2_ef.csv")
rownames(re_h2_ef) <- colnames(re_h2_ef)
Gb <- G.tuneup(G = as.matrix(re_h2_ef), bend = TRUE, eig.tol = 1e-06)$Gb
re_h2_ef_inv<- G.inverse(G = Gb , sparseform = T)
saveRDS(re_h2_ef_inv, file = "./data/relmatrices/re_h2_ef_inv.rds")

# ---------------------calculating relationship matrix for "mw" location---------------------
# ---------------------reading mw blues data---------------------
wavebluesMW_wider<- fread("./output/wavebluesMW_wider.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[2]<- "name2"
wavebluesMW_wider_corrected<- 
  wavebluesMW_wider |>
  left_join(Names_WEST) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) & grepl("PI", Corrected_names), NA, Corrected_names)) %>% 
  filter(!is.na(taxa)) %>% select(-name2, -Corrected_names ) %>% select(taxa, everything())
wavebluesMW_wider_corrected<- wavebluesMW_wider_corrected %>%  filter(wavebluesMW_wider_corrected$taxa %in% selected_lines)
fwrite(wavebluesMW_wider_corrected, "./output/wavebluesMW_wider_corrected.csv")
wavebluesMW_wider_corrected<- fread("./output/wavebluesMW_wider_corrected.csv")
#generating relationship matrix for whole spectra
re_w_mw <- calculate_relationship(wavebluesMW_wider_corrected)
write.csv(re_w_mw,'./data/relmatrices/re_w_MW.csv',row.names= F)
re_w_mw<- read.csv('./data/relmatrices/re_w_MW.csv')
rownames(re_w_mw) <- colnames(re_w_mw)
Gb <- G.tuneup(G = as.matrix(re_w_mw), bend = TRUE, eig.tol = 1e-06)$Gb
re_w_mw_inv<- G.inverse(G = Gb , sparseform = T)
saveRDS(re_w_mw_inv, file = "./data/relmatrices/re_w_mw_inv.rds")

# ---------------------relationship matrix for NIRS wavelength---------------------
wavenirsmw<- wavebluesMW_wider_corrected%>% select(taxa,"wave_860":"wave_1660")
fwrite(wavenirsmw,"./output/intermediate/wavenirsmw.csv")
re_nirs_mw<- calculate_relationship(wavenirsmw)
write.csv(re_nirs_mw,'./data/relmatrices/re_nirs_mw.csv', row.names = F)
re_nirs_mw<- read.csv('./data/relmatrices/re_nirs_mw.csv')
rownames(re_nirs_mw) <- colnames(re_nirs_mw)
Gb <- G.tuneup(G = as.matrix(re_nirs_mw), bend = TRUE, eig.tol = 1e-06)$Gb
nirs_inv_mw<- G.inverse(G = Gb , sparseform = T)
saveRDS(nirs_inv_mw, file = "./data/relmatrices/nirs_inv_mw.rds")
nirs_inv_mw <- readRDS(file = "./data/relmatrices/nirs_inv_mw.rds")
#for selecting genotype with higher heritability
h2MW_binded<-fread('./output/h2MW_binded.csv',data.table = F)#reading h2 data
#visualizing trend in plot with wavelegth vs heritability 
h2MW_binded<- h2MW_binded %>%
  mutate(wavenumber= as.numeric(str_replace(wave,"wave_", "")))
ggplot(data= h2MW_binded, aes(x = wavenumber, y = h2)) + geom_point(color= "purple", size= 0.1) + labs(title= "heritability per wavelength", x= "wavenumber",y= "heritability")
ggplot(data = h2MW_binded[h2MW_binded$h2 > 0.4, ], aes(x = wavenumber, y = h2)) + geom_point(color= "purple", size= 0.1) + labs(title= "heritability per wavelength", x= "wavenumber",y= "heritability")
highh2wave<- h2MW_binded[h2MW_binded$h2 > 0.4, c( "h2","wave")]
wavebluesMW_wider_corrected<- fread("./output/wavebluesMW_wider_corrected.csv")
wavebluesMW_long_corrected<- wavebluesMW_wider_corrected %>% 
  pivot_longer(cols= 2:2152,names_to = "wave", values_to = "predicted.value")
highwaveblues<- wavebluesMW_long_corrected%>% 
  filter(wave %in% highh2wave$wave)
highwaveblues<- highwaveblues %>% 
  pivot_wider(names_from = "wave",values_from = "predicted.value")
# ---------------------high heritable relationship matrix---------------------
re_h2_mw<- calculate_relationship(highwaveblues)
write.csv(re_h2_mw, "./data/relmatrices/re_h2_mw.csv", row.names = F)
re_h2_mw<-read.csv("./data/relmatrices/re_h2_mw.csv")
rownames(re_h2_mw) <- colnames(re_h2_mw)
Gb <- G.tuneup(G = as.matrix(re_h2_mw), bend = TRUE, eig.tol = 1e-06)$Gb
re_h2_mw_inv<- G.inverse(G = Gb , sparseform = T)
saveRDS(re_h2_mw_inv, file = "./data/relmatrices/re_h2_mw_inv.rds")
# ---------------------creating three different relationship matrix for joint location---------------------
wavebluesjoint_wider<- wavejoint_binded %>% pivot_wider(names_from = "wave",values_from = "predicted.value")
fwrite(wavebluesjoint_wider, "./output/wavebluesjoint_wider.csv", row.names = F)
# ---------------------reading the processed wavejoint data---------------------
wavebluesjoint_wider<- fread("./output/wavebluesjoint_wider.csv")
Names_WEST<- fread("./data/Names_WEST_SF.csv", data.table = F)
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[2]<- "name2"
wavebluesjoint_corrected<- 
  wavebluesjoint_wider |>
  left_join(Names_WEST) %>%
  mutate(taxa = ifelse((name2 != Corrected_names) & grepl("PI", Corrected_names), NA, Corrected_names)) %>% 
  filter(!is.na(taxa)) %>% select(-name2, -Corrected_names ) %>% select(taxa, everything())
wavebluesjoint_corrected<- wavebluesjoint_corrected %>%  filter(wavebluesjoint_corrected$taxa %in% selected_lines)

fwrite(wavebluesjoint_corrected, "./output/wavebluesjoint_corrected.csv")
wavebluesjoint_corrected<-fread("./output/wavebluesjoint_corrected.csv")
wavebluesjoint_long<- pivot_longer(wavebluesjoint_corrected,cols= 2:2152,values_to = "predicted.values", names_to = "wave")
# ---------------------calculating joint whole spectra relationship matrix---------------------
re_w_joint<- calculate_relationship(wavebluesjoint_corrected)
write.csv(re_w_joint, "./data/relmatrices/re_w_joint.csv", row.names = F)
re_w_joint<- read.csv("./data/relmatrices/re_w_joint.csv")
rownames(re_w_joint) <- colnames(re_w_joint)
Gb <- G.tuneup(G = as.matrix(re_w_joint), bend = TRUE, eig.tol = 1e-06)$Gb
re_w_joint_inv<- G.inverse(G = Gb , sparseform = T)
saveRDS(re_w_joint_inv, file = "./data/relmatrices/re_w_joint_inv.rds")
# ---------------------calculating NIRS matrix for joint model---------------------
wavenirsjoint<- wavebluesjoint_corrected %>% 
  select(taxa,"wave_860":"wave_1660")
fwrite(wavenirsjoint,"./output/intermediate/wavenirsjoint.csv")
re_nirs_joint<- calculate_relationship(wavenirsjoint)
write.csv(re_nirs_joint,'./data/relmatrices/re_nirs_joint.csv', row.names = F)
re_nirs_joint<- read.csv('./data/relmatrices/re_nirs_joint.csv')
rownames(re_nirs_joint) <- colnames(re_nirs_joint)
Gb <- G.tuneup(G = as.matrix(re_nirs_joint), bend = TRUE, eig.tol = 1e-06)$Gb
nirs_inv<- G.inverse(G = Gb , sparseform = T)
saveRDS(nirs_inv, file = "./data/relmatrices/nirs_inv.rds")
# ---------------------calculating highly h2 matrix---------------------
h2joint<- fread("./output/h2.csv", data.table = FALSE)
# ---------------------visualizing trend in plot with wavelegth vs heritability ---------------------
h2joint<- h2joint %>% mutate(wavenumber= as.numeric(str_replace(wave,"wave_", "")))
ggplot(data= h2joint, aes(x = wavenumber, y = h2)) + geom_point(color= "purple", size= 0.1) + labs(title= "heritability per wavelength", x= "wavenumber",y= "heritability")
ggplot(data = h2joint[h2joint$h2 > 0.4, ], aes(x = wavenumber, y = h2)) + geom_point(color= "purple", size= 0.1) + labs(title= "heritability per wavelength", x= "wavenumber",y= "heritability")
highh2wave_joint<- h2joint[h2joint$h2 > 0.4, c( "h2","wave")]
highwaveblues_joint<- wavebluesjoint_long %>% filter(wave %in% highh2wave_joint$wave)
highwaveblues<- highwaveblues_joint %>% pivot_wider(names_from = "wave",values_from = "predicted.values")
# ---------------------calculating relationship matrix for highly h2 wavelength for joint location---------------------
re_h2_joint<- calculate_relationship(highwaveblues)
write.csv(re_h2_joint, "./data/relmatrices/re_h2_joint.csv", row.names = F)
re_h2_joint<- read.csv("./data/relmatrices/re_h2_joint.csv")
rownames(re_h2_joint) <- colnames(re_h2_joint)
Gb <- G.tuneup(G = as.matrix(re_h2_joint), bend = TRUE, eig.tol = 1e-06)$Gb
re_h2_joint_inv<- G.inverse(G = Gb , sparseform = T)
saveRDS(re_h2_joint_inv, file = "./data/relmatrices/re_h2_joint_inv.rds")

#removing 10%, 25% and 50% individual from kinship matrix
kin<- fread('./data/relmatrices/kinship_additive.txt', data.table= F)
# ---------------------10%---------------------
set.seed(1)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.1))
# --------------------for GBLUP--------------------
G10<- kin[-remove, ]
G10<- G10[, -remove]
G10 <- na.omit(G10)
fwrite(G10, "./data/relmatrices/GBLUP/G10.csv", row.names = F)
rownames(G10) <- colnames(G10)
Gb <- G.tuneup(G = as.matrix(G10), bend = TRUE, eig.tol = 1e-06)$Gb
G10inv<- G.inverse(G = Gb , sparseform = T)
saveRDS(G10inv, file = "./data/relmatrices/GBLUP/G10inv.rds")
# ---------------------25%---------------------
set.seed(2)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.25))
# --------------------for GBLUP--------------------
G25<- kin[-remove, ]
G25<- G25[, -remove]
G25 <- na.omit(G25)
fwrite(G25, "./data/relmatrices/GBLUP/G25.csv", row.names = F)
rownames(G25) <- colnames(G25)
Gb <- G.tuneup(G = as.matrix(G25), bend = TRUE, eig.tol = 1e-06)$Gb
G25inv<- G.inverse(G = Gb , sparseform = T)
saveRDS(G25inv, file = "./data/relmatrices/GBLUP/G25inv.rds")
# ---------------------50%---------------------
set.seed(5)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.50))
# --------------------for GBLUP--------------------
G50<- kin[-remove, ]
G50<- G50[, -remove]
G50 <- na.omit(G50)
fwrite(G50, "./data/relmatrices/GBLUP/G50.csv", row.names = F)
rownames(G50) <- colnames(G50)
Gb <- G.tuneup(G = as.matrix(G50), bend = TRUE, eig.tol = 1e-06)$Gb
G50inv<- G.inverse(G = Gb , sparseform = T)
saveRDS(G50inv, file = "./data/relmatrices/GBLUP/G50inv.rds")
# ----replacing random 10%,25%,50% from kinship and replacing 10% by high heritable matrix for joint location---------------------
kin<- fread('./data/relmatrices/kinship_additive.txt', data.table= F)
re_h2_joint<- read.csv("./data/relmatrices/re_h2_joint.csv")
# ---------------------10%---------------------
set.seed(1)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.1))
Gh2<- kin
# --------------------imputing kin removed info with highh2 matrix---------------------
Gh2[remove, ] <- re_h2_joint[remove, ]
Gh2[,remove ] <- re_h2_joint[,remove]
fwrite(Gh2, "./data/relmatrices/Gh2/joint/10/Gh2.csv",row.names=F)
Gh2<- read.csv("./data/relmatrices/Gh2/joint/10/Gh2.csv")
rownames(Gh2) <- colnames(Gh2)
Gb <- G.tuneup(G = as.matrix(Gh2), bend = TRUE, eig.tol = 1e-06)$Gb
Gh2inv_joint<- G.inverse(G = Gb , sparseform = T)
saveRDS(Gh2inv_joint, file = "./data/relmatrices/Gh2/joint/10/Gh2inv_joint.rds")

# ---------------------25%---------------------
set.seed(2)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.25))
Gh2 <- kin
# --------------------imputing kin removed info with high h2 matrix---------------------
Gh2[remove, ] <- re_h2_joint[remove, ]
Gh2[,remove ] <- re_h2_joint[,remove]
fwrite(Gh2, "./data/relmatrices/Gh2/joint/25/Gh2.csv",row.names=F)
Gh2<- read.csv("./data/relmatrices/Gh2/joint/25/Gh2.csv")
rownames(Gh2) <- colnames(Gh2)
Gb <- G.tuneup(G = as.matrix(Gh2), bend = TRUE, eig.tol = 1e-06)$Gb
Gh2inv_joint<- G.inverse(G = Gb , sparseform = T)
saveRDS(Gh2inv_joint, file = "./data/relmatrices/Gh2/joint/25/Gh2inv_joint.rds")
# ---------------------50%---------------------
set.seed(5)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.50))
Gh2 <- kin
# --------------------imputing kin removed info with high h2 matrix---------------------
Gh2[remove, ] <- re_h2_joint[remove, ]
Gh2[,remove ] <- re_h2_joint[,remove]
fwrite(Gh2, "./data/relmatrices/Gh2/joint/50/Gh2.csv",row.names=F)
Gh2<- read.csv("./data/relmatrices/Gh2/joint/50/Gh2.csv")
rownames(Gh2) <- colnames(Gh2)
Gb <- G.tuneup(G = as.matrix(Gh2), bend = TRUE, eig.tol = 1e-06)$Gb
Gh2inv_joint<- G.inverse(G = Gb , sparseform = T)
saveRDS(Gh2inv_joint, file = "./data/relmatrices/Gh2/joint/50/Gh2inv_joint.rds")
# ----replacing random 10%,25%,50% from kinship and replacing 10% by high heritable matrix for ef location---------------------
kin<- fread('./data/relmatrices/kinship_additive.txt', data.table= F)
re_h2_ef<- read.csv("./data/relmatrices/re_h2_ef.csv")
# ---------------------10%---------------------
set.seed(7)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.1))
# --------------------imputing kin removed info with high h2 matrix---------------------
Gh2<- kin
Gh2[remove, ] <- re_h2_ef[remove, ]
Gh2[,remove ] <- re_h2_ef[,remove]
fwrite(Gh2, "./data/relmatrices/Gh2/ef/10/Gh2.csv",row.names=F)
Gh2<- read.csv("./data/relmatrices/Gh2/ef/10/Gh2.csv")
rownames(Gh2) <- colnames(Gh2)
Gb <- G.tuneup(G = as.matrix(Gh2), bend = TRUE, eig.tol = 1e-06)$Gb
Gh2inv_ef<- G.inverse(G = Gb , sparseform = T)
saveRDS(Gh2inv_ef, file = "./data/relmatrices/Gh2/ef/10/Gh2inv_ef.rds")
# ---------------------25---------------------
set.seed(6)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.25))
Gh2 <- kin
# --------------------imputing kin removed info with high h2 matrix---------------------
Gh2[remove, ] <- re_h2_ef[remove, ]
Gh2[,remove ] <- re_h2_ef[,remove]
fwrite(Gh2, "./data/relmatrices/Gh2/ef/25/Gh2.csv",row.names=F)
Gh2<- read.csv("./data/relmatrices/Gh2/ef/25/Gh2.csv")
rownames(Gh2) <- colnames(Gh2)
Gb <- G.tuneup(G = as.matrix(Gh2), bend = TRUE, eig.tol = 1e-06)$Gb
Gh2inv_ef<- G.inverse(G = Gb , sparseform = T)
saveRDS(Gh2inv_ef, file = "./data/relmatrices/Gh2/ef/25/Gh2inv_ef.rds")

# ---------------------50%---------------------
set.seed(9)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.50))
Gh2 <- kin
# --------------------imputing kin removed info with high h2 matrix---------------------
Gh2[remove, ] <- re_h2_ef[remove, ]
Gh2[,remove ] <- re_h2_ef[,remove]
fwrite(Gh2, "./data/relmatrices/Gh2/ef/50/Gh2.csv",row.names=F)
Gh2<- read.csv("./data/relmatrices/Gh2/ef/50/Gh2.csv")
rownames(Gh2) <- colnames(Gh2)
Gb <- G.tuneup(G = as.matrix(Gh2), bend = TRUE, eig.tol = 1e-06)$Gb
Gh2inv_ef<- G.inverse(G = Gb , sparseform = T)
saveRDS(Gh2inv_ef, file = "./data/relmatrices/Gh2/ef/50/Gh2inv_ef.rds")
# ----replacing random 10%,25%,50% from kinship and replacing 10% by high heritable matrix for mw location---------------------
kin<- fread('./data/relmatrices/kinship_additive.txt', data.table= F)
re_h2_mw<- read.csv("./data/relmatrices/re_h2_mw.csv")
# ---------------------10%---------------------
set.seed(7)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.1))
# --------------------imputing kin removed info with high h2 matrix---------------------
Gh2<- kin
Gh2[remove, ] <- re_h2_mw[remove, ]
Gh2[,remove ] <- re_h2_mw[,remove]
fwrite(Gh2, "./data/relmatrices/Gh2/mw/10/Gh2.csv",row.names=F)
Gh2<- read.csv("./data/relmatrices/Gh2/mw/10/Gh2.csv")
rownames(Gh2) <- colnames(Gh2)
Gb <- G.tuneup(G = as.matrix(Gh2), bend = TRUE, eig.tol = 1e-06)$Gb
Gh2inv_mw<- G.inverse(G = Gb , sparseform = T)
saveRDS(Gh2inv_mw, file = "./data/relmatrices/Gh2/mw/10/Gh2inv_mw.rds")
# ---------------------25%---------------------
set.seed(6)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.25))
Gh2 <- kin
# --------------------imputing kin removed info with high h2 matrix---------------------
Gh2[remove, ] <- re_h2_mw[remove, ]
Gh2[,remove ] <- re_h2_mw[,remove]
fwrite(Gh2, "./data/relmatrices/Gh2/mw/25/Gh2.csv",row.names=F)
Gh2<- read.csv("./data/relmatrices/Gh2/mw/25/Gh2.csv")
rownames(Gh2) <- colnames(Gh2)
Gb <- G.tuneup(G = as.matrix(Gh2), bend = TRUE, eig.tol = 1e-06)$Gb
Gh2inv_mw<- G.inverse(G = Gb , sparseform = T)
saveRDS(Gh2inv_mw, file = "./data/relmatrices/Gh2/mw/25/Gh2inv_mw.rds")
# ---------------------50%---------------------
set.seed(9)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.50))
Gh2 <- kin
# --------------------imputing kin removed info with high h2 matrix---------------------
Gh2[remove, ] <- re_h2_mw[remove, ]
Gh2[,remove ] <- re_h2_mw[,remove]
fwrite(Gh2, "./data/relmatrices/Gh2/mw/50/Gh2.csv",row.names=F)
Gh2<- read.csv("./data/relmatrices/Gh2/mw/50/Gh2.csv")
rownames(Gh2) <- colnames(Gh2)
Gb <- G.tuneup(G = as.matrix(Gh2), bend = TRUE, eig.tol = 1e-06)$Gb
Gh2inv_mw<- G.inverse(G = Gb , sparseform = T)
saveRDS(Gh2inv_mw, file = "./data/relmatrices/Gh2/mw/50/Gh2inv_mw.rds")

# ----replacing random 10%,25%,50% from kinship and replacing 10% by nirs matrix for joint location---------------------
kin<- fread('./data/relmatrices/kinship_additive.txt', data.table= F)
re_nirs_joint<- read.csv("./data/relmatrices/re_nirs_joint.csv")
# ---------------------10%---------------------
set.seed(7)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.1))
Gnirs<- kin
# --------------------imputing kin removed info with nirs matrix---------------------
Gnirs[remove, ] <- re_nirs_joint[remove, ]
Gnirs[,remove ] <- re_nirs_joint[,remove]
fwrite(Gnirs, "./data/relmatrices/Gnirs/joint/10/Gnirs.csv",row.names=F)
Gnirs<- read.csv("./data/relmatrices/Gnirs/joint/10/Gnirs.csv")
rownames(Gnirs) <- colnames(Gnirs)
Gb <- G.tuneup(G = as.matrix(Gnirs), bend = TRUE, eig.tol = 1e-06)$Gb
Gnirsinv_joint<- G.inverse(G = Gb , sparseform = T)
saveRDS(Gnirsinv_joint, file = "./data/relmatrices/Gnirs/joint/10/Gnirsinv_joint.rds")

# ---------------------25%---------------------
set.seed(6)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.25))
Gnirs <- kin
# --------------------imputing kin removed info with nirs matrix---------------------
Gnirs[remove, ] <- re_nirs_joint[remove, ]
Gnirs[,remove ] <- re_nirs_joint[,remove]
fwrite(Gnirs, "./data/relmatrices/Gnirs/joint/25/Gnirs.csv",row.names=F)
Gnirs<- read.csv("./data/relmatrices/Gnirs/joint/25/Gnirs.csv")
rownames(Gnirs) <- colnames(Gnirs)
Gb <- G.tuneup(G = as.matrix(Gnirs), bend = TRUE, eig.tol = 1e-06)$Gb
Gnirsinv_joint<- G.inverse(G = Gb , sparseform = T)
saveRDS(Gnirsinv_joint, file = "./data/relmatrices/Gnirs/joint/25/Gnirsinv_joint.rds")
# ---------------------50%---------------------
set.seed(9)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.50))
Gnirs <- kin
# --------------------imputing kin removed info with nirs matrix---------------------
Gnirs[remove, ] <- re_nirs_joint[remove, ]
Gnirs[,remove ] <- re_nirs_joint[,remove]
fwrite(Gnirs, "./data/relmatrices/Gnirs/joint/50/Gnirs.csv",row.names=F)
Gnirs<- read.csv("./data/relmatrices/Gnirs/joint/50/Gnirs.csv")
rownames(Gnirs) <- colnames(Gnirs)
Gb <- G.tuneup(G = as.matrix(Gnirs), bend = TRUE, eig.tol = 1e-06)$Gb
Gnirsinv_joint<- G.inverse(G = Gb , sparseform = T)
saveRDS(Gnirsinv_joint, file = "./data/relmatrices/Gnirs/joint/50/Gnirsinv_joint.rds")

# ----replacing random 10%,25%,50% from kinship and replacing 10% by nirs matrix for ef loc---------------------
kin<- fread('./data/relmatrices/kinship_additive.txt', data.table= F)
re_nirs_ef<- read.csv("./data/relmatrices/re_nirs_ef.csv")
# ---------------------10%---------------------
set.seed(7)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.1))
Gnirs<- kin
# --------------------imputing kin removed info with nirs matrix--------------------
Gnirs[remove, ] <- re_nirs_ef[remove, ]
Gnirs[,remove ] <- re_nirs_ef[,remove]
fwrite(Gnirs, "./data/relmatrices/Gnirs/ef/10/Gnirs.csv",row.names=F)
Gnirs<- read.csv("./data/relmatrices/Gnirs/ef/10/Gnirs.csv")
rownames(Gnirs) <- colnames(Gnirs)
Gb <- G.tuneup(G = as.matrix(Gnirs), bend = TRUE, eig.tol = 1e-06)$Gb
Gnirsinv_ef<- G.inverse(G = Gb , sparseform = T)
saveRDS(Gnirsinv_ef, file = "./data/relmatrices/Gnirs/ef/10/Gnirsinv_ef.rds")

# ---------------------25%---------------------
set.seed(6)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.25))
Gnirs<- kin
# --------------------imputing kin removed info with nirs matrix--------------------
Gnirs[remove, ] <- re_nirs_ef[remove, ]
Gnirs[,remove ] <- re_nirs_ef[,remove]
fwrite(Gnirs, "./data/relmatrices/Gnirs/ef/25/Gnirs.csv",row.names=F)
Gnirs<- read.csv("./data/relmatrices/Gnirs/ef/25/Gnirs.csv")
rownames(Gnirs) <- colnames(Gnirs)
Gb <- G.tuneup(G = as.matrix(Gnirs), bend = TRUE, eig.tol = 1e-06)$Gb
Gnirsinv_ef<- G.inverse(G = Gb , sparseform = T)
saveRDS(Gnirsinv_ef, file = "./data/relmatrices/Gnirs/ef/25/Gnirsinv_ef.rds")
# ---------------------50%---------------------
set.seed(9)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.50))
Gnirs<- kin
# --------------------imputing kin removed info with nirs matrix--------------------
Gnirs[remove, ] <- re_nirs_ef[remove, ]
Gnirs[,remove ] <- re_nirs_ef[,remove]
fwrite(Gnirs, "./data/relmatrices/Gnirs/ef/50/Gnirs.csv",row.names=F)
Gnirs<- read.csv("./data/relmatrices/Gnirs/ef/50/Gnirs.csv")
rownames(Gnirs) <- colnames(Gnirs)
Gb <- G.tuneup(G = as.matrix(Gnirs), bend = TRUE, eig.tol = 1e-06)$Gb
Gnirsinv_ef<- G.inverse(G = Gb , sparseform = T)
saveRDS(Gnirsinv_ef, file = "./data/relmatrices/Gnirs/ef/50/Gnirsinv_ef.rds")
# ----replacing random 10%,25%,50% from kinship and replacing 10% by nirs matrix for mw loc---------------------
kin<- fread('./data/relmatrices/kinship_additive.txt', data.table= F)
re_nirs_mw<- read.csv("./data/relmatrices/re_nirs_mw.csv")
# ---------------------10%---------------------
set.seed(7)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.1))
Gnirs<- kin
# --------------------imputing kin removed info with nirs matrix--------------------
Gnirs[remove, ] <- re_nirs_mw[remove, ]
Gnirs[,remove ] <- re_nirs_mw[,remove]
fwrite(Gnirs, "./data/relmatrices/Gnirs/mw/10/Gnirs.csv",row.names=F)
Gnirs<- read.csv("./data/relmatrices/Gnirs/mw/10/Gnirs.csv")
rownames(Gnirs) <- colnames(Gnirs)
Gb <- G.tuneup(G = as.matrix(Gnirs), bend = TRUE, eig.tol = 1e-06)$Gb
Gnirsinv_mw<- G.inverse(G = Gb , sparseform = T)
saveRDS(Gnirsinv_mw, file = "./data/relmatrices/Gnirs/mw/10/Gnirsinv_mw.rds")

# ---------------------25%---------------------
set.seed(6)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.25))
Gnirs<- kin
# --------------------imputing kin removed info with nirs matrix--------------------
Gnirs[remove, ] <- re_nirs_mw[remove, ]
Gnirs[,remove ] <- re_nirs_mw[,remove]
fwrite(Gnirs, "./data/relmatrices/Gnirs/mw/25/Gnirs.csv",row.names=F)
Gnirs<- read.csv("./data/relmatrices/Gnirs/mw/25/Gnirs.csv")
rownames(Gnirs) <- colnames(Gnirs)
Gb <- G.tuneup(G = as.matrix(Gnirs), bend = TRUE, eig.tol = 1e-06)$Gb
Gnirsinv_mw<- G.inverse(G = Gb , sparseform = T)
saveRDS(Gnirsinv_mw, file = "./data/relmatrices/Gnirs/mw/25/Gnirsinv_mw.rds")
# ---------------------50%---------------------
set.seed(9)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.50))
Gnirs<- kin
# --------------------imputing kin removed info with nirs matrix--------------------
Gnirs[remove, ] <- re_nirs_mw[remove, ]
Gnirs[,remove ] <- re_nirs_mw[,remove]
fwrite(Gnirs, "./data/relmatrices/Gnirs/mw/50/Gnirs.csv",row.names=F)
Gnirs<- read.csv("./data/relmatrices/Gnirs/mw/50/Gnirs.csv")
rownames(Gnirs) <- colnames(Gnirs)
Gb <- G.tuneup(G = as.matrix(Gnirs), bend = TRUE, eig.tol = 1e-06)$Gb
Gnirsinv_mw<- G.inverse(G = Gb , sparseform = T)
saveRDS(Gnirsinv_mw, file = "./data/relmatrices/Gnirs/mw/50/Gnirsinv_mw.rds")

# ----replacing random 10%,25%,50% from kinship and replacing 10% by whole wave matrix for joint location---------------------
kin<- fread('./data/relmatrices/kinship_additive.txt', data.table= F)
re_w_joint<- read.csv("./data/relmatrices/re_w_joint.csv")
# ---------------------10%---------------------
set.seed(7)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.1))
GWW<- kin
# --------------------imputing kin removed info with WW matrix---------------------
GWW[remove, ] <- re_w_joint[remove, ]
GWW[,remove ] <- re_w_joint[,remove]
fwrite(GWW, "./data/relmatrices/GWW/joint/10/GWW.csv",row.names=F)
GWW<- read.csv("./data/relmatrices/GWW/joint/10/GWW.csv")
rownames(GWW) <- colnames(GWW)
Gb <- G.tuneup(G = as.matrix(GWW), bend = TRUE, eig.tol = 1e-06)$Gb
GWWinv_joint<- G.inverse(G = Gb , sparseform = T)
saveRDS(GWWinv_joint, file = "./data/relmatrices/GWW/joint/10/GWWinv_joint.rds")

# ---------------------25%---------------------
set.seed(6)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.25))
GWW <- kin
# --------------------imputing kin removed info with WW matrix---------------------
GWW[remove, ] <- re_w_joint[remove, ]
GWW[,remove ] <- re_w_joint[,remove]
fwrite(GWW, "./data/relmatrices/GWW/joint/25/GWW.csv",row.names=F)
GWW<- read.csv("./data/relmatrices/GWW/joint/25/GWW.csv")
rownames(GWW) <- colnames(GWW)
Gb <- G.tuneup(G = as.matrix(GWW), bend = TRUE, eig.tol = 1e-06)$Gb
GWWinv_joint<- G.inverse(G = Gb , sparseform = T)
saveRDS(GWWinv_joint, file = "./data/relmatrices/GWW/joint/25/GWWinv_joint.rds")

# ---------------------50%---------------------
set.seed(56)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.50))
GWW <- kin
# --------------------imputing kin removed info with WW matrix---------------------
GWW[remove, ] <- re_w_joint[remove, ]
GWW[,remove ] <- re_w_joint[,remove]
fwrite(GWW, "./data/relmatrices/GWW/joint/50/GWW.csv",row.names=F)
GWW<- read.csv("./data/relmatrices/GWW/joint/50/GWW.csv")
rownames(GWW) <- colnames(GWW)
Gb <- G.tuneup(G = as.matrix(GWW), bend = TRUE, eig.tol = 1e-06)$Gb
GWWinv_joint<- G.inverse(G = Gb , sparseform = T)
saveRDS(GWWinv_joint, file = "./data/relmatrices/GWW/joint/50/GWWinv_joint.rds")

# ----replacing random 10%,25%,50% from kinship and replacing 10% by wholewave matrix for ef loc---------------------
kin<- fread('./data/relmatrices/kinship_additive.txt', data.table= F)
re_w_ef<- read.csv("./data/relmatrices/re_w_ef.csv")
# ---------------------10%---------------------
set.seed(7)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.1))
# --------------------imputing kin removed info with WW matrix---------------------
GWW<- kin
GWW[remove, ] <- re_w_ef[remove, ]
GWW[,remove ] <- re_w_ef[,remove]
fwrite(GWW, "./data/relmatrices/GWW/ef/10/GWW.csv",row.names=F)
GWW<- read.csv("./data/relmatrices/GWW/ef/10/GWW.csv")
rownames(GWW) <- colnames(GWW)
Gb <- G.tuneup(G = as.matrix(GWW), bend = TRUE, eig.tol = 1e-06)$Gb
GWWinv_ef<- G.inverse(G = Gb , sparseform = T)
saveRDS(GWWinv_ef, file = "./data/relmatrices/GWW/ef/10/GWWinv_ef.rds")

# ---------------------25%---------------------
set.seed(6)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.25))
GWW<- kin
# --------------------imputing kin removed info with WW matrix---------------------
GWW[remove, ] <- re_w_ef[remove, ]
GWW[,remove ] <- re_w_ef[,remove]
fwrite(GWW, "./data/relmatrices/GWW/ef/25/GWW.csv",row.names=F)
GWW<- read.csv("./data/relmatrices/GWW/ef/25/GWW.csv")
rownames(GWW) <- colnames(GWW)
Gb <- G.tuneup(G = as.matrix(GWW), bend = TRUE, eig.tol = 1e-06)$Gb
GWWinv_ef<- G.inverse(G = Gb , sparseform = T)
saveRDS(GWWinv_ef, file = "./data/relmatrices/GWW/ef/25/GWWinv_ef.rds")
# ---------------------50%---------------------
set.seed(9)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.50))
GWW<- kin
# --------------------imputing kin removed info with WW matrix---------------------
GWW[remove, ] <- re_w_ef[remove, ]
GWW[,remove ] <- re_w_ef[,remove]
fwrite(GWW, "./data/relmatrices/GWW/ef/50/GWW.csv",row.names=F)
GWW<- read.csv("./data/relmatrices/GWW/ef/50/GWW.csv")
rownames(GWW) <- colnames(GWW)
Gb <- G.tuneup(G = as.matrix(GWW), bend = TRUE, eig.tol = 1e-06)$Gb
GWWinv_ef<- G.inverse(G = Gb , sparseform = T)
saveRDS(GWWinv_ef, file = "./data/relmatrices/GWW/ef/50/GWWinv_ef.rds")

# ----replacing random 10%,25%,50% from kinship and replacing 10% by wholewave matrix for mw loc---------------------
kin<- fread('./data/relmatrices/kinship_additive.txt', data.table= F)
re_w_mw<- read.csv("./data/relmatrices/re_w_MW.csv")
# ---------------------10%---------------------
set.seed(7)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.1))
GWW<- kin
# --------------------imputing kin removed info with WW matrix---------------------
GWW[remove, ] <- re_w_mw[remove, ]
GWW[,remove ] <- re_w_mw[,remove]
fwrite(GWW, "./data/relmatrices/GWW/mw/10/GWW.csv",row.names=F)
GWW<- read.csv("./data/relmatrices/GWW/mw/10/GWW.csv")
rownames(GWW) <- colnames(GWW)
Gb <- G.tuneup(G = as.matrix(GWW), bend = TRUE, eig.tol = 1e-06)$Gb
GWWinv_mw<- G.inverse(G = Gb , sparseform = T)
saveRDS(GWWinv_mw, file = "./data/relmatrices/GWW/mw/10/GWWinv_mw.rds")

# ---------------------25%---------------------
set.seed(6)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.25))
GWW<- kin
# --------------------imputing kin removed info with WW matrix---------------------
GWW[remove, ] <- re_w_mw[remove, ]
GWW[,remove ] <- re_w_mw[,remove]
fwrite(GWW, "./data/relmatrices/GWW/mw/25/GWW.csv",row.names=F)
GWW<- read.csv("./data/relmatrices/GWW/mw/25/GWW.csv")
rownames(GWW) <- colnames(GWW)
Gb <- G.tuneup(G = as.matrix(GWW), bend = TRUE, eig.tol = 1e-06)$Gb
GWWinv_mw<- G.inverse(G = Gb , sparseform = T)
saveRDS(GWWinv_mw, file = "./data/relmatrices/GWW/mw/25/GWWinv_mw.rds")
# ---------------------50%---------------------
set.seed(9)
N <- nrow(kin)
remove <- sample(1:N, size = round(N*0.50))
GWW<- kin
# --------------------imputing kin removed info with WW matrix---------------------
GWW[remove, ] <- re_w_mw[remove, ]
GWW[,remove ] <- re_w_mw[,remove]
fwrite(GWW, "./data/relmatrices/GWW/mw/50/GWW.csv",row.names=F)
GWW<- read.csv("./data/relmatrices/GWW/mw/50/GWW.csv")
rownames(GWW) <- colnames(GWW)
Gb <- G.tuneup(G = as.matrix(GWW), bend = TRUE, eig.tol = 1e-06)$Gb
GWWinv_mw<- G.inverse(G = Gb , sparseform = T)
saveRDS(GWWinv_mw, file = "./data/relmatrices/GWW/mw/50/GWWinv_mw.rds")


#replication 2-creating Gh2, GWW and Gnirs 10%,25% and 50% combination matrix
# ----replacing random 10%,25%,50% from kinship,imputing by high heritable matrix for joint location---------------------
kin<- fread('./data/relmatrices/kinship_additive.txt', data.table= F)
re_h2_joint<- read.csv("./data/relmatrices/re_h2_joint.csv")

Gh2_2<- process_kinship(kin_matrix = kin, wave_matrix = re_h2_joint,remove_percent = 0.1,
                        seed=345,rep=2)

Gh2_2<- process_kinship(kin_matrix = kin, wave_matrix = re_h2_joint,remove_percent = 0.25,
                        seed=345,rep=2)
Gh2_2<- process_kinship(kin_matrix = kin, wave_matrix = re_h2_joint,remove_percent = 0.50,
                        seed=345,rep=2)
# ----replacing random 10%,25%,50% from kinship,imputing by high heritable matrix for ef location---------------------
kin<- fread('./data/relmatrices/kinship_additive.txt', data.table= F)
re_h2_ef<- read.csv("./data/relmatrices/re_h2_ef.csv")

Gh2_2<- process_kinship(kin_matrix = kin, wave_matrix = re_h2_ef,remove_percent = 0.1,
                        seed=345,rep=2)

Gh2_2<- process_kinship(kin_matrix = kin, wave_matrix = re_h2_ef,remove_percent = 0.25,
                        seed=345,rep=2)
Gh2_2<- process_kinship(kin_matrix = kin, wave_matrix = re_h2_ef,remove_percent = 0.50,
                        seed=345,rep=2)
# ----replacing random 10%,25%,50% from kinship,imputing by high heritable matrix for mw location---------------------
kin<- fread('./data/relmatrices/kinship_additive.txt', data.table= F)
re_h2_mw<- read.csv("./data/relmatrices/re_h2_mw.csv")

Gh2_2<- process_kinship(kin_matrix = kin, wave_matrix = re_h2_mw,remove_percent = 0.1,
                        seed=345,rep=2)

Gh2_2<- process_kinship(kin_matrix = kin, wave_matrix = re_h2_mw,remove_percent = 0.25,
                        seed=345,rep=2)
Gh2_2<- process_kinship(kin_matrix = kin, wave_matrix = re_h2_mw,remove_percent = 0.50,
                        seed=345,rep=2)
# ----replacing random 10%,25%,50% from kinship,imputing by wholewave matrix for joint location---------------------
kin<- fread('./data/relmatrices/kinship_additive.txt', data.table= F)
re_w_joint<- read.csv("./data/relmatrices/re_w_joint.csv")

GWW_2<- process_kinship(kin_matrix = kin, wave_matrix = re_w_joint,remove_percent = 0.1,
                        seed=345,rep=2, type= "WW")

GWW_2<- process_kinship(kin_matrix = kin, wave_matrix = re_w_joint,remove_percent = 0.25,
                        seed=345,rep=2, type= "WW")
GWW_2<- process_kinship(kin_matrix = kin, wave_matrix = re_w_joint,remove_percent = 0.50,
                        seed=345,rep=2, type = "WW")
# ----replacing random 10%,25%,50% from kinship,imputing by high heritable matrix for ef location---------------------
kin<- fread('./data/relmatrices/kinship_additive.txt', data.table= F)
re_w_ef<- read.csv("./data/relmatrices/re_w_ef.csv")

GWW_2<- process_kinship(kin_matrix = kin, wave_matrix = re_w_ef,remove_percent = 0.1,
                        seed=345,rep=2,type="WW")

GWW_2<- process_kinship(kin_matrix = kin, wave_matrix = re_w_ef,remove_percent = 0.25,
                        seed=345,rep=2, type = "WW")
GWW_2<- process_kinship(kin_matrix = kin, wave_matrix = re_w_ef,remove_percent = 0.50,
                        seed=345,rep=2, type= "WW")
# ----replacing random 10%,25%,50% from kinship,imputing by whole wave matrix for mw location---------------------
kin<- fread('./data/relmatrices/kinship_additive.txt', data.table= F)
re_w_mw<- read.csv("./data/relmatrices/re_w_mw.csv")

GWW_2<- process_kinship(kin_matrix = kin, wave_matrix = re_w_mw,remove_percent = 0.1,
                        seed=345,rep=2, type= "WW")

GWW_2<- process_kinship(kin_matrix = kin, wave_matrix = re_w_mw,remove_percent = 0.25,
                        seed=345,rep=2, type="WW")
GWW_2<- process_kinship(kin_matrix = kin, wave_matrix = re_w_mw,remove_percent = 0.50,
                        seed=345,rep=2, type= "WW")
# ----replacing random 10%,25%,50% from kinship,imputing by nirs matrix for joint location---------------------
kin<- fread('./data/relmatrices/kinship_additive.txt', data.table= F)
re_nirs_joint<- read.csv("./data/relmatrices/re_nirs_joint.csv")

Gnirs_2<- process_kinship(kin_matrix = kin, wave_matrix = re_nirs_joint,remove_percent = 0.1,
                        seed=345,rep=2)
Gnirs_2<- process_kinship(kin_matrix = kin, wave_matrix = re_nirs_joint,remove_percent = 0.25,
                          seed=345,rep=2)
Gnirs_2<- process_kinship(kin_matrix = kin, wave_matrix = re_nirs_joint,remove_percent = 0.50,
                          seed=345,rep=2)

# ----replacing random 10%,25%,50% from kinship,imputing by nirs matrix for ef location---------------------
kin<- fread('./data/relmatrices/kinship_additive.txt', data.table= F)
re_nirs_ef<- read.csv("./data/relmatrices/re_nirs_ef.csv")
Gnirs_2<- process_kinship(kin_matrix = kin, wave_matrix = re_nirs_ef,remove_percent = 0.1,
                          seed=345,rep=2)
Gnirs_2<- process_kinship(kin_matrix = kin, wave_matrix = re_nirs_ef,remove_percent = 0.25,
                          seed=345,rep=2)
Gnirs_2<- process_kinship(kin_matrix = kin, wave_matrix = re_nirs_ef,remove_percent = 0.50,
                          seed=345,rep=2)
# ----replacing random 10%,25%,50% from kinship,imputing b nirs matrix for mw location---------------------
re_nirs_mw<- read.csv("./data/relmatrices/re_nirs_mw.csv")
Gnirs_2<- process_kinship(kin_matrix = kin, wave_matrix = re_nirs_mw,remove_percent = 0.1,
                          seed=345,rep=2)
Gnirs_2<- process_kinship(kin_matrix = kin, wave_matrix = re_nirs_mw,remove_percent = 0.25,
                          seed=345,rep=2)
Gnirs_2<- process_kinship(kin_matrix = kin, wave_matrix = re_nirs_mw,remove_percent = 0.50,
                          seed=345,rep=2)


#rep-3


