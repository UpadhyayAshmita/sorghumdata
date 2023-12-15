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
Gb <- G.tuneup(G = as.matrix(kin), bend = TRUE, eig.tol = 1e-06)$Gb
GINV <- G.inverse(G = Gb , sparseform = T)
saveRDS(GINV, file = "./data/relmatrices/GINV.rds")
GINV <- readRDS(file = "./data/relmatrices/GINV.rds")
# ---------------------reading the data for wavelength with blues and wavelength that had na(which didn't converged)---------------------
wavebluesEF_wider<- fread("./output/wavebluesEF_wider.csv")
# ---------------------left joining corrected name with name2 from wavelength blues to filter individual that are also present in genotypic data ---------------------
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
# ---------------------creating relationship matrix for nirs data---------------------
wavenirs<-wavebluesEF_wider_corrected %>%select(taxa,"wave_860":"wave_1660")
re_nirs_ef<- calculate_relationship(wavenirs)
write.csv(re_nirs_ef,'./data/relmatrices/re_nirs_ef.csv', row.names= F)
re_nirs_ef<- read.csv('./data/relmatrices/re_nirs_ef.csv')
# ---------------------for selecting genotype with higher heritability---------------------
# ---------------------reading h2 for ef location---------------------
h2ef<-fread('./output/h2EF.csv',data.table = F) #reading h2 data
#h2ef<- h2ef[-2152,]
#visualizing trend in plot with wavelegth vs heritability 
h2ef<- h2ef %>% mutate(wavenumber= as.numeric(str_replace(wave,"wave_", "")))
ggplot(data= h2ef, aes(x = wavenumber, y = h2)) + geom_point(color= "purple", size= 0.1) + labs(title= "heritability per wavelength", x= "wavenumber",y= "heritability")
ggplot(data = h2ef[h2ef$h2 > 0.3, ], aes(x = wavenumber, y = h2)) + geom_point(color= "purple", size= 0.1) + labs(title= "heritability per wavelength", x= "wavenumber",y= "heritability")
# ---------------------filtering wavelength that are highly heritable according to the trend and the threshold---------------------
highh2wave<- h2ef[h2ef$h2 > 0.3, c( "h2","wave")]
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
# ---------------------relationship matrix for NIRS wavelength---------------------
wavenirsmw<- wavebluesMW_wider_corrected%>% select(taxa,"wave_860":"wave_1660")
fwrite(wavenirsmw,"./output/intermediate/wavenirsmw.csv")
re_nirs_mw<- calculate_relationship(wavenirsmw)
write.csv(re_nirs_mw,'./data/relmatrices/re_nirs_mw.csv', row.names = F)
re_nirs_mw<- read.csv('./data/relmatrices/re_nirs_mw.csv')
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
# ---------------------calculating NIRS matrix for joint model---------------------
wavenirsjoint<- wavebluesjoint_corrected %>% 
  select(taxa,"wave_860":"wave_1660")
fwrite(wavenirsjoint,"./output/intermediate/wavenirsjoint.csv")
re_nirs_joint<- calculate_relationship(wavenirsjoint)
write.csv(re_nirs_joint,'./data/relmatrices/re_nirs_joint.csv', row.names = F)
re_nirs_joint<- read.csv('./data/relmatrices/re_nirs_joint.csv')
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


