# ---------------------loading packages---------------------
library(data.table)
library(tidyverse)
library(ASRgenomics)
library(asreml)
source("./script/aux_functions.R")
# ---------------------reading the data for wavelength with blues and wavelength that had na(which didn't converged)---------------------
wavebluesef_wider<- wavebluesEF_binded %>% pivot_wider(names_from = "wave",values_from = "predicted.value")
wavebluesef_wider<- wavebluesef_wider[-855,]
fwrite(wavebluesef_wider,"./output/wavebluesef_wider.csv")
# ---------------------left joining corrected name with name2 from wavelength blues to filter individual that are also present in genotypic data ---------------------
Names_WEST<- fread("./data/Names_WEST_SF.csv", data.table = F)
wavebluesef_wider<- fread("./output/wavebluesef_wider.csv")
Names_WEST$Name2 <- gsub(" ", "", Names_WEST$Name2)
colnames(Names_WEST)[colnames(Names_WEST)== "Name2"]<- "name2"
wavebluesef_wider <- wavebluesef_wider |> left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>% mutate(taxa = ifelse((name2 != Corrected_names) & grepl("PI", Corrected_names), NA, Corrected_names))
wavebluesef_wider<- subset(wavebluesef_wider, !is.na(taxa))
length(unique(wavebluesef_wider$taxa))
wavebluesef_wider <- wavebluesef_wider %>% select(-name2, -Corrected_names) 
wavebluesef_wider<- relocate(wavebluesef_wider, taxa,.before = "wave_350")
colnames(wavebluesef_wider)[colnames(wavebluesef_wider)== "taxa"]<- "name2"
#fwrite(wavebluesef_wider_corrected, "./output/wavebluesef_wider_corrected.csv")
#reading wavebluesef data
wavebluesef_wider<- fread("./output/wavebluesef_wider_corrected.csv")
# ---------------------matrix for whole spectra---------------------
re_w_ef<- calculate_relationship(wavebluesef_wider_corrected)
# ---------------------saving relationship matrix for whole location---------------------
write.csv(re_w_ef,'./relmatrices/re_w_ef.csv',row.names= F)
re_w_ef<- read.csv('./relmatrices/re_w_ef.csv')
# ---------------------creating relationship matrix for nirs data---------------------
wavenirs<-wavebluesef_wider_corrected %>%select(name2,"wave_860":"wave_1660")
re_nirs_ef<- calculate_relationship(wavenirs)
#write.csv(re_nirs_ef,'./relmatrices/re_nirs_ef.csv', row.names= F)
re_nirs_ef<- read.csv('./relmatrices/re_nirs_ef.csv')
# ---------------------for selecting genotype with higher heritability---------------------
# ---------------------reading h2 for ef location---------------------
h2ef<-fread('./output/h2EF.csv',data.table = F) #reading h2 data
h2ef<- h2ef[-2152,]
#visualizing trend in plot with wavelegth vs heritability 
h2ef<- h2ef %>% mutate(wavenumber= as.numeric(str_replace(wave,"wave_", "")))
ggplot(data= h2ef, aes(x = wavenumber, y = h2)) + geom_point(color= "purple", size= 0.1) + labs(title= "heritability per wavelength", x= "wavenumber",y= "heritability")
ggplot(data = h2ef[h2ef$h2 > 0.3, ], aes(x = wavenumber, y = h2)) + geom_point(color= "purple", size= 0.1) + labs(title= "heritability per wavelength", x= "wavenumber",y= "heritability")
# ---------------------filtering wavelength that are highly heritable according to the trend and the threshold---------------------
highh2wave<- h2ef[h2ef$h2 > 0.3, c("wavenumber", "h2")]
high2wavedata<- right_join(h2ef,highh2wave, by= "wavenumber")
high2wavedata<- high2wavedata[,-3]
high2wavedata<-high2wavedata[,-3]
wave_high2<-high2wavedata$wave
# ---------------------filtering bluesef data for wavelength that are highly heritable---------------------
wavebluesef_corrected<- fread("./output/wavebluesef_corrected.csv")
wavebluesef_corrected<- pivot_longer(names_to = "wave", values_to = "predicted.value")
highwaveblues<- wavebluesef_corrected %>%  filter(wave %in% wave_high2)
highwaveblues<- highwaveblues %>% pivot_wider(names_from = "wave",values_from = "predicted.value")
# ---------------------saving the high heritable waveblues ---------------------
fwrite(highwaveblues, "./intermediate/highwaveblues.csv", row.names = F)
# ---------------------high heritable relationship matrix---------------------
highwaveblues<- fread("./intermediate/highwaveblues.csv")
highwaveblues_filtered<- highwaveblues%>%
  filter(name2 %in% wavebluesef_wider_corrected$name2)
re_h2_ef<- calculate_relationship(highwaveblues_filtered)
write.csv(re_h2_ef, "./relmatrices/re_h2_ef.csv", row.names = F)
re_h2_ef<-read.csv("./relmatrices/re_h2_ef.csv")

# ---------------------calculating relationship matrix for "mw" location---------------------
# wavebluesmw<- fread("./output/wavebluesmw.csv")
# wavebluesmw_wider <- wavebluesmw %>% pivot_wider(names_from = "wave",values_from = "predicted.value")
# wavebluesmw_wider<- wavebluesmw_wider[-855,]
# fwrite(wavebluesmw_wider, "./output/wavebluesmw_wider.csv",row.names = F)

# ---------------------reading mw blues data---------------------

wavebluesmw_wider<- fread("./output/wavebluesmw_wider.csv")
wavebluesmw_wider_corrected<- wavebluesmw_wider |> left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>% mutate(taxa = ifelse((name2 != Corrected_names) & grepl("PI", Corrected_names), NA, Corrected_names))
wavebluesmw_corrected<- subset(wavebluesmw_corrected, !is.na(taxa))
wavebluesmw_wider_corrected<- wavebluesmw_wider_corrected %>% select(-name2, -Corrected_names) 
wavebluesmw_wider<- relocate(wavebluesmw_wider, taxa,.before = "wave_350")
colnames(wavebluesmw_wider_corrected)[colnames(wavebluesmw_wider_corrected)== "taxa"]<- "name2"
#fwrite(wavebluesmw_wider_corrected, "./output/wavebluesmw_wider_corrected.csv", row.names = F)
#generating relationship matrix for whole spectra

re_w_mw <- calculate_relationship(wavebluesmw_wider_corrected)
write.csv(re_w_mw,'./relmatrices/re_w_mw.csv',row.names= F)
re_w_mw<- read.csv('./relmatrices/re_w_mw.csv')
# ---------------------relationship matrix for NIRS wavelength---------------------

wavenirsmw<- wavebluesmw_wider_corrected%>% select(name2,"wave_860":"wave_1660")
fwrite(wavenirsmw,"./intermediate/wavenirsmw.csv")

re_nirs_mw<- calculate_relationship(wavenirsmw)

write.csv(re_nirs_mw,'./relmatrices/re_nirs_mw.csv', row.names = F)
re_nirs_mw<- read.csv('./relmatrices/re_nirs_mw.csv')

#for selecting genotype with higher heritability
h2mw<-fread('./output/h2MW.csv',data.table = F)#reading h2 data
#visualizing trend in plot with wavelegth vs heritability 
h2mw<- h2mw %>% mutate(wavenumber= as.numeric(str_replace(wave,"wave_", "")))
ggplot(data= h2mw, aes(x = wavenumber, y = h2)) + geom_point(color= "purple", size= 0.1) + labs(title= "heritability per wavelength", x= "wavenumber",y= "heritability")
ggplot(data = h2mw[h2mw$h2 > 0.4, ], aes(x = wavenumber, y = h2)) + geom_point(color= "purple", size= 0.1) + labs(title= "heritability per wavelength", x= "wavenumber",y= "heritability")

highh2wave<- h2mw[h2mw$h2 > 0.4, c( "h2","wave")]

wavebluesmw<- fread("./output/wavebluesmw.csv")
wave_high2<-highh2wave$wave
highwaveblues<- wavebluesmw%>% filter(wave %in% wave_high2)

highwaveblues<- highwaveblues %>% pivot_wider(names_from = "wave",values_from = "predicted.value")
highwaveblues<- highwaveblues[-855, ]

# ---------------------high heritable relationship matrix---------------------
re_h2_mw<- calculate_relationship(highwaveblues)

write.csv(re_h2_mw, "./relmatrices/re_h2_mw.csv", row.names = F)
re_h2_mw<-read.csv("./relmatrices/re_h2_mw.csv")


# ---------------------creating three different relationship matrix for joint location---------------------
# ---------------------preprocessing of wavejoint data---------------------
wavejoint<- fread("./output/waveblues_sbf.csv")
wavejoint<- wavejoint[-855,]
wavejoint_long<- pivot_longer(wavejoint, cols = 2:2152,names_to = "wave",values_to = "predicted.value") %>% 
  mutate(temp = as.integer(str_sub(wave , start = 6)))%>% 
  arrange(temp, name2) %>% 
  select(-temp)
na_wave <- subset(wavejoint_long, is.na(wavejoint_long$predicted.value)) 
# ---------------------binding na wave ---------------------
nawaveblues<- fread("./intermediate/nawaveblues.csv")

nawaveblues_wide<- pivot_wider(nawaveblues, values_from = "predicted.value", names_from = "wave")

wavejoint_long <- drop_na(wavejoint_long)

wavebluesjoint<- rbind(wavejoint_long,nawaveblues) %>% mutate(temp = as.integer(str_sub(wave , start = 6)))%>% 
  arrange(temp, name2) %>% 
  select(-temp)

wavebluesjoint_wider<- wavebluesjoint %>% pivot_wider(names_from = "wave",values_from = "predicted.value")
fwrite(wavebluesjoint_wider, "./output/wavebluesjoint.csv", row.names = F)
# ---------------------reading the processed wavejoint data---------------------
wavebluesjoint<- fread("./output/wavebluesjoint.csv")
wavebluesjoint_corrected<- wavebluesjoint|> left_join(Names_WEST %>% dplyr::select(name2, Corrected_names)) %>% mutate(taxa = ifelse((name2 != Corrected_names) & grepl("PI", Corrected_names), NA, Corrected_names))
wavebluesjoint_corrected<- subset(wavebluesjoint_corrected, !is.na(taxa))
length(unique(wavebluesjoint_corrected$taxa))
wavebluesjoint_corrected <- wavebluesjoint_corrected %>% select(-name2, -Corrected_names) 
wavebluesjoint_corrected<- relocate(wavebluesjoint_corrected, taxa,.before = "wave_350")
colnames(wavebluesjoint_corrected)[colnames(wavebluesjoint_corrected)== "taxa"]<- "name2"
fwrite(wavebluesjoint_corrected, "./output/wavebluesjoint_corrected.csv", row.names = F)

wavebluesjoint_long<- pivot_longer(wavebluesjoint_corrected,cols= 2:2152,values_to = "predicted.values", names_to = "wave")

# ---------------------calculating joint whole spectra relationship matrix---------------------

re_w_joint<- calculate_relationship(wavebluesjoint_corrected)
write.csv(re_w_joint, "./relmatrices/re_w_joint.csv", row.names = F)
re_w_joint<- read.csv("./relmatrices/re_w_joint.csv")



# ---------------------calculating NIRS matrix for joint model---------------------

wavenirsjoint<- wavebluesjoint_corrected %>% select(name2,"wave_860":"wave_1660")
fwrite(wavenirsjoint,"./intermediate/wavenirsjoint.csv")

re_nirs_joint<- calculate_relationship(wavenirsjoint)

write.csv(re_nirs_joint,'./relmatrices/re_nirs_joint.csv', row.names = F)
re_nirs_joint<- read.csv('./relmatrices/re_nirs_joint.csv')
# ---------------------calculating highly h2 matrix---------------------
h2joint<- fread("./output/h2.csv", data.table = FALSE)

# ---------------------visualizing trend in plot with wavelegth vs heritability ---------------------
h2joint<- h2joint %>% mutate(wavenumber= as.numeric(str_replace(wave,"wave_", "")))
ggplot(data= h2joint, aes(x = wavenumber, y = h2)) + geom_point(color= "purple", size= 0.1) + labs(title= "heritability per wavelength", x= "wavenumber",y= "heritability")
ggplot(data = h2joint[h2joint$h2 > 0.4, ], aes(x = wavenumber, y = h2)) + geom_point(color= "purple", size= 0.1) + labs(title= "heritability per wavelength", x= "wavenumber",y= "heritability")

highh2wave_joint<- h2joint[h2joint$h2 > 0.4, c( "h2","wave")]
wave_high2<-highh2wave_joint$wave
highwaveblues_joint<- wavebluesjoint_long %>% filter(wave %in% wave_high2)
highwaveblues<- highwaveblues_joint %>% pivot_wider(names_from = "wave",values_from = "predicted.values")

# ---------------------calculating relationship matrix for highly h2 wavelength for joint location---------------------

re_h2_joint<- calculate_relationship(highwaveblues)

write.csv(re_h2_joint, "./relmatrices/re_h2_joint.csv", row.names = F)
re_h2_joint<- read.csv("./relmatrices/re_h2_joint.csv")


