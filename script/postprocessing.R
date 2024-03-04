#loading apckages
library(dplyr)
library(RColorBrewer)
source("./script/aux_functions.R")
#selecting topk ind from loc ef and mw by reading blues from first stage
#top 20% indv for slablues ef and their mean
temp_data <- slablues_ef_correct %>% arrange(desc(sla)) %>% na.omit()
num_indv<- 165
n1 <- nrow(temp_data)  # total number of genotypes in data 
topk1 <- min(n1, num_indv)  # Select minimum of 165 and the number of genotypes
taxa_blues <- temp_data[1:topk1, ]
mean_value <- mean(taxa_blues[["sla"]], na.rm = TRUE)#309.0880
top_taxa_names <- taxa_blues$taxa
write.csv(top_taxa_names, "./output/selected/top_slablues_ef.csv")

#top 20% indv for slablues mw and their mean
temp_data <- slablues_mw_correct %>% arrange(desc(sla)) %>% na.omit()
num_indv<- 165
n1 <- nrow(temp_data)  # total number of genotypes in data 
topk1 <- min(n1, num_indv)  # Select minimum of 165 and the number of genotypes
taxa_blues <- temp_data[1:topk1, ]
mean_value <- mean(taxa_blues[["sla"]], na.rm = TRUE)#312.0316
top_taxa_names <- taxa_blues$taxa
write.csv(top_taxa_names, "./output/selected/top_slablues_mw.csv")

#top 20% indv for nareablues ef and their mean
temp_data <- nareablues_ef_correct %>% arrange(desc(narea)) %>% na.omit()
num_indv<- 165
n1 <- nrow(temp_data)  # total number of genotypes in data 
topk1 <- min(n1, num_indv)  # Select minimum of 165 and the number of genotypes
taxa_blues <- temp_data[1:topk1, ]
mean_value <- mean(taxa_blues[["narea"]], na.rm = TRUE)#1.5611
top_taxa_names <- taxa_blues$taxa
write.csv(top_taxa_names, "./output/selected/top_nareablues_ef.csv")

#top 20% indv for nareablues mw and their mean
temp_data <- nareablues_mw_correct %>% arrange(desc(narea)) %>% na.omit()
num_indv<- 165
n1 <- nrow(temp_data)  # total number of genotypes in data 
topk1 <- min(n1, num_indv)  # Select minimum of 165 and the number of genotypes
taxa_blues <- temp_data[1:topk1, ]
mean_value <- mean(taxa_blues[["narea"]], na.rm = TRUE)#1.56841
top_taxa_names <- taxa_blues$taxa
write.csv(top_taxa_names, "./output/selected/top_nareablues_mw.csv")

#selecting top20%indv and their means for narea ef from GBLUP_efmw
GBLUP_efmw_narea<- read.csv("./output/GBLUP/narea_efmw.csv") %>% na.omit()
GBLUP_efmw_narea<- top_selected(data= GBLUP_efmw_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_efmw_narea$result_top_taxa), each = length(GBLUP_efmw_narea$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_efmw_narea$result_top_taxa))

GBLUP_topn_efmw <- merge(selected_df,nareablues_ef_correct, by ='taxa', all.x = "TRUE")
write.csv(GBLUP_topn_efmw,"./output/selected/GBLUP_topn_efmw.csv", row.names=F)

#selecting top20%indv and their means for narea mw from GBLUP_mwef
GBLUP_mwef_narea<- read.csv("./output/GBLUP/narea_mwef.csv")
GBLUP_mwef_narea<- top_selected(data= GBLUP_mwef_narea)
#extracting mean and selected individual in a data frame per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_mwef_narea$result_top_taxa), each = length(GBLUP_mwef_narea$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_mwef_narea$result_top_taxa))
GBLUP_topn_mwef <- merge(selected_df,nareablues_mw_correct, by ='taxa', all.x = "TRUE")
write.csv(GBLUP_topn_mwef,"./output/selected/GBLUP_topn_mwef.csv", row.names=F)

#selecting top20%indv and their means for sla ef from GBLUP_efmw
GBLUP_efmw_sla<- read.csv("./output/GBLUP/sla_efmw.csv")
GBLUP_efmw_sla<- top_selected(data= GBLUP_efmw_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_efmw_sla$result_top_taxa), each = length(GBLUP_efmw_sla$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_efmw_sla$result_top_taxa))
GBLUP_tops_efmw<- merge(selected_df,slablues_ef_correct, by ='taxa', all.x = "TRUE")
write.csv(GBLUP_tops_efmw,"./output/selected/GBLUP_tops_efmw.csv", row.names=F)

#selecting top20%indv and their means for sla mw from GBLUP_mwef
GBLUP_mwef_sla<- read.csv("./output/GBLUP/sla_mwef.csv")
GBLUP_mwef_sla<- top_selected(data= GBLUP_mwef_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_mwef_sla$result_top_taxa), each = length(GBLUP_mwef_sla$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_mwef_sla$result_top_taxa))
GBLUP_tops_mwef <- merge(selected_df,slablues_mw_correct, by ='taxa', all.x = "TRUE")
write.csv(GBLUP_tops_mwef,"./output/selected/GBLUP_tops_mwef.csv", row.names=F)

#selecting top20%indv and their means for sla ef from Gh2_efmw
Gh2_efmw10_sla<- read.csv("./output/Gh2_efmw_sla/sla_Gh2_10_efmw.csv")
Gh2_efmw10_sla<- top_selected(data= Gh2_efmw10_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gh2_efmw10_sla$result_top_taxa), each = length(Gh2_efmw10_sla$result_top_taxa[[1]])),
                          taxa = unlist(Gh2_efmw10_sla$result_top_taxa))
Gh2_tops10_efmw <- merge(selected_df,slablues_ef_correct, by ='taxa', all.x = "TRUE")
write.csv(Gh2_tops10_efmw,"./output/selected/Gh2_tops10_efmw.csv", row.names=F)

#selecting top20%indv and their means for sla mw from Gh2_mwef
Gh2_mwef10_sla<- read.csv("./output/Gh2_mwef_sla/sla_Gh2_10_mwef.csv")
Gh2_mwef10_sla<- top_selected(data= Gh2_mwef10_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gh2_mwef10_sla$result_top_taxa), each = length(Gh2_mwef10_sla$result_top_taxa[[1]])),
                          taxa = unlist(Gh2_mwef10_sla$result_top_taxa))
Gh2_tops10_mwef<- merge(selected_df,slablues_mw_correct, by ='taxa', all.x = "TRUE")

write.csv(Gh2_tops10_mwef,"./output/selected/Gh2_tops10_mwef.csv", row.names=F)

#selecting top20%indv and their means for sla ef from Gh2_efmw
Gh2_efmw25_sla<- read.csv("./output/Gh2_efmw_sla/sla_Gh2_25_efmw.csv")
Gh2_efmw25_sla<- top_selected(data= Gh2_efmw25_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gh2_efmw25_sla$result_top_taxa), each = length(Gh2_efmw25_sla$result_top_taxa[[1]])),
                          taxa = unlist(Gh2_efmw25_sla$result_top_taxa))

Gh2_tops25_efmw<- merge(selected_df,slablues_ef_correct, by ='taxa', all.x = "TRUE")
write.csv(Gh2_tops25_efmw,"./output/selected/Gh2_tops25_efmw.csv", row.names=F)

#selecting top20%indv and their means for sla mw from Gh2_mwef
Gh2_mwef25_sla<- read.csv("./output/Gh2_mwef_sla/sla_Gh2_25_mwef.csv")
Gh2_mwef25_sla<- top_selected(data= Gh2_mwef25_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gh2_mwef25_sla$result_top_taxa), each = length(Gh2_mwef25_sla$result_top_taxa[[1]])),
                          taxa = unlist(Gh2_mwef25_sla$result_top_taxa))
Gh2_tops25_mwef <- merge(selected_df, slablues_mw_correct, by ='taxa', all.x = "TRUE")
write.csv(Gh2_tops25_mwef,"./output/selected/Gh2_tops25_mwef.csv", row.names=F)

#selecting top20%indv and their means for sla ef from Gh2_efmw
Gh2_efmw50_sla<- read.csv("./output/Gh2_efmw_sla/sla_Gh2_50_efmw.csv")
Gh2_efmw50_sla<- top_selected(data= Gh2_efmw50_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gh2_efmw50_sla$result_top_taxa), each = length(Gh2_efmw50_sla$result_top_taxa[[1]])),
                          taxa = unlist(Gh2_efmw50_sla$result_top_taxa))
Gh2_tops50_efmw <- merge(selected_df,slablues_ef_correct, by ='taxa', all.x = "TRUE" )
write.csv(Gh2_tops50_efmw,"./output/selected/Gh2_tops50_efmw.csv", row.names=F)

#selecting top20%indv and their means for sla mw from Gh2_mwef
Gh2_mwef50_sla<- read.csv("./output/Gh2_mwef_sla/sla_Gh2_50_mwef.csv")
Gh2_mwef50_sla<- top_selected(data= Gh2_mwef50_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gh2_mwef50_sla$result_top_taxa), each = length(Gh2_mwef50_sla$result_top_taxa[[1]])),
                          taxa = unlist(Gh2_mwef50_sla$result_top_taxa))
Gh2_tops50_mwef <- merge(selected_df, slablues_mw_correct, by ='taxa', all.x = "TRUE" )
write.csv(Gh2_tops50_mwef,"./output/selected/Gh2_tops50_mwef.csv", row.names=F)

#selecting top20%indv and their means for narea ef from Gh2_efmw
Gh2_efmw10_narea<- read.csv("./output/Gh2_efmw_narea/narea_Gh2_10_efmw.csv")
Gh2_efmw10_narea<- top_selected(data= Gh2_efmw10_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gh2_efmw10_narea$result_top_taxa), each = length(Gh2_efmw10_narea$result_top_taxa[[1]])),
                          taxa = unlist(Gh2_efmw10_narea$result_top_taxa))
Gh2_topn10_efmw <- merge(selected_df,nareablues_ef_correct, by ='taxa', all.x = "TRUE")
write.csv(Gh2_topn10_efmw,"./output/selected/Gh2_topn10_efmw.csv", row.names=F)

#selecting top20%indv and their means for narea mw from Gh2_mwef
Gh2_mwef10_narea<- read.csv("./output/Gh2_mwef_narea/narea_Gh2_10_mwef.csv")
Gh2_mwef10_narea<- top_selected(data= Gh2_mwef10_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gh2_mwef10_narea$result_top_taxa), each = length(Gh2_mwef10_narea$result_top_taxa[[1]])),
                          taxa = unlist(Gh2_mwef10_narea$result_top_taxa))
Gh2_topn10_mwef <- merge(selected_df, nareablues_mw_correct, by ='taxa', all.x = "TRUE" )
write.csv(Gh2_topn10_mwef,"./output/selected/Gh2_topn10_mwef.csv", row.names=F)



#selecting top20%indv and their means for narea ef from Gh2_efmw
Gh2_efmw25_narea<- read.csv("./output/Gh2_efmw_narea/narea_Gh2_25_efmw.csv")
Gh2_efmw25_narea<- top_selected(data= Gh2_efmw25_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gh2_efmw25_narea$result_top_taxa), each = length(Gh2_efmw25_narea$result_top_taxa[[1]])),
                          taxa = unlist(Gh2_efmw25_narea$result_top_taxa))
Gh2_topn25_efmw <- merge(selected_df, nareablues_ef_correct, by ='taxa', all.x = "TRUE" )
write.csv(Gh2_topn25_efmw,"./output/selected/Gh2_topn25_efmw.csv", row.names=F)

#selecting top20%indv and their means for narea mw from Gh2_mwef
Gh2_mwef25_narea<- read.csv("./output/Gh2_mwef_narea/narea_Gh2_25_mwef.csv")
Gh2_mwef25_narea<- top_selected(data= Gh2_mwef25_narea)
#extracting  and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gh2_mwef25_narea$result_top_taxa), each = length(Gh2_mwef25_narea$result_top_taxa[[1]])),
                          taxa = unlist(Gh2_mwef25_narea$result_top_taxa))
Gh2_topn25_mwef <- merge(selected_df, nareablues_mw_correct, by ='taxa', all.x = "TRUE" )
write.csv(Gh2_topn25_mwef,"./output/selected/Gh2_topn25_mwef.csv", row.names=F)
#selecting top20%indv and their means for narea ef from Gh2_efmw
Gh2_efmw50_narea<- read.csv("./output/Gh2_efmw_narea/narea_Gh2_50_efmw.csv")
Gh2_efmw50_narea<- top_selected(data= Gh2_efmw50_narea)
#extracting  and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gh2_efmw50_narea$result_top_taxa), each = length(Gh2_efmw50_narea$result_top_taxa[[1]])),
                          taxa = unlist(Gh2_efmw50_narea$result_top_taxa))
Gh2_topn50_efmw <- merge(selected_df, nareablues_ef_correct, by ='taxa', all.x = "TRUE" )
write.csv(Gh2_topn50_efmw,"./output/selected/Gh2_topn50_efmw.csv", row.names=F)

#selecting top20%indv and their means for narea mw from Gh2_mwef
Gh2_mwef50_narea<- read.csv("./output/Gh2_mwef_narea/narea_Gh2_50_mwef.csv")
Gh2_mwef50_narea<- top_selected(data= Gh2_mwef50_narea)
#extracting  and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gh2_mwef50_narea$result_top_taxa), each = length(Gh2_mwef50_narea$result_top_taxa[[1]])),
                          taxa = unlist(Gh2_mwef50_narea$result_top_taxa))
Gh2_topn50_mwef <- merge(selected_df, nareablues_mw_correct, by ='taxa', all.x = "TRUE" )
write.csv(Gh2_topn50_mwef,"./output/selected/Gh2_topn50_mwef.csv", row.names=F)




#GWW selected individual
#selecting top20%indv and their means for sla ef from GWW efmw
GWW_efmw10_sla<- read.csv("./output/GWW_efmw_sla/sla_GWW_10_efmw.csv")
GWW_efmw10_sla<- top_selected(data= GWW_efmw10_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GWW_efmw10_sla$result_top_taxa), each = length(GWW_efmw10_sla$result_top_taxa[[1]])),
                          taxa = unlist(GWW_efmw10_sla$result_top_taxa))
GWW_tops10_efmw <- merge(selected_df, slablues_ef_correct, by ='taxa', all.x = "TRUE" )
write.csv(GWW_tops10_efmw,"./output/selected/GWW_tops10_efmw.csv", row.names=F)

#selecting top20%indv and their means for sla mw from GWW_mwef
GWW_mwef10_sla<- read.csv("./output/GWW_mwef_sla/sla_GWW_10_mwef.csv")
GWW_mwef10_sla<- top_selected(data= GWW_mwef10_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GWW_mwef10_sla$result_top_taxa), each = length(GWW_mwef10_sla$result_top_taxa[[1]])),
                          taxa = unlist(GWW_mwef10_sla$result_top_taxa))
GWW_tops10_mwef <- merge(selected_df, slablues_mw_correct, by ='taxa', all.x = "TRUE" )
write.csv(GWW_tops10_mwef,"./output/selected/GWW_tops10_mwef.csv", row.names=F)

#selecting top20%indv and their means for sla ef from GWW_efmw
GWW_efmw25_sla<- read.csv("./output/GWW_efmw_sla/sla_GWW_25_efmw.csv")
GWW_efmw25_sla<- top_selected(data= GWW_efmw25_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GWW_efmw25_sla$result_top_taxa), each = length(GWW_efmw25_sla$result_top_taxa[[1]])),
                          taxa = unlist(GWW_efmw25_sla$result_top_taxa))
GWW_tops25_efmw <- merge(selected_df, slablues_ef_correct, by ='taxa', all.x = "TRUE" )
write.csv(GWW_tops25_efmw,"./output/selected/GWW_tops25_efmw.csv", row.names=F)

#selecting top20%indv and their means for sla mw from GWW_mwef
GWW_mwef25_sla<- read.csv("./output/GWW_mwef_sla/sla_GWW_25_mwef.csv")
GWW_mwef25_sla<- top_selected(data= GWW_mwef25_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GWW_mwef25_sla$result_top_taxa), each = length(GWW_mwef25_sla$result_top_taxa[[1]])),
                          taxa = unlist(GWW_mwef25_sla$result_top_taxa))
GWW_tops25_mwef <- merge(selected_df, slablues_mw_correct, by ='taxa', all.x = "TRUE" )
write.csv(GWW_tops25_mwef,"./output/selected/GWW_tops25_mwef.csv", row.names=F)

#selecting top20%indv and their means for sla ef from GWW_efmw
GWW_efmw50_sla<- read.csv("./output/GWW_efmw_sla/sla_GWW_50_efmw.csv")
GWW_efmw50_sla<- top_selected(data= GWW_efmw50_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GWW_efmw50_sla$result_top_taxa), each = length(GWW_efmw50_sla$result_top_taxa[[1]])),
                          taxa = unlist(GWW_efmw50_sla$result_top_taxa))
GWW_tops50_efmw <- merge(selected_df, slablues_ef_correct, by ='taxa', all.x = "TRUE" )
write.csv(GWW_tops50_efmw,"./output/selected/GWW_tops50_efmw.csv", row.names=F)

#selecting top20%indv and their means for sla mw from GWW_mwef
GWW_mwef50_sla<- read.csv("./output/GWW_mwef_sla/sla_GWW_50_mwef.csv")
GWW_mwef50_sla<- top_selected(data= GWW_mwef50_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GWW_mwef50_sla$result_top_taxa), each = length(GWW_mwef50_sla$result_top_taxa[[1]])),
                          taxa = unlist(GWW_mwef50_sla$result_top_taxa))
GWW_tops50_mwef <- merge(selected_df, slablues_mw_correct, by ='taxa', all.x = "TRUE" )
write.csv(GWW_tops50_mwef,"./output/selected/GWW_tops50_mwef.csv", row.names=F)

#selecting top20%indv and their means for narea ef from GWW_efmw
GWW_efmw10_narea<- read.csv("./output/GWW_efmw_narea/narea_GWW_10_efmw.csv")
GWW_efmw10_narea<- top_selected(data= GWW_efmw10_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GWW_efmw10_narea$result_top_taxa), each = length(GWW_efmw10_narea$result_top_taxa[[1]])),
                          taxa = unlist(GWW_efmw10_narea$result_top_taxa))
GWW_topn10_efmw <- merge(selected_df, nareablues_ef_correct, by ='taxa', all.x = "TRUE" )
write.csv(GWW_topn10_efmw,"./output/selected/GWW_topn10_efmw.csv", row.names=F)

#selecting top20%indv and their means for narea mw from GWW_mwef
GWW_mwef10_narea<- read.csv("./output/GWW_mwef_narea/narea_GWW_10_mwef.csv")
GWW_mwef10_narea<- top_selected(data= GWW_mwef10_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GWW_mwef10_narea$result_top_taxa), each = length(GWW_mwef10_narea$result_top_taxa[[1]])),
                          taxa = unlist(GWW_mwef10_narea$result_top_taxa))
GWW_topn10_mwef <- merge(selected_df, nareablues_mw_correct, by ='taxa', all.x = "TRUE")
write.csv(GWW_topn10_mwef,"./output/selected/GWW_topn10_mwef.csv", row.names=F)

#selecting top20%indv and their means for narea mw from GWW_mwef
GWW_mwef10_narea<- read.csv("./output/GWW_mwef_narea/narea_GWW_10_mwef.csv")
GWW_mwef10_narea<- top_selected(data= GWW_mwef10_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GWW_mwef10_narea$result_top_taxa), each = length(GWW_mwef10_narea$result_top_taxa[[1]])),
                          taxa = unlist(GWW_mwef10_narea$result_top_taxa))
GWW_topn10_mwef <- merge(selected_df, nareablues_mw_correct, by ='taxa', all.x = "TRUE")
write.csv(GWW_topn10_mwef,"./output/selected/GWW_topn10_mwef.csv", row.names=F)

#selecting top20%indv and their means for narea ef from GWW_efmw
GWW_efmw25_narea<- read.csv("./output/GWW_efmw_narea/narea_GWW_25_efmw.csv")
GWW_efmw25_narea<- top_selected(data= GWW_efmw25_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GWW_efmw25_narea$result_top_taxa), each = length(GWW_efmw25_narea$result_top_taxa[[1]])),
                          taxa = unlist(GWW_efmw25_narea$result_top_taxa))
GWW_topn25_efmw <- merge(selected_df, nareablues_ef_correct, by ='taxa', all.x = "TRUE")
write.csv(GWW_topn25_efmw,"./output/selected/GWW_topn25_efmw.csv", row.names=F)

#selecting top20%indv and their means for narea mw from GWW_mwef
GWW_mwef25_narea<- read.csv("./output/GWW_mwef_narea/narea_GWW_25_mwef.csv")
GWW_mwef25_narea<- top_selected(data= GWW_mwef25_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GWW_mwef25_narea$result_top_taxa), each = length(GWW_mwef25_narea$result_top_taxa[[1]])),
                          taxa = unlist(GWW_mwef25_narea$result_top_taxa))
GWW_topn25_mwef <- merge(selected_df, nareablues_mw_correct, by ='taxa', all.x = "TRUE")
write.csv(GWW_topn25_mwef,"./output/selected/GWW_topn25_mwef.csv", row.names=F)

#selecting top20%indv and their means for narea ef from GWW_efmw
GWW_efmw50_narea<- read.csv("./output/GWW_efmw_narea/narea_GWW_50_efmw.csv")
GWW_efmw50_narea<- top_selected(data= GWW_efmw50_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GWW_efmw50_narea$result_top_taxa), each = length(GWW_efmw50_narea$result_top_taxa[[1]])),
                          taxa = unlist(GWW_efmw50_narea$result_top_taxa))
GWW_topn50_efmw <- merge(selected_df, nareablues_ef_correct, by ='taxa', all.x = "TRUE")
write.csv(GWW_topn50_efmw,"./output/selected/GWW_topn50_efmw.csv", row.names=F)

#selecting top20%indv and their means for narea mw from GWW_mwef
GWW_mwef50_narea<- read.csv("./output/GWW_mwef_narea/narea_GWW_50_mwef.csv")
GWW_mwef50_narea<- top_selected(data= GWW_mwef50_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(GWW_mwef50_narea$result_top_taxa), each = length(GWW_mwef50_narea$result_top_taxa[[1]])),
                          taxa = unlist(GWW_mwef50_narea$result_top_taxa))
GWW_topn50_mwef <- merge(selected_df, nareablues_mw_correct, by ='taxa', all.x = "TRUE")
write.csv(GWW_topn50_mwef,"./output/selected/GWW_topn50_mwef.csv", row.names=F)

#Gnirs individual selection
#selecting top20%indv and their means for sla ef from Gnirs efmw
Gnirs_efmw10_sla<- read.csv("./output/Gnirs_efmw_sla/sla_Gnirs_10_efmw.csv")
Gnirs_efmw10_sla<- top_selected(data= Gnirs_efmw10_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gnirs_efmw10_sla$result_top_taxa), each = length(Gnirs_efmw10_sla$result_top_taxa[[1]])),
                          taxa = unlist(Gnirs_efmw10_sla$result_top_taxa))
Gnirs_tops10_efmw <- merge(selected_df, slablues_ef_correct, by ='taxa', all.x = "TRUE")
write.csv(Gnirs_tops10_efmw,"./output/selected/Gnirs_tops10_efmw.csv", row.names=F)

#selecting top20%indv and their means for sla mw from Gnirs_mwef
Gnirs_mwef10_sla<- read.csv("./output/Gnirs_mwef_sla/sla_Gnirs_10_mwef.csv")
Gnirs_mwef10_sla<- top_selected(data= Gnirs_mwef10_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gnirs_mwef10_sla$result_top_taxa), each = length(Gnirs_mwef10_sla$result_top_taxa[[1]])),
                          taxa = unlist(Gnirs_mwef10_sla$result_top_taxa))
Gnirs_tops10_mwef <- merge(selected_df, slablues_mw_correct, by ='taxa', all.x = "TRUE")
write.csv(Gnirs_tops10_mwef,"./output/selected/Gnirs_tops10_mwef.csv", row.names=F)

#selecting top20%indv and their means for sla ef from Gnirs_efmw
Gnirs_efmw25_sla<- read.csv("./output/Gnirs_efmw_sla/sla_Gnirs_25_efmw.csv")
Gnirs_efmw25_sla<- top_selected(data= Gnirs_efmw25_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gnirs_efmw25_sla$result_top_taxa), each = length(Gnirs_efmw25_sla$result_top_taxa[[1]])),
                          taxa = unlist(Gnirs_efmw25_sla$result_top_taxa))
Gnirs_tops25_efmw <- merge(selected_df, slablues_ef_correct, by ='taxa', all.x = "TRUE")
write.csv(Gnirs_tops25_efmw,"./output/selected/Gnirs_tops25_efmw.csv", row.names=F)

#selecting top20%indv and their means for sla mw from Gnirs_mwef
Gnirs_mwef25_sla<- read.csv("./output/Gnirs_mwef_sla/sla_Gnirs_25_mwef.csv")
Gnirs_mwef25_sla<- top_selected(data= Gnirs_mwef25_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gnirs_mwef25_sla$result_top_taxa), each = length(Gnirs_mwef25_sla$result_top_taxa[[1]])),
                          taxa = unlist(Gnirs_mwef25_sla$result_top_taxa))
Gnirs_tops25_mwef <- merge(selected_df, slablues_mw_correct, by ='taxa', all.x = "TRUE")
write.csv(Gnirs_tops25_mwef,"./output/selected/Gnirs_tops25_mwef.csv", row.names=F)

#selecting top20%indv and their means for sla ef from Gnirs_efmw
Gnirs_efmw50_sla<- read.csv("./output/Gnirs_efmw_sla/sla_Gnirs_50_efmw.csv")
Gnirs_efmw50_sla<- top_selected(data= Gnirs_efmw50_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gnirs_efmw50_sla$result_top_taxa), each = length(Gnirs_efmw50_sla$result_top_taxa[[1]])),
                          taxa = unlist(Gnirs_efmw50_sla$result_top_taxa))
Gnirs_tops50_efmw <- merge(selected_df, slablues_ef_correct, by ='taxa', all.x = "TRUE")
write.csv(Gnirs_tops50_efmw,"./output/selected/Gnirs_tops50_efmw.csv", row.names=F)

#selecting top20%indv and their means for sla mw from Gnirs_mwef
Gnirs_mwef50_sla<- read.csv("./output/Gnirs_mwef_sla/sla_Gnirs_50_mwef.csv")
Gnirs_mwef50_sla<- top_selected(data= Gnirs_mwef50_sla)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gnirs_mwef50_sla$result_top_taxa), each = length(Gnirs_mwef50_sla$result_top_taxa[[1]])),
                          taxa = unlist(Gnirs_mwef50_sla$result_top_taxa))
Gnirs_tops50_mwef <- merge(selected_df, slablues_mw_correct, by ='taxa', all.x = "TRUE")
write.csv(Gnirs_tops50_mwef,"./output/selected/Gnirs_tops50_mwef.csv", row.names=F)

#selecting top20%indv and their means for narea ef from Gnirs_efmw
Gnirs_efmw10_narea<- read.csv("./output/Gnirs_efmw_narea/narea_Gnirs_10_efmw.csv")
Gnirs_efmw10_narea<- top_selected(data= Gnirs_efmw10_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gnirs_efmw10_narea$result_top_taxa), each = length(Gnirs_efmw10_narea$result_top_taxa[[1]])),
                          taxa = unlist(Gnirs_efmw10_narea$result_top_taxa))
Gnirs_topn10_efmw <- merge(selected_df, nareablues_ef_correct, by ='taxa', all.x = "TRUE")
write.csv(Gnirs_topn10_efmw,"./output/selected/Gnirs_topn10_efmw.csv", row.names=F)

#selecting top20%indv and their means for narea mw from Gnirs_mwef
Gnirs_mwef10_narea<- read.csv("./output/Gnirs_mwef_narea/narea_Gnirs_10_mwef.csv")
Gnirs_mwef10_narea<- top_selected(data= Gnirs_mwef10_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gnirs_mwef10_narea$result_top_taxa), each = length(Gnirs_mwef10_narea$result_top_taxa[[1]])),
                          taxa = unlist(Gnirs_mwef10_narea$result_top_taxa))
Gnirs_topn10_mwef <- merge(selected_df, nareablues_mw_correct, by ='taxa', all.x = "TRUE")
write.csv(Gnirs_topn10_mwef,"./output/selected/Gnirs_topn10_mwef.csv", row.names=F)

#selecting top20%indv and their means for narea mw from Gnirs_mwef
Gnirs_mwef10_narea<- read.csv("./output/Gnirs_mwef_narea/narea_Gnirs_10_mwef.csv")
Gnirs_mwef10_narea<- top_selected(data= Gnirs_mwef10_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gnirs_mwef10_narea$result_top_taxa), each = length(Gnirs_mwef10_narea$result_top_taxa[[1]])),
                          taxa = unlist(Gnirs_mwef10_narea$result_top_taxa))
Gnirs_topn10_mwef <- merge(selected_df, nareablues_mw_correct, by ='taxa', all.x = "TRUE")
write.csv(Gnirs_topn10_mwef,"./output/selected/Gnirs_topn10_mwef.csv", row.names=F)

#selecting top20%indv and their means for narea ef from Gnirs_efmw
Gnirs_efmw25_narea<- read.csv("./output/Gnirs_efmw_narea/narea_Gnirs_25_efmw.csv")
Gnirs_efmw25_narea<- top_selected(data= Gnirs_efmw25_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gnirs_efmw25_narea$result_top_taxa), each = length(Gnirs_efmw25_narea$result_top_taxa[[1]])),
                          taxa = unlist(Gnirs_efmw25_narea$result_top_taxa))
Gnirs_topn25_efmw <- merge(selected_df, nareablues_ef_correct, by ='taxa', all.x = "TRUE")
write.csv(Gnirs_topn25_efmw,"./output/selected/Gnirs_topn25_efmw.csv", row.names=F)

#selecting top 165 indv and their means for narea mw from Gnirs_mwef
Gnirs_mwef25_narea<- read.csv("./output/Gnirs_mwef_narea/narea_Gnirs_25_mwef.csv")
Gnirs_mwef25_narea<- top_selected(data= Gnirs_mwef25_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gnirs_mwef25_narea$result_top_taxa), each = length(Gnirs_mwef25_narea$result_top_taxa[[1]])),
                          taxa = unlist(Gnirs_mwef25_narea$result_top_taxa))
Gnirs_topn25_mwef <- merge(selected_df,nareablues_mw_correct, by ='taxa', all.x = "TRUE")
write.csv(Gnirs_topn25_mwef,"./output/selected/Gnirs_topn25_mwef.csv", row.names=F)

#
#selecting top20%indv and their means for narea ef from Gnirs_efmw
Gnirs_efmw50_narea<- read.csv("./output/Gnirs_efmw_narea/narea_Gnirs_50_efmw.csv")
Gnirs_efmw50_narea<- top_selected(data= Gnirs_efmw50_narea)
#extracting mean and selected individual in a dataframe per replication
selected_df <- data.frame(rep = rep(seq_along(Gnirs_efmw50_narea$result_top_taxa), each = length(Gnirs_efmw50_narea$result_top_taxa[[1]])),
                          taxa = unlist(Gnirs_efmw50_narea$result_top_taxa))
Gnirs_topn50_efmw <- merge(selected_df, nareablues_ef_correct, by ='taxa', all.x = "TRUE")
write.csv(Gnirs_topn50_efmw,"./output/selected/Gnirs_topn50_efmw.csv", row.names=F)

#selecting top20%indv and their means for narea mw from Gnirs_mwef
Gnirs_mwef50_narea<- read.csv("./output/Gnirs_mwef_narea/narea_Gnirs_50_mwef.csv")
Gnirs_mwef50_narea<- top_selected(data= Gnirs_mwef50_narea)
#extracting mean and selected individual in a data frame per replication
selected_df <- data.frame(rep = rep(seq_along(Gnirs_mwef50_narea$result_top_taxa), each = length(Gnirs_mwef50_narea$result_top_taxa[[1]])),
                          taxa = unlist(Gnirs_mwef50_narea$result_top_taxa))
Gnirs_topn50_mwef <- merge(selected_df, nareablues_mw_correct, by ='taxa', all.x = "TRUE")
write.csv(Gnirs_topn50_mwef,"./output/selected/Gnirs_topn50_mwef.csv", row.names=F)
#selecting top 165 indv from GBLUP reduced 
GBLUP_redefmw10_narea<- read.csv("./output/GBLUP_reduced/narea_efmw10.csv") %>% na.omit()
GBLUP_redefmw10_narea<- top_selected(data= GBLUP_redefmw10_narea)
#extracting mean and selected individual in a data frame per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_redefmw10_narea$result_top_taxa), each = length(GBLUP_redefmw10_narea$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_redefmw10_narea$result_top_taxa))
GBLUPred_topn10_efmw <- merge(selected_df, nareablues_ef_correct, by ='taxa', all.x = "TRUE")
write.csv(GBLUPred_topn10_efmw,"./output/selected/GBLUPred_topn10_efmw.csv", row.names=F)
#selecting top 165 indv from GBLUP reduced for 25 and 50%
GBLUP_redefmw25_narea<- read.csv("./output/GBLUP_reduced/narea_efmw25.csv") %>% na.omit()
GBLUP_redefmw25_narea<- top_selected(data= GBLUP_redefmw25_narea)
#extracting mean and selected individual in a data frame per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_redefmw25_narea$result_top_taxa), each = length(GBLUP_redefmw25_narea$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_redefmw25_narea$result_top_taxa))
GBLUPred_topn25_efmw <- merge(selected_df, nareablues_ef_correct, by ='taxa', all.x = "TRUE")
write.csv(GBLUPred_topn25_efmw,"./output/selected/GBLUPred_topn25_efmw.csv", row.names=F)
#selecting top 165 indv from GBLUP reduced for ef and 50%
GBLUP_redefmw50_narea<- read.csv("./output/GBLUP_reduced/narea_efmw50.csv") %>% na.omit()
GBLUP_redefmw50_narea<- top_selected(data= GBLUP_redefmw50_narea)
#extracting mean and selected individual in a data frame per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_redefmw50_narea$result_top_taxa), each = length(GBLUP_redefmw50_narea$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_redefmw50_narea$result_top_taxa))
GBLUPred_topn50_efmw <- merge(selected_df, nareablues_ef_correct, by ='taxa', all.x = "TRUE")

write.csv(GBLUPred_topn50_efmw,"./output/selected/GBLUPred_topn50_efmw.csv", row.names=F)

#selecting top 165 indv from GBLUP reduced 
GBLUP_redmwef10_narea<- read.csv("./output/GBLUP_reduced/narea_mwef10.csv") %>% na.omit()
GBLUP_redmwef10_narea<- top_selected(data= GBLUP_redmwef10_narea)
#extracting mean and selected individual in a data frame per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_redmwef10_narea$result_top_taxa), each = length(GBLUP_redmwef10_narea$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_redmwef10_narea$result_top_taxa))
GBLUPred_topn10_mwef <- merge(selected_df, nareablues_mw_correct, by ='taxa', all.x = "TRUE")
write.csv(GBLUPred_topn10_mwef,"./output/selected/GBLUPred_topn10_mwef.csv", row.names=F)
#selecting top 165 indv from GBLUP reduced for 25 and 50%
GBLUP_redmwef25_narea<- read.csv("./output/GBLUP_reduced/narea_mwef25.csv") %>% na.omit()
GBLUP_redmwef25_narea<- top_selected(data= GBLUP_redmwef25_narea)
#extracting mean and selected individual in a data frame per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_redmwef25_narea$result_top_taxa), each = length(GBLUP_redmwef25_narea$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_redmwef25_narea$result_top_taxa))
GBLUPred_topn25_mwef <- merge(selected_df, nareablues_mw_correct, by ='taxa', all.x = "TRUE")
write.csv(GBLUPred_topn25_mwef,"./output/selected/GBLUPred_topn25_mwef.csv", row.names=F)
#selecting top 165 indv from GBLUP reduced for ef and 50%
GBLUP_redmwef50_narea<- read.csv("./output/GBLUP_reduced/narea_mwef50.csv") %>% na.omit()
GBLUP_redmwef50_narea<- top_selected(data= GBLUP_redmwef50_narea)
#extracting mean and selected individual in a data frame per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_redmwef50_narea$result_top_taxa), each = length(GBLUP_redmwef50_narea$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_redmwef50_narea$result_top_taxa))
GBLUPred_topn50_mwef <- merge(selected_df, nareablues_mw_correct, by ='taxa', all.x = "TRUE")
write.csv(GBLUPred_topn50_mwef,"./output/selected/GBLUPred_topn50_mwef.csv", row.names=F)

#selecting top 165 indv from GBLUP reduced 
GBLUP_redmwef10_sla<- read.csv("./output/GBLUP_reduced/sla_mwef10.csv") %>% na.omit()
GBLUP_redmwef10_sla<- top_selected(data= GBLUP_redmwef10_sla)
#extracting mean and selected individual in a data frame per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_redmwef10_sla$result_top_taxa), each = length(GBLUP_redmwef10_sla$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_redmwef10_sla$result_top_taxa))
GBLUPred_tops10_mwef <- merge(selected_df, slablues_mw_correct, by ='taxa', all.x = "TRUE")
write.csv(GBLUPred_tops10_mwef,"./output/selected/GBLUPred_tops10_mwef.csv", row.names=F)
#selecting top 165 indv from GBLUP reduced for 25 and 50%
GBLUP_redmwef25_sla<- read.csv("./output/GBLUP_reduced/sla_mwef25.csv") %>% na.omit()
GBLUP_redmwef25_sla<- top_selected(data= GBLUP_redmwef25_sla)
#extracting mean and selected individual in a data frame per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_redmwef25_sla$result_top_taxa), each = length(GBLUP_redmwef25_sla$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_redmwef25_sla$result_top_taxa))
GBLUPred_tops25_mwef <- merge(selected_df, slablues_mw_correct, by ='taxa', all.x = "TRUE")
write.csv(GBLUPred_tops25_mwef,"./output/selected/GBLUPred_tops25_mwef.csv", row.names=F)
#selecting top 165 indv from GBLUP reduced for ef and 50%
GBLUP_redmwef50_sla<- read.csv("./output/GBLUP_reduced/sla_mwef50.csv") %>% na.omit()
GBLUP_redmwef50_sla<- top_selected(data= GBLUP_redmwef50_sla)
#extracting mean and selected individual in a data frame per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_redmwef50_sla$result_top_taxa), each = length(GBLUP_redmwef50_sla$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_redmwef50_sla$result_top_taxa))
GBLUPred_tops50_mwef <- merge(selected_df, slablues_mw_correct, by ='taxa', all.x = "TRUE")
write.csv(GBLUPred_tops50_mwef,"./output/selected/GBLUPred_tops50_mwef.csv", row.names=F)


#selecting top 165 indv from GBLUP reduced 
GBLUP_redefmw10_sla<- read.csv("./output/GBLUP_reduced/sla_efmw10.csv") %>% na.omit()
GBLUP_redefmw10_sla<- top_selected(data= GBLUP_redefmw10_sla)
#extracting mean and selected individual in a data frame per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_redefmw10_sla$result_top_taxa), each = length(GBLUP_redefmw10_sla$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_redefmw10_sla$result_top_taxa))
GBLUPred_tops10_efmw <- merge(selected_df, slablues_ef_correct, by ='taxa', all.x = "TRUE")
write.csv(GBLUPred_tops10_efmw,"./output/selected/GBLUPred_tops10_efmw.csv", row.names=F)
#selecting top 165 indv from GBLUP reduced for 25 and 50%
GBLUP_redefmw25_sla<- read.csv("./output/GBLUP_reduced/sla_efmw25.csv") %>% na.omit()
GBLUP_redefmw25_sla<- top_selected(data= GBLUP_redefmw25_sla)
#extracting mean and selected individual in a data frame per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_redefmw25_sla$result_top_taxa), each = length(GBLUP_redefmw25_sla$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_redefmw25_sla$result_top_taxa))
GBLUPred_tops25_efmw <- merge(selected_df,slablues_ef_correct, by ='taxa', all.x = "TRUE")
write.csv(GBLUPred_tops25_efmw,"./output/selected/GBLUPred_tops25_efmw.csv", row.names=F)
#selecting top 165 indv from GBLUP reduced for ef and 50%
GBLUP_redefmw50_sla<- read.csv("./output/GBLUP_reduced/sla_efmw50.csv") %>% na.omit()
GBLUP_redefmw50_sla<- top_selected(data= GBLUP_redefmw50_sla)
#extracting mean and selected individual in a data frame per replication
selected_df <- data.frame(rep = rep(seq_along(GBLUP_redefmw50_sla$result_top_taxa), each = length(GBLUP_redefmw50_sla$result_top_taxa[[1]])),
                          taxa = unlist(GBLUP_redefmw50_sla$result_top_taxa))
GBLUPred_tops50_efmw <- merge(selected_df, slablues_ef_correct, by ='taxa', all.x = "TRUE")
write.csv(GBLUPred_tops50_efmw,"./output/selected/GBLUPred_tops50_efmw.csv", row.names=F)


# reading data for narea mw for all model and blues
top_nareablues_mw <- read.csv("./output/selected/top_nareablues_mw.csv")
GBLUP_topn10_mwef <- read.csv("./output/selected/GBLUP_topn_mwef.csv") %>% na.omit()
Gh2_topn10_mwef <- read.csv("./output/selected/Gh2_topn10_mwef.csv")%>% na.omit()
GWW_topn10_mwef <- read.csv("./output/selected/GWW_topn10_mwef.csv")%>% na.omit()
Gnirs_topn10_mwef <- read.csv("./output/selected/Gnirs_topn10_mwef.csv")%>% na.omit()
GBLUPred_topn10_mwef <- read.csv("./output/selected/GBLUPred_topn10_mwef.csv")%>% na.omit()
GBLUP_topn25_mwef <- read.csv("./output/selected/GBLUP_topn_mwef.csv") %>% na.omit()
Gh2_topn25_mwef <- read.csv("./output/selected/Gh2_topn25_mwef.csv")%>% na.omit()
GWW_topn25_mwef <- read.csv("./output/selected/GWW_topn25_mwef.csv")%>% na.omit()
Gnirs_topn25_mwef <- read.csv("./output/selected/Gnirs_topn25_mwef.csv")%>% na.omit()
GBLUPred_topn25_mwef <- read.csv("./output/selected/GBLUPred_topn25_mwef.csv")%>% na.omit()
Gh2_topn50_mwef <- read.csv("./output/selected/Gh2_topn50_mwef.csv")%>% na.omit()
GWW_topn50_mwef <- read.csv("./output/selected/GWW_topn50_mwef.csv")%>% na.omit()
Gnirs_topn50_mwef <- read.csv("./output/selected/Gnirs_topn50_mwef.csv")%>% na.omit()
GBLUPred_topn50_mwef <- read.csv("./output/selected/GBLUPred_topn50_mwef.csv")%>% na.omit()
GBLUP_topn50_mwef <- read.csv("./output/selected/GBLUP_topn_mwef.csv") %>% na.omit()

# Combining the dataset for different schemes
mean_mw_narea <- bind_rows(
  mutate(Gh2_topn10_mwef, Model = "Gh2", Scheme = "10%"),
  mutate(Gnirs_topn10_mwef, Model = "Gnirs", Scheme = "10%"),
  mutate(GWW_topn10_mwef, Model = "WW", Scheme = "10%"),
  mutate(GBLUP_topn10_mwef, Model = "GBLUP", Scheme = "10%"),
  mutate(GBLUPred_topn10_mwef, Model = "GBLUP_red", Scheme = "10%"),
  mutate(Gh2_topn25_mwef, Model = "Gh2", Scheme = "25%"),
  mutate(Gnirs_topn25_mwef, Model = "Gnirs", Scheme = "25%"),
  mutate(GWW_topn25_mwef, Model = "WW", Scheme = "25%"),
  mutate(GBLUP_topn25_mwef, Model = "GBLUP", Scheme = "25%"),
  mutate(GBLUPred_topn25_mwef, Model = "GBLUP_red", Scheme = "25%"),
  mutate(Gh2_topn50_mwef, Model = "Gh2", Scheme = "50%"),
  mutate(GBLUP_topn50_mwef, Model = "GBLUP", Scheme = "50%"),
  mutate(Gnirs_topn50_mwef, Model = "Gnirs", Scheme = "50%"),
  mutate(GWW_topn50_mwef, Model = "WW", Scheme = "50%"),
  mutate(GBLUPred_topn50_mwef, Model = "GBLUP_red", Scheme = "50%")
)
blues_values <- c(1.56841, 1.56841,1.56841)
model_colors_cb <- brewer.pal(5, "Set1")

# Plotting the combined data with facets for different schemes and colored boxplots 
ggplot(mean_mw_narea, aes(x = Model, y = narea, fill = Model)) +
  geom_boxplot() +
  geom_hline(yintercept = blues_values, linetype = "dashed", color = "blue") +
  labs(x = "Model", y = "mean value") +
  ggtitle("Model Performance; BLUES vs model schemes for Narea in loc mw") +
  theme_minimal() +
  scale_fill_manual(values = model_colors_cb) +  
  facet_wrap(~Scheme, scales = "free") 


# reading data for narea mw for all model and blues
top_nareablues_ef <- read.csv("./output/selected/top_nareablues_ef.csv")
GBLUP_topn10_efmw <- read.csv("./output/selected/GBLUP_topn_efmw.csv") %>% na.omit()
Gh2_topn10_efmw <- read.csv("./output/selected/Gh2_topn10_efmw.csv")%>% na.omit()
GWW_topn10_efmw <- read.csv("./output/selected/GWW_topn10_efmw.csv")%>% na.omit()
Gnirs_topn10_efmw <- read.csv("./output/selected/Gnirs_topn10_efmw.csv")%>% na.omit()
GBLUPred_topn10_efmw <- read.csv("./output/selected/GBLUPred_topn10_efmw.csv")%>% na.omit()
GBLUP_topn25_efmw <- read.csv("./output/selected/GBLUP_topn_efmw.csv") %>% na.omit()
Gh2_topn25_efmw <- read.csv("./output/selected/Gh2_topn25_efmw.csv")%>% na.omit()
GWW_topn25_efmw <- read.csv("./output/selected/GWW_topn25_efmw.csv")%>% na.omit()
Gnirs_topn25_efmw <- read.csv("./output/selected/Gnirs_topn25_efmw.csv")%>% na.omit()
GBLUPred_topn25_efmw <- read.csv("./output/selected/GBLUPred_topn25_efmw.csv")%>% na.omit()
Gh2_topn50_efmw <- read.csv("./output/selected/Gh2_topn50_efmw.csv")%>% na.omit()
GWW_topn50_efmw <- read.csv("./output/selected/GWW_topn50_efmw.csv")%>% na.omit()
Gnirs_topn50_efmw <- read.csv("./output/selected/Gnirs_topn50_efmw.csv")%>% na.omit()
GBLUPred_topn50_efmw <- read.csv("./output/selected/GBLUPred_topn50_efmw.csv")%>% na.omit()
GBLUP_topn50_efmw <- read.csv("./output/selected/GBLUP_topn_efmw.csv") %>% na.omit()


# Combining the dataset for different schemes
mean_ef_narea <- bind_rows(
  mutate(Gh2_topn10_efmw, Model = "Gh2", Scheme = "10%"),
  mutate(Gnirs_topn10_efmw, Model = "Gnirs", Scheme = "10%"),
  mutate(GWW_topn10_efmw, Model = "WW", Scheme = "10%"),
  mutate(GBLUP_topn10_efmw, Model = "GBLUP", Scheme = "10%"),
  mutate(GBLUPred_topn10_efmw, Model = "GBLUP_red", Scheme = "10%"),
  mutate(Gh2_topn25_efmw, Model = "Gh2", Scheme = "25%"),
  mutate(Gnirs_topn25_efmw, Model = "Gnirs", Scheme = "25%"),
  mutate(GWW_topn25_efmw, Model = "WW", Scheme = "25%"),
  mutate(GBLUP_topn25_efmw, Model = "GBLUP", Scheme = "25%"),
  mutate(GBLUPred_topn25_efmw, Model = "GBLUP_red", Scheme = "25%"),
  mutate(Gh2_topn50_efmw, Model = "Gh2", Scheme = "50%"),
  mutate(GBLUP_topn50_efmw, Model = "GBLUP", Scheme = "50%"),
  mutate(Gnirs_topn50_efmw, Model = "Gnirs", Scheme = "50%"),
  mutate(GWW_topn50_efmw, Model = "WW", Scheme = "50%"),
  mutate(GBLUPred_topn50_efmw, Model = "GBLUP_red", Scheme = "50%")
)
blues_values <- c(1.5611, 1.5611,1.5611)
model_colors_cb <- brewer.pal(5, "Set1")

# Plotting the combined data with facets for different schemes and colored boxplots 
ggplot(mean_ef_narea, aes(x = Model, y = narea, fill = Model)) +
  geom_boxplot() +
  geom_hline(yintercept = blues_values, linetype = "dashed", color = "blue") +
  labs(x = "Model", y = "mean value") +
  ggtitle("Model Performance; BLUES vs model schemes for Narea in loc ef") +
  theme_minimal() +
  scale_fill_manual(values = model_colors_cb) +  
  facet_wrap(~Scheme,scales = "free")


# reading the data for different schemes for sla ef
top_slablues_ef <- read.csv("./output/selected/top_slablues_ef.csv")
GBLUP_tops10_efmw <- read.csv("./output/selected/GBLUP_tops_efmw.csv") %>% na.omit()
Gh2_tops10_efmw <- read.csv("./output/selected/Gh2_tops10_efmw.csv")%>% na.omit()
GWW_tops10_efmw <- read.csv("./output/selected/GWW_tops10_efmw.csv")%>% na.omit()
Gnirs_tops10_efmw <- read.csv("./output/selected/Gnirs_tops10_efmw.csv")%>% na.omit()
GBLUPred_tops10_efmw <- read.csv("./output/selected/GBLUPred_tops10_efmw.csv")%>% na.omit()
GBLUP_tops25_efmw <- read.csv("./output/selected/GBLUP_tops_efmw.csv") %>% na.omit()
Gh2_tops25_efmw <- read.csv("./output/selected/Gh2_tops25_efmw.csv")%>% na.omit()
GWW_tops25_efmw <- read.csv("./output/selected/GWW_tops25_efmw.csv")%>% na.omit()
Gnirs_tops25_efmw <- read.csv("./output/selected/Gnirs_tops25_efmw.csv")%>% na.omit()
GBLUPred_tops25_efmw <- read.csv("./output/selected/GBLUPred_tops25_efmw.csv")%>% na.omit()
Gh2_tops50_efmw <- read.csv("./output/selected/Gh2_tops50_efmw.csv")%>% na.omit()
GWW_tops50_efmw <- read.csv("./output/selected/GWW_tops50_efmw.csv")%>% na.omit()
Gnirs_tops50_efmw <- read.csv("./output/selected/Gnirs_tops50_efmw.csv")%>% na.omit()
GBLUPred_tops50_efmw <- read.csv("./output/selected/GBLUPred_tops50_efmw.csv")%>% na.omit()
GBLUP_tops50_efmw <- read.csv("./output/selected/GBLUP_tops_efmw.csv") %>% na.omit()


# Combining the dataset for different schemes
# Combining the dataset for different schemes
mean_ef_sla <- bind_rows(
  mutate(Gh2_tops10_efmw, Model = "Gh2", Scheme = "10%"),
  mutate(Gnirs_tops10_efmw, Model = "Gnirs", Scheme = "10%"),
  mutate(GWW_tops10_efmw, Model = "WW", Scheme = "10%"),
  mutate(GBLUP_tops10_efmw, Model = "GBLUP", Scheme = "10%"),
  mutate(GBLUPred_tops10_efmw, Model = "GBLUP_red", Scheme = "10%"),
  mutate(Gh2_tops25_efmw, Model = "Gh2", Scheme = "25%"),
  mutate(Gnirs_tops25_efmw, Model = "Gnirs", Scheme = "25%"),
  mutate(GWW_tops25_efmw, Model = "WW", Scheme = "25%"),
  mutate(GBLUP_tops25_efmw, Model = "GBLUP", Scheme = "25%"),
  mutate(GBLUPred_tops25_efmw, Model = "GBLUP_red", Scheme = "25%"),
  mutate(Gh2_tops50_efmw, Model = "Gh2", Scheme = "50%"),
  mutate(GBLUP_tops50_efmw, Model = "GBLUP", Scheme = "50%"),
  mutate(Gnirs_tops50_efmw, Model = "Gnirs", Scheme = "50%"),
  mutate(GWW_tops50_efmw, Model = "WW", Scheme = "50%"),
  mutate(GBLUPred_tops50_efmw, Model = "GBLUP_red", Scheme = "50%")
)
blues_values <- c(309.008,309.008,309.008 )
model_colors_cb <- brewer.pal(5, "Set1")

# Plotting the combined data with facets for different schemes and colored boxplots 
ggplot(mean_ef_sla, aes(x = Model, y = sla, fill = Model)) +
  geom_boxplot() +
  geom_hline(yintercept = blues_values, linetype = "dashed", color = "blue") +
  labs(x = "Model", y = "Mean Value") +
  ggtitle("Model Performance; BLUES vs model schemes for sla in loc ef") +
  theme_minimal() +
  scale_fill_manual(values = model_colors_cb) +  
  facet_wrap(~Scheme, scales = "free")

# reading the data for different schemes for sla mw
top_slablues_mw <- read.csv("./output/selected/top_slablues_mw.csv")
GBLUP_tops10_mwef <- read.csv("./output/selected/GBLUP_tops_mwef.csv") %>% na.omit()
Gh2_tops10_mwef <- read.csv("./output/selected/Gh2_tops10_mwef.csv")%>% na.omit()
GWW_tops10_mwef <- read.csv("./output/selected/GWW_tops10_mwef.csv")%>% na.omit()
Gnirs_tops10_mwef <- read.csv("./output/selected/Gnirs_tops10_mwef.csv")%>% na.omit()
GBLUPred_tops10_mwef <- read.csv("./output/selected/GBLUPred_tops10_mwef.csv")%>% na.omit()
GBLUP_tops25_mwef <- read.csv("./output/selected/GBLUP_tops_mwef.csv") %>% na.omit()
Gh2_tops25_mwef <- read.csv("./output/selected/Gh2_tops25_mwef.csv")%>% na.omit()
GWW_tops25_mwef <- read.csv("./output/selected/GWW_tops25_mwef.csv")%>% na.omit()
Gnirs_tops25_mwef <- read.csv("./output/selected/Gnirs_tops25_mwef.csv")%>% na.omit()
GBLUPred_tops25_mwef <- read.csv("./output/selected/GBLUPred_tops25_mwef.csv")%>% na.omit()
Gh2_tops50_mwef <- read.csv("./output/selected/Gh2_tops50_mwef.csv")%>% na.omit()
GWW_tops50_mwef <- read.csv("./output/selected/GWW_tops50_mwef.csv")%>% na.omit()
Gnirs_tops50_mwef <- read.csv("./output/selected/Gnirs_tops50_mwef.csv")%>% na.omit()
GBLUPred_tops50_mwef <- read.csv("./output/selected/GBLUPred_tops50_mwef.csv")%>% na.omit()
GBLUP_tops50_mwef <- read.csv("./output/selected/GBLUP_tops_mwef.csv") %>% na.omit()

# Combining dataset for different schemes for sla mw
mean_mw_sla <- bind_rows(
  mutate(Gh2_tops10_mwef, Model = "Gh2", Scheme = "10%"),
  mutate(Gnirs_tops10_mwef, Model = "Gnirs", Scheme = "10%"),
  mutate(GWW_tops10_mwef, Model = "WW", Scheme = "10%"),
  mutate(GBLUP_tops10_mwef, Model = "GBLUP", Scheme = "10%"),
  mutate(GBLUPred_tops10_mwef, Model = "GBLUP_red", Scheme = "10%"),
  mutate(Gh2_tops25_mwef, Model = "Gh2", Scheme = "25%"),
  mutate(Gnirs_tops25_mwef, Model = "Gnirs", Scheme = "25%"),
  mutate(GWW_tops25_mwef, Model = "WW", Scheme = "25%"),
  mutate(GBLUP_tops25_mwef, Model = "GBLUP", Scheme = "25%"),
  mutate(GBLUPred_tops25_mwef, Model = "GBLUP_red", Scheme = "25%"),
  mutate(Gh2_tops50_mwef, Model = "Gh2", Scheme = "50%"),
  mutate(GBLUP_tops50_mwef, Model = "GBLUP", Scheme = "50%"),
  mutate(Gnirs_tops50_mwef, Model = "Gnirs", Scheme = "50%"),
  mutate(GWW_tops50_mwef, Model = "WW", Scheme = "50%"),
  mutate(GBLUPred_tops50_mwef, Model = "GBLUP_red", Scheme = "50%")
)
blues_values <- c(312.0316,312.0316,312.0316 )
model_colors_cb <- brewer.pal(5, "Set1")

# Plotting the combined data with facets for different schemes and colored boxplots 
ggplot(mean_mw_sla, aes(x = Model, y = sla, fill = Model)) +
  geom_boxplot() +
  geom_hline(yintercept = blues_values, linetype = "dashed", color = "blue") +
  labs(x = "Model", y = "Mean Value") +
  ggtitle("Model Performance; BLUES vs model schemes for sla in loc mw") +
  theme_minimal() +
  scale_fill_manual(values = model_colors_cb) +  
  facet_wrap(~Scheme, scales = "free")



