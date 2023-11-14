library(data.table)
library(cvTools)
library(bestNormalize)
library(fs)
library(ASRgenomics)
library(asreml)
library(dplyr)
library(janitor)
source("./script/aux_functions.R")


#-----reading the blues data file-----
blues<- fread("./traitsoutput/blues.csv")




#reading blues from first stage analysis for narea 
# blues for narea trait

Nitroblues <- subset(blues, select = c("taxa", "narea"))
fwrite(Nitroblues, "./traitsoutput/nitroblues.csv",row.names = FALSE)

nitroblues<- fread("./traitsoutput/nitroblues.csv")
# blues for sla trait
slablues <- subset(blues, select = c("taxa", "sla"))
#fwrite(slablues, "./traitsoutput/slablues.csv", row.names = F)

nitroblues<- fread("./traitsoutput/nitroblues.csv")
nareabluesef<- fread("./traitsoutput/nareabluesef.csv")
nareabluesmw<- fread("./traitsoutput/nareabluesmw.csv")

# ---------------------performing cross validation for joint location for narea trait---------------------


# ---------------------creating list with 5 fold and 20 reps----------------------

set.seed(134)
sort<-list()
nl <- length(unique(nitroblues$taxa))
for(a in 1:20){
  folds <- cvFolds(nl,type ="random", K=5, R = 1)
  Sample<-cbind(folds$which,folds$subsets)
  cv<-split(levels(nitroblues$taxa)[Sample[,2]], f=Sample[,1])#a list of subsets (folds) of taxa
  sort[[a]]<-cv
}
rm(a, folds, cv, Sample)


# ---------------------doing cross validation in narea---------------------

sort= sort
train= nitroblues
validation= nitroblues
GINV= GINV
traitname = "narea"

# ---------------------using the function crossv for crossvalidation---------------------


gbevN<- crossv(sort =sort, train= nitroblues, validation = nitroblues, GINV, traitname = "narea")



# ---------------------crossvalidation for joint location for sla trait---------------------




slablues<- fread("./traitsoutput/slablues.csv")
slablues<- transform(slablues,
                     taxa= factor(taxa))

# ---------------------creating list with 5 fold and 20 reps---------------------

set.seed(135)
sort<-list()
nl <- length(unique(slablues$taxa))
for(a in 1:20){
  folds <- cvFolds(nl,type ="random", K=5, R = 1)
  Sample<-cbind(folds$which,folds$subsets)
  cv<-split(levels(slablues$taxa)[Sample[,2]], f=Sample[,1])#a list of subsets (folds) of taxa
  sort[[a]]<-cv
}
rm(a, folds, cv, Sample)

#----- second stage for trait sla --------

sort = sort
train = slablues
validation = slablues
GINV = GINV
traitname ="sla"

gbevsla<- crossv(sort= sort,train = slablues, validation = slablues,GINV = GINV,traitname = "sla")



# ---------------------fitting model and cross validation using wavebluesef and waveblues mw---------------------

wavebluesef<- fread("./output/wavebluesef.csv")
wavebluesef<- transform(wavebluesef,
                     name2= factor(name2))

wavebluesmw<- fread("./output/wavebluesmw.csv")
wavebluesmw<- transform(wavebluesmw,
                        name2= factor(name2))

nareabluesef<-fread("./traitsoutput/nareabluesef.csv")

# ---------------------cross validation set---------------------

set.seed(135)
sort<-list()
nl <- length(unique(wavebluesef$name2))

for(a in 1:20){
  folds <- cvFolds(nl,type ="random", K=5, R = 1)
  Sample<-cbind(folds$which,folds$subsets)
  cv<-split(levels(wavebluesef$name2)[Sample[,2]], f=Sample[,1])#a list of subsets (folds) of taxa
  sort[[a]]<-cv
}
rm(a, folds, cv, Sample)


# ------------------ obtaining inverse of whole spectra relationship matrix---------------------

rel<- fread('./relmatrices/re_w_ef.csv', data.table= F)
rownames(rel) <- colnames(rel)
Gb <- G.tuneup(G = as.matrix(rel), bend = TRUE, eig.tol = 1e-06)$Gb
GINV <- G.inverse(G = Gb , sparseform = T)


# --------------------fitting model---------------------

sort= sort
train= wavebluesef
validation= wavebluesmw
GINV= GINV


gbev_narea<- crossv(sort= sort,train = wavebluesef, validation = wavebluesmw, GINV = GINV,traitname = "narea")

