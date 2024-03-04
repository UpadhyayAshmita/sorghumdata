#loading library-----------------
library(tidyverse)
library(asreml)
source('./script/outlier.R')
library(data.table)
#reading field design and wavelength data and filtering year 16 ("wavelength data available for year 16)

design<- read.csv('./data/design.csv') %>% filter(year!= 16)
averagedspectra <- read.csv('./data/Averaged_Spectra.csv')

averagedspectra <- averagedspectra %>%
  separate_wider_delim(
    Spectra,
    delim = "_",
    names = c("Spectra", "plotnum"),
    cols_remove = T
  ) %>%
  select(-plotnum) %>%
  group_by(Spectra) %>%
  summarise_all(mean) %>%
  ungroup() %>%
  rename(PlotID = Spectra)

designnew<- left_join(design,averagedspectra, by = "PlotID") # Apply left_join dplyr function

