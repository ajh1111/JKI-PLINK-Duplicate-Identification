# This script integrates data from the Julius Kuhn Institute with data from the Jim Dunckley Orchard and Plant & Food Research, to then find duplicates with PLINK.
# Prior to using this script, JKI data was formatted to match PLINK specifications, but still in need of transposing.

#Load packages
library(tibble)
library(igraph)

library(dplyr)
library(tidyr)

#Set wd
setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/JKI PLINK Duplicate Identification/Inputs")

#Load JKI data and transpose
JKI <- read.delim("JKI_Samples_Fixed.txt", header = FALSE, stringsAsFactors = FALSE)
JKI_t <- as.data.frame(t(JKI))

#Add empty columns for PLINK formatting
JKI_t <- add_column(JKI_t, Fa = 0, Mo = 0, Se = 0, Ph = 0, .before = "V2" )
JKI_t <- add_column(JKI_t, Fid = 0, .before = "V1" )

#Remove header row 
JKI_t = JKI_t[-1, ]

#Remove row and column names
colnames(JKI_t) <- NULL
rownames(JKI_t) <- NULL

#Load in PLINK .ped
ped <- read.csv("JD_PFR_PLINK.ped", header = FALSE,sep = "\t")

#Bind JKI data onto PLINK .ped
names(JKI_t) <- names(ped)
combined_ped <- bind_rows(ped, JKI_t)


