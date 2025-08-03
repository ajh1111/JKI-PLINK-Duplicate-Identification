# This script integrates data from the Julius Kuhn Institute with data from the Jim Dunckley Orchard and Plant & Food Research, to then find duplicates with PLINK.
# Prior to using this script, JKI data was formatted to match PLINK specifications, but still in need of transposing.

#Load packages
library(tibble)
library(igraph)

library(dplyr)
library(tidyr)

#Set wd
setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/JKI PLINK Duplicate Identification/Inputs")

#Load JKI data
JKI <- read.delim("JKI_Samples_Fixed.txt", header = FALSE, stringsAsFactors = FALSE)