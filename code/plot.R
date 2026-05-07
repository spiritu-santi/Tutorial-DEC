library(RevGadgets)
library(coda)
library(ggplot2)
library(ggtree)
library(grid)
library(gridExtra)
library(tidyverse)

file <- "output/mi_corrida.1.bg.ase.tre"

# Create the labels vector.
# This is a named vector where names correspond 
# to the computer-readable numbers generated 
# in the biogeographic analysis and the values 
# are character strings of whatever you'd like 
# as labels on the figure. The state.labels.txt
# file produced in the analysis links the 
# computer-readable numbers with presence/ absence
# data for individual ranges.

labs <- read.table("mi_corrida.1.state_labels.txt",sep=",",header=TRUE) %>% 
  as_tibble() %>% select(-range,) %>% deframe()

# pass the labels vector and file name to the processing script
dec_example <- processAncStates(file, state_labels = labs)

plotAncStatesPie(dec_example, cladogenetic = T, tip_labels_offset = 0.2)



