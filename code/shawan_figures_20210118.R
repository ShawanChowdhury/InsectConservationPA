# Filename: shawan_figures.R
# R code to visualise conservation coverage for insect families

# Mixing phylogeny with circular plots

# Load libraries ----------------------------------------------------------
library(treeio)
library(ape)
library(ggtree)
library(tidyverse)

# Shawan's coverage data set ----------------------------------------------------------------
AH <- read.csv("data/AH_interpolation.csv", header = TRUE) %>% tibble()
# AH <- read.csv("data/interpolation_AH.csv", header = TRUE) %>% tibble()

## Calculate the mean coverage for each insect family 
AH_df <- AH %>% 
  dplyr::mutate(Family = str_extract(Family, regex('[A-Za-z]+'))) %>% 
  arrange(Family) %>% 
  group_by(Order, Family) %>% 
  dplyr::mutate(pos = sum(gap > 0),
                neg = sum(gap < 0)) %>% 
  dplyr::summarise(positive = mean(pos),
                   negative = mean(neg)) %>% 
  dplyr::mutate(pa_coverage = ifelse(negative < 1, 100, 
                              ifelse(negative >= 0, (negative / (positive + negative)) * 100, "NA"))) 

### Circular plot
source("code/circular_plotter_function.R")

AH_df %>% 
  group_by(Order) %>% 
  summarise(count = n()) %>% 
  arrange(Order) %>% 
  tail(20)


set_1 <- AH_df %>% filter(Order %in% c("Archaeognatha", "Blattodea", "Coleoptera", "Dernaotera"))
set_2 <- AH_df %>% filter(Order %in% c("Diptera", "Embioptera", "Ephemeroptera", "Grylloblattodea"))
set_3 <- AH_df %>% filter(Order %in% c("Hemiptera", "Hymenoptera", "Isoptera"))
set_4 <- AH_df %>% filter(Order %in% c("Isoptera", "Lepidoptera", "Mantodea", "MantoPhasmatodea", "Mecoptera")) 
set_5 <- AH_df %>% filter(Order %in% c("Megaloptera", "Neuroptera", "Odonata", "Orthoptera", "Phasmotodea", "Trichoptera", "Zygentoma")) 

full_p <- circular_plotter(.df = AH_df)

p1 <- circular_plotter(.df = set_1)
p2 <- circular_plotter(.df = set_2)
p3 <- circular_plotter(.df = set_3)
p4 <- circular_plotter(.df = set_4)
p5 <- circular_plotter(.df = set_5)


ggsave(filename = 'output/AH_circular.pdf', plot = full_p)



library(gridExtra)
print_plot <- list(p1, p2, p3, p4, p5)

pdf("shawan_plot.pdf", onefile = TRUE)

for (i in seq_along(print_plot)) {
  print(print_plot[[i]])
}

dev.off()





# time tree of insect order -----------------------------------------------
# insect family tree from timetree.org. Download newick tree file 
# read newick file with ggtree:read.tree
order_tree <- treeio::read.tree("data/insecta_order.nwk")

# convert tree into tibble so we can join other data later
order_tree_tbl <- order_tree %>% treeio::as_tibble() 

order_tree_tbl_ah <- dplyr::left_join(order_tree_tbl, AH_df, by = c('label' = 'Order')) 

## Dropping tips CH
taxa_to_drop_CH <- order_tree_tbl_ah %>% 
  dplyr::filter(is.na(pa_coverage), # filters those that do not have mean percent_gap_CH data
                stringr::str_detect(label, pattern = regex("a"))) %>% 
  dplyr::select(label)

tip_drop_CH <- taxa_to_drop_CH$label

## Prune tree to only ones we have name for 
sub_tree_CH <- order_tree %>%
  ape::drop.tip(., tip_drop_CH) # keep named tips

## Convert to dataframe for sub_tree_CH 
sub_tree_CH_df <- sub_tree_CH %>% treeio::as_tibble() 


plot(sub_tree_CH, type = 'fan')



### CONTINUE HERE ####

