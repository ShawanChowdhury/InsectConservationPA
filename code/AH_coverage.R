# Filename: Shawan_tree.R
# R code to visualise conservation coverage for insect families

# Visualise CH for insect families on phylogeny

# Load libraries ----------------------------------------------------------
library(treeio)
library(ape)
library(ggtree)
library(tidyverse)

# Shawan's coverage data set ----------------------------------------------------------------
AH <- read.csv("data/interpolation_AH.csv", header = TRUE)
# AH_update <- read.csv("data/AH_interpolation.csv", header = TRUE)

## Calculate the mean coverage for each insect family 
AH_df_pre <- AH %>% 
  dplyr::mutate(Family = str_extract(Family, regex('[A-Za-z]+'))) %>% 
  dplyr::group_by(Family, Order) %>% # group by order and family
  # dplyr::summarise(percent_gap = (mean(Gap))) %>% 
  dplyr::summarise(pos = sum(per.ov.ah - y >0),
                   neg = sum(per.ov.ah - y <0)) %>%
  dplyr::mutate(percent_gap = (neg / (pos+neg) * 100))

percentile_cov <- quantile(AH_df_pre$percent_gap, seq(0, 1, by = 0.01), na.rm = TRUE)
as.numeric(percentile_cov['2%']) # 98th percentile
as.numeric(percentile_cov['10%']) # 90th percentile
as.numeric(percentile_cov['20%']) # 80th percentile
as.numeric(percentile_cov['30%']) # 70th percentile
as.numeric(percentile_cov['50%']) # 50th percentile
as.numeric(percentile_cov['90%']) # 50th percentile

AH_df <- AH_df_pre %>% 
  mutate(percentile = ifelse(percent_gap >= (percentile_cov['0%']) & percent_gap < as.numeric(percentile_cov['2%']), '98th',
                      ifelse(percent_gap >= (percentile_cov['2%']) & percent_gap < as.numeric(percentile_cov['10%']), '90th',
                      ifelse(percent_gap >= as.numeric(percentile_cov['10%']) & percent_gap < as.numeric(percentile_cov['20%']), '80th',
                      ifelse(percent_gap >= as.numeric(percentile_cov['20%']) & percent_gap < as.numeric(percentile_cov['30%']), '70th',
                      ifelse(percent_gap >= as.numeric(percentile_cov['30%']) & percent_gap <= as.numeric(percentile_cov['50%']), '50th', 
                      ifelse(percent_gap > (percentile_cov['50%']) & percent_gap <= as.numeric(percentile_cov['90%']), '<50th', NA)))))))  %>% 
  mutate(percent_cat = ifelse(percent_gap == 0, '<10', 
                       ifelse(percent_gap >= 10 & percent_gap < 60, '10-60',
                       ifelse(percent_gap >= 60 & percent_gap < 70, '60-70', 
                       ifelse(percent_gap >= 70 & percent_gap < 80, '70-80', 
                       ifelse(percent_gap >= 80 & percent_gap < 90, '80-90',               
                       ifelse(percent_gap >= 90 & percent_gap < 95, '90-95',                 
                       ifelse(percent_gap >= 95, '>95', NA))))))))                    


# tree from timetree ------------------------------------------------------
# insect family tree from timetree.org. Download newick tree file 
# read newick file with ggtree:read.tree
full_tree <- treeio::read.tree("data/insecta_family.nwk")

# convert tree into tibble so we can join other data later
tbl_tree <- full_tree %>% treeio::as_tibble() 


# Prune tree data --------------------------------------------------------------
## Figure out what data we have for AH
tbl_tree_AH_df_tmp <- dplyr::left_join(tbl_tree, AH_df, by = c('label' = 'Family')) 

## Dropping tips AH
taxa_to_drop_AH <- tbl_tree_AH_df_tmp %>% 
  dplyr::filter(is.na(percent_gap), # filters those that do not have mean ch.area data
                stringr::str_detect(label, pattern = regex("idae"))) %>% 
  dplyr::select(label)

tip_drop_ah <- taxa_to_drop_AH$label

## Prune tree to only ones we have name for 
sub_tree_ah <- full_tree %>%
  ape::drop.tip(., tip_drop_ah) # keep named tips


## Convert to dataframe for sub_tree_ah 
sub_tree_ah_df <- sub_tree_ah %>% treeio::as_tibble() 

## Join tree with sub tree data
tbl_tree_AH_df <- dplyr::left_join(sub_tree_ah_df, AH_df, by = c('label' = 'Family')) 

# Plotting ----------------------------------------------------------------
## Plotting just data that have data

## Annotate Order on the phylogeny.
# Before doing that, we must get the node 'id' for each order.
node_label_df_AH <- tbl_tree_AH_df %>%
  dplyr::group_by(Order) %>%
  dplyr::summarise(minnode = min(node),
                   maxnode = max(node)) %>%
  ungroup() %>%
  as.data.frame()


# Convert table with phylo + coverage data back into phylo object
y_tree_ch <- ape::as.phylo(tbl_tree_AH_df) # full data set

# Get number of node that represents the insect Order. We use this for annotating different orders later

## MCRA for CH
mcra_ch <- c() # create an empty vector for the for loop. We do this by finding the most recent common ancestors of the family.

for (i in 1:dim(node_label_df_AH)[1]-1) {
  
  mcra_ch[i] <- c(ape::getMRCA(y_tree_ch, tip = c(node_label_df_AH$minnode[i], node_label_df_AH$maxnode[i])))
}

mcra_ch_vector <- c(mcra_ch, NA)

# Append the empty vector as a column in our main dataframe
node_label_df_AH$MRCA_node <- mcra_ch_vector

## Convert back to phylo
y_tree_phylo_ah <- as.treedata(tbl_tree_AH_df)

# PLOTTING ----------------------------------------------------------------

hist(tbl_tree_AH_df$percent_gap)

## AH Protection Gap discrete
colours_percent_cat <- c("<10" = "#0f4601", "10-60" = "#1b7c2b", 
                     "60-70" = "#3cabd3", "70-80" = "lightblue", 
                     "80-90" = "blue", "90-95" = "orange",
                     ">95" = "red")

AH_percentile_discrete <- ggtree(y_tree_phylo_ah, aes(colour = percent_cat), layout = 'circular',
                  ladderize = TRUE, size = 0.3) +
        labs(color = "Protection Gap %") +
        scale_colour_manual(values = colours_percent_cat,
                            breaks = c("<10", "10-60", "60-70", "70-80", "80-90", "90-95", ">95"),
                            labels = c("<10", "10-60", "60-70", "70-80", "80-90", "90-95", ">95"), 
                            na.value = "black") +
  # Insect order labels
  geom_cladelabel(node = node_label_df_AH$MRCA_node[1], label = paste(node_label_df_AH$Order[1]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[2], label = paste(node_label_df_AH$Order[2]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[3], label = paste(node_label_df_AH$Order[3]), barsize = 0.3, colour = 'black', fontsize = 2, hjust =1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[4], label = paste(node_label_df_AH$Order[4]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[5], label = paste(node_label_df_AH$Order[5]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[6], label = paste(node_label_df_AH$Order[6]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[7], label = paste(node_label_df_AH$Order[7]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[8], label = paste(node_label_df_AH$Order[8]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[9], label = paste(node_label_df_AH$Order[9]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[10], label = paste(node_label_df_AH$Order[10]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[11], label = paste(node_label_df_AH$Order[11]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[12], label = paste(node_label_df_AH$Order[12]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[13], label = paste(node_label_df_AH$Order[13]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=-1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[14], label = paste(node_label_df_AH$Order[14]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=-0.75) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[15], label = paste(node_label_df_AH$Order[15]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[16], label = paste(node_label_df_AH$Order[16]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[17], label = paste(node_label_df_AH$Order[17]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[18], label = paste(node_label_df_AH$Order[18]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[19], label = paste(node_label_df_AH$Order[19]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[20], label = paste(node_label_df_AH$Order[20]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[21], label = paste(node_label_df_AH$Order[21]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=-1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[22], label = paste(node_label_df_AH$Order[22]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[23], label = paste(node_label_df_AH$Order[23]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[24], label = paste(node_label_df_AH$Order[24]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[25], label = paste(node_label_df_AH$Order[25]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[26], label = paste(node_label_df_AH$Order[26]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[27], label = paste(node_label_df_AH$Order[27]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[28], label = paste(node_label_df_AH$Order[28]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[29], label = paste(node_label_df_AH$Order[29]), barsize = 0.3, colour = 'black', fontsize = 2)

AH_percentile_discrete

ggsave(filename = 'output/AH_tree_discrete.pdf', plot = AH_percentile_discrete)

## AH Coverage

colours_percent <- c("<50th" = "red", "50th" = "orange", 
                     "70th" = "blue", "80th" = "lightblue", 
                     "90th" = "#28a8b9", "98th" = "darkgreen")

AH_percentile <- ggtree(y_tree_phylo_ah, aes(colour = percentile), layout = 'circular',
                  ladderize = TRUE, size = 0.3) +
        labs(color = "Percentile Coverage") +
        scale_colour_manual(values = colours_percent, na.value = "black") +
  # Insect order labels
  geom_cladelabel(node = node_label_df_AH$MRCA_node[1], label = paste(node_label_df_AH$Order[1]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[2], label = paste(node_label_df_AH$Order[2]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[3], label = paste(node_label_df_AH$Order[3]), barsize = 0.3, colour = 'black', fontsize = 2, hjust =1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[4], label = paste(node_label_df_AH$Order[4]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[5], label = paste(node_label_df_AH$Order[5]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[6], label = paste(node_label_df_AH$Order[6]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[7], label = paste(node_label_df_AH$Order[7]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[8], label = paste(node_label_df_AH$Order[8]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[9], label = paste(node_label_df_AH$Order[9]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[10], label = paste(node_label_df_AH$Order[10]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[11], label = paste(node_label_df_AH$Order[11]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[12], label = paste(node_label_df_AH$Order[12]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[13], label = paste(node_label_df_AH$Order[13]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=-1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[14], label = paste(node_label_df_AH$Order[14]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=-0.75) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[15], label = paste(node_label_df_AH$Order[15]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[16], label = paste(node_label_df_AH$Order[16]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[17], label = paste(node_label_df_AH$Order[17]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[18], label = paste(node_label_df_AH$Order[18]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[19], label = paste(node_label_df_AH$Order[19]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[20], label = paste(node_label_df_AH$Order[20]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[21], label = paste(node_label_df_AH$Order[21]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=-1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[22], label = paste(node_label_df_AH$Order[22]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[23], label = paste(node_label_df_AH$Order[23]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[24], label = paste(node_label_df_AH$Order[24]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[25], label = paste(node_label_df_AH$Order[25]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[26], label = paste(node_label_df_AH$Order[26]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[27], label = paste(node_label_df_AH$Order[27]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[28], label = paste(node_label_df_AH$Order[28]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[29], label = paste(node_label_df_AH$Order[29]), barsize = 0.3, colour = 'black', fontsize = 2)


AH_percentile

ggsave(filename = 'output/AH_tree_percentile.pdf', plot = AH_percentile)

### Gradient colours

library(colorRamps)

pokemonpal <- colorRampPalette(c('red', 'blue', 'green'))
pokemonpal(100)[100:1]

AH_tree_gradient <- ggtree(y_tree_phylo_ah, aes(colour = percent_gap), layout = 'circular',
                  ladderize = TRUE, size = 0.3) +
  scale_colour_gradientn(name = "Protection Gap %", colours = pokemonpal(100)[100:1], na.value = 'grey') +
  # scale_colour_stepsn(name = "% Coverage CH", colours = c("black", "lightblue", "blue", "orange", "darkgreen"), n.breaks = 6, nice.breaks = TRUE) +
  # scale_colour_stepsn(name = "% Coverage CH", colours = c("black", "lightblue", "darkblue", "orange", "#FF0000"), breaks = c(7, 8, 9,10,11)) +
  
  # Insect order labels
  # geom_cladelabel(node = node_label_df_AH$MRCA_node[1], label = paste(node_label_df_AH$Order[1]), barsize = 0.3, colour = 'black', fontsize = 2) +
  # geom_cladelabel(node = node_label_df_AH$MRCA_node[2], label = paste(node_label_df_AH$Order[2]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[3], label = paste(node_label_df_AH$Order[3]), barsize = 0.3, colour = 'black', fontsize = 2, hjust =2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[4], label = paste(node_label_df_AH$Order[4]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[5], label = paste(node_label_df_AH$Order[5]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[6], label = paste(node_label_df_AH$Order[6]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[7], label = paste(node_label_df_AH$Order[7]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[8], label = paste(node_label_df_AH$Order[8]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[9], label = paste(node_label_df_AH$Order[9]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[10], label = paste(node_label_df_AH$Order[10]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[11], label = paste(node_label_df_AH$Order[11]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[12], label = paste(node_label_df_AH$Order[12]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[13], label = paste(node_label_df_AH$Order[13]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=-1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[14], label = paste(node_label_df_AH$Order[14]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=-0.75) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[15], label = paste(node_label_df_AH$Order[15]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[16], label = paste(node_label_df_AH$Order[16]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[17], label = paste(node_label_df_AH$Order[17]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[18], label = paste(node_label_df_AH$Order[18]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[19], label = paste(node_label_df_AH$Order[19]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[20], label = paste(node_label_df_AH$Order[20]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[21], label = paste(node_label_df_AH$Order[21]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=-1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[22], label = paste(node_label_df_AH$Order[22]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[23], label = paste(node_label_df_AH$Order[23]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[24], label = paste(node_label_df_AH$Order[24]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[25], label = paste(node_label_df_AH$Order[25]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[26], label = paste(node_label_df_AH$Order[26]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[27], label = paste(node_label_df_AH$Order[27]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[28], label = paste(node_label_df_AH$Order[28]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[29], label = paste(node_label_df_AH$Order[29]), barsize = 0.3, colour = 'black', fontsize = 2)

AH_tree_gradient

# class(CH_tree)
# ggsave(filename = 'output/CH_tree.pdf', plot = CH_tree)
ggsave(filename = 'output/AH_tree_gradient.pdf', plot = AH_tree)






tbl_tree_AH_df %>% 
  # distinct(label) %>% 
  filter(str_detect(label, pattern = regex("idae"))) %>% 
  glimpse()






# additional tips to add --------------------------------------------------
additional_df <- read.csv("data/additional_family_names.csv", header = TRUE)


## Upon joining, I see there are some families that do not have coverage data
# The following code filters those that do not have mean ch.area data
family_na_data <- tbl_tree_AH_df %>%
  dplyr::filter(is.na(mean.ch.area), # filters those that do not have mean ch.area data
                stringr::str_detect(label, pattern = regex("dae"))) # must be a label (some labels are numbers)

## There are about 150 families that do not have coverage data but this may be because synonyms or slightly different spelling. 
family_na_data 

# I suggest going back through your interpolation data to check if you have misspelled something or are able to fill data for these missing families.


