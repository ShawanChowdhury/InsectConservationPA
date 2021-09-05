# Filename: Shawan_tree.R
# R code to visualise conservation coverage for insect families

# Visualise CH for insect families on phylogeny

# Load libraries ----------------------------------------------------------
library(treeio)
library(ape)
library(ggtree)
library(tidyverse)

# Shawan's coverage data set ----------------------------------------------------------------
CH <- read.csv("data/interpolation_CH.csv", header = TRUE)

## Calculate the mean coverage for each insect family 
CH_df <- CH %>% 
  dplyr::mutate(Family = str_extract(Family, regex('[A-Za-z]+'))) %>% 
  dplyr::group_by(Family, Order) %>% # group by order and family
  dplyr::summarise(pos = sum(per.ov.ch - y >0),
                neg = sum(per.ov.ch - y <0)) %>% 
  dplyr::mutate(coverage_ch = (pos / (pos+neg) * 100))

# tree from timetree ------------------------------------------------------
# insect family tree from timetree.org. Download newick tree file 
# read newick file with ggtree:read.tree
full_tree <- treeio::read.tree("data/insecta_family.nwk")

# convert tree into tibble so we can join other data later
tbl_tree <- full_tree %>% treeio::as_tibble() 

# dplyr::filter(is.na(coverage), # filters those that do not have mean ch.area data
#               stringr::str_detect(label, pattern = regex("dae"))) # must be a label (some labels are numbers)

# prune tree --------------------------------------------------------------
## Figure out what data we have for AH and CH
tbl_tree_CH_df <- dplyr::left_join(tbl_tree, CH_df, by = c('label' = 'Family'))
# tbl_tree_AH_df <- dplyr::left_join(tbl_tree, AH_df, by = c('label' = 'Family'))

## Dropping tips CH
taxa_to_drop_CH <- tbl_tree_CH_df %>% 
  dplyr::filter(is.na(coverage_ch), # filters those that do not have mean ch.area data
                stringr::str_detect(label, pattern = regex("idae"))) %>% 
  dplyr::select(label)

tip_drop_ch <- taxa_to_drop_CH$label

## Prune tree to only ones we have name for 
sub_tree_ch <- full_tree %>%
  ape::drop.tip(., tip_drop_ch) # keep named tips

## Convert to dataframe for sub_tree_ah and sub_tree_ch
sub_tree_ch_df <- sub_tree_ch %>% treeio::as_tibble() 

# JOIN TREE DATA ----------------------------------------------------------
## join coverage data with tree data
tbl_tree_CH_df <- dplyr::left_join(sub_tree_ch_df, CH_df, by = c('label' = 'Family'))

# Plotting ----------------------------------------------------------------
## Plotting just data that have data

## Annotate Order on the phylogeny.
# Before doing that, we must get the node 'id' for each order.
node_label_df_CH <- tbl_tree_CH_df %>%
  dplyr::group_by(Order) %>%
  dplyr::summarise(minnode = min(node),
                   maxnode = max(node)) %>%
  ungroup() %>%
  as.data.frame()


# Convert table with phylo + coverage data back into phylo object
y_tree_ch <- ape::as.phylo(tbl_tree_CH_df) # full data set

# For loop to get number of node that represents the insect Order. We use this for annotating different orders later

## MCRA for CH
mcra_ch <- c() # create an empty vector for the for loop. We do this by finding the most recent common ancestors of the family.
for (i in 1:dim(node_label_df_CH)[1]-1) {
  mcra_ch[i] <- c(ape::getMRCA(y_tree_ch, tip = c(node_label_df_CH$minnode[i], node_label_df_CH$maxnode[i])))
}
mcra_ch_vector <- c(mcra_ch, NA)

# Append the empty vector as a column in our main dataframe
node_label_df_CH$MRCA_node <- mcra_ch_vector


## Convert back to phylo
y_tree_phylo_ch <- as.treedata(tbl_tree_CH_df)

tbl_tree_CH_df %>% 
  # distinct(label) %>% 
  filter(str_detect(label, pattern = regex("idae"))) %>% 
  glimpse()


hist(tbl_tree_CH_df$coverage_ch)

# PLOTTING ----------------------------------------------------------------

## CH Coverage

CH_tree <- ggtree(y_tree_phylo_ch, aes(colour = coverage_ch), layout = 'circular',
                  ladderize = TRUE, size = 0.3) +
  scale_colour_gradientn(name = "CH Protection Gap %", colours = pokemonpal(100), na.value = 'grey') +
  # scale_colour_gradient(name = "CH Protection Gap %", low = "red", high = "black", na.value = 'grey') +
  # scale_colour_stepsn(name = "% Coverage CH", colours = c("black", "lightblue", "blue", "orange", "darkgreen"), n.breaks = 6, nice.breaks = TRUE) +
  # scale_colour_stepsn(name = "% Coverage CH", colours = c("black", "lightblue", "darkblue", "orange", "#FF0000"), breaks = c(7, 8, 9,10,11)) +
  
  # Insect order labels
  # geom_cladelabel(node = node_label_df_CH$MRCA_node[1], label = paste(node_label_df_CH$Order[1]), barsize = 0.3, colour = 'black', fontsize = 2) +
  # geom_cladelabel(node = node_label_df_CH$MRCA_node[2], label = paste(node_label_df_CH$Order[2]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[3], label = paste(node_label_df_CH$Order[3]), barsize = 0.3, colour = 'black', fontsize = 2, hjust =2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[4], label = paste(node_label_df_CH$Order[4]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[5], label = paste(node_label_df_CH$Order[5]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[6], label = paste(node_label_df_CH$Order[6]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[7], label = paste(node_label_df_CH$Order[7]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[8], label = paste(node_label_df_CH$Order[8]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[9], label = paste(node_label_df_CH$Order[9]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[10], label = paste(node_label_df_CH$Order[10]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[11], label = paste(node_label_df_CH$Order[11]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[12], label = paste(node_label_df_CH$Order[12]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[13], label = paste(node_label_df_CH$Order[13]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=-1) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[14], label = paste(node_label_df_CH$Order[14]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[15], label = paste(node_label_df_CH$Order[15]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[16], label = paste(node_label_df_CH$Order[16]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[17], label = paste(node_label_df_CH$Order[17]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[18], label = paste(node_label_df_CH$Order[18]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[19], label = paste(node_label_df_CH$Order[19]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[20], label = paste(node_label_df_CH$Order[20]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[21], label = paste(node_label_df_CH$Order[21]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=-1) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[22], label = paste(node_label_df_CH$Order[22]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[23], label = paste(node_label_df_CH$Order[23]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[24], label = paste(node_label_df_CH$Order[24]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=1) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[25], label = paste(node_label_df_CH$Order[25]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[26], label = paste(node_label_df_CH$Order[26]), barsize = 0.3, colour = 'black', fontsize = 2, hjust=-0.5) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[27], label = paste(node_label_df_CH$Order[27]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[28], label = paste(node_label_df_CH$Order[28]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[29], label = paste(node_label_df_CH$Order[29]), barsize = 0.3, colour = 'black', fontsize = 2)

CH_tree

# class(CH_tree)
# ggsave(filename = 'output/CH_tree.pdf', plot = CH_tree)
ggsave(filename = 'output/CH_tree_gradient.pdf', plot = CH_tree)



# additional tips to add --------------------------------------------------
additional_df <- read.csv("data/additional_family_names.csv", header = TRUE)


## Upon joining, I see there are some families that do not have coverage data
# The following code filters those that do not have mean ch.area data
family_na_data <- tbl_tree_CH_df %>%
  dplyr::filter(is.na(mean.ch.area), # filters those that do not have mean ch.area data
                stringr::str_detect(label, pattern = regex("dae"))) # must be a label (some labels are numbers)

## There are about 150 families that do not have coverage data but this may be because synonyms or slightly different spelling. 
family_na_data 

# I suggest going back through your interpolation data to check if you have misspelled something or are able to fill data for these missing families.


