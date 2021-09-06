# filename: shawan_phylo_figure_20210119.R

# R code to visualise conservation coverage for insect families

# two figures (a circular barplot, and a phylogenetic tree) for both AH (alpha hull) and CH (convex hull). 

# Load libraries ----------------------------------------------------------
library(treeio)
library(ape)
library(ggtree)
library(tidyverse)

# Shawan's coverage data set ----------------------------------------------------------------
AH <- read.csv("data/AH_interpolation.csv", header = TRUE) %>% tibble()

# Correct classifications
# Grylloblattodea and MantoPhasmotodea are lumped under Notoptera (Matthew and Whiting 2005)
AH[which(AH[,'Order'] == 'Grylloblattodea' | AH[,'Order'] == 'MantoPhasmatodea'), 'Order'] <- 'Notoptera'
AH[which(AH[,'Family'] == 'Braconidae'), 'Order'] <- 'Hymenoptera'
AH[which(AH[,'Family'] == 'Sialidae'), 'Order'] <- 'Megaloptera'
AH[which(AH[,'Family'] == 'Stenotritidae'), 'Order'] <- 'Hymenoptera'

## Calculate the mean coverage for each insect family 
AH_df <- AH %>% 
  dplyr::mutate(Family = str_extract(Family, regex('[A-Za-z]+'))) %>% 
  arrange(Family) %>% 
  group_by(Family, Order) %>% 
  dplyr::mutate(pos = sum(gap > 0),
                neg = sum(gap <= 0)) %>% 
  dplyr::summarise(not_meet_target = mean(pos),
                   meet_target = mean(neg)) %>% 
  mutate(prop_meet_target = meet_target / (not_meet_target + meet_target) * 100) 

## How many species representing each family
AH_fam_rep <- AH %>% 
  dplyr::mutate(Family = str_extract(Family, regex('[A-Za-z]+'))) %>% 
  group_by(Family) %>% 
  arrange(Family)

# write.csv(x = AH_fam_rep, file = "output/ah_fam_represented.csv", row.names = FALSE)



## Basic patterns

hist(AH_df$prop_meet_target)

summ_stat <- AH_df %>% 
  group_by(Order) %>% 
  summarise(order_meet_target = sum(meet_target),
            order_not_meet = sum(not_meet_target)) %>% 
  mutate(total = order_meet_target + order_not_meet,
         order_prop_meet = order_meet_target / (order_meet_target + order_not_meet) * 100)

summ_stat %>% 
  arrange(desc(total, order_prop_meet)) %>% 
  filter(total <5000) %>% 
  ggplot(aes(x = total, y = order_prop_meet, label = Order)) +
  geom_text(size = 4)

plot(data = summ_stat, order_prop_meet ~ total)


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
  dplyr::filter(is.na(prop_meet_target), # filters those that do not have mean percent_gap_ah data
                stringr::str_detect(label, pattern = regex("idae"))) %>% 
  dplyr::select(label)

tip_drop_ah <- taxa_to_drop_AH$label

## Prune tree to only ones we have name for 
sub_tree_ah <- full_tree %>%
  ape::drop.tip(., tip_drop_ah) # keep named tips

# plot(sub_tree_ah)

sub_tree_ah$tip.label

## Convert to dataframe for sub_tree_ah 
sub_tree_ah_df <- sub_tree_ah %>% treeio::as_tibble() 

plot(as.phylo(sub_tree_ah_df))

AH_df_filt <- AH_df %>% 
  filter(Family %in% sub_tree_ah$tip.label) %>% 
  group_by(Family) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count))

## Join tree with sub tree data
tbl_tree_AH_df <- dplyr::left_join(sub_tree_ah_df, AH_df, by = c('label' = 'Family')) 

# tbl_tree_AH_df



# pie chart for order -----------------------------------------------------
### Pie chart showing proportion of targets met in each order. 2021-08-16

pie_stat_df <- summ_stat %>% 
  arrange(desc(total, order_prop_meet)) %>% 
  mutate(order_prop_not = 100 - order_prop_meet) %>% 
  select(Order, order_prop_meet, order_prop_not) %>% 
  pivot_longer(!Order, names_to = "Proportion", values_to = "value")
  
pie_charts_order_proportion <- ggplot(data = pie_stat_df, aes(x=" ", y = value, group = Proportion, colour = Proportion, fill = Proportion)) +
  geom_bar(width = 1, stat = "identity", colour = "white") +
  scale_fill_brewer(palette = "Set1") +
  coord_polar("y", start=0) +
  facet_wrap(.~ Order) +theme_void()
  
  
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

ggtree(y_tree_phylo_ah, aes(colour = prop_meet_target), ladderize = TRUE, size = 0.4)


### Gradient colours

library(colorRamps)

library(viridis)
viridis(n = 100, option = "D")


pokemonpal <- colorRampPalette(c('#CC79A7', '#0072B2', '#009E73'))
rev(pokemonpal(100)[100:1])

viridispal <- colorRampPalette(viridis(n = 100, option = "C"))
viridispal(10)[10:1]



# IBM Pallette
ibmpal <- colorRampPalette(c('#648FFF', '#785EF0', '#DC267F', '#FE6100', '#FFB000'))

AH_tree_gradient <- ggtree(y_tree_phylo_ah, aes(colour = prop_meet_target), layout = 'circular',
                           ladderize = TRUE, size = 0.4) +
  scale_colour_gradientn(name = "", colours = ibmpal(100)[100:1], na.value = 'grey') +
  geom_treescale(x= 0, y = 0, width = 100, color = 'red', linesize = 1, offset = 1) + 
  # scale_colour_stepsn(name = "% Coverage CH", colours = c("black", "lightblue", "blue", "orange", "darkgreen"), n.breaks = 6, nice.breaks = TRUE) +
  # scale_colour_stepsn(name = "% Coverage CH", colours = c("black", "lightblue", "darkblue", "orange", "#FF0000"), breaks = c(7, 8, 9,10,11)) +
  
  # Insect order labels
  geom_cladelabel(node = node_label_df_AH$MRCA_node[1], label = paste(node_label_df_AH$Order[1]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[2], label = paste(node_label_df_AH$Order[2]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[3], label = paste(node_label_df_AH$Order[3]), hjust = +1, barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[4], label = paste(node_label_df_AH$Order[4]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[5], label = paste(node_label_df_AH$Order[5]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[6], label = paste(node_label_df_AH$Order[6]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[7], label = paste(node_label_df_AH$Order[7]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[8], label = paste(node_label_df_AH$Order[8]), hjust = +1, barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[9], label = paste(node_label_df_AH$Order[9]), hjust = +1, barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[10], label = paste(node_label_df_AH$Order[10]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[11], label = paste(node_label_df_AH$Order[11]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[12], label = paste(node_label_df_AH$Order[12]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[13], label = paste(node_label_df_AH$Order[13]), hjust = +1, barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[14], label = paste(node_label_df_AH$Order[14]), hjust = +1,  barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[15], label = paste(node_label_df_AH$Order[15]), hjust = +1, barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[16], label = paste(node_label_df_AH$Order[16]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[17], label = paste(node_label_df_AH$Order[17]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[18], label = paste(node_label_df_AH$Order[18]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[19], label = paste(node_label_df_AH$Order[19]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[20], label = paste(node_label_df_AH$Order[20]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[21], label = paste(node_label_df_AH$Order[21]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[22], label = paste(node_label_df_AH$Order[22]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[23], label = paste(node_label_df_AH$Order[23]), hjust = +1, barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[24], label = paste(node_label_df_AH$Order[24]), hjust = +1, barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[25], label = paste(node_label_df_AH$Order[25]), hjust = +1, barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[26], label = paste(node_label_df_AH$Order[26]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[27], label = paste(node_label_df_AH$Order[27]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[28], label = paste(node_label_df_AH$Order[28]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_AH$MRCA_node[29], label = paste(node_label_df_AH$Order[29]), barsize = 0.7, colour = 'black', fontsize = 2)

AH_tree_gradient

# class(CH_tree)
# ggsave(filename = 'output/CH_tree.pdf', plot = AH_tree_gr)
# ggsave(filename = 'output/AH_tree_gradient_1.pdf', plot = AH_tree_gradient)
# ggsave(filename = 'output/AH_tree_gradient.tiff', plot = AH_tree_gradient)


# PIE CHART ---------------------------------------------------------------

pie_df <- summ_stat %>% 
  mutate(order_prop_not_meet = 100 - order_prop_meet) %>% 
  select(Order, order_prop_meet, order_prop_not_meet) %>% 
  pivot_longer(!Order, names_to = "order_prop_meet")



pie_charts_order_proportion <- ggplot(data=pie_df, aes(x=" ", y=value, group=Order, colour=order_prop_meet, fill=order_prop_meet)) +
  geom_bar(width = 1, stat = "identity", colour = "white") +
  scale_fill_manual("Proportions", values=c("#648FFF", "#FFB000"), labels = c("Meet target", "Did not meet target")) +
  coord_polar("y") + 
  facet_wrap(.~ Order) +theme_void()



# PRINT -------------------------------------------------------------------


print_plot <- list(AH_tree_gradient, pie_charts_order_proportion)

dev.off()
# 
pdf("20210610_ah_gradients.pdf")

for (i in seq_along(print_plot)) {
  print(print_plot[[i]])
}

dev.off()