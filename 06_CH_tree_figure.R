# filename: CH_tree_figure.R

# R code to visualise conservation coverage for insect families
# Make sure all packages and dependencies are installed. 
# two figures (a circular barplot, and a phylogenetic tree) for both CH (alpha hull) and CH (convex hull). 

# Load libraries ----------------------------------------------------------
library(treeio); library(ape); library(ggtree); library(tidyverse)

# Shawan's coverage data set ----------------------------------------------------------------
CH <- read.csv("data/CH_interpolation.csv", header = TRUE) %>% tibble()

# Correct classifications
# Grylloblattodea and MantoPhasmotodea are lumped under Notoptera (Matthew and Whiting 2005)
CH[which(CH[,'Order'] == 'Grylloblattodea' | CH[,'Order'] == 'MantoPhasmatodea'), 'Order'] <- 'Notoptera'
CH[which(CH[,'Family'] == 'Braconidae'), 'Order'] <- 'Hymenoptera'
CH[which(CH[,'Family'] == 'Sialidae'), 'Order'] <- 'Megaloptera'
CH[which(CH[,'Family'] == 'Stenotritidae'), 'Order'] <- 'Hymenoptera'

# Isoptera now part of Blattodea
CH[which(CH[,'Order'] == 'Isoptera'), 'Order'] <- 'Blattodea'

## Calculate the mean coverage for each insect family 
CH_df <- CH %>% 
  dplyr::mutate(Family = str_extract(Family, regex('[A-Za-z]+'))) %>% 
  arrange(Family) %>% 
  group_by(Family, Order) %>% 
  dplyr::mutate(pos = sum(gap > 0),
                neg = sum(gap <= 0)) %>% 
  dplyr::summarise(not_meet_target = mean(pos),
                   meet_target = mean(neg)) %>% 
  mutate(prop_meet_target = meet_target / (not_meet_target + meet_target) * 100) 

## How many species representing each family
CH_fam_rep <- CH %>% 
  dplyr::mutate(Family = str_extract(Family, regex('[A-Za-z]+'))) %>% 
  group_by(Family) %>% 
  arrange(Family)

# write.csv(x = CH_fam_rep, file = "output/ah_fam_represented.csv", row.names = FALSE)

## Basic patterns

hist(CH_df$prop_meet_target)

summ_stat <- CH_df %>% 
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

# tree from timetree ------------------------------------------------------
# insect family tree from timetree.org. Download newick tree file 
# read newick file with ggtree:read.tree
full_tree <- treeio::read.tree("data/insecta_family.nwk")

# convert tree into tibble so we can join other data later
tbl_tree <- full_tree %>% treeio::as_tibble() 

# Prune tree data --------------------------------------------------------------
## Figure out what data we have for CH
tbl_tree_CH_df_tmp <- dplyr::left_join(tbl_tree, CH_df, by = c('label' = 'Family')) 

## Dropping tips CH
taxa_to_drop_CH <- tbl_tree_CH_df_tmp %>% 
  dplyr::filter(is.na(prop_meet_target), # filters those that do not have mean percent_gap_ah data
                stringr::str_detect(label, pattern = regex("idae"))) %>% 
  dplyr::select(label)

tip_drop_ah <- taxa_to_drop_CH$label

## Prune tree to only ones we have name for 
sub_tree_ah <- full_tree %>%
  ape::drop.tip(., tip_drop_ah) # keep named tips

# plot(sub_tree_ah)
# nodelabels(text=1:sub_tree_ah$Nnode,node=1:sub_tree_ah$Nnode+Ntip(sub_tree_ah))

sub_tree_ah$tip.label

## Convert to dataframe for sub_tree_ah 
sub_tree_ah_df <- sub_tree_ah %>% treeio::as_tibble() 

# plot(as.phylo(sub_tree_ah_df))

CH_df_filt <- CH_df %>% 
  dplyr::filter(Family %in% sub_tree_ah$tip.label) %>% 
  dplyr::group_by(Family) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::arrange(desc(count))

## Join tree with sub tree data
tbl_tree_CH_df <- dplyr::left_join(sub_tree_ah_df, CH_df, by = c('label' = 'Family')) 

# pie chart for order -----------------------------------------------------
### Pie chart showing proportion of targets met in each order. 2021-08-16

pie_stat_df <- summ_stat %>% 
  arrange(desc(order_prop_meet)) %>% 
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
node_label_df_CH <- tbl_tree_CH_df %>%
  dplyr::group_by(Order) %>%
  dplyr::summarise(minnode = min(node),
                   maxnode = max(node)) %>%
  ungroup() %>%
  as.data.frame()

# Convert table with phylo + coverage data back into phylo object
y_tree_ch <- ape::as.phylo(tbl_tree_CH_df) # full data set

# Get number of node that represents the insect Order. We use this for annotating different orders later

## MCRA for CH
mcra_ch <- c() # create an empty vector for the for loop. We do this by finding the most recent common ancestors of the family.

for (i in 1:dim(node_label_df_CH)[1]-1) {
  
  mcra_ch[i] <- c(ape::getMRCA(y_tree_ch, tip = c(node_label_df_CH$minnode[i], node_label_df_CH$maxnode[i])))
}

mcra_ch_vector <- c(mcra_ch, NA)

# Append the empty vector as a column in our main dataframe
node_label_df_CH$MRCA_node <- mcra_ch_vector


## Convert back to phylo
y_tree_phylo_ah <- as.treedata(tbl_tree_CH_df)

# PLOTTING ----------------------------------------------------------------

ggtree(y_tree_phylo_ah, aes(colour = prop_meet_target), ladderize = TRUE, size = 0.4)


### Gradient colours

# library(colorRamps)
# library(viridis)
# viridis(n = 100, option = "D")
# pokemonpal <- colorRampPalette(c('#CC79A7', '#0072B2', '#009E73'))
# rev(pokemonpal(100)[100:1])
# 
# viridispal <- colorRampPalette(viridis(n = 100, option = "C"))
# viridispal(10)[10:1]

# IBM Pallette
ibmpal <- colorRampPalette(c('#648FFF', '#785EF0', '#DC267F', '#FE6100', '#FFB000'))

CH_tree_gradient <- ggtree(y_tree_phylo_ah, aes(colour = prop_meet_target), layout = 'circular',
                           ladderize = TRUE, size = 0.4) +
  scale_colour_gradientn(name = "", colours = ibmpal(100)[100:1], na.value = 'grey') +
  geom_treescale(x= 0, y = 0, width = 100, color = 'red', linesize = 1, offset = 1) + 
  # scale_colour_stepsn(name = "% Coverage CH", colours = c("black", "lightblue", "blue", "orange", "darkgreen"), n.breaks = 6, nice.breaks = TRUE) +
  # scale_colour_stepsn(name = "% Coverage CH", colours = c("black", "lightblue", "darkblue", "orange", "#FF0000"), breaks = c(7, 8, 9,10,11)) +
  # Insect order labels
  geom_cladelabel(node = node_label_df_CH$MRCA_node[1], label = paste(node_label_df_CH$Order[1]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[2], label = paste(node_label_df_CH$Order[2]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[3], label = paste(node_label_df_CH$Order[3]), hjust = +1, barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[4], label = paste(node_label_df_CH$Order[4]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[5], label = paste(node_label_df_CH$Order[5]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[6], label = paste(node_label_df_CH$Order[6]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[7], label = paste(node_label_df_CH$Order[7]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[8], label = paste(node_label_df_CH$Order[8]), hjust = +1, barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[9], label = paste(node_label_df_CH$Order[9]), hjust = +1, barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[10], label = paste(node_label_df_CH$Order[10]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[11], label = paste(node_label_df_CH$Order[11]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[12], label = paste(node_label_df_CH$Order[12]), hjust = +1, barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[13], label = paste(node_label_df_CH$Order[13]), hjust = +1, barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[14], label = paste(node_label_df_CH$Order[14]), hjust = +1,  barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[15], label = paste(node_label_df_CH$Order[15]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[16], label = paste(node_label_df_CH$Order[16]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[17], label = paste(node_label_df_CH$Order[17]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[18], label = paste(node_label_df_CH$Order[18]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[19], label = paste(node_label_df_CH$Order[19]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[20], label = paste(node_label_df_CH$Order[20]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[21], label = paste(node_label_df_CH$Order[21]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[22], label = paste(node_label_df_CH$Order[22]), hjust = +1, barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[23], label = paste(node_label_df_CH$Order[23]), hjust = +1, barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[24], label = paste(node_label_df_CH$Order[24]), hjust = +1, barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[25], label = paste(node_label_df_CH$Order[25]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[26], label = paste(node_label_df_CH$Order[26]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[27], label = paste(node_label_df_CH$Order[27]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[28], label = paste(node_label_df_CH$Order[28]), barsize = 0.7, colour = 'black', fontsize = 2) 

CH_tree_gradient

# class(CH_tree)
# ggsave(filename = 'output/CH_tree.pdf', plot = CH_tree_gr)
# ggsave(filename = 'output/CH_tree_gradient_1.pdf', plot = CH_tree_gradient)
# ggsave(filename = 'output/CH_tree_gradient.tiff', plot = CH_tree_gradient)


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
# print_plot <- list(CH_tree_gradient, pie_charts_order_proportion)
# 
# dev.off()
# # 
# pdf("20220218_ch_gradients.pdf")
# 
# for (i in seq_along(print_plot)) {
#   print(print_plot[[i]])
# }
# 
# dev.off()
