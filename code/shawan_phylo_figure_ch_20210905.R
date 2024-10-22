# filename: shawan_phylo_figure_ch_20210119.R

# R code to visualise conservation coverage for insect families

# two figures (a circular barplot, and a phylogenetic tree) for CH (convex hull). 

# Load libraries ----------------------------------------------------------
library(treeio)
library(ape)
library(ggtree)
library(tidyverse)

# Shawan's coverage data set ----------------------------------------------------------------
CH <- read.csv("data/CH_interpolation.csv", header = TRUE) %>% tibble()

# Correct classifications
# Grylloblattodea and MantoPhasmotodea are lumped under Notoptera (Matthew and Whiting 2005)
CH[which(CH[,'Family'] == 'Grylloblattidae' | CH[,'Family'] == 'Mantophasmatidae'), 'Order'] <- 'Notoptera'
CH[which(CH[,'Family'] == 'Braconidae'), 'Order'] <- 'Hymenoptera'
CH[which(CH[,'Family'] == 'Sialidae'), 'Order'] <- 'Megaloptera'
CH[which(CH[,'Family'] == 'Stenotritidae'), 'Order'] <- 'Hymenoptera'

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


result <- lm(prop_meet_target ~ Order, data = CH_df)
summary(result)

## Basic model

hist(CH_df$prop_meet_target)

summ_stat <- CH_df %>% 
  group_by(Order) %>% 
  summarise(order_meet_target = sum(meet_target),
            order_not_meet = sum(not_meet_target)) %>% 
  mutate(total = order_meet_target + order_not_meet,
         order_prop_meet = order_meet_target / (order_meet_target + order_not_meet) * 100)
          
CH_df %>% 
  group_by(Order) %>% 
  summarise(mean_meet_target =mean(prop_meet_target),
            median_meet_target = median(prop_meet_target)) %>% 
  arrange(mean_meet_target) %>% 
  filter(mean_meet_target < 40)


summ_stat %>% 
  arrange(desc(median_meet_target))

median(summ_stat$mean_meet_target) 

summ_stat %>% 
  arrange(mean_meet_target) %>% 
  ggplot(aes(x = reorder(Order, -mean_meet_target), y = mean_meet_target)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))
  
  
  
## How many species representing each family
CH_fam_rep <- CH %>% 
  dplyr::mutate(Family = str_extract(Family, regex('[A-Za-z]+'))) %>% 
  group_by(Family) %>% 
  count()

# write.csv(x = CH_fam_rep, file = "output/ch_fam_represented.csv", row.names = FALSE)

# tree from timetree ------------------------------------------------------
# insect family tree from timetree.org. Download newick tree file 
# read newick file with ggtree:read.tree
full_tree <- treeio::read.tree("data/insecta_family.nwk")

# convert tree into tibble so we can join other data later
tbl_tree <- full_tree %>% treeio::as_tibble() 

# Prune tree data --------------------------------------------------------------
## Figure out what data we have for CH
tbl_tree_CH_df_tmp <- dplyr::left_join(tbl_tree, CH_df, by = c('label' = 'Family')) 

tbl_tree_CH_df_tmp %>% filter(stringr::str_detect(label, pattern = regex("idae"))) %>%  distinct(label) 

## Dropping tips CH
taxa_to_drop_CH <- tbl_tree_CH_df_tmp %>% 
  dplyr::filter(is.na(prop_meet_target), # filters those that do not have mean percent_gap_CH data
                stringr::str_detect(label, pattern = regex("idae"))) %>% 
  dplyr::select(label)

tip_drop_CH <- taxa_to_drop_CH$label

## Species lost
CH %>% 
  filter(Family %in% tip_drop_CH)


## Prune tree to only ones we have name for 
sub_tree_CH <- full_tree %>%
  ape::drop.tip(., tip_drop_CH) # keep named tips

# plot(sub_tree_CH)

sub_tree_CH$tip.label

## Convert to dataframe for sub_tree_CH 
sub_tree_CH_df <- sub_tree_CH %>% treeio::as_tibble() 

CH_df_filt <- CH_df %>% 
  filter(Family %in% sub_tree_CH$tip.label) %>% 
  group_by(Family) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count))

## Join tree with sub tree data
tbl_tree_CH_df <- dplyr::left_join(sub_tree_CH_df, CH_df, by = c('label' = 'Family')) 

tbl_tree_CH_df

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
y_tree_phylo_CH <- as.treedata(tbl_tree_CH_df)

# Reverse of prop_meet_target is Percentages of species meeting representation target
y_tree_phylo_CH@data$per_target <- y_tree_phylo_CH@data$prop_meet_target 

# PLOTTING ----------------------------------------------------------------

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

# PLOTTING ----------------------------------------------------------------

tree1 <- ggtree(y_tree_phylo_CH, aes(colour = per_target), layout = 'circular',
                           ladderize = TRUE, size = 0.4) +
  # scale_fill_gradientn(colours=brewer.pal(n=100, name="PuBuGn"))
  scale_colour_gradientn(name = "", colours = ibmpal(100)[100:1], na.value = 'grey') +
  geom_treescale(x= 0, y = 0, width = 100, color = 'red', linesize = 1, offset = 1) + 
  # Insect order labels
  geom_cladelabel(node = node_label_df_CH$MRCA_node[1],  align = T, label = paste(node_label_df_CH$Order[1]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[2],  align = T, label = paste(node_label_df_CH$Order[2]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[3],  align = T, label = paste(node_label_df_CH$Order[3]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[4],  align = T, label = paste(node_label_df_CH$Order[4]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[5],  align = T, label = paste(node_label_df_CH$Order[5]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[6],  align = T, label = paste(node_label_df_CH$Order[6]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[7],  align = T, label = paste(node_label_df_CH$Order[7]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[8],  align = T, label = paste(node_label_df_CH$Order[8]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[9],  align = T, label = paste(node_label_df_CH$Order[9]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[10], align = T, label = paste(node_label_df_CH$Order[10]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[11], align = T, label = paste(node_label_df_CH$Order[11]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[12], align = T, label = paste(node_label_df_CH$Order[12]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[13], align = T, label = paste(node_label_df_CH$Order[13]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[14], align = T, label = paste(node_label_df_CH$Order[14]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[15], align = T, label = paste(node_label_df_CH$Order[15]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[16], align = T, label = paste(node_label_df_CH$Order[16]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[17], align = T, label = paste(node_label_df_CH$Order[17]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[18], align = T, label = paste(node_label_df_CH$Order[18]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[19], align = T, label = paste(node_label_df_CH$Order[19]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[20], align = T, label = paste(node_label_df_CH$Order[20]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[21], align = T, label = paste(node_label_df_CH$Order[21]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[22], align = T, label = paste(node_label_df_CH$Order[22]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[23], align = T, label = paste(node_label_df_CH$Order[23]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[24], align = T, label = paste(node_label_df_CH$Order[24]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[25], align = T, label = paste(node_label_df_CH$Order[25]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[26], align = T, label = paste(node_label_df_CH$Order[26]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[27], align = T, label = paste(node_label_df_CH$Order[27]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[28], align = T, label = paste(node_label_df_CH$Order[28]), barsize = 0.3, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = node_label_df_CH$MRCA_node[29], align = T, label = paste(node_label_df_CH$Order[29]), barsize = 0.3, colour = 'black', fontsize = 2)

tree1


# g2 <- tree_gradienter(colorRampPalette(c('red', 'blue', 'green')))


# PIE CHARTS --------------------------------------------------------------

pie_df_ch <- summ_stat %>% 
  mutate(order_prop_not_meet = 100 - order_prop_meet) %>% 
  select(Order, order_prop_meet, order_prop_not_meet) %>% 
  pivot_longer(!Order, names_to = "order_prop_meet")



pie_charts_order_proportion_ch <- ggplot(data=pie_df_ch, aes(x=" ", y=value, group=Order, colour=order_prop_meet, fill=order_prop_meet)) +
  geom_bar(width = 1, stat = "identity", colour = "white") +
  scale_fill_manual("Proportions CH", values=c("#648FFF", "#FFB000"), labels = c("Meet target", "Did not meet target")) +
  coord_polar("y") + 
  facet_wrap(.~ Order) +theme_void()




# PRINT -------------------------------------------------------------------

library(gridExtra)
print_plot <- list(tree1, pie_charts_order_proportion_ch)

dev.off()
# 
pdf("20210906_ch_gradients.pdf")

for (i in seq_along(print_plot)) {
  print(print_plot[[i]])
}

dev.off()



# class(CH_tree)
ggsave(filename = 'output/CH_tree.pdf', plot = g2)
# ggsave(filename = 'output/CH_tree_gradient.tiff', plot = g2)

# Circular plot -----------------------------------------------------------
source("code/circular_plotter_function.R")

full_p_ch <- circular_plotter(.df = CH_df)

ggsave(filename = 'output/CH_circular_no_lab.pdf', plot = full_p_ch)

