# filename: AH_tree_figure.R

# Make sure all packages and dependencies are installed. 

# Libraries ---------------------------------------------------------------
library(stringr)
library(ape); library(phytools); library(treeio)
library(geiger); library(dplyr); 
library(ggplot2); library(ggtree)
library(colorRamps) # colour palette 


# Read data ---------------------------------------------------------------
# Insect family tree downloaded from timetree.org
tree <- ape::read.tree('data/insecta_family.nwk')

# Check tree by plotting
phytools::plotTree(tree)

# Read in data
dat <- read.csv("data/AH_interpolation.csv", header = TRUE) %>% tibble()

# Data corrections --------------------------------------------------------
# Correct classifications
# Grylloblattodea and MantoPhasmotodea are lumped under Notoptera (Matthew and Whiting 2005)
dat[which(dat[,'Order'] == 'Grylloblattodea' | dat[,'Order'] == 'MantoPhasmatodea'), 'Order'] <- 'Notoptera'
# Isoptera order now part of Blattodea
dat[which(dat[,'Order'] == 'Isoptera'), 'Order'] <- 'Blattodea'
# Phthiraptera and Psocoptera now lumped into Psocodea (de Moya et al 2021. Systematic Biol)
dat[which(dat[,'Order'] == 'Phthiraptera'), 'Order'] <- 'Psocodea'
dat[which(dat[,'Order'] == 'Psocoptera'), 'Order'] <- 'Psocodea'

# Family changes
dat[which(dat[,'Family'] == 'Braconidae'), 'Order'] <- 'Hymenoptera'
dat[which(dat[,'Family'] == 'Sialidae'), 'Order'] <- 'Megaloptera'
dat[which(dat[,'Family'] == 'Stenotritidae'), 'Order'] <- 'Hymenoptera'

## Calculate the mean coverage for each insect family 
dat_df <- dat %>% 
  dplyr::mutate(Family = str_extract(Family, regex('[A-Za-z]+'))) %>% 
  dplyr::arrange(Family) %>% 
  dplyr::group_by(Family, Order) %>% 
  dplyr::mutate(pos = sum(gap > 0),
                neg = sum(gap <= 0)) %>% 
  dplyr::summarise(not_meet_target = mean(pos),
                   meet_target = mean(neg),
                   mean_target_proportion = mean(target_proportion),
                   mean_gap = mean(gap),
                   mean_ah = mean(AH), 
                   mean_area = mean(ah.area),
                   mean_ov = mean(ov.ah)) %>% 
  dplyr::mutate(prop_meet_target = meet_target / (not_meet_target + meet_target) * 100) %>% 
  dplyr::select(Order, Family, mean_area, mean_target_proportion, everything()) %>% 
  dplyr::arrange(Order, Family)

# write.csv(dat_df, file = 'output/summary-overall-ah.csv')

dat_df %>% distinct(Family)


# PRUNE TREE --------------------------------------------------------------

# Cross reference data with tree tips, first need to convet to data frame
dat_df <- as.data.frame(dat_df)
rownames(dat_df) <- dat_df$Family

head(dat_df)

## GEIGER::NAME.CHECK
chk_name <- geiger::name.check(phy = tree, data = dat_df)
# vector of taxa in the tree that does not have any data
chk_name$tree_not_data

## PRUNING THE TREE
sub_tree <- ape::drop.tip(phy = tree, tip = chk_name$tree_not_data) 
plotTree(sub_tree)

# Joining data into tree

## Convert to dataframe for sub_tree 
sub_tree_df <- sub_tree %>% treeio::as_tibble() 

dat_df_filt <- dat_df %>% 
  dplyr::filter(Family %in% sub_tree$tip.label) %>% 
  dplyr::group_by(Family) %>% 
  dplyr::summarise(count = n()) %>% 
  dplyr::arrange(desc(count))

summary_in_tree <- dat_df %>% 
  dplyr::filter(Family %in% sub_tree$tip.label)

# write.csv(summary_in_tree, file = 'output/summary-data-of-tree-ah.csv')

## Join tree with sub tree data
tbl_tree_dat_df <- dplyr::left_join(sub_tree_df, dat_df, by = c('label' = 'Family')) 

# Convert table with phylo + coverage data back into phylo object
phylo_data <- ape::as.phylo(tbl_tree_dat_df) # full data set
plot(phylo_data)


#### Convert back tree object with treeio::as.treedata for plotting with ggtree ####
joined_phylo <- treeio::as.treedata(tbl_tree_dat_df)

#### test plotting with ggtree ####
ggtree(joined_phylo, aes(colour = prop_meet_target), ladderize = TRUE, size = 0.4)

# Node labels -------------------------------------------------------------

node_label_data <- tbl_tree_dat_df %>%
  dplyr::group_by(Order) %>%
  dplyr::summarise(minnode = min(node),
                   maxnode = max(node)) %>%
  ungroup() %>%
  as.data.frame()

# Most recent common ancestor

# MRCA for data
mrca <- c() # create an empty vector for the for loop. We do this by finding the most recent common ancestors of the family.

for (i in 1:dim(node_label_data)[1]-1) {
  
  mrca[i] <- c(ape::getMRCA(phylo_data, tip = c(node_label_data$minnode[i], node_label_data$maxnode[i])))
}

mrca_vector <- c(mrca, NA)

# Append the empty vector as a column in our main dataframe
node_label_data$MRCA_node <- mrca_vector
node_label_data

# GGTREE + COLOURS --------------------------------------------------------
# IBM Pallette
ibmpal <- colorRampPalette(c('#648FFF', '#785EF0', '#DC267F', '#FE6100', '#FFB000'))

# shorten node_label_data to nd for shorter code
nd <- node_label_data
nd

tree_plot <- ggtree(joined_phylo, aes(colour = prop_meet_target), ladderize = TRUE, size = 0.4, layout = 'circular')+
  scale_colour_gradientn(name = "", colours = ibmpal(100)[100:1], na.value = 'grey') +
  geom_treescale(x= 0, y = 0, width = 100, color = 'red', linesize = 1, offset = 1) + 
  geom_cladelabel(node = nd$MRCA_node[1], label = paste(nd$Order[1]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = nd$MRCA_node[2], label = paste(nd$Order[2]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = nd$MRCA_node[3], label = paste(nd$Order[3]), hjust = +1, barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = nd$MRCA_node[4], label = paste(nd$Order[4]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = nd$MRCA_node[5], label = paste(nd$Order[5]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = nd$MRCA_node[6], label = paste(nd$Order[6]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = nd$MRCA_node[7], label = paste(nd$Order[7]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = nd$MRCA_node[8], label = paste(nd$Order[8]), hjust = +1, barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = nd$MRCA_node[9], label = paste(nd$Order[9]), hjust = +1, barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = nd$MRCA_node[10], label = paste(nd$Order[10]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = nd$MRCA_node[11], label = paste(nd$Order[11]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = nd$MRCA_node[12], label = paste(nd$Order[12]), hjust = +1, barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = nd$MRCA_node[13], label = paste(nd$Order[13]), hjust = +1, barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = nd$MRCA_node[14], label = paste(nd$Order[14]), hjust = +1,  barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = nd$MRCA_node[15], label = paste(nd$Order[15]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = nd$MRCA_node[16], label = paste(nd$Order[16]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = nd$MRCA_node[17], label = paste(nd$Order[17]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = nd$MRCA_node[18], label = paste(nd$Order[18]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = nd$MRCA_node[19], label = paste(nd$Order[19]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = nd$MRCA_node[20], label = paste(nd$Order[20]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = nd$MRCA_node[21], label = paste(nd$Order[21]), hjust = +1, barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = nd$MRCA_node[22], label = paste(nd$Order[22]), hjust = +1, barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = nd$MRCA_node[23], label = paste(nd$Order[23]), hjust = +1, barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = nd$MRCA_node[24], label = paste(nd$Order[24]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = nd$MRCA_node[25], label = paste(nd$Order[25]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = nd$MRCA_node[26], label = paste(nd$Order[26]), barsize = 0.7, colour = 'black', fontsize = 2) +
  geom_cladelabel(node = nd$MRCA_node[27], label = paste(nd$Order[27]), barsize = 0.7, colour = 'black', fontsize = 2) +
  theme(legend.position = c(0.56, 0.45), legend.direction = 'horizontal')

tree_plot
# Other data --------------------------------------------------------------

## Basic patterns
hist(dat_df$prop_meet_target)

# Summary statistics
summ_stat <- dat_df %>% 
  dplyr::group_by(Order) %>% 
  dplyr::summarise(order_meet_target = sum(meet_target),
                   order_not_meet = sum(not_meet_target)) %>% 
  dplyr::mutate(total = order_meet_target + order_not_meet,
                order_prop_meet = order_meet_target / (order_meet_target + order_not_meet) * 100)

summ_stat %>% 
  arrange(desc(order_prop_meet)) %>% 
  filter(total <5000) %>% 
  ggplot(aes(x = total, y = order_prop_meet, label = Order)) +
  geom_point() +
  geom_text(size = 4)

# Summary data of the tree
dat_df %>% 
  dplyr::group_by(Family, Order) %>% 
  dplyr::summarise(order_meet_target = sum(meet_target),
                   order_not_meet = sum(not_meet_target)) %>% 
  dplyr::mutate(total = order_meet_target + order_not_meet,
                order_prop_meet = order_meet_target / (order_meet_target + order_not_meet) * 100)


# Pie chart for order -----------------------------------------------------
### Pie chart showing proportion of targets met in each order. 2021-08-16
pie_stat_df <- summ_stat %>% 
  dplyr::arrange(desc(total, order_prop_meet)) %>% 
  dplyr::mutate(order_prop_not = 100 - order_prop_meet) %>% 
  dplyr::select(Order, order_prop_meet, order_prop_not) %>% 
  tidyr::pivot_longer(!Order, names_to = "Proportion", values_to = "value")

pie_charts_order_proportion <- ggplot(data = pie_stat_df, aes(x=" ", y = value, group = Proportion, colour = Proportion, fill = Proportion)) +
  geom_bar(width = 1, stat = "identity", colour = "white") +
  scale_fill_brewer(palette = "Set1") +
  scale_fill_manual("Proportions AH", values=c("#648FFF", "#FFB000"), labels = c("Meet target", "Did not meet target"))+
  coord_polar("y", start=0) +
  facet_wrap(.~ Order) +theme_void() 

# PRINT -------------------------------------------------------------------

# print_plot <- list(tree_plot, pie_charts_order_proportion)

# dev.off()
# # 
# pdf("20220221_ah_gradients.pdf")
# 
# for (i in seq_along(print_plot)) {
#   print(print_plot[[i]])
# }
# 
# dev.off()
