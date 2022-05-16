#Produce a phylogenetic tree of all strains top 1000 ranked by readscore then readcount
#Tree was made in MEGA from file
#ALL_0.1percent_allTRUE_normReadcount_normReadscore_top1000_wStrainName
#.nwk tree was read into R

library(ape)
library(dplyr)
library(TreeTools)
library(ips)
#library(msa)
library(ggtree)
library(tidyr)
library(ggplot2)
library(chisq.posthoc.test)
#library(adegenet)
#library(phangorn)
#library(stats)
#library(ade4)


setwd("~/Box/centromere_ChIP_seq/CENP-A_ChIP/data/readscore/RankReadscoreReadcount")
top1000_all <- read.tree(file = "AllStrains_top1000.nwk")


#preorder for subtree function
preorder_top1000_all <- Preorder(top1000_all)

tibbletree <- as_tibble(preorder_top1000_all)

#From tibble tree, take top length branches and partition

#5337
#308 tips
tree1 <- Subtree(preorder_top1000_all, 5337)

#5652
#171 tips
tree2 <- Subtree(preorder_top1000_all, 5652)

#6835
#1165 tips
tree3 <- Subtree(preorder_top1000_all, 6835)

#6597
#239 tips
tree4 <- Subtree(preorder_top1000_all, 6597)

#4011
#1326 tips
tree5 <- Subtree(preorder_top1000_all, 4011)

#5822
#200 tips
tree6 <- Subtree(preorder_top1000_all, 5822)

nodeslist <- c(5337, 5652, 6835, 6597, 4011, 5882, 4001)

#Make sure the tip names don't overlap
intersect(tree1$tip.label, tree2$tip.label)
intersect(tree1$tip.label, tree3$tip.label)
intersect(tree1$tip.label, tree4$tip.label)
intersect(tree1$tip.label, tree5$tip.label)
intersect(tree1$tip.label, tree6$tip.label)

intersect(tree2$tip.label, tree3$tip.label)
intersect(tree2$tip.label, tree4$tip.label)
intersect(tree2$tip.label, tree5$tip.label)
intersect(tree2$tip.label, tree6$tip.label)

intersect(tree3$tip.label, tree4$tip.label)
intersect(tree3$tip.label, tree5$tip.label)
intersect(tree3$tip.label, tree6$tip.label)

intersect(tree4$tip.label, tree5$tip.label)
intersect(tree4$tip.label, tree6$tip.label)

intersect(tree5$tip.label, tree6$tip.label)


t1tiplabel <- as.data.frame(tree1$tip.label)
colnames(t1tiplabel) <- c('TipLabel')
t2tiplabel <- as.data.frame(tree2$tip.label)
colnames(t2tiplabel) <- c('TipLabel')
t3tiplabel <- as.data.frame(tree3$tip.label)
colnames(t3tiplabel) <- c('TipLabel')
t4tiplabel <- as.data.frame(tree4$tip.label)
colnames(t4tiplabel) <- c('TipLabel')
t5tiplabel <- as.data.frame(tree5$tip.label)
colnames(t5tiplabel) <- c('TipLabel')
t6tiplabel <- as.data.frame(tree6$tip.label)
colnames(t6tiplabel) <- c('TipLabel')

tiplabel_tree <- rbind(t1tiplabel,t2tiplabel,t3tiplabel,t4tiplabel,t5tiplabel,t6tiplabel)
tiplabel_tree <- tiplabel_tree$TipLabel
tiplabel_all <- preorder_top1000_all$tip.label
tiplabels_diff <- setdiff(tiplabel_all,tiplabel_tree)

tiplabels_diff <- as.data.frame(tiplabels_diff)

#separate column name by "_" character and count number of strain values for pie chart

#Node1
t1tiplabel <- separate(t1tiplabel, col = TipLabel, into = c("Strain","Score","NormFreq","NormScore","rank"), sep = "_")

t2tiplabel <- separate(t2tiplabel, col = TipLabel, into = c("Strain","Score","NormFreq","NormScore","rank"), sep = "_")

t3tiplabel <- separate(t3tiplabel, col = TipLabel, into = c("Strain","Score","NormFreq","NormScore","rank"), sep = "_")

t4tiplabel <- separate(t4tiplabel, col = TipLabel, into = c("Strain","Score","NormFreq","NormScore","rank"), sep = "_")

t5tiplabel <- separate(t5tiplabel, col = TipLabel, into = c("Strain","Score","NormFreq","NormScore","rank"), sep = "_")

t6tiplabel <- separate(t6tiplabel, col = TipLabel, into = c("Strain","Score","NormFreq","NormScore","rank"), sep = "_")

tiplabels_diff <- separate(tiplabels_diff, col = tiplabels_diff, into = c("Strain","Score","NormFreq","NormScore","rank"), sep = "_")

##Plot tree (with scale bar)

plot(top1000_all, show.node.label = FALSE)
add.scale.bar()


#######
#donut plot
#Node1
t1tiplabel %>% 
  group_by(Strain) %>%
  dplyr::summarize(count = n()) %>%
  as.data.frame() -> t1tiplabel

# Compute percentages
t1tiplabel$fraction <- t1tiplabel$count / sum(t1tiplabel$count)

# Compute the cumulative percentages (top of each rectangle)
t1tiplabel$ymax <- cumsum(t1tiplabel$fraction)

# Compute the bottom of each rectangle
t1tiplabel$ymin <- c(0, head(t1tiplabel$ymax, n=-1))

# Compute label position
t1tiplabel$labelPosition <- (t1tiplabel$ymax + t1tiplabel$ymin) / 2

# Compute a good label
t1tiplabel$label <- paste0(t1tiplabel$Strain, "\n value: ", t1tiplabel$count)

# Make the plot
ggplot(t1tiplabel, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Strain)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

t1tiplabel$Expected <- sum(t1tiplabel$count)/4
t1tiplabel$OminusE <- t1tiplabel$count-t1tiplabel$Expected
t1tiplabel$OminusEsquared <- (t1tiplabel$OminusE)^2
t1tiplabel$chiSQ <- t1tiplabel$OminusEsquared/t1tiplabel$Expected
1 - pchisq(q=sum(t1tiplabel$chiSQ), df=3)
#[1] 0.07747818

###
t1tiplabel$Strain <- as.factor(t1tiplabel$Strain)
chisq.test(t1tiplabel[,2])
#X-squared = 6.8312, df = 3, p-value = 0.07748
chisq.posthoc.test(t1tiplabel[,2])

#t1tiplabel$node <- 5337


#Node2
t2tiplabel %>% 
  group_by(Strain) %>%
  dplyr::summarize(count = n()) %>%
  as.data.frame() -> t2tiplabel

# Compute percentages
t2tiplabel$fraction <- t2tiplabel$count / sum(t2tiplabel$count)

# Compute the cumulative percentages (top of each rectangle)
t2tiplabel$ymax <- cumsum(t2tiplabel$fraction)

# Compute the bottom of each rectangle
t2tiplabel$ymin <- c(0, head(t2tiplabel$ymax, n=-1))

# Compute label position
t2tiplabel$labelPosition <- (t2tiplabel$ymax + t2tiplabel$ymin) / 2

# Compute a good label
t2tiplabel$label <- paste0(t2tiplabel$Strain, "\n value: ", t2tiplabel$count)

# Make the plot
ggplot(t2tiplabel, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Strain)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

t2tiplabel$Expected <- sum(t2tiplabel$count)/4
t2tiplabel$OminusE <- t2tiplabel$count-t2tiplabel$Expected
t2tiplabel$OminusEsquared <- (t2tiplabel$OminusE)^2
t2tiplabel$chiSQ <- t2tiplabel$OminusEsquared/t2tiplabel$Expected
1 - pchisq(q=sum(t2tiplabel$chiSQ), df=3)
#[1] 0.8238752

#t2tiplabel$node <- 5652


#Node3
t3tiplabel %>% 
  group_by(Strain) %>%
  dplyr::summarize(count = n()) %>%
  as.data.frame() -> t3tiplabel

# Compute percentages
t3tiplabel$fraction <- t3tiplabel$count / sum(t3tiplabel$count)

# Compute the cumulative percentages (top of each rectangle)
t3tiplabel$ymax <- cumsum(t3tiplabel$fraction)

# Compute the bottom of each rectangle
t3tiplabel$ymin <- c(0, head(t3tiplabel$ymax, n=-1))

# Compute label position
t3tiplabel$labelPosition <- (t3tiplabel$ymax + t3tiplabel$ymin) / 2

# Compute a good label
t3tiplabel$label <- paste0(t3tiplabel$Strain, "\n value: ", t3tiplabel$count)

# Make the plot
ggplot(t3tiplabel, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Strain)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

t3tiplabel$Expected <- sum(t3tiplabel$count)/4
t3tiplabel$OminusE <- t3tiplabel$count-t3tiplabel$Expected
t3tiplabel$OminusEsquared <- (t3tiplabel$OminusE)^2
t3tiplabel$chiSQ <- t3tiplabel$OminusEsquared/t3tiplabel$Expected
1 - pchisq(q=sum(t3tiplabel$chiSQ), df=3)
#[1] 0.005489182

#t3tiplabel$node <- 6835


#Node4
t4tiplabel %>% 
  group_by(Strain) %>%
  dplyr::summarize(count = n()) %>%
  as.data.frame() -> t4tiplabel

# Compute percentages
t4tiplabel$fraction <- t4tiplabel$count / sum(t4tiplabel$count)

# Compute the cumulative percentages (top of each rectangle)
t4tiplabel$ymax <- cumsum(t4tiplabel$fraction)

# Compute the bottom of each rectangle
t4tiplabel$ymin <- c(0, head(t4tiplabel$ymax, n=-1))

# Compute label position
t4tiplabel$labelPosition <- (t4tiplabel$ymax + t4tiplabel$ymin) / 2

# Compute a good label
t4tiplabel$label <- paste0(t4tiplabel$Strain, "\n value: ", t4tiplabel$count)

# Make the plot
ggplot(t4tiplabel, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Strain)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

t4tiplabel$Expected <- sum(t4tiplabel$count)/4
t4tiplabel$OminusE <- t4tiplabel$count-t4tiplabel$Expected
t4tiplabel$OminusEsquared <- (t4tiplabel$OminusE)^2
t4tiplabel$chiSQ <- t4tiplabel$OminusEsquared/t4tiplabel$Expected
1 - pchisq(q=sum(t4tiplabel$chiSQ), df=3)
#[1] 0.1326343

#t4tiplabel$node <- 6597

#Node5
t5tiplabel %>% 
  group_by(Strain) %>%
  dplyr::summarize(count = n()) %>%
  as.data.frame() -> t5tiplabel

# Compute percentages
t5tiplabel$fraction <- t5tiplabel$count / sum(t5tiplabel$count)

# Compute the cumulative percentages (top of each rectangle)
t5tiplabel$ymax <- cumsum(t5tiplabel$fraction)

# Compute the bottom of each rectangle
t5tiplabel$ymin <- c(0, head(t5tiplabel$ymax, n=-1))

# Compute label position
t5tiplabel$labelPosition <- (t5tiplabel$ymax + t5tiplabel$ymin) / 2

# Compute a good label
t5tiplabel$label <- paste0(t5tiplabel$Strain, "\n value: ", t5tiplabel$count)

# Make the plot
ggplot(t5tiplabel, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Strain)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

t5tiplabel$Expected <- sum(t5tiplabel$count)/4
t5tiplabel$OminusE <- t5tiplabel$count-t5tiplabel$Expected
t5tiplabel$OminusEsquared <- (t5tiplabel$OminusE)^2
t5tiplabel$chiSQ <- t5tiplabel$OminusEsquared/t5tiplabel$Expected
1 - pchisq(q=sum(t5tiplabel$chiSQ), df=3)
#[1] 0.01419381

#t5tiplabel$node <- 4011

#Node6
t6tiplabel %>% 
  group_by(Strain) %>%
  dplyr::summarize(count = n()) %>%
  as.data.frame() -> t6tiplabel

# Compute percentages
t6tiplabel$fraction <- t6tiplabel$count / sum(t6tiplabel$count)

# Compute the cumulative percentages (top of each rectangle)
t6tiplabel$ymax <- cumsum(t6tiplabel$fraction)

# Compute the bottom of each rectangle
t6tiplabel$ymin <- c(0, head(t6tiplabel$ymax, n=-1))

# Compute label position
t6tiplabel$labelPosition <- (t6tiplabel$ymax + t6tiplabel$ymin) / 2

# Compute a good label
t6tiplabel$label <- paste0(t6tiplabel$Strain, "\n value: ", t6tiplabel$count)

# Make the plot
ggplot(t6tiplabel, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Strain)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

t6tiplabel$Expected <- sum(t6tiplabel$count)/4
t6tiplabel$OminusE <- t6tiplabel$count-t6tiplabel$Expected
t6tiplabel$OminusEsquared <- (t6tiplabel$OminusE)^2
t6tiplabel$chiSQ <- t6tiplabel$OminusEsquared/t6tiplabel$Expected
1 - pchisq(q=sum(t6tiplabel$chiSQ), df=3)
#[1] 6.739227e-06
#t6tiplabel$node <- 5822

#Node7
tiplabels_diff %>% 
  group_by(Strain) %>%
  dplyr::summarize(count = n()) %>%
  as.data.frame() -> tiplabels_diff

# Compute percentages
tiplabels_diff$fraction <- tiplabels_diff$count / sum(tiplabels_diff$count)

# Compute the cumulative percentages (top of each rectangle)
tiplabels_diff$ymax <- cumsum(tiplabels_diff$fraction)

# Compute the bottom of each rectangle
tiplabels_diff$ymin <- c(0, head(tiplabels_diff$ymax, n=-1))

# Compute label position
tiplabels_diff$labelPosition <- (tiplabels_diff$ymax + tiplabels_diff$ymin) / 2

# Compute a good label
tiplabels_diff$label <- paste0(tiplabels_diff$Strain, "\n value: ", tiplabels_diff$count)

# Make the plot
ggplot(tiplabels_diff, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Strain)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=6) +
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

tiplabels_diff$Expected <- sum(tiplabels_diff$count)/4
tiplabels_diff$OminusE <- tiplabels_diff$count-tiplabels_diff$Expected
tiplabels_diff$OminusEsquared <- (tiplabels_diff$OminusE)^2
tiplabels_diff$chiSQ <- tiplabels_diff$OminusEsquared/tiplabels_diff$Expected
1 - pchisq(q=sum(tiplabels_diff$chiSQ), df=3)

#[1] 0.006354228
#tiplabels_diff$node <- 4001

###Make a df with all clades and then run chi square test

all_clades <- cbind(as.data.frame(t1tiplabel[,2]),as.data.frame(t2tiplabel[,2]),as.data.frame(t3tiplabel[,2]),as.data.frame(t4tiplabel[,2]),as.data.frame(t5tiplabel[,2]),as.data.frame(t6tiplabel[,2]),as.data.frame(tiplabels_diff[,2]))
colnames(all_clades) <- c("6","5","1","2","7","4","3")
rownames(all_clades) <- c("B6","CAST","LEWES","PWK")
chisq.test(all_clades)
chisq.posthoc.test(all_clades)

