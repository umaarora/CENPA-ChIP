#install.packages("phangorn")
#install.packages("seqinr")

library(ape)
library(phangorn)
library(seqinr)
library(ggplot2)

#Making phylo bar chart with reads ranked by cumulative sum abundance (across CENP-A and input samples)
setwd("~/Box/centromere_ChIP_seq/CENP-A_ChIP/data/Top100NormAbundance_SplitSat")

minor <- read.table(file = "ChIP_Forward_mapped_normalized_top100sum_minor.txt", header = TRUE)
minor_reads <- data.frame(sequence = character(),
                          stringsAsFactors = FALSE)
minor_prep <- data.frame(sequence = character(),
                         stringsAsFactors = FALSE)

for (i in 1:100) {
  minor_prep[1,1] <- paste(">",i, sep = "")
  minor_prep[2,1] <- as.character(minor[i,1])
  minor_reads <- rbind(minor_reads,minor_prep)
}

write.table(minor_reads , file = "reads_top100sum_minor.fa", row.names = FALSE, col.names = FALSE, quote = FALSE)


major <- read.table(file = "ChIP_Forward_mapped_normalized_top100sum_major.txt", header = TRUE)
major_reads <- data.frame(sequence = character(),
                          stringsAsFactors = FALSE)
major_prep <- data.frame(sequence = character(),
                         stringsAsFactors = FALSE)

for (i in 1:100) {
  major_prep[1,1] <- paste(">",i, sep = "")
  major_prep[2,1] <- as.character(major[i,1])
  major_reads <- rbind(major_reads,major_prep)
}

write.table(major_reads , file = "reads_top100sum_major.fa", row.names = FALSE, col.names = FALSE, quote = FALSE)

Ymin <- read.table(file = "ChIP_Forward_mapped_normalized_top100sum_Ymin.txt", header = TRUE)
Ymin_reads <- data.frame(sequence = character(),
                         stringsAsFactors = FALSE)
Ymin_prep <- data.frame(sequence = character(),
                        stringsAsFactors = FALSE)

for (i in 1:100) {
  Ymin_prep[1,1] <- paste(">",i, sep = "")
  Ymin_prep[2,1] <- as.character(Ymin[i,1])
  Ymin_reads <- rbind(Ymin_reads,Ymin_prep)
}

write.table(Ymin_reads , file = "reads_top100sum_Ymin.fa", row.names = FALSE, col.names = FALSE, quote = FALSE)

####read in sequences for phylogenetic tree
minor <- read.dna("reads_top100sum_minor.fa", format = "fasta")
major <- read.dna("reads_top100sum_major.fa", format = "fasta")
Ymin <- read.dna("reads_top100sum_Ymin.fa", format = "fasta")

#Making phylo bar chart with reads ranked by CENP-A/input variance inndependently

setwd("~/Desktop/centromere_ChIP_seq/CENP-A_ChIP/data/Top100NormAbundance_SplitSat/varOrdered")

#Minor
minor <- read.table(file = "ChIP_Forward_mapped_normalized_vartop100_minor.txt", header = TRUE)
minor_reads <- data.frame(sequence = character(),
                          stringsAsFactors = FALSE)
minor_prep <- data.frame(sequence = character(),
                         stringsAsFactors = FALSE)

for (i in 1:100) {
  minor_prep[1,1] <- paste(">",i, sep = "")
  minor_prep[2,1] <- as.character(minor[i,1])
  minor_reads <- rbind(minor_reads,minor_prep)
}

write.table(minor_reads , file = "reads_top100var_minor.fa", row.names = FALSE, col.names = FALSE, quote = FALSE)

#Major
major <- read.table(file = "ChIP_Forward_mapped_normalized_vartop100_major.txt", header = TRUE)
major_reads <- data.frame(sequence = character(),
                          stringsAsFactors = FALSE)
major_prep <- data.frame(sequence = character(),
                         stringsAsFactors = FALSE)

for (i in 1:100) {
  major_prep[1,1] <- paste(">",i, sep = "")
  major_prep[2,1] <- as.character(major[i,1])
  major_reads <- rbind(major_reads,major_prep)
}

write.table(major_reads , file = "reads_top100var_major.fa", row.names = FALSE, col.names = FALSE, quote = FALSE)

#Y
Ymin <- read.table(file = "ChIP_Forward_mapped_normalized_vartop100_Ymin.txt", header = TRUE)
Ymin_reads <- data.frame(sequence = character(),
                         stringsAsFactors = FALSE)
Ymin_prep <- data.frame(sequence = character(),
                        stringsAsFactors = FALSE)

for (i in 1:100) {
  Ymin_prep[1,1] <- paste(">",i, sep = "")
  Ymin_prep[2,1] <- as.character(Ymin[i,1])
  Ymin_reads <- rbind(Ymin_reads,Ymin_prep)
}

write.table(Ymin_reads , file = "reads_top100var_Ymin.fa", row.names = FALSE, col.names = FALSE, quote = FALSE)

#Build tree from fasta with ClustalW (online on clustal omega)
#Read in 
minor_tree <- read.tree("reads_top100var_minor.dnd")
major_tree <- read.tree("reads_top100var_major.dnd")
Ymin_tree <- read.tree("reads_top100var_Ymin.dnd")

plot(minor_tree, type = "cladogram")
plot(major_tree, type = "cladogram")
plot(Ymin_tree, type = "cladogram")

#Set order for bar chart to align with tip label order
minor_order <- minor_tree$tip.label
major_order <- major_tree$tip.label
Ymin_order <- Ymin_tree$tip.label

#Add row number as a column in df 
minor$row <- row.names(minor)
minor$row <- as.factor(minor$row)

major$row <- row.names(major)
major$row <- as.factor(major$row)

Ymin$row <- row.names(Ymin)
Ymin$row <- as.factor(Ymin$row)

#Order rows with phylogram

minor$row <- ordered(minor$row, levels = minor_order)
major$row <- ordered(major$row, levels = major_order)
Ymin$row <- ordered(Ymin$row, levels = Ymin_order)

#Make bar charts

B6min <- ggplot(minor, aes(x = row, y = B6mean)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin=B6mean, ymax=B6mean+B6sd)) +
  theme_bw() + 
  coord_flip()

CASTmin <- ggplot(minor, aes(x = row, y = CASTmean)) +
  geom_bar(stat = "identity", fill = "forestgreen") +
  geom_errorbar(aes(ymin=CASTmean, ymax=CASTmean+CASTsd)) +
  theme_bw() + 
  coord_flip()

PWKmin <- ggplot(minor, aes(x = row, y = PWKmean)) +
  geom_bar(stat = "identity", fill = "red") +
  geom_errorbar(aes(ymin=PWKmean, ymax=PWKmean+PWKsd)) +
  theme_bw() + 
  coord_flip()

LEWESmin <- ggplot(minor, aes(x = row, y = LEWESmean)) +
  geom_bar(stat = "identity", fill = "plum") +
  geom_errorbar(aes(ymin=LEWESmean, ymax=LEWESmean+LEWESsd)) +
  theme_bw() + 
  coord_flip()


## installing the package
install.packages("ggmsa")
## loading the package
library("ggmsa")


#Visualize alignment in a different way with seqinr
library(seqinr)
minor_aln <- read.alignment(file = "minor.clustal", format = "clustal")

#Visualize alignment with msa package
myClustalWAlignment <- msa(minor, "ClustalW")



#########051921
#Making phylo bar chart with reads ranked by cumulative sum abundance (across CENP-A and input samples)
setwd("~/Box/centromere_ChIP_seq/CENP-A_ChIP/data/Top100NormAbundance_SplitSat/FstatsOrdered")

minor <- read.table(file = "minortop100_ChIP_Forward_mapped_normalized_SatSeparated_Fstats.txt", header = TRUE)
minor_reads <- data.frame(sequence = character(),
                          stringsAsFactors = FALSE)
minor_prep <- data.frame(sequence = character(),
                         stringsAsFactors = FALSE)

for (i in 1:100) {
  minor_prep[1,1] <- paste(">",i, sep = "")
  minor_prep[2,1] <- as.character(minor[i,1])
  minor_reads <- rbind(minor_reads,minor_prep)
}

write.table(minor_reads , file = "reads_top100Fstats_minor.fa", row.names = FALSE, col.names = FALSE, quote = FALSE)


major <- read.table(file = "majortop100_ChIP_Forward_mapped_normalized_SatSeparated_Fstats.txt", header = TRUE)
major_reads <- data.frame(sequence = character(),
                          stringsAsFactors = FALSE)
major_prep <- data.frame(sequence = character(),
                         stringsAsFactors = FALSE)

for (i in 1:100) {
  major_prep[1,1] <- paste(">",i, sep = "")
  major_prep[2,1] <- as.character(major[i,1])
  major_reads <- rbind(major_reads,major_prep)
}

write.table(major_reads , file = "reads_top100Fstats_major.fa", row.names = FALSE, col.names = FALSE, quote = FALSE)

#Build tree from fasta with ClustalW (online on clustal omega)
#Read in 
minor_tree <- read.tree("reads_top100Fstats_minor.dnd")
major_tree <- read.tree("reads_top100Fstats_major.dnd")

plot(minor_tree, type = "cladogram")
plot(major_tree, type = "cladogram")

#Set order for bar chart to align with tip label order
minor_order <- minor_tree$tip.label
major_order <- major_tree$tip.label

#Add row number as a column in df 
minor$row <- row.names(minor)
minor$row <- as.factor(minor$row)

major$row <- row.names(major)
major$row <- as.factor(major$row)

#Order rows with phylogram
minor$row <- ordered(minor$row, levels = minor_order)
major$row <- ordered(major$row, levels = major_order)

#Add standard deviation
minor$B6_CENPA_sd <- rowSds(as.matrix(minor[,3:26]), cols = c(3,6))
minor$CAST_CENPA_sd <- rowSds(as.matrix(minor[,3:26]), cols = c(7,10,11))
minor$LEWES_CENPA_sd <- rowSds(as.matrix(minor[,3:26]), cols = c(12,13,16))
minor$PWK_CENPA_sd <- rowSds(as.matrix(minor[,3:26]), cols = c(17,18,21))

#Make bar charts

B6min <- ggplot(minor, aes(x = row, y = B6_CENPA_mean)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin=B6_CENPA_mean, ymax=B6_CENPA_mean+B6_CENPA_sd)) +
  theme_bw() + 
  ylim(0,200) +
  coord_flip() 

CASTmin <- ggplot(minor, aes(x = row, y = CAST_CENPA_mean)) +
  geom_bar(stat = "identity", fill = "forestgreen") +
  geom_errorbar(aes(ymin=CAST_CENPA_mean, ymax=CAST_CENPA_mean+CAST_CENPA_sd)) +
  theme_bw() + 
  ylim(0,200) +
  coord_flip()

LEWESmin <- ggplot(minor, aes(x = row, y = LEWES_CENPA_mean)) +
  geom_bar(stat = "identity", fill = "plum") +
  geom_errorbar(aes(ymin=LEWES_CENPA_mean, ymax=LEWES_CENPA_mean+LEWES_CENPA_sd)) +
  theme_bw() + 
  ylim(0,200) +
  coord_flip()

PWKmin <- ggplot(minor, aes(x = row, y = PWK_CENPA_mean)) +
  geom_bar(stat = "identity", fill = "red") +
  geom_errorbar(aes(ymin=PWK_CENPA_mean, ymax=PWK_CENPA_mean+PWK_CENPA_sd)) +
  theme_bw() + 
  ylim(0,200) +
  coord_flip()


