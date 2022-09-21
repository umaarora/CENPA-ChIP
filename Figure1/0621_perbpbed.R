library(tidyverse)
library(matrixStats)
library(corrplot)

setwd("~/Box/centromere_ChIP_seq/CENP-A_ChIP/data/0621_perbpbed")

##########set path for where files are located###############
path = "~/Box/centromere_ChIP_seq/CENP-A_ChIP/data/0621_perbpbed"
fs = list.files(path, pattern = glob2rx("*.bed")) 

##########read in all files###############

for(f in fs){      
  fname = file.path(path,f)
  n = gsub("", "", f)
  df = read.table(fname, header = FALSE)
  assign(n,df)
  rm(df)
}

##Our data
#Merge dataframes by satellite consensus mapped to and position

all <- full_join(fastp_B6_F_092920_CENPA_sorted.bam.bed_perbp.bed,
                 fastp_B6_F_092920_input_sorted.bam.bed_perbp.bed, 
                 by = c("V1","V2"))

all <- full_join(all,
                 fastp_B6_M_CENPA_sorted.bam.bed_perbp.bed, 
                 by = c("V1","V2"))

all <- full_join(all,
                 fastp_B6_M_input_sorted.bam.bed_perbp.bed, 
                 by = c("V1","V2"))

all <- full_join(all,
                 fastp_CAST_F_011921_CENPA_sorted.bam.bed_perbp.bed, 
                 by = c("V1","V2"))

all <- full_join(all,
                 fastp_CAST_F_011921_input_sorted.bam.bed_perbp.bed, 
                 by = c("V1","V2"))

all <- full_join(all,
                 fastp_CAST_M_011921_CENPA_sorted.bam.bed_perbp.bed, 
                 by = c("V1","V2"))

all <- full_join(all,
                 fastp_CAST_M_011921_input_sorted.bam.bed_perbp.bed, 
                 by = c("V1","V2"))

all <- full_join(all,
                 fastp_CAST_M_101620_CENPA_sorted.bam.bed_perbp.bed, 
                 by = c("V1","V2"))

all <- full_join(all,
                 fastp_CAST_M_101620_input_sorted.bam.bed_perbp.bed, 
                 by = c("V1","V2"))

all <- full_join(all,
                 fastp_LEWES_F_011921_CENPA_sorted.bam.bed_perbp.bed, 
                 by = c("V1","V2"))

all <- full_join(all,
                 fastp_LEWES_F_011921_input_sorted.bam.bed_perbp.bed, 
                 by = c("V1","V2"))

all <- full_join(all,
                 fastp_LEWES_F_012621_CENPA_sorted.bam.bed_perbp.bed, 
                 by = c("V1","V2"))

all <- full_join(all,
                 fastp_LEWES_F_012621_input_sorted.bam.bed_perbp.bed, 
                 by = c("V1","V2"))

all <- full_join(all,
                 fastp_LEWES_M_011921_CENPA_sorted.bam.bed_perbp.bed, 
                 by = c("V1","V2"))

all <- full_join(all,
                 fastp_LEWES_M_011921_input_sorted.bam.bed_perbp.bed, 
                 by = c("V1","V2"))

all <- full_join(all,
                 fastp_PWK_F_012621_CENPA_sorted.bam.bed_perbp.bed, 
                 by = c("V1","V2"))

all <- full_join(all,
                 fastp_PWK_F_012621_input_sorted.bam.bed_perbp.bed, 
                 by = c("V1","V2"))

all <- full_join(all,
                 fastp_PWK_F_092920_CENPA_sorted.bam.bed_perbp.bed, 
                 by = c("V1","V2"))

all <- full_join(all,
                 fastp_PWK_F_092920_input_sorted.bam.bed_perbp.bed, 
                 by = c("V1","V2"))

all <- full_join(all,
                 fastp_PWK_M_101620_CENPA_sorted.bam.bed_perbp.bed, 
                 by = c("V1","V2"))

all <- full_join(all,
                 fastp_PWK_M_101620_input_sorted.bam.bed_perbp.bed, 
                 by = c("V1","V2"))


#change column names
all_columns <- c("Satellite Consensus",
                 "Position",
                 "B6_F_092920_CENPA",
                 "B6_F_092920_input",
                 "B6_M_CENPA", 
                 "B6_M_input",
                 "CAST_F_011921_CENPA",
                 "CAST_F_011921_input",
                 "CAST_M_011921_CENPA",
                 "CAST_M_011921_input",
                 "CAST_M_101620_CENPA",
                 "CAST_M_101620_input",
                 "LEWES_F_011921_CENPA",
                 "LEWES_F_011921_input",
                 "LEWES_F_012621_CENPA",
                 "LEWES_F_012621_input",
                 "LEWES_M_011921_CENPA",
                 "LEWES_M_011921_input",
                 "PWK_F_012621_CENPA",
                 "PWK_F_012621_input",
                 "PWK_F_092920_CENPA",
                 "PWK_F_092920_input",
                 "PWK_M_101620_CENPA",
                 "PWK_M_101620_input")

colnames(all) <- all_columns

##Separate by consensus sequence
Major <- all[all$`Satellite Consensus` == 'MajorSatelliteConsensusWongandRattner1988x3',]
Minor <- all[all$`Satellite Consensus` == 'MinorSatelliteConsensusWongandRattner1988x3',]
Y <- all[all$`Satellite Consensus` == 'eco12_Ym841',]

##Normalize columns by number of reads that mapped to the satellite consensus
#Major
Major$B6_F_092920_CENPA <- (Major$B6_F_092920_CENPA / 1148517)
Major$B6_F_092920_input <- (Major$B6_F_092920_input / 2470606)
Major$B6_M_CENPA <- (Major$B6_M_CENPA/897924)
Major$B6_M_input <- (Major$B6_M_input/6758158)
Major$CAST_F_011921_CENPA <- (Major$CAST_F_011921_CENPA/3730742)
Major$CAST_F_011921_input <- (Major$CAST_F_011921_input/3322422)
Major$CAST_M_011921_CENPA <- (Major$CAST_M_011921_CENPA/794096)
Major$CAST_M_011921_input <- (Major$CAST_M_011921_input/5499293)
Major$CAST_M_101620_CENPA <- (Major$CAST_M_101620_CENPA/4892367)
Major$CAST_M_101620_input <- (Major$CAST_M_101620_input/7092234)
Major$LEWES_F_011921_CENPA <- (Major$LEWES_F_011921_CENPA/2869564)
Major$LEWES_F_011921_input <- (Major$LEWES_F_011921_input/4871780)
Major$LEWES_F_012621_CENPA <- (Major$LEWES_F_012621_CENPA/3475430)
Major$LEWES_F_012621_input <- (Major$LEWES_F_012621_input/2063226)
Major$LEWES_M_011921_CENPA <- (Major$LEWES_M_011921_CENPA/2680739)
Major$LEWES_M_011921_input <- (Major$LEWES_M_011921_input/4025499)
Major$PWK_F_012621_CENPA <- (Major$PWK_F_012621_CENPA/1161893)
Major$PWK_F_012621_input <- (Major$PWK_F_012621_input/3570179)
Major$PWK_F_092920_CENPA <- (Major$PWK_F_092920_CENPA/378688)
Major$PWK_F_092920_input <- (Major$PWK_F_092920_input/1877087)
Major$PWK_M_101620_CENPA <- (Major$PWK_M_101620_CENPA/857868)
Major$PWK_M_101620_input <- (Major$PWK_M_101620_input/4789314)

##ChIP over input ratios
#B6
Major$B6_F_092920_CENPAoverInput <- Major$B6_F_092920_CENPA/Major$B6_F_092920_input
Major$B6_M_CENPAoverInput <- Major$B6_M_CENPA/Major$B6_M_input
Major$B6mean <- rowMeans(Major[,25:26])
Major$B6sd <- rowSds(as.matrix(Major[,3:27]), cols = c(23,24))

#CAST
Major$CAST_F_011921_CENPAoverInput <- Major$CAST_F_011921_CENPA/Major$CAST_F_011921_input
Major$CAST_M_011921_CENPAoverInput <- Major$CAST_M_011921_CENPA/Major$CAST_M_011921_input
Major$CAST_M_101620_CENPAoverInput <- Major$CAST_M_101620_CENPA/Major$CAST_M_101620_input
Major$CASTmean <- rowMeans(Major[,29:31])
Major$CASTsd <- rowSds(as.matrix(Major[,3:32]), cols = c(27,28,29))

#LEWES
Major$LEWES_F_011921_CENPAoverInput <- Major$LEWES_F_011921_CENPA/Major$LEWES_F_011921_input
Major$LEWES_F_012621_CENPAoverInput <- Major$LEWES_F_012621_CENPA/Major$LEWES_F_012621_input
Major$LEWES_M_011921_CENPAoverInput <- Major$LEWES_M_011921_CENPA/Major$LEWES_M_011921_input
Major$LEWESmean <- rowMeans(Major[,34:36])
Major$LEWESsd <- rowSds(as.matrix(Major[,3:37]), cols = c(32,33,34))

#PWK
Major$PWK_F_012621_CENPAoverInput <- Major$PWK_F_012621_CENPA/Major$PWK_F_012621_input
Major$PWK_F_092920_CENPAoverInput <- Major$PWK_F_092920_CENPA/Major$PWK_F_092920_input
Major$PWK_M_101620_CENPAoverInput <- Major$PWK_M_101620_CENPA/Major$PWK_M_101620_input
Major$PWKmean <- rowMeans(Major[,39:41])
Major$PWKsd <- rowSds(as.matrix(Major[,3:42]), cols = c(37,38,39))

##Df for ggplot

Major_B6 <- Major[,c(1,2,27,28)]
colnames(Major_B6) <- c("Satellite Consensus","Position","Mean","SD")
Major_B6$Strain <- c("C57BL/6J")

Major_CAST <- Major[,c(1,2,32,33)]
colnames(Major_CAST) <- c("Satellite Consensus","Position","Mean","SD")
Major_CAST$Strain <- c("CAST/EiJ")

Major_LEWES <- Major[,c(1,2,37,38)]
colnames(Major_LEWES) <- c("Satellite Consensus","Position","Mean","SD")
Major_LEWES$Strain <- c("LEWES/EiJ")

Major_PWK <- Major[,c(1,2,42,43)]
colnames(Major_PWK) <- c("Satellite Consensus","Position","Mean","SD")
Major_PWK$Strain <- c("PWK/PhJ")

Major_ggplot <- rbind(Major_B6,Major_CAST,Major_LEWES,Major_PWK)

ggplot(Major_ggplot, aes(x = Position, y = Mean, group = Strain)) +
  geom_line(aes(color = Strain)) +
  geom_ribbon(aes(y = Mean, ymin = Mean - SD, ymax = Mean + SD, fill = Strain), alpha = .1) +
  scale_fill_manual(values = c("plum", "forestgreen", "purple","red")) +
  scale_color_manual(values = c("plum", "forestgreen", "purple","red")) +
  theme_classic()

#Minor

Minor$B6_F_092920_CENPA <- (Minor$B6_F_092920_CENPA / 9603511)
Minor$B6_F_092920_input <- (Minor$B6_F_092920_input / 194174)
Minor$B6_M_CENPA <- (Minor$B6_M_CENPA/16231650)
Minor$B6_M_input <- (Minor$B6_M_input/383006)
Minor$CAST_F_011921_CENPA <- (Minor$CAST_F_011921_CENPA/7939934)
Minor$CAST_F_011921_input <- (Minor$CAST_F_011921_input/205561)
Minor$CAST_M_011921_CENPA <- (Minor$CAST_M_011921_CENPA/7644849)
Minor$CAST_M_011921_input <- (Minor$CAST_M_011921_input/355884)
Minor$CAST_M_101620_CENPA <- (Minor$CAST_M_101620_CENPA/4751560)
Minor$CAST_M_101620_input <- (Minor$CAST_M_101620_input/449453)
Minor$LEWES_F_011921_CENPA <- (Minor$LEWES_F_011921_CENPA/6467193)
Minor$LEWES_F_011921_input <- (Minor$LEWES_F_011921_input/469558)
Minor$LEWES_F_012621_CENPA <- (Minor$LEWES_F_012621_CENPA/7201272)
Minor$LEWES_F_012621_input <- (Minor$LEWES_F_012621_input/230072)
Minor$LEWES_M_011921_CENPA <- (Minor$LEWES_M_011921_CENPA/15839180)
Minor$LEWES_M_011921_input <- (Minor$LEWES_M_011921_input/366551)
Minor$PWK_F_012621_CENPA <- (Minor$PWK_F_012621_CENPA/10144218)
Minor$PWK_F_012621_input <- (Minor$PWK_F_012621_input/1124228)
Minor$PWK_F_092920_CENPA <- (Minor$PWK_F_092920_CENPA/4860899)
Minor$PWK_F_092920_input <- (Minor$PWK_F_092920_input/588714)
Minor$PWK_M_101620_CENPA <- (Minor$PWK_M_101620_CENPA/8984407)
Minor$PWK_M_101620_input <- (Minor$PWK_M_101620_input/1063132)

##ChIP over input ratios
#B6
Minor$B6_F_092920_CENPAoverInput <- Minor$B6_F_092920_CENPA/Minor$B6_F_092920_input
Minor$B6_M_CENPAoverInput <- Minor$B6_M_CENPA/Minor$B6_M_input
Minor$B6mean <- rowMeans(Minor[,25:26])
Minor$B6sd <- rowSds(as.matrix(Minor[,3:27]), cols = c(23,24))

#CAST
Minor$CAST_F_011921_CENPAoverInput <- Minor$CAST_F_011921_CENPA/Minor$CAST_F_011921_input
Minor$CAST_M_011921_CENPAoverInput <- Minor$CAST_M_011921_CENPA/Minor$CAST_M_011921_input
Minor$CAST_M_101620_CENPAoverInput <- Minor$CAST_M_101620_CENPA/Minor$CAST_M_101620_input
Minor$CASTmean <- rowMeans(Minor[,29:31])
Minor$CASTsd <- rowSds(as.matrix(Minor[,3:32]), cols = c(27,28,29))

#LEWES
Minor$LEWES_F_011921_CENPAoverInput <- Minor$LEWES_F_011921_CENPA/Minor$LEWES_F_011921_input
Minor$LEWES_F_012621_CENPAoverInput <- Minor$LEWES_F_012621_CENPA/Minor$LEWES_F_012621_input
Minor$LEWES_M_011921_CENPAoverInput <- Minor$LEWES_M_011921_CENPA/Minor$LEWES_M_011921_input
Minor$LEWESmean <- rowMeans(Minor[,34:36])
Minor$LEWESsd <- rowSds(as.matrix(Minor[,3:37]), cols = c(32,33,34))

#PWK
Minor$PWK_F_012621_CENPAoverInput <- Minor$PWK_F_012621_CENPA/Minor$PWK_F_012621_input
Minor$PWK_F_092920_CENPAoverInput <- Minor$PWK_F_092920_CENPA/Minor$PWK_F_092920_input
Minor$PWK_M_101620_CENPAoverInput <- Minor$PWK_M_101620_CENPA/Minor$PWK_M_101620_input
Minor$PWKmean <- rowMeans(Minor[,39:41])
Minor$PWKsd <- rowSds(as.matrix(Minor[,3:42]), cols = c(37,38,39))

##Df for ggplot

Minor_B6 <- Minor[,c(1,2,27,28)]
colnames(Minor_B6) <- c("Satellite Consensus","Position","Mean","SD")
Minor_B6$Strain <- c("C57BL/6J")

Minor_CAST <- Minor[,c(1,2,32,33)]
colnames(Minor_CAST) <- c("Satellite Consensus","Position","Mean","SD")
Minor_CAST$Strain <- c("CAST/EiJ")

Minor_LEWES <- Minor[,c(1,2,37,38)]
colnames(Minor_LEWES) <- c("Satellite Consensus","Position","Mean","SD")
Minor_LEWES$Strain <- c("LEWES/EiJ")

Minor_PWK <- Minor[,c(1,2,42,43)]
colnames(Minor_PWK) <- c("Satellite Consensus","Position","Mean","SD")
Minor_PWK$Strain <- c("PWK/PhJ")

Minor_ggplot <- rbind(Minor_B6,Minor_CAST,Minor_LEWES,Minor_PWK)

ggplot(Minor_ggplot, aes(x = Position, y = Mean, group = Strain)) +
  geom_line(aes(color = Strain)) +
  geom_ribbon(aes(y = Mean, ymin = Mean - SD, ymax = Mean + SD, fill = Strain), alpha = .1) +
  scale_fill_manual(values = c("plum", "forestgreen", "purple","red")) +
  scale_color_manual(values = c("plum", "forestgreen", "purple","red")) +
  theme_classic()


##Iwata Otsubo B6 data vs our B6 data
all <- full_join(fastp_B6_F_092920_CENPA_sorted.bam.bed_perbp.bed,
                 fastp_B6_F_092920_input_sorted.bam.bed_perbp.bed, 
                 by = c("V1","V2"))

all <- full_join(all,
                 fastp_B6_M_CENPA_sorted.bam.bed_perbp.bed, 
                 by = c("V1","V2"))

all <- full_join(all,
                 fastp_B6_M_input_sorted.bam.bed_perbp.bed, 
                 by = c("V1","V2"))

all <- full_join(all,
                 `fastp_BL6-1-CENPA_sorted.bam.bed_perbp.bed`, 
                 by = c("V1","V2"))

all <- full_join(all,
                 `fastp_BL6-1-INP_sorted.bam.bed_perbp.bed`, 
                 by = c("V1","V2"))

all <- full_join(all,
                 `fastp_BL6-2-CENPA_sorted.bam.bed_perbp.bed`, 
                 by = c("V1","V2"))

all <- full_join(all,
                 `fastp_BL6-2-INP_sorted.bam.bed_perbp.bed`, 
                 by = c("V1","V2"))

all <- full_join(all,
                 `fastp_BL6-3-CENPA_sorted.bam.bed_perbp.bed`, 
                 by = c("V1","V2"))

all <- full_join(all,
                 `fastp_BL6-3-INP_sorted.bam.bed_perbp.bed`, 
                 by = c("V1","V2"))

all_columns <- c("Satellite Consensus",
                 "Position",
                 "B6_F_092920_CENPA",
                 "B6_F_092920_input",
                 "B6_M_CENPA", 
                 "B6_M_input",
                 "BL6-1-CENPA",
                 "BL6-1-INP",
                 "BL6-2-CENPA",
                 "BL6-2-INP",
                 "BL6-3-CENPA",
                 "BL6-3-INP")

colnames(all) <- all_columns

##Separate by consensus sequence
Major <- all[all$`Satellite Consensus` == 'MajorSatelliteConsensusWongandRattner1988x3',]
Minor <- all[all$`Satellite Consensus` == 'MinorSatelliteConsensusWongandRattner1988x3',]
Y <- all[all$`Satellite Consensus` == 'eco12_Ym841',]

##Normalize columns by number of reads that mapped to the satellite consensus
#Major
Major$B6_F_092920_CENPA <- (Major$B6_F_092920_CENPA / 1148517)
Major$B6_F_092920_input <- (Major$B6_F_092920_input / 2470606)
Major$B6_M_CENPA <- (Major$B6_M_CENPA/897924)
Major$B6_M_input <- (Major$B6_M_input/6758158)
Major$`BL6-1-CENPA` <- (Major$`BL6-1-CENPA`/1551991)
Major$`BL6-1-INP` <- (Major$`BL6-1-INP`/1617881)
Major$`BL6-2-CENPA` <- (Major$`BL6-2-CENPA`/124810)
Major$`BL6-2-INP` <- (Major$`BL6-2-INP`/1526742)
Major$`BL6-3-CENPA` <- (Major$`BL6-3-CENPA`/88788)
Major$`BL6-3-INP` <- (Major$`BL6-3-INP`/1408327)

##ChIP over input ratios
#B6
Major$B6_F_092920_CENPAoverInput <- Major$B6_F_092920_CENPA/Major$B6_F_092920_input
Major$B6_M_CENPAoverInput <- Major$B6_M_CENPA/Major$B6_M_input
Major$B6mean <- rowMeans(Major[,13:14])
Major$B6sd <- rowSds(as.matrix(Major[,3:15]), cols = c(11,12))

#B6 Iwata Otsubo
Major$BL61_CENPAoverInput <- Major$`BL6-1-CENPA`/Major$`BL6-1-INP`
Major$BL62_CENPAoverInput <- Major$`BL6-2-CENPA`/Major$`BL6-2-INP`
Major$BL63_CENPAoverInput <- Major$`BL6-3-CENPA`/Major$`BL6-3-INP`
Major$B6IwataMean <- rowMeans(Major[,17:19])
Major$B6Iwatasd <- rowSds(as.matrix(Major[,3:20]), cols = c(15,16,17))

##Df for ggplot

Major_B6 <- Major[,c(1,2,15,16)]
colnames(Major_B6) <- c("Satellite Consensus","Position","Mean","SD")
Major_B6$Strain <- c("C57BL/6J")

Major_B6Iwata <- Major[,c(1,2,20,21)]
colnames(Major_B6Iwata) <- c("Satellite Consensus","Position","Mean","SD")
Major_B6Iwata$Strain <- c("C57BL/6J Iwata")

Major_ggplot <- rbind(Major_B6,Major_B6Iwata)

ggplot(Major_ggplot, aes(x = Position, y = Mean, group = Strain)) +
  geom_line(aes(color = Strain)) +
  geom_ribbon(aes(y = Mean, ymin = Mean - SD, ymax = Mean + SD, fill = Strain), alpha = .1) +
  scale_fill_manual(values = c("plum", "grey")) +
  scale_color_manual(values = c("plum", "grey")) +
  theme_classic()

#Minor
Minor$B6_F_092920_CENPA <- (Minor$B6_F_092920_CENPA / 9603511)
Minor$B6_F_092920_input <- (Minor$B6_F_092920_input / 194174)
Minor$B6_M_CENPA <- (Minor$B6_M_CENPA/16231650)
Minor$B6_M_input <- (Minor$B6_M_input/383006)
Minor$`BL6-1-CENPA` <- (Minor$`BL6-1-CENPA`/17107374)
Minor$`BL6-1-INP` <- (Minor$`BL6-1-INP`/135300)
Minor$`BL6-2-CENPA` <- (Minor$`BL6-2-CENPA`/784431)
Minor$`BL6-2-INP` <- (Minor$`BL6-2-INP`/88039)
Minor$`BL6-3-CENPA` <- (Minor$`BL6-3-CENPA`/742549)
Minor$`BL6-3-INP` <- (Minor$`BL6-3-INP`/74505)

##ChIP over input ratios
#B6
Minor$B6_F_092920_CENPAoverInput <- Minor$B6_F_092920_CENPA/Minor$B6_F_092920_input
Minor$B6_M_CENPAoverInput <- Minor$B6_M_CENPA/Minor$B6_M_input
Minor$B6mean <- rowMeans(Minor[,13:14])
Minor$B6sd <- rowSds(as.matrix(Minor[,3:15]), cols = c(11,12))

#B6 Iwata Otsubo
Minor$BL61_CENPAoverInput <- Minor$`BL6-1-CENPA`/Minor$`BL6-1-INP`
Minor$BL62_CENPAoverInput <- Minor$`BL6-2-CENPA`/Minor$`BL6-2-INP`
Minor$BL63_CENPAoverInput <- Minor$`BL6-3-CENPA`/Minor$`BL6-3-INP`
Minor$B6IwataMean <- rowMeans(Minor[,17:19])
Minor$B6Iwatasd <- rowSds(as.matrix(Minor[,3:20]), cols = c(15,16,17))

##Df for ggplot

Minor_B6 <- Minor[,c(1,2,15,16)]
colnames(Minor_B6) <- c("Satellite Consensus","Position","Mean","SD")
Minor_B6$Strain <- c("C57BL/6J")

Minor_B6Iwata <- Minor[,c(1,2,20,21)]
colnames(Minor_B6Iwata) <- c("Satellite Consensus","Position","Mean","SD")
Minor_B6Iwata$Strain <- c("C57BL/6J Iwata")

Minor_ggplot <- rbind(Minor_B6,Minor_B6Iwata)

ggplot(Minor_ggplot, aes(x = Position, y = Mean, group = Strain)) +
  geom_line(aes(color = Strain)) +
  geom_ribbon(aes(y = Mean, ymin = Mean - SD, ymax = Mean + SD, fill = Strain), alpha = .1) +
  scale_fill_manual(values = c("plum", "grey")) +
  scale_color_manual(values = c("plum", "grey")) +
  theme_classic()



##Y sequence #Our data
##Normalize columns by number of reads that mapped to the satellite consensus
#Y

Y$B6_M_CENPA <- (Y$B6_M_CENPA/171300)
Y$B6_M_input <- (Y$B6_M_input/374)
Y$CAST_M_011921_CENPA <- (Y$CAST_M_011921_CENPA/97459)
Y$CAST_M_011921_input <- (Y$CAST_M_011921_input/1683)
Y$CAST_M_101620_CENPA <- (Y$CAST_M_101620_CENPA/50900)
Y$CAST_M_101620_input <- (Y$CAST_M_101620_input/2294)
Y$LEWES_M_011921_CENPA <- (Y$LEWES_M_011921_CENPA/233458)
Y$LEWES_M_011921_input <- (Y$LEWES_M_011921_input/6886)
Y$PWK_M_101620_CENPA <- (Y$PWK_M_101620_CENPA/73449)
Y$PWK_M_101620_input <- (Y$PWK_M_101620_input/1957)

##ChIP over input ratios
#B6
Y$B6_M_CENPAoverInput <- Y$B6_M_CENPA/Y$B6_M_input

#CAST
Y$CAST_M_011921_CENPAoverInput <- Y$CAST_M_011921_CENPA/Y$CAST_M_011921_input
Y$CAST_M_101620_CENPAoverInput <- Y$CAST_M_101620_CENPA/Y$CAST_M_101620_input
Y$CASTmean <- rowMeans(Y[,26:27])
Y$CASTsd <- rowSds(as.matrix(Y[,3:27]), cols = c(24,25))

#LEWES
Y$LEWES_M_011921_CENPAoverInput <- Y$LEWES_M_011921_CENPA/Y$LEWES_M_011921_input

#PWK
Y$PWK_M_101620_CENPAoverInput <- Y$PWK_M_101620_CENPA/Y$PWK_M_101620_input


#Df for ggplot

Y_B6 <- Y[,c(1,2,25)]
colnames(Y_B6) <- c("Satellite Consensus","Position","Mean")
Y_B6$Strain <- c("C57BL/6J")

Y_CAST <- Y[,c(1,2,28)]
colnames(Y_CAST) <- c("Satellite Consensus","Position","Mean")
Y_CAST$Strain <- c("CAST/EiJ")

Y_LEWES <- Y[,c(1,2,30)]
colnames(Y_LEWES) <- c("Satellite Consensus","Position","Mean")
Y_LEWES$Strain <- c("LEWES/EiJ")

Y_PWK <- Y[,c(1,2,31)]
colnames(Y_PWK) <- c("Satellite Consensus","Position","Mean")
Y_PWK$Strain <- c("PWK/PhJ")

Y_ggplot <- rbind(Y_B6,Y_CAST,Y_LEWES,Y_PWK)

ggplot(Y_ggplot, aes(x = Position, y = Mean, group = Strain)) +
  geom_line(aes(color = Strain)) +
  scale_fill_manual(values = c("plum", "forestgreen", "purple","red")) +
  scale_color_manual(values = c("plum", "forestgreen", "purple","red")) +
  theme_classic()

##########################################################################################
#################081921 can we consolidate to one repeat?##############################
##########################################################################################

setwd("~/Box/centromere_ChIP_seq/CENP-A_ChIP/data/0621_perbpbed")

##########set path for where files are located###############
path = "~/Box/centromere_ChIP_seq/CENP-A_ChIP/data/0621_perbpbed"
fs = list.files(path, pattern = glob2rx("*.bed")) 

##########read in all files###############

for(f in fs){      
  fname = file.path(path,f)
  n = gsub("", "", f)
  df = read.table(fname, header = FALSE)
  
  dfmajor <- df[df$V1 == 'MajorSatelliteConsensusWongandRattner1988x3', ]
  dfnewmajor <- dfmajor[1:234,3] + dfmajor[235:468,3] + dfmajor[469:702,3]
  dfnewmajor <- as.data.frame(dfnewmajor)
  nmajor <- paste(n,"_major",sep = "")
  assign(nmajor,dfnewmajor)
  
  dfminor <- df[df$V1 == 'MinorSatelliteConsensusWongandRattner1988x3', ]
  dfnewminor <- dfminor[1:120,3] + dfminor[121:240,3] + dfminor[241:360,3]
  dfnewminor <- as.data.frame(dfnewminor)
  nminor <- paste(n,"_minor",sep = "")
  assign(nminor,dfnewminor)

  rm(df)
  rm(dfmajor)
  rm(dfnewmajor)
  rm(dfminor)
  rm(dfnewminor)
}


#######MAJOR SATELLITE
Major <- cbind(fastp_B6_F_092920_CENPA_sorted.bam.bed_perbp.bed_major,
               fastp_B6_F_092920_input_sorted.bam.bed_perbp.bed_major,
               fastp_B6_M_CENPA_sorted.bam.bed_perbp.bed_major,
               fastp_B6_M_input_sorted.bam.bed_perbp.bed_major,
               fastp_CAST_F_011921_CENPA_sorted.bam.bed_perbp.bed_major,
               fastp_CAST_F_011921_input_sorted.bam.bed_perbp.bed_major,
               fastp_CAST_M_011921_CENPA_sorted.bam.bed_perbp.bed_major,
               fastp_CAST_M_011921_input_sorted.bam.bed_perbp.bed_major,
               fastp_CAST_M_101620_CENPA_sorted.bam.bed_perbp.bed_major,
               fastp_CAST_M_101620_input_sorted.bam.bed_perbp.bed_major,
               fastp_LEWES_F_011921_CENPA_sorted.bam.bed_perbp.bed_major,
               fastp_LEWES_F_011921_input_sorted.bam.bed_perbp.bed_major,
               fastp_LEWES_F_012621_CENPA_sorted.bam.bed_perbp.bed_major,
               fastp_LEWES_F_012621_input_sorted.bam.bed_perbp.bed_major,
               fastp_LEWES_M_011921_CENPA_sorted.bam.bed_perbp.bed_major,
               fastp_LEWES_M_011921_input_sorted.bam.bed_perbp.bed_major,
               fastp_PWK_F_012621_CENPA_sorted.bam.bed_perbp.bed_major,
               fastp_PWK_F_012621_input_sorted.bam.bed_perbp.bed_major,
               fastp_PWK_F_092920_CENPA_sorted.bam.bed_perbp.bed_major,
               fastp_PWK_F_092920_input_sorted.bam.bed_perbp.bed_major,
               fastp_PWK_M_101620_CENPA_sorted.bam.bed_perbp.bed_major,
               fastp_PWK_M_101620_input_sorted.bam.bed_perbp.bed_major)

###Change column names

all_columns <- c("B6_F_092920_CENPA",
                 "B6_F_092920_input",
                 "B6_M_CENPA", 
                 "B6_M_input",
                 "CAST_F_011921_CENPA",
                 "CAST_F_011921_input",
                 "CAST_M_011921_CENPA",
                 "CAST_M_011921_input",
                 "CAST_M_101620_CENPA",
                 "CAST_M_101620_input",
                 "LEWES_F_011921_CENPA",
                 "LEWES_F_011921_input",
                 "LEWES_F_012621_CENPA",
                 "LEWES_F_012621_input",
                 "LEWES_M_011921_CENPA",
                 "LEWES_M_011921_input",
                 "PWK_F_012621_CENPA",
                 "PWK_F_012621_input",
                 "PWK_F_092920_CENPA",
                 "PWK_F_092920_input",
                 "PWK_M_101620_CENPA",
                 "PWK_M_101620_input")


colnames(Major) <- all_columns

#Major
Major$B6_F_092920_CENPA <- (Major$B6_F_092920_CENPA / 1148517)
Major$B6_F_092920_input <- (Major$B6_F_092920_input / 2470606)
Major$B6_M_CENPA <- (Major$B6_M_CENPA/897924)
Major$B6_M_input <- (Major$B6_M_input/6758158)
Major$CAST_F_011921_CENPA <- (Major$CAST_F_011921_CENPA/3730742)
Major$CAST_F_011921_input <- (Major$CAST_F_011921_input/3322422)
Major$CAST_M_011921_CENPA <- (Major$CAST_M_011921_CENPA/794096)
Major$CAST_M_011921_input <- (Major$CAST_M_011921_input/5499293)
Major$CAST_M_101620_CENPA <- (Major$CAST_M_101620_CENPA/4892367)
Major$CAST_M_101620_input <- (Major$CAST_M_101620_input/7092234)
Major$LEWES_F_011921_CENPA <- (Major$LEWES_F_011921_CENPA/2869564)
Major$LEWES_F_011921_input <- (Major$LEWES_F_011921_input/4871780)
Major$LEWES_F_012621_CENPA <- (Major$LEWES_F_012621_CENPA/3475430)
Major$LEWES_F_012621_input <- (Major$LEWES_F_012621_input/2063226)
Major$LEWES_M_011921_CENPA <- (Major$LEWES_M_011921_CENPA/2680739)
Major$LEWES_M_011921_input <- (Major$LEWES_M_011921_input/4025499)
Major$PWK_F_012621_CENPA <- (Major$PWK_F_012621_CENPA/1161893)
Major$PWK_F_012621_input <- (Major$PWK_F_012621_input/3570179)
Major$PWK_F_092920_CENPA <- (Major$PWK_F_092920_CENPA/378688)
Major$PWK_F_092920_input <- (Major$PWK_F_092920_input/1877087)
Major$PWK_M_101620_CENPA <- (Major$PWK_M_101620_CENPA/857868)
Major$PWK_M_101620_input <- (Major$PWK_M_101620_input/4789314)

##ChIP over input ratios
#B6
Major$B6_F_092920_CENPAoverInput <- Major$B6_F_092920_CENPA/Major$B6_F_092920_input
Major$B6_M_CENPAoverInput <- Major$B6_M_CENPA/Major$B6_M_input
Major$B6mean <- rowMeans(Major[,23:24])
Major$B6sd <- rowSds(as.matrix(Major), cols = c(23,24))

#CAST
Major$CAST_F_011921_CENPAoverInput <- Major$CAST_F_011921_CENPA/Major$CAST_F_011921_input
Major$CAST_M_011921_CENPAoverInput <- Major$CAST_M_011921_CENPA/Major$CAST_M_011921_input
Major$CAST_M_101620_CENPAoverInput <- Major$CAST_M_101620_CENPA/Major$CAST_M_101620_input
Major$CASTmean <- rowMeans(Major[,27:29])
Major$CASTsd <- rowSds(as.matrix(Major), cols = c(27,28,29))

#LEWES
Major$LEWES_F_011921_CENPAoverInput <- Major$LEWES_F_011921_CENPA/Major$LEWES_F_011921_input
Major$LEWES_F_012621_CENPAoverInput <- Major$LEWES_F_012621_CENPA/Major$LEWES_F_012621_input
Major$LEWES_M_011921_CENPAoverInput <- Major$LEWES_M_011921_CENPA/Major$LEWES_M_011921_input
Major$LEWESmean <- rowMeans(Major[,32:34])
Major$LEWESsd <- rowSds(as.matrix(Major), cols = c(32,33,34))

#PWK
Major$PWK_F_012621_CENPAoverInput <- Major$PWK_F_012621_CENPA/Major$PWK_F_012621_input
Major$PWK_F_092920_CENPAoverInput <- Major$PWK_F_092920_CENPA/Major$PWK_F_092920_input
Major$PWK_M_101620_CENPAoverInput <- Major$PWK_M_101620_CENPA/Major$PWK_M_101620_input
Major$PWKmean <- rowMeans(Major[,37:39])
Major$PWKsd <- rowSds(as.matrix(Major), cols = c(37,38,39))


Major_B6 <- Major[,c(25,26)]
colnames(Major_B6) <- c("Mean","SD")
Major_B6$Strain <- c("C57BL/6J")
Major_B6$Position <- 1:234

Major_CAST <- Major[,c(30,31)]
colnames(Major_CAST) <- c("Mean","SD")
Major_CAST$Strain <- c("CAST/EiJ")
Major_CAST$Position <- 1:234


Major_LEWES <- Major[,c(35,36)]
colnames(Major_LEWES) <- c("Mean","SD")
Major_LEWES$Strain <- c("LEWES/EiJ")
Major_LEWES$Position <- 1:234

Major_PWK <- Major[,c(40,41)]
colnames(Major_PWK) <- c("Mean","SD")
Major_PWK$Strain <- c("PWK/PhJ")
Major_PWK$Position <- 1:234


Major_ggplot <- rbind(Major_B6,Major_CAST,Major_LEWES,Major_PWK)

ggplot(Major_ggplot, aes(x = Position, y = Mean, group = Strain)) +
  geom_line(aes(color = Strain)) +
  geom_ribbon(aes(y = Mean, ymin = Mean - SD, ymax = Mean + SD, fill = Strain), alpha = .1) +
  scale_fill_manual(values = c("plum", "forestgreen", "purple","red")) +
  scale_color_manual(values = c("plum", "forestgreen", "purple","red")) +
  theme_classic()


####Correlation test of values
cor(Major$B6mean, Major$CASTmean)
#[1] -0.4162681
cor(Major$B6mean, Major$LEWESmean)
#[1] -0.181959
cor(Major$B6mean, Major$PWKmean)
#[1] 0.9566204
cor(Major$LEWESmean, Major$CASTmean)
#[1] 0.6365139
cor(Major$PWKmean, Major$CASTmean)
#[1] -0.3982783
cor(Major$LEWESmean, Major$PWKmean)
#[1] -0.139759

#Kruskal Wallis
major_strain_kw <- kruskal.test(Major_ggplot$Mean ~ Major_ggplot$Strain)
#data:  Major_ggplot$Mean by Major_ggplot$Strain
#Kruskal-Wallis chi-squared = 9.6153, df = 3, p-value = 0.02214

####Perform a sliding window analysis to identify specific windows that
####are significantly different between strains

###First, create a dataframe with each strain being a column
df <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("B6", "CAST", "LEWES","PWK"))

for (i in 1:229) {
  start <- i
  end <- i + 5
  B6 <- mean(Major[start:end,25])
  CAST <- mean(Major[start:end,30])
  LEWES <- mean(Major[start:end,35])
  PWK <- mean(Major[start:end,40])
  vector <- c(B6, CAST, LEWES, PWK)
  df <- rbind(df, vector)
  colnames(df) <- c("B6", "CAST", "LEWES","PWK")
}

df$kw_PVALUE <- 0

for (i in 1:nrow(df)) {
  calc <- as.data.frame(t(df[i,1:4]))
  calc[,2] <- c("B6", "CAST", "LEWES","PWK")
  kw <- kruskal.test(calc[,1] ~ calc[,2])
  df[i,5] <- kw$p.value
}


######MINOR SATELLITE

Minor <- cbind(fastp_B6_F_092920_CENPA_sorted.bam.bed_perbp.bed_minor,
               fastp_B6_F_092920_input_sorted.bam.bed_perbp.bed_minor,
               fastp_B6_M_CENPA_sorted.bam.bed_perbp.bed_minor,
               fastp_B6_M_input_sorted.bam.bed_perbp.bed_minor,
               fastp_CAST_F_011921_CENPA_sorted.bam.bed_perbp.bed_minor,
               fastp_CAST_F_011921_input_sorted.bam.bed_perbp.bed_minor,
               fastp_CAST_M_011921_CENPA_sorted.bam.bed_perbp.bed_minor,
               fastp_CAST_M_011921_input_sorted.bam.bed_perbp.bed_minor,
               fastp_CAST_M_101620_CENPA_sorted.bam.bed_perbp.bed_minor,
               fastp_CAST_M_101620_input_sorted.bam.bed_perbp.bed_minor,
               fastp_LEWES_F_011921_CENPA_sorted.bam.bed_perbp.bed_minor,
               fastp_LEWES_F_011921_input_sorted.bam.bed_perbp.bed_minor,
               fastp_LEWES_F_012621_CENPA_sorted.bam.bed_perbp.bed_minor,
               fastp_LEWES_F_012621_input_sorted.bam.bed_perbp.bed_minor,
               fastp_LEWES_M_011921_CENPA_sorted.bam.bed_perbp.bed_minor,
               fastp_LEWES_M_011921_input_sorted.bam.bed_perbp.bed_minor,
               fastp_PWK_F_012621_CENPA_sorted.bam.bed_perbp.bed_minor,
               fastp_PWK_F_012621_input_sorted.bam.bed_perbp.bed_minor,
               fastp_PWK_F_092920_CENPA_sorted.bam.bed_perbp.bed_minor,
               fastp_PWK_F_092920_input_sorted.bam.bed_perbp.bed_minor,
               fastp_PWK_M_101620_CENPA_sorted.bam.bed_perbp.bed_minor,
               fastp_PWK_M_101620_input_sorted.bam.bed_perbp.bed_minor)

colnames(Minor) <- all_columns

Minor$B6_F_092920_CENPA <- (Minor$B6_F_092920_CENPA / 9603511)
Minor$B6_F_092920_input <- (Minor$B6_F_092920_input / 194174)
Minor$B6_M_CENPA <- (Minor$B6_M_CENPA/16231650)
Minor$B6_M_input <- (Minor$B6_M_input/383006)
Minor$CAST_F_011921_CENPA <- (Minor$CAST_F_011921_CENPA/7939934)
Minor$CAST_F_011921_input <- (Minor$CAST_F_011921_input/205561)
Minor$CAST_M_011921_CENPA <- (Minor$CAST_M_011921_CENPA/7644849)
Minor$CAST_M_011921_input <- (Minor$CAST_M_011921_input/355884)
Minor$CAST_M_101620_CENPA <- (Minor$CAST_M_101620_CENPA/4751560)
Minor$CAST_M_101620_input <- (Minor$CAST_M_101620_input/449453)
Minor$LEWES_F_011921_CENPA <- (Minor$LEWES_F_011921_CENPA/6467193)
Minor$LEWES_F_011921_input <- (Minor$LEWES_F_011921_input/469558)
Minor$LEWES_F_012621_CENPA <- (Minor$LEWES_F_012621_CENPA/7201272)
Minor$LEWES_F_012621_input <- (Minor$LEWES_F_012621_input/230072)
Minor$LEWES_M_011921_CENPA <- (Minor$LEWES_M_011921_CENPA/15839180)
Minor$LEWES_M_011921_input <- (Minor$LEWES_M_011921_input/366551)
Minor$PWK_F_012621_CENPA <- (Minor$PWK_F_012621_CENPA/10144218)
Minor$PWK_F_012621_input <- (Minor$PWK_F_012621_input/1124228)
Minor$PWK_F_092920_CENPA <- (Minor$PWK_F_092920_CENPA/4860899)
Minor$PWK_F_092920_input <- (Minor$PWK_F_092920_input/588714)
Minor$PWK_M_101620_CENPA <- (Minor$PWK_M_101620_CENPA/8984407)
Minor$PWK_M_101620_input <- (Minor$PWK_M_101620_input/1063132)

#B6
Minor$B6_F_092920_CENPAoverInput <- Minor$B6_F_092920_CENPA/Minor$B6_F_092920_input
Minor$B6_M_CENPAoverInput <- Minor$B6_M_CENPA/Minor$B6_M_input
Minor$B6mean <- rowMeans(Minor[,23:24])
Minor$B6sd <- rowSds(as.matrix(Minor), cols = c(23,24))

#CAST
Minor$CAST_F_011921_CENPAoverInput <- Minor$CAST_F_011921_CENPA/Minor$CAST_F_011921_input
Minor$CAST_M_011921_CENPAoverInput <- Minor$CAST_M_011921_CENPA/Minor$CAST_M_011921_input
Minor$CAST_M_101620_CENPAoverInput <- Minor$CAST_M_101620_CENPA/Minor$CAST_M_101620_input
Minor$CASTmean <- rowMeans(Minor[,27:29])
Minor$CASTsd <- rowSds(as.matrix(Minor), cols = c(27,28,29))

#LEWES
Minor$LEWES_F_011921_CENPAoverInput <- Minor$LEWES_F_011921_CENPA/Minor$LEWES_F_011921_input
Minor$LEWES_F_012621_CENPAoverInput <- Minor$LEWES_F_012621_CENPA/Minor$LEWES_F_012621_input
Minor$LEWES_M_011921_CENPAoverInput <- Minor$LEWES_M_011921_CENPA/Minor$LEWES_M_011921_input
Minor$LEWESmean <- rowMeans(Minor[,32:34])
Minor$LEWESsd <- rowSds(as.matrix(Minor), cols = c(32,33,34))

#PWK
Minor$PWK_F_012621_CENPAoverInput <- Minor$PWK_F_012621_CENPA/Minor$PWK_F_012621_input
Minor$PWK_F_092920_CENPAoverInput <- Minor$PWK_F_092920_CENPA/Minor$PWK_F_092920_input
Minor$PWK_M_101620_CENPAoverInput <- Minor$PWK_M_101620_CENPA/Minor$PWK_M_101620_input
Minor$PWKmean <- rowMeans(Minor[,37:39])
Minor$PWKsd <- rowSds(as.matrix(Minor), cols = c(37,38,39))


Minor_B6 <- Minor[,c(25,26)]
colnames(Minor_B6) <- c("Mean","SD")
Minor_B6$Strain <- c("C57BL/6J")
Minor_B6$Position <- 1:120

Minor_CAST <- Minor[,c(30,31)]
colnames(Minor_CAST) <- c("Mean","SD")
Minor_CAST$Strain <- c("CAST/EiJ")
Minor_CAST$Position <- 1:120


Minor_LEWES <- Minor[,c(35,36)]
colnames(Minor_LEWES) <- c("Mean","SD")
Minor_LEWES$Strain <- c("LEWES/EiJ")
Minor_LEWES$Position <- 1:120

Minor_PWK <- Minor[,c(40,41)]
colnames(Minor_PWK) <- c("Mean","SD")
Minor_PWK$Strain <- c("PWK/PhJ")
Minor_PWK$Position <- 1:120


Minor_ggplot <- rbind(Minor_B6,Minor_CAST,Minor_LEWES,Minor_PWK)


ggplot(Minor_ggplot, aes(x = Position, y = Mean, group = Strain)) +
  geom_line(aes(color = Strain)) +
  geom_hline(yintercept = 1) +
  geom_ribbon(aes(y = Mean, ymin = Mean - SD, ymax = Mean + SD, fill = Strain), alpha = .1) +
  scale_fill_manual(values = c("plum", "forestgreen", "purple","red")) +
  scale_color_manual(values = c("plum", "forestgreen", "purple","red")) +
  theme_classic()

####Correlation test of values
#Pearson's product-moment correlation

cor.test(Minor$B6mean, Minor$CASTmean)
#cor 0.4876965
#p-value = 1.609e-08
#95 percent confidence interval:
#0.3380017 0.6133239

cor.test(Minor$B6mean, Minor$LEWESmean)
#cor 0.5213508
#p-value = 1.026e-09
#95 percent confidence interval:
#0.3773742 0.6407186

cor.test(Minor$B6mean, Minor$PWKmean)
#cor 0.869829
#p-value < 2.2e-16
#95 percent confidence interval:
#0.8181440 0.9075715

cor.test(Minor$LEWESmean, Minor$CASTmean)
#cor 0.9624259
#p-value < 2.2e-16
#95 percent confidence interval:
#0.9464541 0.9736980

cor.test(Minor$PWKmean, Minor$CASTmean)
#cor 0.8109182
#p-value < 2.2e-16
#95 percent confidence interval:
#0.7391059 0.8645038

cor.test(Minor$LEWESmean, Minor$PWKmean)
#cor 0.8515345
#p-value < 2.2e-16
#95 percent confidence interval:
#0.7933881 0.8942815

hist(Minor$B6mean)
hist(Minor$CASTmean)
hist(Minor$LEWESmean)
hist(Minor$PWKmean)

##ANOVA

minor_strain_aov <- aov(Minor_ggplot$Mean ~ Minor_ggplot$Strain)
summary(minor_strain_aov)

minor_strain_positions_aov <- aov(Minor_ggplot$Mean ~ Minor_ggplot$Strain + Minor_ggplot$Position)
summary(minor_strain_positions_aov)

#Kruskal Wallis
minor_strain_kw <- kruskal.test(Minor_ggplot$Mean ~ Minor_ggplot$Strain)
#data:  Minor_ggplot$Mean by Minor_ggplot$Strain
#Kruskal-Wallis chi-squared = 14.913, df = 3, p-value = 0.001893

###############################################################
#######################Correlation heatmap#####################
###############################################################
#CENP-A ChIP only
MinorCorCENPAChIP <- Minor[,c(1,3,5,7,9,11,13,15,17,19,21)]
colnames(MinorCorCENPAChIP) <- c("B61","B62","CAST1","CAST2","CAST3","LEWES1","LEWES2","LEWES3","PWK1","PWK2","PWK3")
corMinorCorCENPAChIP <- cor(MinorCorCENPAChIP)

corrplot(corMinorCorCENPAChIP, diag = TRUE,type = 'lower') %>%
  corrRect(name = c('B61', 'CAST1', 'LEWES1', 'PWK1','PWK3'))

#ChIP over input relative abundance
MinorCor <- Minor[,c(23,24,27,28,29,32,33,34,37,38,39)]
colnames(MinorCor) <- c("B61","B62","CAST1","CAST2","CAST3","LEWES1","LEWES2","LEWES3","PWK1","PWK2","PWK3")
corMinorCor <- cor(MinorCor)

corrplot(corMinorCor, diag = TRUE,type = 'lower',col = COL1('Blues', 10), addCoef.col = 'grey') %>%
  corrRect(name = c('B61', 'CAST1', 'LEWES1', 'PWK1','PWK3')) 


corrplot(corMinorCor, diag = TRUE,order = 'hclust',col = COL1('Blues', 10), addCoef.col = 'grey')


########################################################
#####################IWATA OTSUBO#######################
########################################################

Minor <- cbind(fastp_B6_F_092920_CENPA_sorted.bam.bed_perbp.bed_minor,
               fastp_B6_F_092920_input_sorted.bam.bed_perbp.bed_minor,
               fastp_B6_M_CENPA_sorted.bam.bed_perbp.bed_minor,
               fastp_B6_M_input_sorted.bam.bed_perbp.bed_minor,
               `fastp_BL6-1-CENPA_sorted.bam.bed_perbp.bed_minor`,
               `fastp_BL6-1-INP_sorted.bam.bed_perbp.bed_minor`,
               `fastp_BL6-2-CENPA_sorted.bam.bed_perbp.bed_minor`,
               `fastp_BL6-2-INP_sorted.bam.bed_perbp.bed_minor`,
               `fastp_BL6-3-CENPA_sorted.bam.bed_perbp.bed_minor`,
               `fastp_BL6-3-INP_sorted.bam.bed_perbp.bed_minor`,
               `fastp_CF1-1-CENPA_sorted.bam.bed_perbp.bed_minor`,
               `fastp_CF1-1-INP_sorted.bam.bed_perbp.bed_minor`,
               `fastp_CHPO-1-CENPA_sorted.bam.bed_perbp.bed_minor`,
               `fastp_CHPO-1-INP_sorted.bam.bed_perbp.bed_minor`,
               `fastp_CHPO-2-CENPA_sorted.bam.bed_perbp.bed_minor`,
               `fastp_CHPO-2-INP_sorted.bam.bed_perbp.bed_minor`,
               `fastp_CHPO-3-CENPA_sorted.bam.bed_perbp.bed_minor`,
               `fastp_CHPO-3-INP_sorted.bam.bed_perbp.bed_minor`)

all_columns <- c("B6_F_092920_CENPA",
                 "B6_F_092920_input",
                 "B6_M_CENPA", 
                 "B6_M_input",
                 "BL6-1-CENPA",
                 "BL6-1-INP",
                 "BL6-2-CENPA",
                 "BL6-2-INP",
                 "BL6-3-CENPA",
                 "BL6-3-INP",
                 "CF1-1-CENPA",
                 "CF1-1-INP",
                 "CHPO-1-CENPA",
                 "CHPO-1-INP",
                 "CHPO-2-CENPA",
                 "CHPO-2-INP",
                 "CHPO-3-CENPA",
                 "CHPO-3-INP")

colnames(Minor) <- all_columns

Minor$B6_F_092920_CENPA <- (Minor$B6_F_092920_CENPA / 9603511)
Minor$B6_F_092920_input <- (Minor$B6_F_092920_input / 194174)
Minor$B6_M_CENPA <- (Minor$B6_M_CENPA/16231650)
Minor$B6_M_input <- (Minor$B6_M_input/383006)
Minor$`BL6-1-CENPA` <- (Minor$`BL6-1-CENPA`/17107374)
Minor$`BL6-1-INP` <- (Minor$`BL6-1-INP`/135300)
Minor$`BL6-2-CENPA` <- (Minor$`BL6-2-CENPA`/784431)
Minor$`BL6-2-INP` <- (Minor$`BL6-2-INP`/88039)
Minor$`BL6-3-CENPA` <- (Minor$`BL6-3-CENPA`/742549)
Minor$`BL6-3-INP` <- (Minor$`BL6-3-INP`/74505)
Minor$`CF1-1-CENPA` <- (Minor$`CF1-1-CENPA`/11886612)
Minor$`CF1-1-INP` <- (Minor$`CF1-1-INP`/353530)
Minor$`CHPO-1-CENPA` <- (Minor$`CHPO-1-CENPA`/9158789)
Minor$`CHPO-1-INP` <- (Minor$`CHPO-1-INP`/17495)
Minor$`CHPO-2-CENPA` <- (Minor$`CHPO-2-CENPA`/313520)
Minor$`CHPO-2-INP` <- (Minor$`CHPO-2-INP`/11595)
Minor$`CHPO-3-CENPA` <- (Minor$`CHPO-3-CENPA`/412817)
Minor$`CHPO-3-INP` <- (Minor$`CHPO-3-INP`/10318)

#B6
Minor$B6_F_092920_CENPAoverInput <- Minor$B6_F_092920_CENPA/Minor$B6_F_092920_input
Minor$B6_M_CENPAoverInput <- Minor$B6_M_CENPA/Minor$B6_M_input
Minor$B6mean <- rowMeans(Minor[,19:20])
Minor$B6sd <- rowSds(as.matrix(Minor), cols = c(19,20))

#B6 Iwata Otsubo
Minor$BL61_CENPAoverInput <- Minor$`BL6-1-CENPA`/Minor$`BL6-1-INP`
Minor$BL62_CENPAoverInput <- Minor$`BL6-2-CENPA`/Minor$`BL6-2-INP`
Minor$BL63_CENPAoverInput <- Minor$`BL6-3-CENPA`/Minor$`BL6-3-INP`
Minor$B6IwataMean <- rowMeans(Minor[,23:25])
Minor$B6Iwatasd <- rowSds(as.matrix(Minor), cols = c(23,24,25))

#CF1 Iwata Otsubo
Minor$CF1_CENPAoverInput <- Minor$`CF1-1-CENPA`/Minor$`CF1-1-INP`

#CHPO
Minor$`CHPO-1-CENPAoverInput` <- Minor$`CHPO-1-CENPA`/Minor$`CHPO-1-INP`
Minor$`CHPO-2-CENPAoverInput` <- Minor$`CHPO-2-CENPA`/Minor$`CHPO-2-INP`
Minor$`CHPO-3-CENPAoverInput` <- Minor$`CHPO-3-CENPA`/Minor$`CHPO-3-INP`
Minor$CHPOMean <- rowMeans(Minor[,29:31])
Minor$CHPOsd <- rowSds(as.matrix(Minor), cols = c(29,30,31))


##
Minor_B6 <- Minor[,c(21,22)]
colnames(Minor_B6) <- c("Mean","SD")
Minor_B6$Strain <- c("C57BL/6J")
Minor_B6$Position <- 1:120

Minor_B6Iwata <- Minor[,c(26,27)]
colnames(Minor_B6Iwata) <- c("Mean","SD")
Minor_B6Iwata$Strain <- c("C57BL/6J Iwata")
Minor_B6Iwata$Position <- 1:120


Minor_CF1 <- Minor[,c(28)]
Minor_CF1 <- as.data.frame(Minor_CF1)
Minor_CF1$SD <- 0
colnames(Minor_CF1) <- c("Mean","SD")
Minor_CF1$Strain <- c("CF1")
Minor_CF1$Position <- 1:120

Minor_CHPO <- Minor[,c(32,33)]
colnames(Minor_CHPO) <- c("Mean","SD")
Minor_CHPO$Strain <- c("ZALENDE/EiJ")
Minor_CHPO$Position <- 1:120

Minor_ggplot <- rbind(Minor_B6,Minor_B6Iwata,Minor_CF1,Minor_CHPO)

ggplot(Minor_ggplot, aes(x = Position, y = Mean, group = Strain)) +
  geom_line(aes(color = Strain)) +
  geom_hline(yintercept = 1) +
  geom_ribbon(aes(y = Mean, ymin = Mean - SD, ymax = Mean + SD, fill = Strain), alpha = .1) +
  scale_fill_manual(values = c("plum", "purple1", "purple","green")) +
  scale_color_manual(values = c("plum", "purple1", "purple","green")) +
  theme_classic()



