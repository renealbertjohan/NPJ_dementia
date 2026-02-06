## 2.1 BiTAge: A transcriptome-based aging clock near the theoretical limit of accuracy
# Load libraries
library(tidyverse)
library(ggpubr) # For function stat_compare_means
library(ggh4x) # For function facet_grid2
library(ggtext) # For adding text to figure, if wanted/needed
library(rstatix) # For tukeyhsd() function in figures
library(car) # qqPlot function for Two-way ANOVA

# Load phenoHEROES_RNAAge_FPKM data frame obtained from RNAAgeCalc.R script
phenoHEROES_RNAAge_FPKM <- read.delim("NPJ_dementia/1_Datasets/PhenoData/phenoHEROES_RNAAge_FPKM.csv", sep = ",", dec = ".", header = TRUE)
row.names(phenoHEROES_RNAAge_FPKM) <- phenoHEROES_RNAAge_FPKM$Tube_code
phenoHEROES_RNAAge_FPKM$X <- NULL
phenoHEROES_RNAAge_FPKM

# Load BiT Age results perform through Python script
BiTage <- read.delim("NPJ_dementia/2_Scripts/BiTage/BiT_age.txt", sep = " ", dec = ".", header = FALSE)
colnames(BiTage)[1] <- "Tube_code"
colnames(BiTage)[2] <- "BiT"
BiTage <- BiTage[ order(match(BiTage$Tube_code, phenoHEROES_RNAAge_FPKM$Tube_code)), ]
rownames(BiTage) <- BiTage$Tube_code

# Check if phenoHEROES_RNAAge_FPKM$Tube_code and phenoHEROES$Tube_code are identical, including order
identical(BiTage$Tube_code, phenoHEROES_RNAAge_FPKM$Tube_code)

AgeDif_BiT <- BiTage$BiT - phenoHEROES_RNAAge_FPKM$ChronAge # Add also the calculated age difference (i.e., RNAAge - ChronAge)
phenoHEROES_RNAAge_BiTage <- as.data.frame(append(phenoHEROES_RNAAge_FPKM, list(BiT_Age = BiTage$BiT, AgeDif_BiT = AgeDif_BiT), after = 5))
row.names(phenoHEROES_RNAAge_BiTage) <- phenoHEROES_RNAAge_BiTage$Tube_code
write.csv(phenoHEROES_RNAAge_BiTage, "NPJ_dementia/1_Datasets/PhenoData/phenoHEROES_RNAAge_BiTage.csv")

# Visualization of data and statistical testing
ChronAge <- as.data.frame(cbind(phenoHEROES_RNAAge_BiTage$Status, phenoHEROES_RNAAge_BiTage$Tube_code, phenoHEROES_RNAAge_BiTage$Condition, phenoHEROES_RNAAge_BiTage$ChronAge))
ChronAge$Timepoint <- rep("ChronAge", 20)
BITAge <- as.data.frame(cbind(phenoHEROES_RNAAge_BiTage$Status, phenoHEROES_RNAAge_BiTage$Tube_code, phenoHEROES_RNAAge_BiTage$Condition, phenoHEROES_RNAAge_BiTage$BiT_Age))
BITAge$Timepoint <- rep("BiTAge", 20)
RNAAGECAL <- rbind(ChronAge, BITAge)

colnames(RNAAGECAL)[1] <- "Status"
colnames(RNAAGECAL)[2] <- "Tube_code"
colnames(RNAAGECAL)[3] <- "Condition"
colnames(RNAAGECAL)[4] <- "Age"
RNAAGECAL$Age <- as.numeric(RNAAGECAL$Age)
RNAAGECAL$Status <- factor(RNAAGECAL$Status, levels = c("Control", "PDD", "AD", "DSD"))
factor(RNAAGECAL$Timempoint, levels = c("ChronAge", "BiTAge"))

# Calculations
# Control
mean(RNAAGECAL[RNAAGECAL$Status == "Control" & RNAAGECAL$Timepoint == "BiTAge", "Age"])
sd(RNAAGECAL[RNAAGECAL$Status == "Control" & RNAAGECAL$Timepoint == "BiTAge", "Age"])

# PDD
mean(RNAAGECAL[RNAAGECAL$Status == "PDD" & RNAAGECAL$Timepoint == "BiTAge", "Age"])
sd(RNAAGECAL[RNAAGECAL$Status == "PDD" & RNAAGECAL$Timepoint == "BiTAge", "Age"])

# AD
mean(RNAAGECAL[RNAAGECAL$Status == "AD" & RNAAGECAL$Timepoint == "BiTAge", "Age"])
sd(RNAAGECAL[RNAAGECAL$Status == "AD" & RNAAGECAL$Timepoint == "BiTAge", "Age"])

# DSD
mean(RNAAGECAL[RNAAGECAL$Status == "DSD" & RNAAGECAL$Timepoint == "BiTAge", "Age"])
sd(RNAAGECAL[RNAAGECAL$Status == "DSD" & RNAAGECAL$Timepoint == "BiTAge", "Age"])

## Comparison plot --> Chronological age versus Transcriptional Age
stat.test_RNAAGECAL <- aov(Age ~ Status * Timepoint, data = RNAAGECAL) %>%
  tukey_hsd()
stat.test_RNAAGECAL <- stat.test_RNAAGECAL[c(11,18,24,29),] 

RNAAGECAL$Timepoint <- factor(RNAAGECAL$Timepoint,
                       levels = c('ChronAge','BiTAge'), ordered = TRUE)

png(filename="NPJ_dementia/3_Figures/2_RNAAgeCalc/BiTAge/BiT_Age_Calculation.png",
    width     = 50,
    height    = 15,
    units     = "cm",
    res       = 1200,
    pointsize = 4)
RNAAGECAL
p <- ggplot(RNAAGECAL, aes(x = Timepoint, y = Age, fill = Timepoint)) +
  stat_boxplot(geom = "errorbar",
               width = 0.4,
               color = "black") +
  geom_boxplot(width = 0.8, alpha = 0.7) +
  theme_classic() +
  scale_fill_manual(values = c("seagreen2", "firebrick2"), name = "Age") + 
  scale_y_continuous(limits = c(50, 125), expand = c(0, 0)) +
  theme(axis.text = element_text(color = "black")) +
  geom_point(color = "gray42", size = 5, alpha = 0.7) +
  geom_line(aes(group = Tube_code), color = "gray42", linetype = "41") +
  facet_grid2(. ~ Status, scales = "free", space = "free",
              strip = strip_themed(background_x = list(element_rect(fill = "limegreen"),
                                                       element_rect(fill = "navyblue"),
                                                       element_rect(fill = "gold"),
                                                       element_rect(fill = "darkred")))) + 
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 22, face = "bold"),
        plot.title = element_markdown(hjust = 0.5, size = 30, face = "bold"),
        panel.grid = element_blank(),
        legend.position = "none",
        legend.title = element_text(size = 20, face = "bold", hjust = 0.1),
        legend.text = element_text(size = 10),
        legend.key.size = unit(2, "cm"),
        strip.text = element_text(color = "white", size = 30)) +
  labs(title = "Chronological age *versus* Transcriptional age") +
  ylab(expression(bold("Age (years)")))
anno1 <- annotate("text", label = c("p = 0.427"), size = 8, x = c(1.5), y = 120)
anno2 <- annotate("text", label = c("*"), size = 12, x = c(1.5), y = 120)
anno3 <- annotate("text", label = c("**"), size = 12, x = c(1.5), y = 120)
anno4 <- annotate("text", label = c("****"), size = 12, x = c(1.5), y = 120)
p + at_panel(anno4, PANEL == 4)

dev.off()

# Factorize
phenoHEROES_RNAAge_BiTage$Status <- factor(phenoHEROES_RNAAge_BiTage$Status, levels = c("Control", "PDD", "AD", "DSD"))

## Chronological Age
stat.test_ChronAge <- aov(ChronAge ~ Status, data = phenoHEROES_RNAAge_BiTage) %>%
  tukey_hsd()
stat.test_ChronAge

png(filename="NPJ_dementia/3_Figures/2_RNAAgeCalc/BiTAge/Chron_BiT_Age.png",
    width     = 20,
    height    = 15,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(data = phenoHEROES_RNAAge_BiTage, aes(x = Status, y = ChronAge, fill = factor(Status))) +
  stat_boxplot(geom = "errorbar",
               width = 0.4,
               color = "black") + 
  geom_boxplot(outlier.shape = NA)  +
  geom_jitter(width = 0.2, size = 3, color = "gray42") +
  scale_fill_discrete(breaks = c("Control", "PDD", "AD", "DSD"), type = c("limegreen", "navyblue", "gold", "darkred")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 18, face = "bold"),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
        axis.title.y = element_text(size = 26, face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.caption = element_text(hjust = 0.5, size = 12)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 125), breaks = c(0, 20, 40, 60, 80, 100, 120)) +
  ylab(expression(bold("Age (years)"))) +
  ggtitle(paste("Chronological age")) + 
  stat_pvalue_manual(stat.test_ChronAge, hide.ns = TRUE, y.position = c(95, 105, 115), size = 12)

dev.off()

## Transcriptional BiT Age
stat.test_RNAAge <- aov(BiT_Age ~ Status, data = phenoHEROES_RNAAge_BiTage) %>%
  tukey_hsd()
stat.test_RNAAge

png(filename="NPJ_dementia/3_Figures/2_RNAAgeCalc/BiTAge/BiT_Age.png",
    width     = 20,
    height    = 15,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(data = phenoHEROES_RNAAge_BiTage, aes(x = Status, y = BiT_Age, fill = factor(Status))) +
  stat_boxplot(geom = "errorbar",
               width = 0.4,
               color = "black") + 
  geom_boxplot(outlier.shape = NA)  +
  geom_jitter(width = 0.2, size = 3, color = "gray42") +
  scale_fill_discrete(breaks = c("Control", "PDD", "AD", "DSD"), type = c("limegreen", "navyblue", "gold", "darkred")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 18, face = "bold"),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
        axis.title.y = element_text(size = 26, face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.caption = element_text(hjust = 0.5, size = 12)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 125), breaks = c(0, 20, 40, 60, 80, 100, 120)) +
  ylab(expression(bold("Age (years)"))) +
  ggtitle(paste("Transcriptional age"))

dev.off()

# Age difference comparison between conditions
# To evaluate if the age difference (Diff) was significantly increased between demented individuals with Down sydnrome versus controls.
# I performed one-way ANOVA with Tukey’s multiple comparison with R function ‘aov’ in ggpubr package for the locus coeruleus.
stat.test_AgeDif <- aov(AgeDif_BiT ~ Status, data = phenoHEROES_RNAAge_BiTage) %>%
  tukey_hsd()
stat.test_AgeDif
phenoHEROES_RNAAge_BiTage
png(filename="NPJ_dementia/3_Figures/2_RNAAgeCalc/BiTAge/BiT_Age_Dif.png",
    width     = 20,
    height    = 15,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(data = phenoHEROES_RNAAge_BiTage, aes(x = Status, y = AgeDif_BiT, fill = factor(Status))) +
  stat_boxplot(geom = "errorbar",
               width = 0.4,
               color = "black") + 
  geom_boxplot(outlier.shape = NA)  +
  geom_jitter(width = 0.2, size = 3, color = "gray42") +
  scale_fill_discrete(breaks = c("Control", "PDD", "AD", "DSD"), type = c("limegreen", "navyblue", "gold", "darkred")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 18, face = "bold"),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
        axis.title.y = element_text(size = 26, face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.caption = element_text(hjust = 0.5, size = 12)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-10, 60)) +
  ylab(expression(bold("Age acceleration (years)"))) +
  ggtitle(paste("Differential age")) + 
  stat_pvalue_manual(stat.test_AgeDif, hide.ns = TRUE, y.position = c(40, 45, 50), size = 12)

dev.off()
