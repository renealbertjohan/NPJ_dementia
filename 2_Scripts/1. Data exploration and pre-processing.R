# 1. Data exploration and pre-processing
# Load libraries
library(DT)
library(webshot)
library(tidyverse)
library(dplyr)
library(lmtest) # bptest, testing for Homoscedasticity assumption
library(stats) # kruskal.test function
library(FSA) # dunnTest function
library(ggfortify) # 2D Principal Component Analysis (2D PCA)
library(scatterplot3d) # 3D Principal Component Analysis (3D PCA)
library(ggplot2) # Figures with ggplot
library(ggpubr) # For function stat_compare_means
library(limma) # Required for egdeR package
library(edgeR) # TMM normalization
library(factoextra) # Dendrogram, Hierarchical Clustering
library(rstatix) # For tukeyhsd() function in figures

# PHENOTYPE DATA ----------------------------------------------------------
# Step 1: Load the phenotype data of the HEROES project -------------------
phenoHEROES <- read.delim("NPJ_dementia/1_Datasets/PhenoData.csv", sep = ",", dec = ",", header = TRUE)
row.names(phenoHEROES) <- phenoHEROES$Tube_code
phenoHEROES$Condition <- paste(phenoHEROES$Brain_region, phenoHEROES$Status, sep="_")

# Factorize Status
phenoHEROES$Status <- factor(phenoHEROES$Status, levels = c("Control", "PDD", "AD", "DSD"))

# Clinical diagnosis/Origin of each group
phenoHEROES[phenoHEROES$Status == "Control", "BrainBank"]
phenoHEROES[phenoHEROES$Status == "PDD", "BrainBank"]
phenoHEROES[phenoHEROES$Status == "AD", "BrainBank"]
phenoHEROES[phenoHEROES$Status == "DSD", "BrainBank"]

# Step 2: Create and save datatable -----------------------------------------
# Show relevant data in a table
pheno <- phenoHEROES[, c(2, 1, 3, 4, 5, 11, 12, 13, 14, 10, 8, 7, 15)]
colnames(pheno)[6] <- "AD Braak staging"
colnames(pheno)[7] <- "Thal phasing"
colnames(pheno)[8] <- "CERAD scoring"
colnames(pheno)[9] <- "PD Braak staging"
colnames(pheno)[10] <- "PMI (h)"
colnames(pheno)[11] <- "Source"
colnames(pheno)[12] <- "Identifier (ID)"
colnames(pheno)[13] <- "Brain region"

dtable <- datatable(pheno[,] %>% arrange('Brain region', factor(Status, levels = c("Control", "PDD", "AD", "DSD")), Sex, Age), options = list(pageLength = -1))
html <- "dtable.html"
saveWidget(dtable, html)
webshot(html, "NPJ_dementia/3_Figures/1_Data_Exploration/datatable.png")

# Step 3: Phenotype analyses ----------------------------------------------
# Data for the summarized description in Crans et al., 2026
## Sample sizes
sum(phenoHEROES$Status == "Control") # Sample number of control (non-demented)
sum(phenoHEROES$Status == "PDD") # Sample number of Parkinson´s disease dementia (PDD)
sum(phenoHEROES$Status == "AD") # Sample number of Alzheimer´s disease (AD)
sum(phenoHEROES$Status == "DSD") # Sample number of Down syndrome dementia (DSD)

## Chronological age
png(filename="NPJ_dementia/3_Figures/1_Data_Exploration/Chronological_Age.png",
    width     = 25,
    height    = 15,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(data = phenoHEROES, aes(x = Age, fill = factor(Status))) +
  geom_histogram(binwidth = 5, color = "black") +
  scale_fill_discrete(name = "Status", breaks = c("Control", "PDD", "AD", "DSD"), type = c("limegreen", "navyblue", "gold", "darkred")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, face = "bold"),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, units = "cm")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 5)) +
  ylab(expression(bold("# of individuals"))) +
  xlab(expression(bold("Age (years)"))) + 
  ggtitle(paste("Age distribution"))

dev.off()

# Distribution of the chronological age at death (years), mean, sd, range
# Control (non-demented)
mean(phenoHEROES[phenoHEROES$Status == "Control",]$Age)
sd(phenoHEROES[phenoHEROES$Status == "Control",]$Age)
range(phenoHEROES[phenoHEROES$Status == "Control",]$Age)
# Parkinson´s disease dementia
mean(phenoHEROES[phenoHEROES$Status == "PDD",]$Age)
sd(phenoHEROES[phenoHEROES$Status == "PDD",]$Age)
range(phenoHEROES[phenoHEROES$Status == "PDD",]$Age)
# Alzheimer´s disease
mean(phenoHEROES[phenoHEROES$Status == "AD",]$Age)
sd(phenoHEROES[phenoHEROES$Status == "AD",]$Age)
range(phenoHEROES[phenoHEROES$Status == "AD",]$Age)
# Down syndrome dementia
mean(phenoHEROES[phenoHEROES$Status == "DSD",]$Age)
sd(phenoHEROES[phenoHEROES$Status == "DSD",]$Age)
range(phenoHEROES[phenoHEROES$Status == "DSD",]$Age)

# Test for normal distribution
shapiro.test(residuals(lm(phenoHEROES$Age ~ phenoHEROES$Status)))
# Test for Homoscedasticity
bptest(lm(phenoHEROES$Age ~ phenoHEROES$Status), studentize = FALSE)
# ANOVA with a post-hoc test based on the studentized range distribution (Tukey's Honest Significant Difference test)
anova_Age <- aov(phenoHEROES$Age ~ phenoHEROES$Status)
summary(anova_Age)
TukeyHSD(anova_Age) # Age of Down syndrome dementia is significantly from Control (non-demented), Alzheimer´s disease, and Parkinson´s disease dementia

## Sex
table(factor(phenoHEROES$Sex, levels = c("Male", "Female")), factor(phenoHEROES$Status, levels = c("Control", "PDD", "AD", "DSD")))
# Fisher´s Exact Test for Count Data, a statistical significance test used in the analysis of contingency tables
fisher.test(table(phenoHEROES$Sex, factor(phenoHEROES$Status, levels = c("Control", "PDD", "AD", "DSD"))))
fisher.test(table(phenoHEROES$Sex, factor(phenoHEROES$Status, levels = c("Control", "PDD")))) # Control vs. PDD
fisher.test(table(phenoHEROES$Sex, factor(phenoHEROES$Status, levels = c("Control", "AD")))) # Control vs. AD
fisher.test(table(phenoHEROES$Sex, factor(phenoHEROES$Status, levels = c("Control", "DSD")))) # Control vs. DSD
fisher.test(table(phenoHEROES$Sex, factor(phenoHEROES$Status, levels = c("PDD", "AD")))) # PDD vs. AD
fisher.test(table(phenoHEROES$Sex, factor(phenoHEROES$Status, levels = c("PDD", "DSD")))) # PDD vs. DSD
fisher.test(table(phenoHEROES$Sex, factor(phenoHEROES$Status, levels = c("AD", "DSD")))) # AD vs. DSD

## ApoE genotype
table(phenoHEROES$APOE4, factor(phenoHEROES$Status, levels = c("Control", "PDD", "AD", "DSD")))
# Fisher´s Exact Test for Count Data, a statistical significance test used in the analysis of contingency tables
fisher.test(table(phenoHEROES$APOE4, factor(phenoHEROES$Status, levels = c("Control", "PDD", "AD", "DSD"))))
fisher.test(table(phenoHEROES$APOE4, factor(phenoHEROES$Status, levels = c("Control", "PDD")))) # Control vs. PDD
fisher.test(table(phenoHEROES$APOE4, factor(phenoHEROES$Status, levels = c("Control", "AD")))) # Control vs. AD
fisher.test(table(phenoHEROES$APOE4, factor(phenoHEROES$Status, levels = c("Control", "DSD")))) # Control vs. DSD
fisher.test(table(phenoHEROES$APOE4, factor(phenoHEROES$Status, levels = c("PDD", "AD")))) # PDD vs. AD
fisher.test(table(phenoHEROES$APOE4, factor(phenoHEROES$Status, levels = c("PDD", "DSD")))) # PDD vs. DSD
fisher.test(table(phenoHEROES$APOE4, factor(phenoHEROES$Status, levels = c("AD", "DSD")))) # AD vs. DSD

## AD Braak staging
phenoHEROES[phenoHEROES$Status == "Control",]$AD_Braak
phenoHEROES[phenoHEROES$Status == "PDD",]$AD_Braak
phenoHEROES[phenoHEROES$Status == "AD",]$AD_Braak
phenoHEROES[phenoHEROES$Status == "DSD",]$AD_Braak

## Thal phasing
phenoHEROES[phenoHEROES$Status == "Control",]$Thal_phase
phenoHEROES[phenoHEROES$Status == "PDD",]$Thal_phase
phenoHEROES[phenoHEROES$Status == "AD",]$Thal_phase
phenoHEROES[phenoHEROES$Status == "DSD",]$Thal_phase

## CERAD score
phenoHEROES[phenoHEROES$Status == "Control",]$CERAD_score
phenoHEROES[phenoHEROES$Status == "PDD",]$CERAD_score
phenoHEROES[phenoHEROES$Status == "AD",]$CERAD_score
phenoHEROES[phenoHEROES$Status == "DSD",]$CERAD_score

## PD Braak staging
phenoHEROES[phenoHEROES$Status == "Control",]$PD_Braak
phenoHEROES[phenoHEROES$Status == "PDD",]$PD_Braak
phenoHEROES[phenoHEROES$Status == "AD",]$PD_Braak
phenoHEROES[phenoHEROES$Status == "DSD",]$PD_Braak

## Post-mortem interval (PMI), refers to the time elapsed between death and the examination of a body
stat.test_HEROES_PMI <- aov(PMD ~ Status, data = phenoHEROES) %>%
  tukey_hsd()
stat.test_HEROES_PMI

png(filename="NPJ_dementia/3_Figures/1_Data_Exploration/PMI.png",
    width     = 20,
    height    = 15,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(data = phenoHEROES, aes(x = Status, y = PMD, fill = factor(Status))) +
  stat_boxplot(geom = "errorbar",
               width = 0.4,
               color = "black") + 
  geom_boxplot(outlier.shape = NA)  +
  geom_jitter(width = 0.2, size = 3, color = "gray42") +
  scale_fill_discrete(breaks = c("Control", "PDD", "AD", "DSD"), type = c("limegreen", "navyblue", "gold", "darkred")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, face = "bold"),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 30)) +
  ylab(expression(bold("Post-mortem interval (hrs)"))) +
  ggtitle(paste("Post-mortem interval"))

dev.off()

# Distribution of post-mortem interval (mean, sd, range)
# Control (non-demented)
mean(phenoHEROES[phenoHEROES$Status == "Control",]$PMD)
sd(phenoHEROES[phenoHEROES$Status == "Control",]$PMD)
range(phenoHEROES[phenoHEROES$Status == "Control",]$PMD)
# Parkinson´s disease dementia
mean(phenoHEROES[phenoHEROES$Status == "PDD",]$PMD)
sd(phenoHEROES[phenoHEROES$Status == "PDD",]$PMD)
range(phenoHEROES[phenoHEROES$Status == "PDD",]$PMD)
# Alzheimer´s disease
mean(phenoHEROES[phenoHEROES$Status == "AD",]$PMD)
sd(phenoHEROES[phenoHEROES$Status == "AD",]$PMD)
range(phenoHEROES[phenoHEROES$Status == "AD",]$PMD)
# Down syndrome dementia
mean(phenoHEROES[phenoHEROES$Status == "DSD",]$PMD)
sd(phenoHEROES[phenoHEROES$Status == "DSD",]$PMD)
range(phenoHEROES[phenoHEROES$Status == "DSD",]$PMD)

# Test for normal distribution
shapiro.test(residuals(lm(phenoHEROES$PMD ~ phenoHEROES$Status)))
# Test for Homoscedasticity
bptest(lm(phenoHEROES$PMD ~ phenoHEROES$Status), studentize = FALSE)
# ANOVA with a post-hoc test based on the studentized range distribution (Tukey's Honest Significant Difference test)
anova_PMI <- aov(phenoHEROES$PMD ~ phenoHEROES$Status)
summary(anova_PMI)
TukeyHSD(anova_PMI)

## RNA Integrity Number (RIN)
stat.test_HEROES_RIN <- aov(RIN ~ Status, data = phenoHEROES) %>%
  tukey_hsd()
stat.test_HEROES_RIN

png(filename="NPJ_dementia/3_Figures/1_Data_Exploration/RIN.png",
    width     = 20,
    height    = 15,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(data = phenoHEROES, aes(x = Status, y = RIN, fill = factor(Status))) +
  stat_boxplot(geom = "errorbar",
               width = 0.4,
               color = "black") + 
  geom_boxplot(outlier.shape = NA)  +
  geom_jitter(width = 0.2, size = 3, color = "gray42") +
  scale_fill_discrete(breaks = c("Control", "PDD", "AD", "DSD"), type = c("limegreen", "navyblue", "gold", "darkred")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, face = "bold"),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 10), breaks = c(0, 2, 4, 6, 8, 10)) +
  ylab(expression(bold("RNA Integrity Number (RIN)"))) +
  ggtitle(paste("RNA Integrity Number"))

dev.off()

# Distribution of RIN (mean, sd, range)
# Control (non-demented)
mean(phenoHEROES[phenoHEROES$Status == "Control",]$RIN)
sd(phenoHEROES[phenoHEROES$Status == "Control",]$RIN)
range(phenoHEROES[phenoHEROES$Status == "Control",]$RIN)
# Parkinson´s disease dementia
mean(phenoHEROES[phenoHEROES$Status == "PDD",]$RIN)
sd(phenoHEROES[phenoHEROES$Status == "PDD",]$RIN)
range(phenoHEROES[phenoHEROES$Status == "PDD",]$RIN)
# Alzheimer´s disease
mean(phenoHEROES[phenoHEROES$Status == "AD",]$RIN)
sd(phenoHEROES[phenoHEROES$Status == "AD",]$RIN)
range(phenoHEROES[phenoHEROES$Status == "AD",]$RIN)
# Down syndrome dementia
mean(phenoHEROES[phenoHEROES$Status == "DSD",]$RIN)
sd(phenoHEROES[phenoHEROES$Status == "DSD",]$RIN)
range(phenoHEROES[phenoHEROES$Status == "DSD",]$RIN)

# Test for normal distribution
shapiro.test(residuals(lm(phenoHEROES$RIN ~ phenoHEROES$Status)))
# Test for Homoscedasticity
bptest(lm(phenoHEROES$RIN ~ phenoHEROES$Status), studentize = FALSE)
# ANOVA with a post-hoc test based on the studentized range distribution (Tukey's Honest Significant Difference test)
anova_RIN <- aov(phenoHEROES$RIN ~ phenoHEROES$Status)
summary(anova_RIN)
TukeyHSD(anova_RIN)

## RNA Concentration
stat.test_HEROES_Conc <- aov(Concentration ~ Status, data = phenoHEROES) %>%
  tukey_hsd()
stat.test_HEROES_Conc

png(filename="NPJ_dementia/3_Figures/1_Data_Exploration/RNA_concentration.png",
    width     = 20,
    height    = 15,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(data = phenoHEROES, aes(x = Status, y = Concentration, fill = factor(Status))) +
  stat_boxplot(geom = "errorbar",
               width = 0.4,
               color = "black") + 
  geom_boxplot(outlier.shape = NA)  +
  geom_jitter(width = 0.2, size = 3, color = "gray42") +
  scale_fill_discrete(breaks = c("Control", "PDD", "AD", "DSD"), type = c("limegreen", "navyblue", "gold", "darkred")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 14, face = "bold"),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 40)) +
  ylab(expression(bold("Concentration (ng/µl)"))) +
  ggtitle(paste("RNA concentration"))

dev.off()

# Distribution of RNA concentration (mean, sd, range)
# Control (non-demented)
mean(phenoHEROES[phenoHEROES$Status == "Control",]$Concentration)
sd(phenoHEROES[phenoHEROES$Status == "Control",]$Concentration)
range(phenoHEROES[phenoHEROES$Status == "Control",]$Concentration)
# Parkinson´s disease dementia
mean(phenoHEROES[phenoHEROES$Status == "PDD",]$Concentration)
sd(phenoHEROES[phenoHEROES$Status == "PDD",]$Concentration)
range(phenoHEROES[phenoHEROES$Status == "PDD",]$Concentration)
# Alzheimer´s disease
mean(phenoHEROES[phenoHEROES$Status == "AD",]$Concentration)
sd(phenoHEROES[phenoHEROES$Status == "AD",]$Concentration)
range(phenoHEROES[phenoHEROES$Status == "AD",]$Concentration)
# Down syndrome dementia
mean(phenoHEROES[phenoHEROES$Status == "DSD",]$Concentration)
sd(phenoHEROES[phenoHEROES$Status == "DSD",]$Concentration)
range(phenoHEROES[phenoHEROES$Status == "DSD",]$Concentration)

# Test for normal distribution
shapiro.test(residuals(lm(phenoHEROES$Concentration ~ phenoHEROES$Status)))
# Test for Homoscedasticity
bptest(lm(phenoHEROES$Concentration ~ phenoHEROES$Status), studentize = FALSE)
# ANOVA with a post-hoc test based on the studentized range distribution (Tukey's Honest Significant Difference test)
anova_Conc <- aov(phenoHEROES$Concentration ~ phenoHEROES$Status)
summary(anova_Conc)
TukeyHSD(anova_Conc)
# Based on Homoscedasticity, use the Kruskal Wallis test instead with post-hoc Dunn´s test
kruskal.test(phenoHEROES$Concentration ~ phenoHEROES$Status)
dunnTest(Concentration ~ Status, phenoHEROES)

# COUNTS MATRIX -----------------------------------------------------------
# Step 1: Load Counts Matrix ----------------------------------------------
countsHERO <- read.delim("NPJ_dementia/1_Datasets/Counts.matrix.csv", sep = ",", header = TRUE)
row.names(countsHERO) <- countsHERO[,1]
countsHEROES <- countsHERO[,2:21]

# Step 2: Exploration of Counts Data --------------------------------------
# Calculate the reads per million for each sample
sampleHEROES <- apply(countsHEROES, 2, sum)/10^6

# Create a data frame to make a plot
sample_plot <- data.frame(Sample = names(sampleHEROES),
                          Case = phenoHEROES$Case,
                          Total = sampleHEROES,
                          Group = phenoHEROES$Status)

# Factorize Group
sample_plot$Group <- factor(sample_plot$Group, levels = c("Control", "PDD", "AD", "DSD"))

# Total counts (save)
png(filename="NPJ_dementia/3_Figures/1_Data_Exploration/Total_counts.png",
    width     = 30,
    height    = 15,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

sample_plot %>%
  mutate(Sample = fct_relevel(Case, 
                            "#1", "#2", "#3", "#4", "#5", 
                            "#6", "#7", "#8", "#9", 
                            "#10", "#11", "#12", "#13", "#14", "#15",
                            "#16", "#17", "#18", "#19", "#20")) %>%
  ggplot( aes(x = Sample, y = Total, fill = factor(Group))) + 
  geom_bar(stat = "identity", color = "black") +
  scale_fill_discrete(name = "Status", breaks = c("Control", "PDD", "AD", "DSD"), type = c("limegreen", "navyblue", "gold", "darkred")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12, face = "bold"),
        axis.text.y = element_text(vjust = 0.5, hjust=1, size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, units = "cm")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,80)) +
  ylab(expression(bold("# of reads (million)"))) +
  xlab(expression(bold("Case"))) +
  ggtitle("Read counts per sample")

dev.off()

# Step 3: Filtering Counts Matrix -----------------------------------------
# Remove low expressed genes from the Counts Matrix following a specific criteria.
# I will filter by using CPM values rather than Counts, because they account for differences in sequencing depth between samples!!!

# General Counts Matrix (Filtered)
dim(countsHEROES) # 60645 of genes in 20 samples
countsCPM <- apply(countsHEROES, 2, function(x) (x/sum(x)) * 1000000)
keep <- rowSums(countsCPM > (15/min(sampleHEROES))) >= 6 # At least 6 samples (biggest group AD with n = 6) have 15 (i.e., 0.33) reads per gene (as minimal library size = 30.42 Million)
countsFilter <- countsHEROES[keep, ]
dim(countsFilter) # 24186 of genes, 20 samples

# DATA EXPLORATION --------------------------------------------------------
# Step 1: TMM Normalization -----------------------------------------------
# There are situations in which some genes can accumulate high rates of reads. 
# To correct for these imbalance in the counts composition there are methods such as the Trimmed Mean of M-values (TMM), included in the package edgeR. 
# This normalization is suitable for comparing among the samples, for instance when performing sample aggregations.

## General Counts Matrix (Filtered) Normalization using TMM (edgeR package) 
DGEnes <- DGEList(counts = countsFilter)
NormFactor <- calcNormFactors(DGEnes, method = "TMM")
countsTMM <- edgeR::cpm(NormFactor, log = TRUE)
NoLOGcountsTMM <- edgeR::cpm(NormFactor, log = FALSE)

# Quick distribution check with histogram (change number for sample 1-20)
hist(countsTMM[,20], xlab="log2-ratio", main="TMM", col = "orange")

# Counts and Normalized General Counts Matrix (save)
write.csv(countsFilter, "NPJ_dementia/1_Datasets/Counts_Matrices/countsFilter.csv")
write.csv(countsTMM, "NPJ_dementia/1_Datasets/Counts_Matrices/countsTMM.csv")
write.csv(NoLOGcountsTMM, "NPJ_dementia/1_Datasets/Counts_Matrices/NoLOGcountsTMM")

# Step 2: Hierarchical Clustering -----------------------------------------
genesTTM <- countsTMM
row_genesTTM <- t(countsTMM)
dmatrix_genes <- dist(row_genesTTM, method = "euclidean")
cluster_genes <- hclust(dmatrix_genes, method = "ward.D2")
colored <- factor(phenoHEROES$Status[cluster_genes$order], labels = c("limegreen", "navyblue", "gold", "darkred"), level = c("Control", "PDD", "AD", "DSD"))

# Hierarchical Clustering
png(filename="NPJ_dementia/3_Figures/1_Data_Exploration/Dendrogram.png",
    width     = 30,
    height    = 15,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

fviz_dend(cluster_genes, cex = 0.9, cex.axis = 1, cex.lab = 1.1, k = 4,
                            rect = TRUE, rect_lty = 2,
                            rect_fill = FALSE,
                            rect_border = "lancet",
                            lwd = 0.5, 
                            main = "",
                            xlab = "Samples",
                            ylab = "Euclidean",
                            sub = "",
                            k_colors = "lancet",
                            label_cols = colored,
                            color_labels_by_k = FALSE,
                            labels_track_height = 60,
          ggtheme = theme_classic()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(vjust = 0.5, hjust = 1, size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_blank(),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12)) +
  scale_y_continuous() +
  ylab(expression(bold("Euclidean distance"))) +
  ggtitle("Hierarchical clustering")

dev.off()

# Step 3: Principal Component Analysis ------------------------------------
pca <- prcomp(t(genesTTM), scale = TRUE)
scores <- as.matrix(pca$x[, c("PC1", "PC2", "PC3")])
phenoHEROES <- cbind(phenoHEROES, scores)

# Write phenotype data with only the samples to be analyzed and first three principal components (save)
write.csv(phenoHEROES, "NPJ_dementia/1_Datasets/PhenoData/phenoHEROES.csv")

## 2D PCA
png(filename="NPJ_dementia/3_Figures/1_Data_Exploration/PCA_2D.png",
    width     = 20,
    height    = 20,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(data = phenoHEROES, aes(x = PC1, y = PC2, color = factor(Status), label = Case)) + 
  geom_point(alpha = I(0.7), size = 4) + 
  geom_text(hjust = -0.05, vjust = -0.05, size = 7) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlab(paste("PC1 (", round(summary(pca)$importance[2,1], 2) * 100, "%)")) +
  ylab(paste("PC2 (", round(summary(pca)$importance[2,2], 2) * 100, "%)")) +
  scale_color_manual(name = "Status", breaks = c("Control", "PDD", "AD", "DSD"), values=c("limegreen", "navyblue", "gold", "darkred")) +
  stat_ellipse(type = "euclid", level = 10, size = 1.5) + 
  ylim(-150, 150) +
  xlim(-250, 200) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text = element_text(size = 12), 
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        legend.title = element_text(color = "black", size = 14, face = "bold"),
        legend.text = element_text(size = 10)) +
  labs(color = "Group") +
  ggtitle(paste("2D Principal component analysis")) + 
  xlab(expression(bold("PC1 (34.0%)"))) +
  ylab(expression(bold("PC2 (9.9%)")))

dev.off()

## 3D PCA
# Organize colors and legends
Var <- factor(phenoHEROES$Condition, labels = c("limegreen", "navyblue", "gold", "darkred"), 
              level = c("HC_Control", "HC_PDD", "HC_AD", "HC_DSD"))

Legend <- factor(phenoHEROES$Status, level = c("Control", "PDD", "AD", "DS"))

png(filename="NPJ_dementia/3_Figures/1_Data_Exploration/PCA_3D.png",
    width     = 12,
    height    = 10,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

pca3d <- scatterplot3d(x = phenoHEROES$PC1, y = phenoHEROES$PC2, z = phenoHEROES$PC3,
                       cex.axis = 1.2, cex.lab = 1.5,
                       xlim = c(-200, 200),
                       ylim = c(-150, 100),
                       zlim = c(-100, 150),
                       type = "h",
                       xlab = expression(bold("PC1 (34.0%)")), 
                       ylab = expression(bold("PC2 (9.9%)")), 
                       zlab = expression(bold("PC3 (7.9%)")), 
                       main = "3D Principal component analysis", pch = 20, cex.symbols = 8, cex.main = 2, color = Var, col.grid = "grey")
text(pca3d$xyz.convert(pca$x + 8), labels = phenoHEROES$Case, font = 1, cex = 1.2)
legend("topright", legend = levels(Legend),
       col =  c("limegreen", "navyblue", "gold", "darkred"), pch = 20, pt.cex = 4)

dev.off()

## Multidimensional scaling
VarCol <- factor(phenoHEROES$Condition, labels = 1:4, 
                 level = c("HC_Control", "HC_PDD", "HC_AD", "HC_DS"))
varCol <- gsub("1", "limegreen",Var)
varCol <- gsub("2","navyblue",varCol)
varCol <- gsub("3","gold",varCol)
varCol <- gsub("4","darkred",varCol)

png(filename="NPJ_dementia/3_Figures/1_Data_Exploration/MDS.png",
    width     = 10,
    height    = 10,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

plotMDS(NormFactor, labels = phenoHEROES$Case, col = varCol, xlim = c(-3,2), ylim = c(-2, 2), 
        cex.axis = 1, cex.lab = 1.1, cex.main = 3, cex = 1.5,
        frame = FALSE, main = "Multidimensional scaling plot", 
        xlab = "Leading logFC dim1",
        ylab = "Leading logFC dim2")

dev.off()
