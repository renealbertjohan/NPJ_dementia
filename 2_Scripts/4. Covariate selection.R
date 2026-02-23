# 4. Covariate analyses and selection for DESeq2
# Load libraries
library(tidyverse)
library(DT)
library(ppcor)
library(ggcorrplot)
library(calibrate)
library(ggplot2)
library(corrplot)
library(ggfortify) # 2D Principal Component Analysis (2D PCA)

# PHENOTYPE DATA ----------------------------------------------------------
# Step 1: Load the phenotype data of the HEROES project -------------------
phenoHEROES <- read.delim("NPJ_dementia/1_Datasets/PhenoData/phenoHEROES_RNAAge_FPKM.csv", sep = ",", dec = ",", header = TRUE)
row.names(phenoHEROES) <- phenoHEROES$Tube_code
phenoHEROES$X <- NULL

# Factorize Status
phenoHEROES$Status <- factor(phenoHEROES$Status, levels = c("Control", "PDD", "AD", "DSD"))

# Number of Males and Females
table(phenoHEROES$Sex, factor(phenoHEROES$Status, levels = c("Control", "PDD", "AD", "DSD")))

# Step 2: General ---------------------------------------------------------
datatable(phenoHEROES %>% arrange(factor(Status, levels = c("Control", "PDD", "AD", "DSD")), 
                              Sex, ChronAge), options = list(pageLength = -1))

# Make numeric vector for certain column in the phenotype data
phenoHEROES$ChronAge <- as.numeric(phenoHEROES$ChronAge)
phenoHEROES$RNAAge <- as.numeric(phenoHEROES$RNAAge)
phenoHEROES$PMD <- as.numeric(phenoHEROES$PMD)
phenoHEROES$RIN <- as.numeric(phenoHEROES$RIN)
phenoHEROES$Concentration <- as.numeric(phenoHEROES$Concentration)
phenoHEROES$Batch <- as.numeric(phenoHEROES$Batch)
phenoHEROES$PC1 <- as.numeric(phenoHEROES$PC1)
phenoHEROES$PC2 <- as.numeric(phenoHEROES$PC2)
phenoHEROES$PC3 <- as.numeric(phenoHEROES$PC3)

# Select continuous variables as possible covariates
phenoHEROES_select <- phenoHEROES[,c(2, 3, 4, 20, 21, 23, 24, 25)]
phenoHEROES_select <- as.data.frame(append(phenoHEROES_select, list(Sex = ifelse(phenoHEROES$Sex == "Male", 1, 0), 
                                                            APOE4 = ifelse(phenoHEROES$APOE4 == "E4", 1, 0), 
                                                            'Braak stage' = as.numeric(phenoHEROES$AD_Braak),
                                                            'Brain bank' = as.numeric(phenoHEROES$Batch),
                                                            PMI = phenoHEROES$PMD), after = 2))

rownames(phenoHEROES_select) <- phenoHEROES_select[,1]
phenoHEROES_select <- phenoHEROES_select[,2:13]
phenoHEROES_select
# Checking for normal distribution
shapiro.test(residuals(lm(phenoHEROES$ChronAge ~ phenoHEROES$Status)))
shapiro.test(residuals(lm(phenoHEROES$RNAAge ~ phenoHEROES$Status)))
shapiro.test(residuals(lm(phenoHEROES$APOE4 ~ phenoHEROES$Status)))
shapiro.test(residuals(lm(phenoHEROES$AD_Braak ~ phenoHEROES$Status)))
shapiro.test(residuals(lm(phenoHEROES$PMD ~ phenoHEROES$Status)))
shapiro.test(residuals(lm(phenoHEROES$RIN ~ phenoHEROES$Status)))
shapiro.test(residuals(lm(phenoHEROES$Concentration ~ phenoHEROES$Status)))
shapiro.test(residuals(lm(phenoHEROES$Batch ~ phenoHEROES$Status)))
shapiro.test(residuals(lm(phenoHEROES$PC1 ~ phenoHEROES$Status)))
shapiro.test(residuals(lm(phenoHEROES$PC2 ~ phenoHEROES$Status)))
shapiro.test(residuals(lm(phenoHEROES$PC3 ~ phenoHEROES$Status)))

## Compute a correlation matrix (PEARSON)
corr <- cor(phenoHEROES_select, method = "pearson", use = "pairwise.complete.obs")

# Extract the matrix of correlation p-values (PEARSON)
p_mat <- cor_pmat(phenoHEROES_select, method = "pearson", use = "pairwise.complete.obs")

# Correlation plot PEARSON
png(filename = "NPJ_dementia/3_Figures/4_Covariate_selection/Corrplot_Pearson.png",
    width     = 15,
    height    = 15,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

corrplot(corr,
         method = "circle",
         type = "lower",
         p.mat = p_mat,
         insig = "blank",
         cl.cex = 2,
         col= rev(COL2("RdBu", 200)),
         tl.cex = 2,
         sig.level = 0.05,
         title = "Correlation plot",
         cex.main = 5,
         mar=c(0, 0, 9, 0))

dev.off()

## Compute a correlation matrix (KENDALL)
corr <- cor(phenoHEROES_select, method = "kendall", use = "pairwise.complete.obs")

# Extract the matrix of correlation p-values (KENDALL)
p_mat <- cor_pmat(phenoHEROES_select, method = "kendall", use = "pairwise.complete.obs")

# Correlation plot KENDALL
png(filename = "NPJ_dementia/3_Figures/4_Covariate_selection/Corrplot_Kendall.png",
    width     = 15,
    height    = 15,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

corrplot(corr,
         method = "circle",
         type = "lower",
         p.mat = p_mat,
         insig = "blank",
         cl.cex = 2,
         col= rev(COL2("RdBu", 200)),
         tl.cex = 2,
         sig.level = 0.05,
         title = "Correlation plot",
         cex.main = 5,
         mar=c(0, 0, 9, 0))

dev.off()

## Compute a correlation matrix (SPEARMAN)
corr <- cor(phenoHEROES_select, method = "spearman", use = "pairwise.complete.obs")

# Extract the matrix of correlation p-values (SPEARMAN)
p_mat <- cor_pmat(phenoHEROES_select, method = "spearman", use = "pairwise.complete.obs")

# Correlation plot SPEARMAN
png(filename="NPJ_dementia/3_Figures/4_Covariate_selection/Corrplot_Spearman.png",
    width     = 15,
    height    = 15,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

corrplot(corr,
         method = "circle",
         type = "lower",
         p.mat = p_mat,
         insig = "blank",
         cl.cex = 2,
         col= rev(COL2("RdBu", 200)),
         tl.cex = 2,
         sig.level = 0.05,
         title = "Correlation plot",
         cex.main = 5,
         mar=c(0, 0, 9, 0))

dev.off()

## PCA with RIN
countsTMM <- read.delim("NPJ_dementia/1_Datasets/Counts_Matrices/countsTMM.csv", sep = ",", header = TRUE)
row.names(countsTMM ) <- countsTMM [,1]
countsTMM  <- countsTMM [,2:21]
pca <- prcomp(t(countsTMM), scale = TRUE)

phenoHEROES$RINExp <- scale(phenoHEROES$RIN^2)
lims_RINExp <- range(phenoHEROES$RINExp)

png(filename="NPJ_dementia/3_Figures/4_Covariate_selection/PCA_RIN.png",
    width     = 20,
    height    = 20,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(data = phenoHEROES, aes(x = PC1, y = PC2, color = RINExp, label = Case)) + 
  geom_point(alpha = I(0.7), size = 4) + 
  geom_hline(yintercept = 0) +
  geom_text(hjust = -0.05, vjust = -0.05, size = 7) +
  geom_vline(xintercept = 0) +
  xlab(paste("PC1 (", round(summary(pca)$importance[2,1], 2) * 100, "%)")) +
  ylab(paste("PC2 (", round(summary(pca)$importance[2,2], 2) * 100, "%)")) +
  scale_colour_gradient2(limits = lims_RINExp, low = "blue",
                         mid = "lightgrey",
                         high = "red",
                         midpoint = 0,
                         space = "Lab",
                         na.value = "grey50",
                         transform = "identity",
                         guide = "colourbar",
                         aesthetics = "colour") + 
  ylim(-160, 160) +
  xlim(-220, 220) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text = element_text(size = 12), 
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.size = unit(1, units = "cm")) +
  ggtitle(paste("2D Principal component analysis", "\nRNA Intergrity Number")) +
  xlab(expression(bold("PC1 (34.0%)"))) +
  ylab(expression(bold("PC2 (9.9%)")))

dev.off()

## PCA with RNAAge
phenoHEROES$RNAAgeExp <- scale(phenoHEROES$RNAAge^2)
lims_RNAAgeExp <- range(phenoHEROES$RNAAgeExp)

png(filename="NPJ_dementia/3_Figures/4_Covariate_selection/PCA_RNAAge.png",
    width     = 20,
    height    = 20,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(data = phenoHEROES, aes(x = PC1, y = PC2, color = RNAAgeExp, label = Case)) + 
  geom_point(alpha = I(0.7), size = 4) + 
  geom_hline(yintercept = 0) +
  geom_text(hjust = -0.05, vjust = -0.05, size = 7) +
  geom_vline(xintercept = 0) +
  xlab(paste("PC1 (", round(summary(pca)$importance[2,1], 2) * 100, "%)")) +
  ylab(paste("PC2 (", round(summary(pca)$importance[2,2], 2) * 100, "%)")) +
  scale_colour_gradient2(limits = lims_RNAAgeExp, low = "blue",
                         mid = "lightgrey",
                         high = "red",
                         midpoint = 0,
                         space = "Lab",
                         na.value = "grey50",
                         transform = "identity",
                         guide = "colourbar",
                         aesthetics = "colour") +
  ylim(-160, 160) +
  xlim(-220, 220) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text = element_text(size = 12), 
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.size = unit(1, units = "cm")) +
  ggtitle(paste("2D Principal component analysis", "\nTranscriptional age")) + 
  xlab(expression(bold("PC1 (34.0%)"))) +
  ylab(expression(bold("PC2 (9.9%)")))

dev.off()

## PCA with ChronAge
phenoHEROES$ChronAgeExp <- scale(phenoHEROES$ChronAge^2)
lims_ChronAgeExp <- range(phenoHEROES$ChronAgeExp)

png(filename="NPJ_dementia/3_Figures/4_Covariate_selection/PCA_ChronAge.png",
    width     = 20,
    height    = 20,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(data = phenoHEROES, aes(x = PC1, y = PC2, color = ChronAgeExp, label = Case)) + 
  geom_point(alpha = I(0.7), size = 4) + 
  geom_hline(yintercept = 0) +
  geom_text(hjust = -0.05, vjust = -0.05, size = 7) +
  geom_vline(xintercept = 0) +
  xlab(paste("PC1 (", round(summary(pca)$importance[2,1], 2) * 100, "%)")) +
  ylab(paste("PC2 (", round(summary(pca)$importance[2,2], 2) * 100, "%)")) +
  scale_colour_gradient2(limits = lims_ChronAgeExp, low = "blue",
                         mid = "lightgrey",
                         high = "red",
                         midpoint = 0,
                         space = "Lab",
                         na.value = "grey50",
                         transform = "identity",
                         guide = "colourbar",
                         aesthetics = "colour") + 
  ylim(-160, 160) +
  xlim(-220, 220) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text = element_text(size = 12), 
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.size = unit(1, units = "cm")) +
  ggtitle(paste("2D Principal component analysis", "\nChronological age")) + 
  xlab(expression(bold("PC1 (34.0%)"))) +
  ylab(expression(bold("PC2 (9.9%)")))

dev.off()

## PCA with PMD
phenoHEROES$PMDExp <- scale(phenoHEROES$PMD^2)
lims_PMDExp <- range(phenoHEROES$PMDExp)

png(filename="NPJ_dementia/3_Figures/4_Covariate_selection/PCA_PMD.png",
    width     = 20,
    height    = 20,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(data = phenoHEROES, aes(x = PC1, y = PC2, color = PMDExp, label = Case)) + 
  geom_point(alpha = I(0.7), size = 4) + 
  geom_hline(yintercept = 0) +
  geom_text(hjust = -0.05, vjust = -0.05, size = 7) +
  geom_vline(xintercept = 0) +
  xlab(paste("PC1 (", round(summary(pca)$importance[2,1], 2) * 100, "%)")) +
  ylab(paste("PC2 (", round(summary(pca)$importance[2,2], 2) * 100, "%)")) +
  scale_colour_gradient2(limits = lims_PMDExp, low = "blue",
                         mid = "lightgrey",
                         high = "red",
                         midpoint = 0,
                         space = "Lab",
                         na.value = "grey50",
                         transform = "identity",
                         guide = "colourbar",
                         aesthetics = "colour") +
  ylim(-160, 160) +
  xlim(-220, 220) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text = element_text(size = 12), 
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.size = unit(1, units = "cm")) +
  ggtitle(paste("2D Principal component analysis", "\nPost-mortem interval")) + 
  xlab(expression(bold("PC1 (34.0%)"))) +
  ylab(expression(bold("PC2 (9.9%)")))

dev.off()

## PCA with Concentration
phenoHEROES$ConcExp <- scale(phenoHEROES$Concentration^2)
lims_ConcExp <- range(phenoHEROES$ConcExp)

png(filename="NPJ_dementia/3_Figures/4_Covariate_selection/PCA_Concentration.png",
    width     = 20,
    height    = 20,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(data = phenoHEROES, aes(x = PC1, y = PC2, color = ConcExp, label = Case)) + 
  geom_point(alpha = I(0.7), size = 4) + 
  geom_hline(yintercept = 0) +
  geom_text(hjust = -0.05, vjust = -0.05, size = 7) +
  geom_vline(xintercept = 0) +
  xlab(paste("PC1 (", round(summary(pca)$importance[2,1], 2) * 100, "%)")) +
  ylab(paste("PC2 (", round(summary(pca)$importance[2,2], 2) * 100, "%)")) +
  scale_colour_gradient2(limits = lims_ConcExp, low = "blue",
                         mid = "lightgrey",
                         high = "red",
                         midpoint = 0,
                         space = "Lab",
                         na.value = "grey50",
                         transform = "identity",
                         guide = "colourbar",
                         aesthetics = "colour") +
  ylim(-160, 160) +
  xlim(-220, 220) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text = element_text(size = 12), 
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.size = unit(1, units = "cm")) +
  ggtitle(paste("2D Principal component analysis", "\nRNA Concentration")) + 
  xlab(expression(bold("PC1 (34.0%)"))) +
  ylab(expression(bold("PC2 (9.9%)")))

dev.off()

# PCA with APOE4
phenoHEROES$APOE4_1 <- ifelse(phenoHEROES$APOE4 == "E4", "darkred", "navyblue")

png(filename="NPJ_dementia/3_Figures/4_Covariate_selection/PCA_APOE4.png",
    width     = 20,
    height    = 20,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(data = phenoHEROES, aes(x = PC1, y = PC2, color = APOE4, label = Case)) + 
  geom_point(alpha = I(0.7), size = 4) + 
  geom_hline(yintercept = 0) +
  geom_text(hjust = -0.05, vjust = -0.05, size = 7) +
  geom_vline(xintercept = 0) +
  xlab(paste("PC1 (", round(summary(pca)$importance[2,1], 2) * 100, "%)")) +
  ylab(paste("PC2 (", round(summary(pca)$importance[2,2], 2) * 100, "%)")) +
  scale_color_manual(values = c("darkred", "navyblue")) +
  stat_ellipse(type = "euclid", level = 10, size = 1.5) + 
  ylim(-160, 160) +
  xlim(-220, 220) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text = element_text(size = 12), 
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.size = unit(1, units = "cm")) +
  ggtitle(paste("2D Principal component analysis", "\nAPOE Allele")) + 
  xlab(expression(bold("PC1 (34.0%)"))) +
  ylab(expression(bold("PC2 (9.9%)")))

dev.off()
