# 5. Differential expressed genes (DESeq2)
# Load libraries
library(DT) # Create a nice data table
library(tidyverse) # Arrange function
library(DESeq2) # Statistical methods to normalize and analyze RNA-seq data
library(car) # Function vif (Variance Inflation Factor)
library(ggrepel) # Function geom_text_repel() in volcano plot
library(pheatmap) # Creates pretty heatmaps

# Step 1: Load phenotype data
phenoHEROES <- read.delim("NPJ_dementia/1_Datasets/PhenoData/phenoHEROES_RNAAge_FPKM.csv", sep = ",", header = TRUE)
row.names(phenoHEROES) <- phenoHEROES[,1]
phenoHEROES$X <- NULL

# Step 2: Load previously filtered Counts Matrix (see script pre-processing)
countsHEROES <- read.delim("NPJ_dementia/1_Datasets/Counts_Matrices/countsFilter.csv", sep = ",", header = TRUE)
row.names(countsHEROES) <- countsHEROES[,1]
countsHEROES$X <- NULL
dim(countsHEROES) # 24186 of genes, 20 samples

# Check if rownames(phenoHC) and colnames(countsHEROES_HC) are identical, including order
identical(colnames(countsHEROES), rownames(phenoHEROES))

# Step 3: Factorizing and placing covariates in seperate objects
head(phenoHEROES)
class(phenoHEROES)
str(phenoHEROES)
dim(phenoHEROES)

# Make numeric vector for certain column in the phenotype data
phenoHEROES$Sex <- as.factor(phenoHEROES$Sex)
phenoHEROES$RNAAge <- as.numeric(scale(phenoHEROES$RNAAge))
phenoHEROES$ChronAge <- as.numeric(phenoHEROES$ChronAge)
phenoHEROES$RIN <- as.numeric(scale(phenoHEROES$RIN))
phenoHEROES$PMD <- as.numeric(phenoHEROES$PMD)
phenoHEROES$Concentration <- as.numeric(phenoHEROES$Concentration)
phenoHEROES$APOE4 <- as.factor(phenoHEROES$APOE4)
phenoHEROES$Condition <- as.factor(phenoHEROES$Condition)
phenoHEROES$Batch <- as.numeric(phenoHEROES$Batch)
phenoHEROES$PC1 <- as.numeric(phenoHEROES$PC1)
phenoHEROES$PC2 <- as.numeric(phenoHEROES$PC2)
phenoHEROES$PC3 <- as.numeric(phenoHEROES$PC3)

# Select continuous variables as possible covariates
phenoHEROES_select <- phenoHEROES[,c(3, 4, 20, 21, 23, 24, 25)]
phenoHEROES_select <- as.data.frame(append(phenoHEROES_select, list(Sex = ifelse(phenoHEROES$Sex == "Male", 1, 0), 
                                                            APOE4 = ifelse(phenoHEROES$APOE4 == "E4", 1, 0),
                                                            Condition = phenoHEROES$Condition,
                                                            'Brain bank' = as.numeric(phenoHEROES$Batch),
                                                            PMI = phenoHEROES$PMD), after = 2))

phenoHEROES_select$Condition <- as.numeric(phenoHEROES_select$Condition, labels = c(0,1,2,3))
phenoHEROES_select$Sex <- as.numeric(phenoHEROES_select$Sex)
str(phenoHEROES_select)

# Complete model
model <- lm(Condition ~ ., data = phenoHEROES_select)
vif_values <- vif(model) 
print(vif_values) # APOE4 collinearity (i.e., > 5), thus not taken into the model as covariate

# Step 4: DESeq2
dds <- DESeqDataSetFromMatrix(countData = countsHEROES,
                              colData = phenoHEROES,
                              design = ~ Sex + RNAAge + RIN + Condition)

DESEQ_HEROES <- DESeq(dds)
DESEQ_HEROES <- DESEQ_HEROES[which(mcols(DESEQ_HEROES)$betaConv),]
resultsNames(DESEQ_HEROES)

# Step 5: Results
HC_PDD_vs_Control <- as.data.frame(results(DESEQ_HEROES, pAdjustMethod = "fdr", independentFiltering = FALSE, cooksCutoff=FALSE, alpha = 0.05, contrast=c("Condition", "HC_PDD", "HC_Control")))
colnames(HC_PDD_vs_Control)<- c("baseMean","logFC","lfcSE", "stat", "P.Value", "adj.P.Val")
HC_AD_vs_Control <- as.data.frame(results(DESEQ_HEROES, pAdjustMethod = "fdr", independentFiltering = FALSE, cooksCutoff=FALSE, alpha = 0.05, contrast=c("Condition", "HC_AD", "HC_Control")))
colnames(HC_AD_vs_Control)<- c("baseMean","logFC","lfcSE", "stat", "P.Value", "adj.P.Val")
HC_DSD_vs_Control <- as.data.frame(results(DESEQ_HEROES, pAdjustMethod = "fdr", independentFiltering = FALSE, cooksCutoff=FALSE, alpha = 0.05, contrast=c("Condition", "HC_DSD", "HC_Control")))
colnames(HC_DSD_vs_Control)<- c("baseMean","logFC","lfcSE", "stat", "P.Value", "adj.P.Val")

# Organize the TopTables
HC_PDD_vs_Control <- HC_PDD_vs_Control[order(HC_PDD_vs_Control$logFC, decreasing = TRUE), ]
HC_AD_vs_Control <- HC_AD_vs_Control[order(HC_AD_vs_Control$logFC, decreasing = TRUE), ]
HC_DSD_vs_Control <- HC_DSD_vs_Control[order(HC_DSD_vs_Control$logFC, decreasing = TRUE), ]

# Write the results from DESeq2 (save)
write.csv(HC_PDD_vs_Control, "NPJ_dementia/1_Datasets/DEG_results/HC_PDD_vs_Control.csv")
write.csv(HC_AD_vs_Control, "NPJ_dementia/1_Datasets/DEG_results/HC_AD_vs_Control.csv")
write.csv(HC_DSD_vs_Control, "NPJ_dementia/1_Datasets/DEG_results/HC_DSD_vs_Control.csv")

# Step 6: Extract information per dementia --------------------------------
## Parkinson's disease dementia
HC_PDD_vs_Control$adj.P.Val < 0.05
HC_PDD_vs_Control$diffexpressed <- "NO"
HC_PDD_vs_Control[HC_PDD_vs_Control$P.Value <= 0.05 & HC_PDD_vs_Control$logFC >= 0.3,]$diffexpressed <- "UP"
HC_PDD_vs_Control[HC_PDD_vs_Control$P.Value <= 0.05 & HC_PDD_vs_Control$logFC <= -0.3,]$diffexpressed <- "DOWN"
HC_PDD_vs_Control$ID <- rownames(HC_PDD_vs_Control)
HC_PDD_vs_Control_adj_p_log <- HC_PDD_vs_Control[HC_PDD_vs_Control$adj.P.Val <= 0.05 & abs(HC_PDD_vs_Control$logFC) >= 0.3 & !is.na(HC_PDD_vs_Control$adj.P.Val),]
dim(HC_PDD_vs_Control_adj_p_log) # 2897 DEGs

# Upregulated --> 1296 genes
dim(HC_PDD_vs_Control[HC_PDD_vs_Control$adj.P.Val <= 0.05 & HC_PDD_vs_Control$logFC >= 0.3 & !is.na(HC_PDD_vs_Control$adj.P.Val),])
# Downregulated --> 1601 genes
dim(HC_PDD_vs_Control[HC_PDD_vs_Control$adj.P.Val <= 0.05 & HC_PDD_vs_Control$logFC <= 0.3 & !is.na(HC_PDD_vs_Control$adj.P.Val),])

# Extract gene data from counts_length object (data obtained in RNAAgeCalc scripts through biomaRt search)
genedata <- read.delim("NPJ_dementia/1_Datasets/Counts_Matrices/counts_length.csv", sep = ",", header = TRUE)
genedata$X.1 <- NULL
colnames(genedata)[1] <- "Geneid" 
genedata <- genedata[,-c(2:21)] 

# Extract row of genedata if it is in row HC_PDD_vs_Healthy_adj_p_log
PDD_HC_Table <- genedata[genedata$Geneid %in% HC_PDD_vs_Control_adj_p_log$ID,]
rownames(PDD_HC_Table) <- PDD_HC_Table$Geneid
PDD_HC_Table <- PDD_HC_Table[HC_PDD_vs_Control_adj_p_log$ID,]

colnames(HC_PDD_vs_Control_adj_p_log)
PDD_HC <- cbind(PDD_HC_Table[,c(4)], HC_PDD_vs_Control_adj_p_log[, c(1,2,3,4,5,6,7)], PDD_HC_Table[,c(5,6,7,8,10,11,9)])

colnames(PDD_HC)[1] <- "Chr"
colnames(PDD_HC)[8] <- "DifExpressed"
colnames(PDD_HC)[15] <- "Description"
datatable(PDD_HC %>% arrange(factor(Chr, levels = c("1", "2", "3", "4", "5", "6", "7", "8",
                                                    "9", "10", "11", "12", "13", "14", "15",
                                                    "16", "17", "18", "19", "20", "21", "22",
                                                    "X", "Y", "MT"))), options = list(pageLength = -1))

PDD_vs_control_HC <- PDD_HC %>% arrange(factor(Chr, levels = c("1", "2", "3", "4", "5", "6", "7", "8",
                                                               "9", "10", "11", "12", "13", "14", "15",
                                                               "16", "17", "18", "19", "20", "21", "22",
                                                               "X", "Y", "MT")))

# Parkinson's disease demented DEGs with genedata 
write.csv(PDD_vs_control_HC, "NPJ_dementia/1_Datasets/DEG_results/PDD_vs_control_with_Genedata.csv")

## Alzheimer's disease
HC_AD_vs_Control$adj.P.Val < 0.05
HC_AD_vs_Control$diffexpressed <- "NO"
HC_AD_vs_Control[HC_AD_vs_Control$P.Value <= 0.05 & HC_AD_vs_Control$logFC >= 0.3,]$diffexpressed <- "UP"
HC_AD_vs_Control[HC_AD_vs_Control$P.Value <= 0.05 & HC_AD_vs_Control$logFC <= -0.3,]$diffexpressed <- "DOWN"
HC_AD_vs_Control$ID <- rownames(HC_AD_vs_Control)
HC_AD_vs_Control_adj_p_log <- HC_AD_vs_Control[HC_AD_vs_Control$adj.P.Val <= 0.05 & abs(HC_AD_vs_Control$logFC) >= 0.3 & !is.na(HC_AD_vs_Control$adj.P.Val),]
dim(HC_AD_vs_Control_adj_p_log) #657 DEGs

# Upregulated --> 357 genes
dim(HC_AD_vs_Control[HC_AD_vs_Control$adj.P.Val <= 0.05 & HC_AD_vs_Control$logFC >= 0.3 & !is.na(HC_AD_vs_Control$adj.P.Val),])
# Downregulated --> 300 genes
dim(HC_AD_vs_Control[HC_AD_vs_Control$adj.P.Val <= 0.05 & HC_AD_vs_Control$logFC <= 0.3 & !is.na(HC_AD_vs_Control$adj.P.Val),])

# Extract row of genedata if it is in row HC_AD_vs_Healthy_adj_p_log
AD_HC_Table <- genedata[genedata$Geneid %in% HC_AD_vs_Control_adj_p_log$ID,]
rownames(AD_HC_Table) <- AD_HC_Table$Geneid
AD_HC_Table <- AD_HC_Table[HC_AD_vs_Control_adj_p_log$ID,]
AD_HC <- cbind(AD_HC_Table[,c(4)], HC_AD_vs_Control_adj_p_log[, c(1,2,3,4,5,6,7)], AD_HC_Table[,c(5,6,7,8,10,11,9)])
colnames(AD_HC)[1] <- "Chr"
colnames(AD_HC)[8] <- "DifExpressed"
colnames(AD_HC)[15] <- "Description"
datatable(AD_HC %>% arrange(factor(Chr, levels = c("1", "2", "3", "4", "5", "6", "7", "8",
                                                   "9", "10", "11", "12", "13", "14", "15",
                                                   "16", "17", "18", "19", "20", "21", "22",
                                                   "X", "Y", "MT"))), options = list(pageLength = -1))

AD_vs_control_HC <- AD_HC %>% arrange(factor(Chr, levels = c("1", "2", "3", "4", "5", "6", "7", "8",
                                                             "9", "10", "11", "12", "13", "14", "15",
                                                             "16", "17", "18", "19", "20", "21", "22",
                                                             "X", "Y", "MT")))

write.csv(AD_vs_control_HC, "NPJ_dementia/1_Datasets/DEG_results/AD_vs_control_with_Genedata.csv")

## Down syndrome dementia
HC_DSD_vs_Control$adj.P.Val < 0.05
HC_DSD_vs_Control$diffexpressed <- "NO"
HC_DSD_vs_Control[HC_DSD_vs_Control$P.Value <= 0.05 & HC_DSD_vs_Control$logFC >= 0.3,]$diffexpressed <- "UP"
HC_DSD_vs_Control[HC_DSD_vs_Control$P.Value <= 0.05 & HC_DSD_vs_Control$logFC <= -0.3,]$diffexpressed <- "DOWN"
HC_DSD_vs_Control$ID <- rownames(HC_DSD_vs_Control)
HC_DSD_vs_Control_adj_p_log <- HC_DSD_vs_Control[HC_DSD_vs_Control$adj.P.Val <= 0.05 & abs(HC_DSD_vs_Control$logFC) >= 0.3 & !is.na(HC_DSD_vs_Control$adj.P.Val),]
dim(HC_DSD_vs_Control_adj_p_log) # 196 DEGs

# Upregulated --> 133 genes
dim(HC_DSD_vs_Control[HC_DSD_vs_Control$adj.P.Val <= 0.05 & HC_DSD_vs_Control$logFC >= 0.3 & !is.na(HC_DSD_vs_Control$adj.P.Val),])
# Downregulated --> 63 genes
dim(HC_DSD_vs_Control[HC_DSD_vs_Control$adj.P.Val <= 0.05 & HC_DSD_vs_Control$logFC <= 0.3 & !is.na(HC_DSD_vs_Control$adj.P.Val),])

# Extract row of genedata if it is in row HC_AD_vs_Healthy_adj_p_log
DSD_HC_Table <- genedata[genedata$Geneid %in% HC_DSD_vs_Control_adj_p_log$ID,]
rownames(DSD_HC_Table) <- DSD_HC_Table$Geneid
DSD_HC_Table <- DSD_HC_Table[HC_DSD_vs_Control_adj_p_log$ID,]
DSD_HC <- cbind(DSD_HC_Table[,c(4)], HC_DSD_vs_Control_adj_p_log[, c(1,2,3,4,5,6,7)], DSD_HC_Table[,c(5,6,7,8,10,11,9)])
colnames(DSD_HC)[1] <- "Chr"
colnames(DSD_HC)[8] <- "DifExpressed"
colnames(DSD_HC)[15] <- "Description"
datatable(DSD_HC %>% arrange(factor(Chr, levels = c("1", "2", "3", "4", "5", "6", "7", "8",
                                                    "9", "10", "11", "12", "13", "14", "15",
                                                    "16", "17", "18", "19", "20", "21", "22",
                                                    "X", "Y", "MT"))), options = list(pageLength = -1))

DSD_vs_control_HC <- DSD_HC %>% arrange(factor(Chr, levels = c("1", "2", "3", "4", "5", "6", "7", "8",
                                                               "9", "10", "11", "12", "13", "14", "15",
                                                               "16", "17", "18", "19", "20", "21", "22",
                                                               "X", "Y", "MT")))

write.csv(DSD_vs_control_HC, "NPJ_dementia/1_Datasets/DEG_results/DSD_vs_control_with_Genedata.csv")

# Step 7: Volcano plots ---------------------------------------------------
## Volcano plot: Parkinson's disease demented
png(filename = "NPJ_dementia/3_Figures/5_DEG/Volcano_PDD_vs_Control.png",
    width     = 25,
    height    = 20,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(data = HC_PDD_vs_Control, aes(x = logFC, y = -log10(P.Value), col = diffexpressed, label = ID)) + 
  geom_vline(xintercept = c(-0.3, 0.3), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 1.5) +
  theme_classic() +
  ylim(0, 15) +
  xlim(-8, 8) +
  theme(strip.text = element_text(face = "bold", size = 10, color = "black"),
        strip.background = element_rect(fill = "lightgray", colour = "black", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
        axis.text = element_text(size = 18), 
        axis.title.y = element_text(size = 30, face = "bold"),
        axis.title.x = element_text(size = 30, face = "bold"),
        legend.position = "none", 
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16)) +
  labs(color = "Gene expression", x = expression(bold("log"[2]*"FC")), 
       y = expression(bold("-log"[10]*"(p-value)"))) +
  scale_color_manual(values = c("blue", "lightgrey", "red"), 
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  ggtitle(expression(bold(paste("PDD ",bolditalic("versus"), " Control")))) +
  geom_text_repel(data = filter(HC_PDD_vs_Control, logFC > 0.3 & adj.P.Val < 0.00015 | logFC > 4 & adj.P.Val < 0.05), aes(label = ID), size = 4.5, col = "black", max.overlaps = 200, box.padding = 0.5) +
  geom_text_repel(data = filter(HC_PDD_vs_Control, logFC < -0.3 & adj.P.Val < 0.0008 | logFC < -3 & adj.P.Val < 0.2), aes(label = ID), size = 4.5, col = "black", max.overlaps = 200, box.padding = 0.5)

dev.off()

## Volcano plot: Alzheimer's disease
png(filename = "NPJ_dementia/3_Figures/5_DEG/Volcano_AD_vs_Control.png",
    width     = 25,
    height    = 20,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(data = HC_AD_vs_Control, aes(x = logFC, y = -log10(P.Value), col = diffexpressed, label = ID)) + 
  geom_vline(xintercept = c(-0.3, 0.3), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 1.5) +
  theme_classic() +
  ylim(0, 15) +
  xlim(-8, 8) +
  theme(strip.text = element_text(face = "bold", size = 10, color = "black"),
        strip.background = element_rect(fill = "lightgray", colour = "black", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
        axis.text = element_text(size = 18), 
        axis.title.y = element_text(size = 30, face = "bold"),
        axis.title.x = element_text(size = 30, face = "bold"),
        legend.position = "none", 
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16)) +
  labs(color = "Gene expression", x = expression(bold("log"[2]*"FC")), 
       y = expression(bold("-log"[10]*"(p-value)"))) +
  scale_color_manual(values = c("blue", "lightgrey", "red"), 
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  ggtitle(expression(bold(paste("AD ",bolditalic("versus"), " Control")))) +
  geom_text_repel(data = filter(HC_AD_vs_Control, logFC > 0.3 & adj.P.Val < 0.002 | logFC > 3 & adj.P.Val < 0.05), aes(label = ID), size = 4.5, col = "black", max.overlaps = 200, box.padding = 0.5) +
  geom_text_repel(data = filter(HC_AD_vs_Control, logFC < -0.3 & adj.P.Val < 0.0025| logFC < -3 & adj.P.Val < 0.05), aes(label = ID), size = 4.5, col = "black", max.overlaps = 200, box.padding = 0.5)

dev.off()

## Volcano plot: Down syndrome demented
png(filename = "NPJ_dementia/3_Figures/5_DEG/Volcano_DSD_vs_Control.png",
    width     = 25,
    height    = 20,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(data = HC_DSD_vs_Control, aes(x = logFC, y = -log10(P.Value), col = diffexpressed, label = ID)) + 
  geom_vline(xintercept = c(-0.3, 0.3), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 1.5) +
  theme_classic() +
  ylim(0, 15) +
  xlim(-8, 8) +
  theme(strip.text = element_text(face = "bold", size = 10, color = "black"),
        strip.background = element_rect(fill = "lightgray", colour = "black", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
        axis.text = element_text(size = 18), 
        axis.title.y = element_text(size = 30, face = "bold"),
        axis.title.x = element_text(size = 30, face = "bold"),
        legend.position = "none", 
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16)) +
  labs(color = "Gene expression", x = expression(bold("log"[2]*"FC")), 
       y = expression(bold("-log"[10]*"(p-value)"))) +
  scale_color_manual(values = c("blue", "lightgrey", "red"), 
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  ggtitle(expression(bold(paste("DSD ",bolditalic("versus"), " Control")))) +
  geom_text_repel(data = filter(HC_DSD_vs_Control, logFC > 0.3 & adj.P.Val < 0.003 | logFC > 3 & adj.P.Val < 0.05), aes(label = ID), size = 5, col = "black", max.overlaps = 200, box.padding = 0.5) +
  geom_text_repel(data = filter(HC_DSD_vs_Control, logFC < -0.3 & adj.P.Val < 0.003 | logFC < -3 & adj.P.Val < 0.05), aes(label = ID), size = 5, col = "black", max.overlaps = 200, box.padding = 0.5)

dev.off()

# Step 8: Heatmaps --------------------------------------------------------
# Load countsTMM
countsTMM <- read.delim("NPJ_dementia/1_Datasets/Counts_Matrices/countsTMM.csv", sep = ",", header = TRUE)
row.names(countsTMM) <- countsTMM[,1]
countsTMM$X <- NULL

## Heatmap: Parkinson's disease demented
clus <- phenoHEROES[phenoHEROES$Status == "PDD" | phenoHEROES$Status == "Control", c("Tube_code", "Case", "Condition", "Sex")]
clusz <- clus[clus$Condition == "HC_PDD" | clus$Condition == "HC_Control",]
clusG <- clusz$Sex
clusz <- clusz$Tube_code
HC_PDD_vs_Control_adj_p_log
data.clus <- countsTMM[rownames(HC_PDD_vs_Control_adj_p_log), grep("HC_", colnames(countsTMM))]
data.clus <- data.clus[, colnames(data.clus) %in% clusz]

Condition <- as.factor(phenoHEROES$Condition)
Condition <- Condition[grep("HC_", Condition)]
cond.df <- as.data.frame(Condition)
Condition <- cond.df[cond.df$Condition == "HC_PDD" | cond.df$Condition == "HC_Control",]
Condition <- sub("HC_", "", Condition)

annotation <- data.frame(Condition = factor(Condition), Sex = as.factor(clusG))
annotation_colors = list(
  Condition = c(Control = "limegreen", PDD = "navyblue"),
  Sex = c(Male = "lightblue", Female = "pink"))
rownames(annotation) <- clus$Case
colnames(data.clus) <- clus$Case

## Parkinson's disease demented
png(filename="NPJ_dementia/3_Figures/5_DEG/Heatmap_PDD.png",
    width     = 20,
    height    = 30,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

pheatmap(data.clus, scale = "row", show_rownames = FALSE, clustering_distance_cols = "euclidean", clustering_method = "ward.D2",
         main = expression(bold(paste("PDD ",bolditalic("versus"), " Control"))),
         treeheight_col = 20, treeheight_row = 20, border_color = NA,
         fontsize = 20, fontsize_col = 20, annotation = annotation, 
         annotation_colors = annotation_colors)

dev.off()

## Heatmap: Alzheimer's disease
clusAD <- phenoHEROES[phenoHEROES$Status == "AD" | phenoHEROES$Status == "Control", c("Tube_code", "Case", "Condition", "Sex")]
clusAD_2 <- clusAD[clusAD$Condition == "HC_AD" | clusAD$Condition == "HC_Control",] # clusz
clusAD_3 <- clusAD_2$Sex # clusG
clusAD_2 <- clusAD_2$Tube_code # clusz

data_clus_AD <- countsTMM[rownames(HC_AD_vs_Control_adj_p_log), grep("HC_", colnames(countsTMM))]
data_clus_AD <- data_clus_AD[, colnames(data_clus_AD) %in% clusAD_2]

Condition_AD <- as.factor(phenoHEROES$Condition)
Condition_AD <- Condition_AD[grep("HC_", Condition_AD)]
cond_df_AD <- as.data.frame(Condition_AD)
Condition_AD <- cond_df_AD[cond_df_AD$Condition_AD == "HC_AD" | cond_df_AD$Condition_AD == "HC_Control",]
Condition_AD <- sub("HC_", "", Condition_AD)

annotation_AD <- data.frame(Condition = factor(Condition_AD), Sex = as.factor(clusAD_3))
annotation_AD
annotation_colors_AD = list(
  Condition = c(Control="limegreen", AD="gold"),
  Sex = c(Male="lightblue", Female="pink"))
rownames(annotation_AD) <- clusAD$Case
colnames(data_clus_AD) <- clusAD$Case

## Alzheimer's disease
png(filename="NPJ_dementia/3_Figures/5_DEG/Heatmap_AD.png",
    width     = 20,
    height    = 30,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

pheatmap(data_clus_AD, scale = "row", show_rownames = FALSE, clustering_distance_cols = "euclidean", clustering_method = "ward.D2",
                          main = expression(bold(paste("AD ", bolditalic("versus"), " Control"))),
                          treeheight_col = 20, treeheight_row = 20, border_color = NA,
                          fontsize = 20, fontsize_col = 20, annotation = annotation_AD, annotation_colors = annotation_colors_AD)

dev.off()

## Heatmap: Down syndrome demented
clusDSD <- phenoHEROES[phenoHEROES$Status == "DSD" | phenoHEROES$Status == "Control", c("Tube_code", "Case", "Condition", "Sex")]
clusDSD_2 <- clusDSD[clusDSD$Condition == "HC_DSD" | clusDSD$Condition == "HC_Control",]
clusDSD_3 <- clusDSD_2$Sex
clusDSD_2 <- clusDSD_2$Tube_code

data.clus_DSD <- countsTMM[rownames(HC_DSD_vs_Control_adj_p_log), grep("HC_", colnames(countsTMM))]
data.clus_DSD <- data.clus_DSD[, colnames(data.clus_DSD) %in% clusDSD_2]

Condition_DSD <- as.factor(phenoHEROES$Condition)
Condition_DSD <- Condition_DSD[grep("HC_", Condition_DSD)]
cond.df_DSD <- as.data.frame(Condition_DSD)
Condition_DSD <- cond.df_DSD[cond.df_DSD$Condition == "HC_DSD" | cond.df_DSD$Condition == "HC_Control",]
Condition_DSD <- sub("HC_", "", Condition_DSD)

annotation_DSD <- data.frame(Condition = factor(Condition_DSD), Sex = as.factor(clusDSD_3))
annotation_colors_DSD = list(
  Condition = c(Control = "limegreen", DSD = "darkred"),
  Sex = c(Male = "lightblue", Female = "pink"))
rownames(annotation_DSD) <- clusDSD$Case
colnames(data.clus_DSD) <- clusDSD$Case

## Down syndrome demented 
png(filename="NPJ_dementia/3_Figures/5_DEG/Heatmap_DSD.png",
    width     = 20,
    height    = 30,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

pheatmap(data.clus_DSD, scale = "row", show_rownames = FALSE, clustering_distance_cols = "euclidean", clustering_method = "ward.D2",
         main = expression(bold(paste("DSD ", bolditalic("versus"), " Control"))), 
         treeheight_col = 20, treeheight_row = 20, border_color = NA,
         fontsize = 20, fontsize_col = 20, annotation = annotation_DSD, annotation_colors = annotation_colors_DSD)


dev.off()

# Step 9: Number of gene on chromosomes -----------------------------------
## Chromosome: Parkinson's disease demented
# Already done, but extracted again to be sure about the conditions (p-value adjusted < 0.5 & log2FC > 0.3 & < 0.3)
HC_PDD_vs_Control_adj_p_log <- HC_PDD_vs_Control[HC_PDD_vs_Control$adj.P.Val <= 0.05 & abs(HC_PDD_vs_Control$logFC) >= 0.3,]

# Extract up-regulated genes and gene information data (i.e., chromosomal location)
UP_HC_PDD <- HC_PDD_vs_Control_adj_p_log[HC_PDD_vs_Control_adj_p_log$diffexpressed == "UP",]
UP_HC_PDD <- genedata[genedata$Geneid %in% UP_HC_PDD$ID,]
UP <- UP_HC_PDD$chromosome_name
UP <- as.data.frame(UP)
UP$category <- rep("Up", rep = length(UP))

# Extract down-regulated genes and gene information data (i.e., chromosomal location)
DOWN_HC_PDD <- HC_PDD_vs_Control_adj_p_log[HC_PDD_vs_Control_adj_p_log$diffexpressed == "DOWN",]
DOWN_HC_PDD <- genedata[genedata$Geneid %in% DOWN_HC_PDD$ID,]
DOWN <- DOWN_HC_PDD$chromosome_name
DOWN <- as.data.frame(DOWN)
DOWN$category <- rep("Down", rep = length(DOWN))

# Merge extracted information in one table for ggplot
UP <- as.data.frame(table(UP))
DOWN <- as.data.frame(table(DOWN))
tble <- merge(UP, DOWN, by = 1, all = TRUE)
colnames(tble) <- c("Chr", "Exp", "Freq", "Exp", "Freq")

tble <- bind_rows(tble[,1:3], tble[,c(1,4,5) ])
tble <- tble[order(as.numeric(as.character(tble$Chr))), ]
tble <- tble[!is.na(tble$Exp), ]

# Chromosome Parkinson's disease demented
png(filename="NPJ_dementia/3_Figures/5_DEG/Chromosome_PDD.png",
    width     = 50,
    height    = 15,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(tble, aes(fill = Exp, y = Freq, x = Chr), alpha = 0.8, size = 1.2) +
  geom_bar(position = "stack", stat = "identity", width = 0.9, color = "black") +
  scale_x_discrete(limits=c("1","2","3","4","5","6","7","8","9","10","11","12",
                            "13","14","15","16","17","18","19","20", "21","22", "X", "Y", "MT")) +
  scale_fill_manual(values = c("red", "blue"), labels = c("Upregulated", "Downregulated"),
                    name = "Gene expression") +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, units = "cm")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 300)) +
  labs(y = expression(bold("# DE genes")), 
       x = expression(bold("Chromosome"))) +
  ggtitle(expression(bold(paste("Parkinson's disease demented ",bolditalic("versus"), " Controls")))) 

dev.off()

## Chromosome: Alzheimer's disease
# Already done, but extracted again to be sure about the conditions (p-value adjusted < 0.5 & log2FC > 0.3 & < 0.3)
UP_HC_AD <- HC_AD_vs_Control_adj_p_log[HC_AD_vs_Control_adj_p_log$diffexpressed == "UP",]

# Extract up-regulated genes and gene information data (i.e., chromosomal location)
UP_HC_AD <- genedata[genedata$Geneid %in% UP_HC_AD$ID,]
UP <- UP_HC_AD$chromosome_name
UP <- as.data.frame(UP)
UP$category <- rep("Up", rep = length(UP))

# Extract down-regulated genes and gene information data (i.e., chromosomal location)
DOWN_HC_AD <- HC_AD_vs_Control_adj_p_log[HC_AD_vs_Control_adj_p_log$diffexpressed == "DOWN",]
DOWN_HC_AD <- genedata[genedata$Geneid %in% DOWN_HC_AD$ID,]
DOWN <- DOWN_HC_AD$chromosome_name
DOWN <- as.data.frame(DOWN)
DOWN$category <- rep("Down", rep = length(DOWN))

# Merge extracted information in one table for ggplot
UP <- as.data.frame(table(UP))
DOWN <- as.data.frame(table(DOWN))
tble <- merge(UP, DOWN, by = 1, all = TRUE)
colnames(tble) <- c("Chr", "Exp", "Freq", "Exp", "Freq")

tble <- bind_rows(tble[,1:3], tble[,c(1,4,5) ])
tble <- tble[!is.na(tble$Exp), ]
tble <- tble[order(as.numeric(as.character(tble$Chr))), ]

# Chromosome Alzheimer's disease
png(filename="NPJ_dementia/3_Figures/5_DEG/Chromosome_AD.png",
    width     = 50,
    height    = 15,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(tble, aes(fill = Exp, y = Freq, x = Chr), alpha = 0.8, size = 1.2) +
  geom_bar(position = "stack", stat = "identity", width = 0.9, color = "black") +
  scale_x_discrete(limits=c("1","2","3","4","5","6","7","8","9","10","11","12",
                            "13","14","15","16","17","18","19","20", "21","22", "X", "Y", "MT")) +
  scale_fill_manual(values = c("red", "blue"), labels = c("Upregulated", "Downregulated"),
                    name = "Gene expression") +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, units = "cm")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 300)) +
  labs(y = expression(bold("# DE genes")), 
       x = expression(bold("Chromosome"))) +
  ggtitle(expression(bold(paste("Alzheimer's disease ",bolditalic("versus"), " Controls"))))

dev.off()

## Chromosome: Down syndrome demented
# Already done, but extracted again to be sure about the conditions (p-value adjusted < 0.5 & log2FC > 0.3 & < 0.3)
HC_DSD_vs_Control_adj_p_log <- HC_DSD_vs_Control[HC_DSD_vs_Control$adj.P.Val <= 0.05 & abs(HC_DSD_vs_Control$logFC) >= 0.3,]

# Extract up-regulated genes and gene information data (i.e., chromosomal location)
UP_HC_DSD <- HC_DSD_vs_Control_adj_p_log[HC_DSD_vs_Control_adj_p_log$diffexpressed == "UP",]
UP_HC_DSD <- genedata[genedata$Geneid %in% UP_HC_DSD$ID,]
UP <- UP_HC_DSD$chromosome_name
UP <- as.data.frame(UP)
UP$category <- rep("Up", rep = length(UP))

# Extract down-regulated genes and gene information data (i.e., chromosomal location)
DOWN_HC_DSD <- HC_DSD_vs_Control_adj_p_log[HC_DSD_vs_Control_adj_p_log$diffexpressed == "DOWN",]
DOWN_HC_DSD <- genedata[genedata$Geneid %in% DOWN_HC_DSD$ID,]
DOWN <- DOWN_HC_DSD$chromosome_name
DOWN <- as.data.frame(DOWN)
DOWN$category <- rep("Down", rep = length(DOWN))

# Merge extracted information in one table for ggplot
UP <- as.data.frame(table(UP))
DOWN <- as.data.frame(table(DOWN))
tble <- merge(UP, DOWN, by = 1, all = TRUE)
colnames(tble) <- c("Chr", "Exp", "Freq", "Exp", "Freq")

tble <- bind_rows(tble[,1:3], tble[,c(1,4,5) ])
tble <- tble[order(as.numeric(as.character(tble$Chr))), ]
tble <- tble[!is.na(tble$Exp), ]

# Chromosome Down syndrome demented
png(filename="NPJ_dementia/3_Figures/5_DEG/Chromosome_DSD.png",
    width     = 50,
    height    = 15,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(tble, aes(fill = Exp, y = Freq, x = Chr), alpha = 0.8, size = 1.2) +
  geom_bar(position = "stack", stat = "identity", width = 0.9, color = "black") +
  scale_x_discrete(limits=c("1","2","3","4","5","6","7","8","9","10","11","12",
                            "13","14","15","16","17","18","19", "20", "21","22", "X", "Y", "MT")) +
  scale_fill_manual(values = c("red", "blue"), labels = c("Upregulated", "Downregulated"),
                    name = "Gene expression") +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, units = "cm")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 300)) +
  labs(y = expression(bold("# DE genes")), 
       x = expression(bold("Chromosome"))) +
  ggtitle(expression(bold(paste("Down syndrome demented ", bolditalic("versus"), " Controls")))) 

dev.off()
