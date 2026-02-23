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
class(phenoHEROES_select$Condition)
str(phenoHEROES_select)

# Complete model
model <- lm(Condition ~ ., data = phenoHEROES_select)
vif_values <- vif(model) 
print(vif_values) # APOE4 collinearity (i.e., > 5), thus not taken into the model as covariate

# Step 4: DESeq2
# Create demented and non-demented groups first
phenoHEROES <- phenoHEROES %>%
  mutate(Dementia = ifelse(grepl("Control", Status), "Non_Demented", "Demented"))
phenoHEROES$Dementia <- as.factor(phenoHEROES$Dementia)

dds <- DESeqDataSetFromMatrix(countData = countsHEROES,
                              colData = phenoHEROES,
                              design = ~ Sex + RNAAge + RIN + Dementia)

DESEQ_HEROES <- DESeq(dds)
DESEQ_HEROES <- DESEQ_HEROES[which(mcols(DESEQ_HEROES)$betaConv),]
resultsNames(DESEQ_HEROES)

# Step 5: Results
Demented_vs_Non_Demented <- as.data.frame(results(DESEQ_HEROES, pAdjustMethod = "fdr", independentFiltering = FALSE, cooksCutoff=FALSE, alpha = 0.05, contrast=c("Dementia", "Demented", "Non_Demented")))
colnames(Demented_vs_Non_Demented)<- c("baseMean","logFC","lfcSE", "stat", "P.Value", "adj.P.Val")

# Organize the TopTables
Demented_vs_Non_Demented <- Demented_vs_Non_Demented[order(Demented_vs_Non_Demented$logFC, decreasing = TRUE), ]
Demented_vs_Non_Demented

# Write the results from DESeq2 (save)
write.csv(Demented_vs_Non_Demented, "NPJ_dementia/1_Datasets/DEG_results/Demented_vs_Non_Demented.csv")

# Step 6: Extract information per dementia --------------------------------
# Parkinson's disease dementia
Demented_vs_Non_Demented$adj.P.Val < 0.05
Demented_vs_Non_Demented$diffexpressed <- "NO"
Demented_vs_Non_Demented[Demented_vs_Non_Demented$P.Value <= 0.05 & Demented_vs_Non_Demented$logFC >= 0.3,]$diffexpressed <- "UP"
Demented_vs_Non_Demented[Demented_vs_Non_Demented$P.Value <= 0.05 & Demented_vs_Non_Demented$logFC <= -0.3,]$diffexpressed <- "DOWN"
Demented_vs_Non_Demented$ID <- rownames(Demented_vs_Non_Demented)
Demented_vs_Non_Demented_adj_p_log <- Demented_vs_Non_Demented[Demented_vs_Non_Demented$adj.P.Val <= 0.05 & abs(Demented_vs_Non_Demented$logFC) >= 0.3 & !is.na(Demented_vs_Non_Demented$adj.P.Val),]
dim(Demented_vs_Non_Demented_adj_p_log) # 894 DEGs

# Upregulated --> 505 genes
dim(Demented_vs_Non_Demented[Demented_vs_Non_Demented$adj.P.Val <= 0.05 & Demented_vs_Non_Demented$logFC >= 0.3 & !is.na(Demented_vs_Non_Demented$adj.P.Val),])
# Downregulated --> 395 genes
dim(Demented_vs_Non_Demented[Demented_vs_Non_Demented$adj.P.Val <= 0.05 & Demented_vs_Non_Demented$logFC <= 0.3 & !is.na(Demented_vs_Non_Demented$adj.P.Val),])

# Extract gene data from counts_length object (data obtained in RNAAgeCalc scripts through biomaRt search)
genedata <- read.delim("NPJ_dementia/1_Datasets/Counts_Matrices/counts_length.csv", sep = ",", header = TRUE)
genedata$X.1 <- NULL
colnames(genedata)[1] <- "Geneid" 
genedata <- genedata[,-c(2:21)] 

# Extract row of genedata if it is in row Demented_vs_Non_Demented_adj_p_log
Demented_vs_Non_Demented_Table <- genedata[genedata$Geneid %in% Demented_vs_Non_Demented_adj_p_log$ID,]
rownames(Demented_vs_Non_Demented_Table) <- Demented_vs_Non_Demented_Table$Geneid
Demented_vs_Non_Demented_Table <- Demented_vs_Non_Demented_Table[Demented_vs_Non_Demented_adj_p_log$ID,]

colnames(Demented_vs_Non_Demented_adj_p_log)
Demented_vs_Non_Demented_HC <- cbind(Demented_vs_Non_Demented_Table[,c(4)], Demented_vs_Non_Demented_adj_p_log[, c(1,2,3,4,5,6,7)], Demented_vs_Non_Demented_Table[,c(5,6,7,8,10,11,9)])

colnames(Demented_vs_Non_Demented_HC)[1] <- "Chr"
colnames(Demented_vs_Non_Demented_HC)[8] <- "DifExpressed"
colnames(Demented_vs_Non_Demented_HC)[15] <- "Description"
datatable(Demented_vs_Non_Demented_HC %>% arrange(factor(Chr, levels = c("1", "2", "3", "4", "5", "6", "7", "8",
                                                    "9", "10", "11", "12", "13", "14", "15",
                                                    "16", "17", "18", "19", "20", "21", "22",
                                                    "X", "Y", "MT"))), options = list(pageLength = -1))

Demented_vs_Non_Demented_HC <- Demented_vs_Non_Demented_HC %>% arrange(factor(Chr, levels = c("1", "2", "3", "4", "5", "6", "7", "8",
                                                               "9", "10", "11", "12", "13", "14", "15",
                                                               "16", "17", "18", "19", "20", "21", "22",
                                                               "X", "Y", "MT")))

# Demented versus Non-Demented DEGs with genedata 
write.csv(Demented_vs_Non_Demented_HC, "NPJ_dementia/1_Datasets/DEG_results/Demented_vs_Non_Demented_with_Genedata.csv")

# Step 7: Volcano plots ---------------------------------------------------
# Volcano plot: Parkinson´s disease demented
png(filename = "NPJ_dementia/3_Figures/5_DEG/Demented_vs_Non_Demented.png",
    width     = 25,
    height    = 20,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(data = Demented_vs_Non_Demented, aes(x = logFC, y = -log10(P.Value), col = diffexpressed, label = ID)) + 
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
  ggtitle(expression(bold(paste("Demented ",bolditalic("versus"), " Non-demented")))) +
  geom_text_repel(data = filter(Demented_vs_Non_Demented, logFC > 0.3 & adj.P.Val < 0.0007 | logFC > 3.5 & adj.P.Val < 0.1), aes(label = ID), size = 4.5, col = "black", max.overlaps = 200, box.padding = 0.5) +
  geom_text_repel(data = filter(Demented_vs_Non_Demented, logFC < -0.3 & adj.P.Val < 0.0004 | logFC < -3 & adj.P.Val < 0.2), aes(label = ID), size = 4.5, col = "black", max.overlaps = 200, box.padding = 0.5)

dev.off()

# Step 8: Heatmap --------------------------------------------------------
# Load countsTMM
countsTMM <- read.delim("NPJ_dementia/1_Datasets/Counts_Matrices/countsTMM.csv", sep = ",", header = TRUE)
row.names(countsTMM) <- countsTMM[,1]
countsTMM$X <- NULL

# Heatmap: Demented versus Non-demented
clus <- phenoHEROES[phenoHEROES$Dementia == "Demented" | phenoHEROES$Dementia == "Non_Demented", c("Tube_code", "Dementia", "Sex")]
clusz <- clus[clus$Dementia == "Demented" | clus$Dementia == "Non_Demented",]

clusG <- clusz$Sex
clusz <- clusz$Tube_code
Demented_vs_Non_Demented_adj_p_log
data.clus <- countsTMM[rownames(Demented_vs_Non_Demented_adj_p_log),]
data.clus <- data.clus[, colnames(data.clus) %in% clusz]

Condition <- as.factor(phenoHEROES$Dementia)
cond.df <- as.data.frame(Condition)
Condition <- cond.df[cond.df$Condition == "Demented" | cond.df$Condition == "Non_Demented",]
Condition <- sub("_D", "-d", Condition)

annotation <- data.frame(Condition = factor(Condition), Sex = as.factor(clusG))
annotation_colors = list(
  Condition = c("Non-demented" = "limegreen", Demented = "brown4"),
  Sex = c(Male = "lightblue", Female = "pink"))
rownames(annotation) <- phenoHEROES$Case
colnames(data.clus) <- phenoHEROES$Case

# Demented versus Non-demented
png(filename="NPJ_dementia/3_Figures/5_DEG/Heatmap_Demented_versus_Non_Demented.png",
    width     = 25,
    height    = 35,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

pheatmap(data.clus, scale = "row", show_rownames = FALSE, clustering_distance_cols = "euclidean", clustering_method = "ward.D2",
         main = expression(bold(paste("Demented ",bolditalic("versus"), " Non-demented"))),
         treeheight_col = 20, treeheight_row = 20, border_color = NA,
         fontsize = 20, fontsize_col = 18, annotation = annotation, 
         annotation_colors = annotation_colors)

dev.off()

# Step 9: Number of gene on chromosomes -----------------------------------
# Chromosome: Demented versus Non-demented
# Already done, but extracted again to be sure about the conditions (p-value adjusted < 0.5 & log2FC > 0.3 & < 0.3)
Demented_vs_Non_Demented_adj_p_log <- Demented_vs_Non_Demented[Demented_vs_Non_Demented$adj.P.Val <= 0.05 & abs(Demented_vs_Non_Demented$logFC) >= 0.3,]

# Extract up-regulated genes and gene information data (i.e., chromosomal location)
UP_HC_Demented_vs_Non_Demented <- Demented_vs_Non_Demented_adj_p_log[Demented_vs_Non_Demented_adj_p_log$diffexpressed == "UP",]
UP_HC_Demented_vs_Non_Demented <- genedata[genedata$Geneid %in% UP_HC_Demented_vs_Non_Demented$ID,]
UP <- UP_HC_Demented_vs_Non_Demented$chromosome_name
UP <- as.data.frame(UP)
UP$category <- rep("Up", rep = length(UP))

# Extract down-regulated genes and gene information data (i.e., chromosomal location)
DOWN_HC_Demented_vs_Non_Demented <- Demented_vs_Non_Demented_adj_p_log[Demented_vs_Non_Demented_adj_p_log$diffexpressed == "DOWN",]
DOWN_HC_Demented_vs_Non_Demented <- genedata[genedata$Geneid %in% DOWN_HC_Demented_vs_Non_Demented$ID,]
DOWN <- DOWN_HC_Demented_vs_Non_Demented$chromosome_name
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
tble
# Chromosome Demented
png(filename="NPJ_dementia/3_Figures/5_DEG/Chromosome_Dementia.png",
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
  ggtitle(expression(bold(paste("Demented ",bolditalic("versus"), " Non-demented")))) 

dev.off()
