# 8. Weighted Gene Co-expression Network Analysis (WGCNA)
# Load libraries
library(WGCNA) # Weighted Gene Co-expression Network Analysis
library(edgeR)
library(DESeq2) # For filtering and normalizing your genes
library(reshape2)
library(tidyverse)
library(CorLevelPlot) # For plotting the Module eigenes versus traits (e.g. Dementia, PDD, etc.)
library(gridExtra)
library(ggplot2)
library(flashClust) # Some fancy clustering, but not used at the end
library(uchardet) # To load RCy3
library(RCy3) # To export file into cytoscape, open cytoscape first on your desktop 

# Setting string not as factor
options(stringsAsFactors = FALSE)

#Enable multithread --> enableWGCNAThreads(2)
WGCNAnThreads()
allowWGCNAThreads()

## Step 1: Load phenotype data
phenoHEROES <- read.delim("NPJ_dementia/1_Datasets/PhenoData/phenoHEROES_RNAAge_FPKM.csv", sep = ",", header = TRUE, , stringsAsFactors = FALSE)
row.names(phenoHEROES) <- phenoHEROES[,1]
phenoHEROES$X <- NULL
colnames(phenoHEROES)[1] <- "Disease"
phenoHEROES$Status <- ifelse(phenoHEROES$Disease == "Control", "Non_Demented", "Demented")
phenoHEROES <- phenoHEROES[, -c(23,24,25)] # Remove PC values

## Step 2: Load original Counts Matrix
countsHEROES <- read.delim("NPJ_dementia/1_Datasets/Counts.matrix.csv", sep = ",", header = TRUE)
countsHEROES <- cbind(countsHEROES$X, countsHEROES[, phenoHEROES$Tube_code])
row.names(countsHEROES) <- countsHEROES[,1]
colnames(countsHEROES)[1] <- "GeneID"
countsHEROES <- countsHEROES[,2:21]

# Metadata
head(phenoHEROES)
# Counts data
dim(countsHEROES) # 60645 genes, 20 samples

## Step 3: QC Outlier detection
# Detect outliers, with goodSampleGenes function from the WGCNA package
gsg <- goodSamplesGenes(t(countsHEROES))

summary(gsg)
gsg$allOK
table(gsg$goodGenes)
table(gsg$goodSamples)

# Remove the genes that are detected as outliers
countsHEROES <- countsHEROES[gsg$goodGenes == TRUE,]
dim(countsHEROES) # 55463 genes, 20 samples

# Detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(countsHEROES), method = "manhattan"), method = "ward")
plot(htree)

# Principal component analysis - method 2
pca <- prcomp(t(countsHEROES))

pca_data <- pca$x[, c("PC1", "PC2", "PC3")]
pca_var <- pca$sdev^2
pca_var_percent <- round(pca_var/sum(pca_var)*100, digits = 2)

pca_data <- as.data.frame(pca_data)
phenoHEROES <- cbind(phenoHEROES, pca_data)
phenoHEROES

ggplot(phenoHEROES, aes(PC1, PC2, colour = Disease)) +
  geom_point() +
  geom_text(label = rownames(pca_data)) +
  labs(x = paste0("PC1: ", pca_var_percent[1], " %"),
       y = paste0("PC2: ", pca_var_percent[2], " %")) +
  scale_color_manual(breaks = c("Control", "PDD", "AD", "DSD"), values=c("limegreen", "navyblue", "gold", "darkred")) +
  theme_classic()
  
## Step 4: Create a DESeq2 dataset
all(rownames(phenoHEROES %in% colnames(countsHEROES)))
all(rownames(phenoHEROES) == colnames(countsHEROES))

# Create dds
dds <- DESeqDataSetFromMatrix(countData = countsHEROES,
                              colData = phenoHEROES,
                              design = ~ 1) # Not specifying a model

# Remove all genes with counts <15 in more than 75% of the samples (20*0.75 = 15)
# Suggested by WGCNA on RNAseq FAQ
dds75 <- dds[rowSums(counts(dds) >= 15) >= 15, ]
nrow(dds75) # 23150 genes

# Perform variance stabilization 
dds_norm <- vst(dds75)

# Get normalized counts
norm_counts <- assay(dds_norm) %>%
  t()

head(norm_counts)[1:5, 1:5]

## Step 5: Network construction
# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
powers

# Call network topology analysis function
sft <- pickSoftThreshold(norm_counts,
                         powerVector = powers,
                         networkType = "signed",
                         verbose = 5)
sft_data <- sft$fitIndices
sft$powerEstimate

# Visualization to pick power
# Scale independence
png(filename="NPJ_dementia/3_Figures/8_WGCNA/Scale_Independence.png",
    width     = 15,
    height    = 15,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(sft_data, aes(Power, SFT.R.sq, label = Power)) +
  geom_text(color = "red") +
  geom_hline(yintercept = 0.85, color = "red") +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2), limits = c(0,1)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 12), 
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold")) +
  ggtitle("Scale independence") + 
  labs(x = "Soft threshold (power)", y = expression(bold(Scale~free~topology~signed~R^{2})))

dev.off()

# Mean connectivity
png(filename="NPJ_dementia/3_Figures/8_WGCNA/Mean_Connectivity.png",
    width     = 15,
    height    = 15,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(sft_data, aes(Power, mean.k., label = Power)) +
  geom_text(color = "red") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text = element_text(size = 12), 
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold")) +
  ggtitle("Mean connectivity") + 
  labs(x = "Soft threshold (power)", y = "Mean connectivity")

dev.off()

# Convert matrix to numeric, she is using the blockwisemodule, less computational
# Maybe change here, the rest is pretty good
norm_counts[] <- sapply(norm_counts, as.numeric)

temp_cor <- cor
cor <- WGCNA::cor

# unsigned -> nodes with positive & negative correlation are treated equally 
# signed -> nodes with negative correlation are considered *unconnected*, treated as zero
# Create the Similarity Matrix
sim <- cor(norm_counts)
dim(sim)

# Calculate the Adjacency Matrix
softPower <- 7
adjacency <- adjacency.fromSimilarity(sim, 
                                      power = softPower, 
                                      type = "signed")

rm(sim) # Free some space

# Calculate the Topological Matrix (TOM)
TOM <- TOMsimilarity(adjacency,
                     TOMType = "signed",
                     verbose = TRUE) # Calculates the similarity

TOM_diss <- 1-TOMsimilarity(adjacency,
                            TOMType = "signed",
                            verbose = TRUE) # Calculates the dissimilarity

rm(adjacency) # Free some space
cor <- temp_cor

## Step 6: Module detection and visualization
geneTree <- hclust(as.dist(TOM_diss), method = "average") 

# Plotting the dendrogram
# sizeGrWindow(12,9)
plot(geneTree, 
     xlab="", 
     sub="", 
     main = "Gene clustering on TOM-based dissimilarity", 
     labels = FALSE, 
     hang = 0.04)

Modules <- cutreeDynamic(dendro = geneTree, 
                         distM = TOM_diss,
                         method = "tree",
                         deepSplit = TRUE, 
                         pamRespectsDendro = FALSE, 
                         minClusterSize = 50)

table(Modules)

# Assigns each module number a color
ModuleColors <- labels2colors(Modules)

# Returns the counts for each color (i.e., the number of genes within each module)
table(ModuleColors) 

# Plots the gene dendrogram with the module colors
plotDendroAndColors(geneTree, 
                    ModuleColors,
                    "Module",
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Module Eigengene Identification
MElist <- moduleEigengenes(norm_counts, colors = ModuleColors) 
MEs <- MElist$eigengenes 
head(MEs)

# Calculating the module dissimilarity eigengenes
MEDiss <- 1-cor(MEs)

# Clustering the eigengenes modules
METree <- hclust(as.dist(MEDiss), method = "average")

# Plotting the result
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# Plotting a cut-off line
abline(h = 0.25, col = "red")

merge <- mergeCloseModules(norm_counts, ModuleColors, cutHeight = 0.25)

# Grouping module colors
mergedColors <- merge$colors

# Eigengenes of new grouped modules
mergedMEs <- merge$newMEs
ncol(mergedMEs)

# Cluster dendrogram
png(filename="NPJ_dementia/3_Figures/8_WGCNA/Cluster_Dendrogram.png",
    width     = 30,
    height    = 15,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

plotDendroAndColors(geneTree, cbind(ModuleColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged Dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Cluster dendrogram",
                    autoColorHeight = FALSE,
                    colorHeight = 0.2,
                    cex.main = 5, cex.colorLabels = 3, 
                    cex.rowText = 1, cex.axis = 3, cex.lab = 3,
                    ylab = "",
                    marAll = c(1, 16, 4, 1))
dev.off()

## Step 7: Module to Trait Relationship
# Create traits file - binarize categorical variables
traits1 <- phenoHEROES %>%
  mutate(Control = ifelse(grepl("Non_Demented", Status), 1, 0)) %>% 
  dplyr::select(27)

traits2 <- phenoHEROES %>%
  mutate(Demented = ifelse(grepl("Non_Demented", Status), 0, 1)) %>%
  dplyr::select(27) # Optional and is coded otherway around

# Binarize categorical variables 
phenoHEROES$Disease <- factor(phenoHEROES$Disease, levels = c("Control", "PDD", "AD", "DSD"))

disease_out <- binarizeCategoricalColumns(phenoHEROES$Disease,
                           includePairwise = FALSE,
                           includeLevelVsAll = TRUE)

colnames(disease_out) <- c("PDD", "AD", "DSD")
traits3 <- cbind(traits1, disease_out)
traits3 # Original

traits4 <- cbind(traits1, traits2)
traits4

traits_PDD <- traits3[traits3$Control == 1 | traits3$PDD == 1, c(1,2)] 
traits_PDD  

traits_AD <- traits3[traits3$Control == 1 | traits3$AD == 1, c(1,3)] 
traits_AD

traits_DSD <- traits3[traits3$Control == 1 | traits3$DSD == 1, c(1,4)] 
traits_DSD

# Define number of genes and samples
nSamples <- nrow(norm_counts)
nGenes <- ncol(norm_counts)

# Demented vs Non-demented (post-hoc test)
module_trait_corr_1 <- cor(mergedMEs, traits4, use = "p")
module_trait_corr_pvals_1 <- corPvalueStudent(module_trait_corr_1, nSamples)
module_trait_corr_pvals_1

# GS.adj <- p.adjust(GS.pvalue, method = "BH"), post-hoc test
pvals_fdr_1 <- matrix(
  p.adjust(module_trait_corr_pvals_1, method = "BH"),
  nrow = nrow(module_trait_corr_pvals_1),
  dimnames = dimnames(module_trait_corr_pvals_1)
)
pvals_fdr_1

# PDD versus Non-demented (post-hoc test)
mergedMEs_PDD <- mergedMEs[rownames(mergedMEs) %in% rownames(traits_PDD),]
module_trait_corr_PDD <- cor(mergedMEs_PDD, traits_PDD, use = "p")
module_trait_corr_pvals_PDD <- corPvalueStudent(module_trait_corr_PDD, 9)
module_trait_corr_pvals_PDD
# GS.adj <- p.adjust(GS.pvalue, method = "BH"), post-hoc test
pvals_fdr_PDD <- matrix(
  p.adjust(module_trait_corr_pvals_PDD, method = "BH"),
  nrow = nrow(module_trait_corr_pvals_PDD),
  dimnames = dimnames(module_trait_corr_pvals_PDD)
)
pvals_fdr_PDD

# AD versus Non-demented (post-hoc test)
mergedMEs_AD <- mergedMEs[rownames(mergedMEs) %in% rownames(traits_AD),]
module_trait_corr_AD <- cor(mergedMEs_AD, traits_AD, use = "p")
module_trait_corr_pvals_AD <- corPvalueStudent(module_trait_corr_AD, 11)

# GS.adj <- p.adjust(GS.pvalue, method = "BH"), post-hoc test
pvals_fdr_AD <- matrix(
  p.adjust(module_trait_corr_pvals_AD, method = "BH"),
  nrow = nrow(module_trait_corr_pvals_AD),
  dimnames = dimnames(module_trait_corr_pvals_AD)
)
pvals_fdr_AD

# DSD versus Non-demented (post-hoc test)
mergedMEs_DSD <- mergedMEs[rownames(mergedMEs) %in% rownames(traits_DSD),]
module_trait_corr_DSD <- cor(mergedMEs_DSD, traits_DSD, use = "p")
module_trait_corr_pvals_DSD <- corPvalueStudent(module_trait_corr_DSD, 10)

# GS.adj <- p.adjust(GS.pvalue, method = "BH"), post-hoc test
pvals_fdr_DSD <- matrix(
  p.adjust(module_trait_corr_pvals_DSD, method = "BH"),
  nrow = nrow(module_trait_corr_pvals_DSD),
  dimnames = dimnames(module_trait_corr_pvals_DSD)
)
pvals_fdr_DSD

# Visualize module-trait association as a heatmap
heatmap_data <- merge(mergedMEs, traits3, by = "row.names")
heatmap_data 
heatmap_data <- heatmap_data %>%
 column_to_rownames(var = "Row.names")

heatmap_control <- merge(mergedMEs, traits4, by = "row.names")
heatmap_control 
heatmap_control <- heatmap_control %>%
  column_to_rownames(var = "Row.names")
data_control <- cor(heatmap_control, method = "pearson", use = "pairwise.complete.obs")
Control <- data_control[-c(23,24), 23]

heatmap_data_PDD <- merge(mergedMEs, traits_PDD, by = "row.names")
heatmap_data_PDD 
heatmap_data_PDD <- heatmap_data_PDD %>%
  column_to_rownames(var = "Row.names")
data_PDD <- cor(heatmap_data_PDD, method = "pearson", use = "pairwise.complete.obs")
PDD <- data_PDD[-c(23,24), 24]

heatmap_data_AD <- merge(mergedMEs, traits_AD, by = "row.names")
heatmap_data_AD 
heatmap_data_AD <- heatmap_data_AD %>%
  column_to_rownames(var = "Row.names")
data_AD <- cor(heatmap_data_AD, method = "pearson", use = "pairwise.complete.obs")
AD <- data_AD[-c(23,24), 24]

heatmap_data_DSD <- merge(mergedMEs, traits_DSD, by = "row.names")
heatmap_data_DSD 
heatmap_data_DSD <- heatmap_data_DSD %>%
  column_to_rownames(var = "Row.names")
data_DSD <- cor(heatmap_data_DSD, method = "pearson", use = "pairwise.complete.obs")
DSD <- data_DSD[-c(23,24), 24]

# Create biological useful data frame
complete <- cbind(data_control[-c(23,24), 23], data_PDD[-c(23,24), 24], data_AD[-c(23,24), 24], data_DSD[-c(23,24), 24])
colnames(complete) <- c("Control", "PDD", "AD", "DSD")
heatmap_df$Status
heatmap_df <- complete %>%
  as.data.frame() %>%
  rownames_to_column("Group") %>%
  pivot_longer(-Group,
               names_to = "Status",
               values_to = "Correlation")

heatmap_df <- as.data.frame(heatmap_df)
heatmap_df$Status <-  factor(heatmap_df$Status, levels = c("Control", "PDD", "AD", "DSD"))

heatmap_df$stars <- cut(
  abs(heatmap_df$Correlation),
  breaks = c(0, 0.74, 0.76, 0.84, Inf),
  labels = c("", "*", "", "*"),
  right = FALSE
)
library(ggtext)
# Check for gene modules that are significantly associated with a status (trait)
# Module eigengenes are representative gene expression profiles of a cluster/module 
# Is there a difference in gene expression profile between both these groups?
png(filename="NPJ_dementia/3_Figures/8_WGCNA/Module_trait_relationship.png",
    width     = 20,
    height    = 30,
    units     = "cm",
    res       = 1200,
    pointsize = 4) 

ggplot(heatmap_df, aes(Status, Group, fill = Correlation)) +
  geom_tile() +
  coord_fixed(ratio = 1/2.5) + 
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0,
    limits = c(-1, 1), name = "") +
  labs(title = "Module-trait relationship") +
  scale_x_discrete(labels = c("Control" = "Control vs PDD, AD & DSD",
                              "PDD" = "PDD vs Control",
                              "AD" = "AD vs Control",
                              "DSD" = "DSD vs Control")) +
  geom_text(label = paste0(round(heatmap_df$Correlation, 2), heatmap_df$stars), color = "black", size = 5.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 16, face = "bold", angle = -30, hjust = 0),
        axis.text.y = element_text(size = 16, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_line(),
        axis.ticks.length = unit(2, "mm"),
        plot.title = element_markdown(hjust = 0.5, size = 30, face = "bold"),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.key.height = unit(4.8, "cm"),
        legend.key.width  = unit(0.5, "cm"),
        strip.text = element_text(color = "white", size = 30)) +
  annotate('rect', xmin = 0.4, xmax = 4.6, ymin = 0.4, ymax = 22.6,
           fill = NA, color = 'black', size = 0.3) +
  guides(fill = guide_colourbar(label = TRUE,
                                ticks = FALSE))

dev.off ()

## Step 8: Generate a list of genes
# Darkseagreen4 and maroon, compare this against the previously common DEGs 
DEGs <- read.csv("NPJ_dementia/3_Figures/6_Common expressed_genes/DEGs_Genes.csv", header = TRUE)
myListOfGenes <- as.character(DEGs$X)

colnames(norm_counts)
module_gene_mapping <- as.data.frame(mergedColors)
rownames(module_gene_mapping) <- colnames(norm_counts)
colnames(module_gene_mapping) <- "Color"
head(module_gene_mapping)

darkseagreen <- module_gene_mapping %>%
  filter(Color == "darkseagreen4") %>%
  rownames()

length(darkseagreen)
print(intersect(darkseagreen, myListOfGenes)) # 5 Common DEGs in darkseagreen module!!

## Step 9: Intramodular analysis: Identifying driver genes 
# Calculate the module membership and the associated p-values
# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene
# This quantifies the similarity of all genes on the array to every module
module_membership_measure <- cor(mergedMEs, norm_counts, use = "p")
module_membership_measure_pvals <- corPvalueStudent(module_membership_measure, nSamples)

module_membership_measure_pvals[1:10, 1:10]

# Calculate the gene significance and the associated p-value

gene_sign_corr <- cor(norm_counts, traits1$Control, use = "p")
gene_sign_corr_pvals <- corPvalueStudent(gene_sign_corr, nSamples)

gene_sign_corr_pvals %>%
  as.data.frame() %>%
  arrange(V1) %>%
  head(25)

# Using the gene significance you can identify genes that have a high significance for weight
# Using the module membership measures you can identify genes with high module membership in interesting modules
# gene significance (GS) and modular membership (MM)

## Step 10: Rpubs - Associating Modules and Phenotypes
# The following code allows us to access the genes inside a module, 
# order them by importance order and stores this information in a csv file

#Defining the variable Control containing the column Control of traits
Control <- as.data.frame(traits1$Control)
names(Control) <- "Control"

# names (colors) of the modules
modNames <- substring(names(mergedMEs), 3)
modNames

geneModuleMembership <- as.data.frame(cor(norm_counts, mergedMEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

geneTraitSignificance <- as.data.frame(cor(norm_counts, traits1$Control, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) <- paste("GS.", names(Control), sep="")
names(GSPvalue) <- paste("p.GS.", names(Control), sep="")

# Identifying most important genes for one determined characteristic inside of the cluster
geneInfo0 <- data.frame(EST = colnames(norm_counts),
                       moduleColor = mergedColors,
                       geneTraitSignificance,
                       GSPvalue)

geneInfo0[geneInfo0$moduleColor == "darkseagreen4",]

modOrder <- order(-abs(cor(mergedMEs, traits1$Control, use = "p")))

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames <- names(geneInfo0)
  geneInfo0 <- data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) <- c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
head(geneInfo0)

geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Control))
geneInfo <- geneInfo0[geneOrder, ]
geneInfo2 <- geneInfo[geneInfo$moduleColor == "darkseagreen4",]
rownames(geneInfo2)
rownames(geneInfo2[abs(geneInfo2$MM.darkseagreen4) > 0.80 & abs(geneInfo2$GS.Control) > 0.70,])
geneInfo0$GS.Control

write.csv(geneInfo2, file = "NPJ_dementia/3_Figures/8_WGCNA/geneInfo2.csv")

## Step 11: RPubs - Saving the Analysis
# Last step is to export and save the network. 
# Then you can import it in a software for network visualization as Cytoscape, for example.

# Exporting the network to a cytoscape format
# Recalculating topological overlap, if necessary
# TOM = TOMsimilarityFromExpr(expression0, power = 10);
# Select the modules
modules = "darkseagreen4" #chose modules that u want to export
# Select the gene modules
genes <- colnames(norm_counts)
genes

# if you want export specific colors, substitute the second module colors by above modules
inModule <- is.finite(match(mergedColors, modules))
sum(inModule)
modGenes <- genes[inModule]
modGenes

# Select the corresponding topologic overlap 
modTOM <- TOM[inModule, inModule]
modTOM
dimnames(modTOM) <- list(modGenes, modGenes)
modTOMSignificantes <- which(modTOM > 0.4)

##### warnings()
cyt <- exportNetworkToCytoscape(modTOM,
                                edgeFile = "CytoscapeEdgeFile.txt",
                                nodeFile = "CytoscapeNodeFile.txt",
                                weighted = TRUE,
                                threshold = 0.1,
                                nodeNames = modGenes,
                                nodeAttr = ModuleColors[inModule])

# Step 12: Export to Cytoscape --------------------------------------------
cytoscapePing () # make sure cytoscape is open
cytoscapeVersionInfo ()

################ Darkseagreen module of the merged data (newMEs) ######################
edge <- read.delim("CytoscapeEdgeFile.txt")
colnames(edge)
colnames(edge) <- c("source", "target", "weight", "direction", "fromAltName", "toAltName")

node <- read.delim("CytoscapeNodeFile.txt")
colnames(node)  
colnames(node) <- c("id", "altName", "node_attributes") 

createNetworkFromDataFrames(node, edge, title = "Darkseagreen", collection = "Dementia")

################ customise the network visualization ##################################
# use other pre-set visual style
setVisualStyle('Marquee')
?setVisualStyle

# Set up my own style
style.name = "myStyle"
defaults <- list(NODE_SHAPE = "diamond",
                 NODE_SIZE = 50,
                 EDGE_TRANSPARENCY = 100,
                 NODE_LABEL_POSITION="W,E,c,0.00,0.00")
nodeLabels <- mapVisualProperty('node label','id','p')
nodeFills <- mapVisualProperty('node fill color','node_attributes','d',c("A","B"), c("#FF9900","#66AAAA"))
arrowShapes <- mapVisualProperty('Edge Target Arrow Shape','interaction','d',c("activates","inhibits","interacts"),c("Arrow","T","None"))
edgeWidth <- mapVisualProperty('edge width','weight','p')

createVisualStyle(style.name, defaults, list(nodeLabels,nodeFills,arrowShapes,edgeWidth))
setVisualStyle(style.name)

# Check hub genes distribution over samples/groups
# Load phenotype data
phenoHEROES <- read.delim("NPJ_dementia/1_Datasets/PhenoData/phenoHEROES.csv", sep = ",", header = TRUE)
row.names(phenoHEROES) <- phenoHEROES$Tube_code
phenoHEROES$X <- NULL
head(phenoHEROES)

# Load normalized count data (TTM)
NoLOGcountsTMM <- read.csv("NPJ_dementia/1_Datasets/Counts_Matrices/NoLOGcountsTMM", header = TRUE, row.names = 1)

# Check LMNB2 and EHMT2 in matrix
LMNB2 <- t(NoLOGcountsTMM["LMNB2",])
EHMT2 <- t(NoLOGcountsTMM["EHMT2",])

# Merge normalized data with phenotype data
phenoHEROES_HUB <- cbind(phenoHEROES, LMNB2, EHMT2)
phenoHEROES_HUB$Status <- factor(phenoHEROES_HUB$Status, levels = c("Control", "PDD", "AD", "DSD"))

# Save data
write.csv(phenoHEROES_HUB, "NPJ_dementia/1_Datasets/PhenoData/phenoHEROES_HUB.csv")

# LMNB2 normalized expression in countsMatrix 
stat.test_LMNB2 <- aov(LMNB2 ~ Status, data = phenoHEROES_HUB) %>%
  tukey_hsd()
stat.test_LMNB2

png(filename="NPJ_dementia/3_Figures/8_WGCNA/LMNB2.png",
    width     = 20,
    height    = 15,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(data = phenoHEROES_HUB, aes(x = Status, y = LMNB2, fill = factor(Status))) +
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
  scale_y_continuous(expand = c(0, 0), limits = c(0, 80), breaks = c(0, 20, 40, 60, 80, 100)) +
  ylab(expression(bold("LMNB2 transcript (TTM)"))) +
  ggtitle(paste("LMNB2 expression", "\nHippocampus")) + 
  stat_pvalue_manual(stat.test_LMNB2, hide.ns = TRUE, y.position = c(60, 65, 70), size = 12)

dev.off()

# EHMT2 normalized expression in countsMatrix 
stat.test_EHMT2 <- aov(EHMT2 ~ Status, data = phenoHEROES_HUB) %>%
  tukey_hsd()
stat.test_EHMT2

png(filename="NPJ_dementia/3_Figures/8_WGCNA/EHMT2.png",
    width     = 20,
    height    = 15,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(data = phenoHEROES_HUB, aes(x = Status, y = EHMT2, fill = factor(Status))) +
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
  scale_y_continuous(expand = c(0, 0), limits = c(0, 140), breaks = c(0, 20, 40, 60, 80, 100, 120, 140)) +
  ylab(expression(bold("EHMT2 transcript (TTM)"))) +
  ggtitle(paste("EHMT2 expression", "\nHippocampus")) + 
  stat_pvalue_manual(stat.test_EHMT2, hide.ns = TRUE, y.position = c(105, 115, 125), size = 12)

dev.off()
