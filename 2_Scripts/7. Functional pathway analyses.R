# 7. Functional pathway analyses
# Load libraries
library(biomaRt)
library(ggplot2)
library(dplyr)
library(DT)
library(pheatmap)

# Gene sets
library(org.Hs.eg.db) # keytypes(org.Hs.eg.db)
library(msigdbr)
library(ExperimentHub)
library(msigdb)
library(DOSE)
library(disgenetplus2r) # DisGeNet Plus
library(devtools)
install_gitlab("medbio/disgenet2r")
library(disgenet2r)
API_KEY = "83dba455-8715-4e00-bf1a-309f395cddf2"
Sys.setenv(DISGENET_API_KEY= API_KEY)
sessionInfo()
# Gene ontology
library(pathview)
library(GOstats)
library(clusterProfiler) # gseGO function
library(enrichplot)
library(devtools)
library(GSVA) # Gene Set Variation Analysis
library(GSEABase)
library(GSVAdata)
library(GOSemSim)
library(STRINGdb)
library(poweRlaw)

# Load previous created data to obtain ENSEMBL IDs for HUGO symbols in data sets
GeneInfo <- read.csv("NPJ_dementia/1_Datasets/Counts_Matrices/counts_length.csv", header = TRUE, row.names = 1)
colnames(GeneInfo)[1] <- "ID"
GeneInfo <- GeneInfo[,-c(2:21)]

# Parkinson's disease demented GO, KEGG and DO ----------------------------
# Load result from DESeq2 (see script 5. DEGs.R)
PDD_HC <- read.csv("NPJ_dementia/1_Datasets/DEG_results/HC_PDD_vs_Control.csv", header = TRUE)
colnames(PDD_HC)[1] <- "ID"
PDD_HC <- merge(PDD_HC, GeneInfo, by.x = "ID", by.y = "ID")
PDD_HC <- PDD_HC[order(PDD_HC$logFC, decreasing = TRUE),]

## Make an order in the gene list (ranking) 
## Rank the genes of the sample data results, combining the sign of the logFC and the p-value
# Create a metric to order the gene
sign <- sign(PDD_HC$logFC)
logP <- -log2(PDD_HC$P.Val)
metric <- logP*sign

# Add to the data frame
PDD_HC$metric <- metric

# Extract the metric to create a gene list (logFC)
gene_list_PDD <- PDD_HC$metric
names(gene_list_PDD) <- PDD_HC$ensembl_gene_id
gene_list_PDD <- sort(gene_list_PDD, decreasing = TRUE)

# Check distribution of gene_list_PDD
plot(gene_list_PDD)
barplot(sort(gene_list_PDD, decreasing = TRUE), xaxt = "n")
max(gene_list_PDD)
min(gene_list_PDD)

## GSEA with Gene Ontology: Biological Processes (BP)
gsePDD_BP <- gseGO(geneList = gene_list_PDD, 
                   by = "fgsea",
                   ont = "BP",
                   nPermSimple = 1000,
                   keyType = "ENSEMBL",
                   minGSSize = 3, 
                   maxGSSize = 10000, 
                   eps = 0,
                   pvalueCutoff = 0.05, 
                   verbose = TRUE, 
                   OrgDb = "org.Hs.eg.db", 
                   pAdjustMethod = "fdr")

# Dotplot: Gene Ontology: Biological Processes (BP)
png(filename = "NPJ_dementia/3_Figures/7_Functional_pathway_analyses/PDD/GO_BP_PDD_dotplot.png",
    width     = 30,
    height    = 50,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

dotplot(gsePDD_BP, showCategory = 10, split=".sign") + 
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 20, color = "black"),
        strip.background = element_rect(fill = "lightgray", colour = "black", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.5, units = "cm")) +
  ggtitle("GO:BP - Parkinson's disease demented") + 
  facet_grid(.~.sign, labeller = labeller(.sign = 
                                            c(activated = "Upregulated",
                                              suppressed = "Downregulated")))

dev.off()

# Enrichment Map: Gene Ontology: Biological Processes (BP)
GO_BP_PDD_dotplot2 <- dotplot(gsePDD_BP, showCategory = 30, split = ".sign") + facet_grid(.~.sign)
terms <- as.character(GO_BP_PDD_dotplot2$data$Description)
GO_BP_PDD_dotplot <- pairwise_termsim(gsePDD_BP, showCategory = dim(gsePDD_BP)[1])

# Enrichment Map: UPREGULATED
png(filename = "NPJ_dementia/3_Figures/7_Functional_pathway_analyses/PDD/GO_BP_PDD_UP_emapplot.png",
    width     = 40,
    height    = 30,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

emapplot(GO_BP_PDD_dotplot, showCategory = terms[1:20], size_category = 1, 
         layout = "kk", color = "p.adjust") + 
  theme(strip.text = element_text(face = "bold", size = 10, color = "black"),
        strip.background = element_rect(fill = "lightgray", colour = "black", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text = element_blank(),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.5, units = "cm")) +
  ggtitle("Upregulated GO:BP - Enrichment map Parkinson's disease demented") 

dev.off()

# Enrichment Map: DOWNREGULATED
png(filename = "NPJ_dementia/3_Figures/7_Functional_pathway_analyses/PDD/GO_BP_PDD_DOWN_emapplot.png",
    width     = 40,
    height    = 30,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

emapplot(GO_BP_PDD_dotplot, showCategory = terms[21:50], size_category = 1, 
         layout = "kk", color = "p.adjust") +
  theme(strip.text = element_text(face = "bold", size = 10, color = "black"),
        strip.background = element_rect(fill = "lightgray", colour = "black", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text = element_blank(), 
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.5, units = "cm")) +
  ggtitle("Downregulated GO:BP - Enrichment map Parkinson's disease demented") 

dev.off()

## GSEA with Gene Ontology: Molecular Function (MF)
gsePDD_MF <- gseGO(geneList = gene_list_PDD, 
                   by = "fgsea",
                   ont = "MF", 
                   keyType = "ENSEMBL",
                   nPermSimple = 1000,
                   minGSSize = 3, 
                   maxGSSize = 10000, 
                   eps = 0,
                   pvalueCutoff = 0.05, 
                   verbose = TRUE, 
                   OrgDb = "org.Hs.eg.db", 
                   pAdjustMethod = "fdr")

# Dotplot: Gene Ontology: Molecular Function (MF)
png(filename = "NPJ_dementia/3_Figures/7_Functional_pathway_analyses/PDD/GO_MF_PDD_dotplot.png",
    width     = 30,
    height    = 50,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

dotplot(gsePDD_MF, showCategory = 10, split = ".sign") + 
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 20, color = "black"),
        strip.background = element_rect(fill = "lightgray", colour = "black", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.5, units = "cm")) +
  ggtitle("GO:MF - Parkinson's disease demented") +  
  facet_grid(.~.sign, labeller = labeller(.sign = 
                                            c(activated = "Upregulated",
                                              suppressed = "Downregulated")))

dev.off()

## GSEA with Kyoto Encyclopedia of Genes and Genomes (KEGG)
# Extract entrez gene IDs for KEGG
PDD2 <- PDD_HC[!is.na(PDD_HC$entrezgene_id),]
PDD2 <- PDD2[!duplicated(PDD2$entrezgene_id),]

# Create a vector with the metric
kegg_gene_list_PDD <- PDD2$metric

# Name vector with ENTREZ ids
names(kegg_gene_list_PDD) <- PDD2$entrezgene_id

# omit any NA values 
kegg_gene_list_PDD <- na.omit(kegg_gene_list_PDD)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list_PDD <- sort(kegg_gene_list_PDD, decreasing = TRUE)

gseKEGG_PDD <- gseKEGG(geneList = kegg_gene_list_PDD,
                       by = "fgsea",
                       organism = "hsa",
                       nPermSimple = 10000,
                       minGSSize = 40,
                       maxGSSize = 1000,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "fdr",
                       keyType = "kegg",
                       eps = 0,
                       verbose = TRUE)

# Dotplot: Kyoto Encyclopedia of Genes and Genomes (KEGG)
png(filename = "NPJ_dementia/3_Figures/7_Functional_pathway_analyses/PDD/KEGG_PDD_dotplot.png",
    width     = 30,
    height    = 50,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

dotplot(gseKEGG_PDD, showCategory = 16, split = ".sign") + 
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 20, color = "black"),
        strip.background = element_rect(fill = "lightgray", colour = "black", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.5, units = "cm")) +
  ggtitle("KEGG - Parkinson's disease demented") + 
  facet_grid(.~.sign, labeller = labeller(.sign = 
                                            c(activated = "Upregulated",
                                              suppressed = "Downregulated")))

dev.off()

## GSEA with Disease Ontology (DO)
gseDO_PDD <- gseDO(kegg_gene_list_PDD,
                   by = "fgsea",
                   nPermSimple = 10000,
                   minGSSize = 40,
                   maxGSSize = 1000,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "fdr",
                   eps = 0,
                   verbose = TRUE)

# Dotplot: Disease Ontology (DO)
png(filename = "NPJ_dementia/3_Figures/7_Functional_pathway_analyses/PDD/DO_PDD_dotplot.png",
    width     = 30,
    height    = 50,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

dotplot(gseDO_PDD, showCategory = 19, split = ".sign") + 
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 20, color = "black"),
        strip.background = element_rect(fill = "lightgray", colour = "black", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.5, units = "cm")) +
  ggtitle("DO - Parkinson's disease demented") + 
  facet_grid(.~.sign, labeller = labeller(.sign = 
                                            c(activated = "Upregulated",
                                              suppressed = "Downregulated")))

dev.off()

# Alzheimer's disease GO, KEGG and DO -------------------------------------
# Load result from DESeq2 (see script 5. DEGs.R)
AD_HC <- read.csv("NPJ_dementia/1_Datasets/DEG_results/HC_AD_vs_Control.csv", header = TRUE)
colnames(AD_HC)[1] <- "ID"
AD_HC <- merge(AD_HC, GeneInfo, by.x = "ID", by.y = "ID")
AD_HC <- AD_HC[order(AD_HC$logFC, decreasing = TRUE),]

## Make an order in the gene list (ranking) 
## Rank the genes of the sample data results, combining the sign of the logFC and the p-value
# Create a metric to order the gene
sign <- sign(AD_HC$logFC)
logP <- -log2(AD_HC$P.Val)
metric <- logP*sign

# Add to the data frame
AD_HC$metric <- metric

# Extract the metric to create a gene list (logFC)
gene_list_AD <- AD_HC$metric
names(gene_list_AD) <- AD_HC$ensembl_gene_id
gene_list_AD <- sort(gene_list_AD, decreasing = TRUE)

# Check distribution of gene_list_AD
plot(gene_list_AD)
barplot(sort(gene_list_AD, decreasing = TRUE), xaxt = "n")
max(gene_list_AD)
min(gene_list_AD)

## GSEA with Gene Ontology: Biological Processes (BP)
gseAD_BP <- gseGO(geneList = gene_list_AD, 
                  by = "fgsea",
                  ont = "BP",
                  nPermSimple = 1000,
                  keyType = "ENSEMBL",
                  minGSSize = 3, 
                  maxGSSize = 10000, 
                  eps = 0,
                  pvalueCutoff = 0.05, 
                  verbose = TRUE, 
                  OrgDb = "org.Hs.eg.db", 
                  pAdjustMethod = "fdr")

# Dotplot: Gene Ontology: Biological Processes (BP)
png(filename = "NPJ_dementia/3_Figures/7_Functional_pathway_analyses/AD/GO_BP_AD_dotplot.png",
    width     = 30,
    height    = 50,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

dotplot(gseAD_BP, showCategory = 10, split = ".sign") + 
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 20, color = "black"),
        strip.background = element_rect(fill = "lightgray", colour = "black", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.5, units = "cm")) +
  ggtitle("GO:BP - Alzheimer's disease") + 
  facet_grid(.~.sign, labeller = labeller(.sign = 
                                            c(activated = "Upregulated",
                                              suppressed = "Downregulated")))

dev.off()

# Enrichment Map: Gene Ontology: Biological Processes (BP)
GO_BP_AD_dotplot2 <- dotplot(gseAD_BP, showCategory = 30, split = ".sign") + facet_grid(.~.sign)
terms <- as.character(GO_BP_AD_dotplot2$data$Description)
GO_BP_AD_dotplot <- pairwise_termsim(gseAD_BP, showCategory = dim(gseAD_BP)[1])

# Enrichment Map: UPREGULATED
png(filename = "NPJ_dementia/3_Figures/7_Functional_pathway_analyses/AD/GO_BP_AD_UP_emapplot.png",
    width     = 40,
    height    = 30,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

emapplot(GO_BP_AD_dotplot, showCategory = terms[1:30], size_category = 1, 
         layout = "kk", color = "p.adjust") + 
  theme(strip.text = element_text(face = "bold", size = 10, color = "black"),
        strip.background = element_rect(fill = "lightgray", colour = "black", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text = element_blank(),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.5, units = "cm")) +
  ggtitle("Upregulated GO:BP - Enrichment map Alzheimer's disease") 

dev.off()

# Enrichment Map: DOWNREGULATED
png(filename = "NPJ_dementia/3_Figures/7_Functional_pathway_analyses/AD/GO_BP_AD_DOWN_emapplot.png",
    width     = 40,
    height    = 30,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

emapplot(GO_BP_AD_dotplot, showCategory = terms[31:48], size_category = 1, 
         layout = "kk", color = "p.adjust") +
  theme(strip.text = element_text(face = "bold", size = 10, color = "black"),
        strip.background = element_rect(fill = "lightgray", colour = "black", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text = element_blank(), 
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.5, units = "cm")) +
  ggtitle("Downregulated GO:BP - Enrichment map Alzheimer's disease") 

dev.off()

## GSEA with Gene Ontology: Molecular Function (MF)
gseAD_MF <- gseGO(geneList = gene_list_AD, 
                  by = "fgsea",
                  ont = "MF", 
                  keyType = "ENSEMBL",
                  nPermSimple = 1000,
                  minGSSize = 3, 
                  maxGSSize = 10000, 
                  eps = 0,
                  pvalueCutoff = 0.05, 
                  verbose = TRUE, 
                  OrgDb = "org.Hs.eg.db", 
                  pAdjustMethod = "fdr")

# Dotplot: Gene Ontology: Molecular Function (MF)
png(filename = "NPJ_dementia/3_Figures/7_Functional_pathway_analyses/AD/GO_MF_AD_dotplot.png",
    width     = 30,
    height    = 50,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

dotplot(gseAD_MF, showCategory = 10, split = ".sign") + 
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 20, color = "black"),
        strip.background = element_rect(fill = "lightgray", colour = "black", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.5, units = "cm")) +
  ggtitle("GO:MF - Alzheimer's disease") +  
  facet_grid(.~.sign, labeller = labeller(.sign = 
                                            c(activated = "Upregulated",
                                              suppressed = "Downregulated")))

dev.off()

## GSEA with Kyoto Encyclopedia of Genes and Genomes (KEGG)
# Extract entrez gene IDs for KEGG
AD2 <- AD_HC[!is.na(AD_HC$entrezgene_id),]
AD2 <- AD2[!duplicated(AD2$entrezgene_id),]

# Create a vector with the metric
kegg_gene_list_AD <- AD2$metric

# Name vector with ENTREZ ids
names(kegg_gene_list_AD) <- AD2$entrezgene_id

# omit any NA values 
kegg_gene_list_AD <- na.omit(kegg_gene_list_AD)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list_AD <- sort(kegg_gene_list_AD, decreasing = TRUE)

gseKEGG_AD <- gseKEGG(geneList = kegg_gene_list_AD,
                      by = "fgsea",
                      organism = "hsa",
                      nPermSimple = 10000,
                      minGSSize = 40,
                      maxGSSize = 1000,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "fdr",
                      keyType = "kegg",
                      eps = 0,
                      verbose = TRUE)

# Dotplot: Kyoto Encyclopedia of Genes and Genomes (KEGG)
png(filename = "NPJ_dementia/3_Figures/7_Functional_pathway_analyses/AD/KEGG_AD_dotplot.png",
    width     = 30,
    height    = 50,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

dotplot(gseKEGG_AD, showCategory = 10, split = ".sign") + 
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 20, color = "black"),
        strip.background = element_rect(fill = "lightgray", colour = "black", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.5, units = "cm")) +
  ggtitle("KEGG - Alzheimer's disease") + 
  facet_grid(.~.sign, labeller = labeller(.sign = 
                                            c(activated = "Upregulated",
                                              suppressed = "Downregulated")))

dev.off()

## GSEA with Disease Ontology (DO)
gseDO_AD <- gseDO(kegg_gene_list_AD,
                  by = "fgsea",
                  nPermSimple = 10000,
                  minGSSize = 40,
                  maxGSSize = 1000,
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "fdr",
                  eps = 0,
                  verbose = TRUE)

# Dotplot: Disease Ontology (DO)
png(filename = "NPJ_dementia/3_Figures/7_Functional_pathway_analyses/AD/DO_AD_dotplot.png",
    width     = 30,
    height    = 50,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

dotplot(gseDO_AD, showCategory = 20, split = ".sign") + 
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 20, color = "black"),
        strip.background = element_rect(fill = "lightgray", colour = "black", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.5, units = "cm")) +
  ggtitle("DO - Alzheimer's disease") + 
  facet_grid(.~.sign, labeller = labeller(.sign = 
                                            c(activated = "Upregulated",
                                              suppressed = "Downregulated")))

dev.off()

# Down syndrome demented GO, KEGG and DO ----------------------------------
# Load result from DESeq2 (see script 5. DEGs.R)
DSD_HC <- read.csv("NPJ_dementia/1_Datasets/DEG_results/HC_DSD_vs_Control.csv", header = TRUE)
colnames(DSD_HC)[1] <- "ID"
DSD_HC <- merge(DSD_HC, GeneInfo, by.x = "ID", by.y = "ID")
DSD_HC <- DSD_HC[order(DSD_HC$logFC, decreasing = TRUE),]

## Make an order in the gene list (ranking) 
## Rank the genes of the sample data results, combining the sign of the logFC and the p-value
# Create a metric to order the gene
sign <- sign(DSD_HC$logFC)
logP <- -log2(DSD_HC$P.Val)
metric <- logP*sign

# Add to the data frame
DSD_HC$metric <- metric

# Extract the metric to create a gene list (logFC)
gene_list_DSD <- DSD_HC$metric
names(gene_list_DSD) <- DSD_HC$ensembl_gene_id
gene_list_DSD <- sort(gene_list_DSD, decreasing = TRUE)

# Check distribution of gene_list_DSD
plot(gene_list_DSD)
barplot(sort(gene_list_DSD, decreasing = TRUE), xaxt = "n")
max(gene_list_DSD)
min(gene_list_DSD)

## GSEA with Gene Ontology: Biological Processes (BP)
gseDSD_BP <- gseGO(geneList = gene_list_DSD, 
                     by = "fgsea",
                     ont = "BP",
                     nPermSimple = 1000,
                     keyType = "ENSEMBL",
                     minGSSize = 3, 
                     maxGSSize = 10000, 
                     eps = 0,
                     pvalueCutoff = 0.05, 
                     verbose = TRUE, 
                     OrgDb = "org.Hs.eg.db", 
                     pAdjustMethod = "fdr")

# Dotplot: Gene Ontology: Biological Processes (BP)
png(filename = "NPJ_dementia/3_Figures/7_Functional_pathway_analyses/DSD/GO_BP_DSD_dotplot.png",
    width     = 30,
    height    = 50,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

dotplot(gseDSD_BP, showCategory = 10, split = ".sign") + 
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 20, color = "black"),
        strip.background = element_rect(fill = "lightgray", colour = "black", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.5, units = "cm")) +
  ggtitle("GO:BP - Down syndrome demented") + 
  facet_grid(.~.sign, labeller = labeller(.sign = 
                                            c(activated = "Upregulated",
                                              suppressed = "Downregulated")))

dev.off()

# Enrichment Map: Gene Ontology: Biological Processes (BP)
GO_BP_DSD_dotplot2 <- dotplot(gseDSD_BP, showCategory = 30, split = ".sign") + facet_grid(.~.sign)
terms <- as.character(GO_BP_DSD_dotplot2$data$Description)
GO_BP_DSD_dotplot <- pairwise_termsim(gseDSD_BP, showCategory = dim(gseDSD_BP)[1] )

# Enrichment Map: UPREGULATED
png(filename = "NPJ_dementia/3_Figures/7_Functional_pathway_analyses/DSD/GO_BP_DSD_UP_emapplot.png",
    width     = 40,
    height    = 30,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

emapplot(GO_BP_DSD_dotplot, showCategory = terms[1:30], size_category = 1, 
         layout = "kk", color = "p.adjust") + 
  theme(strip.text = element_text(face = "bold", size = 10, color = "black"),
        strip.background = element_rect(fill = "lightgray", colour = "black", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text = element_blank(),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.5, units = "cm")) +
  ggtitle("Upregulated GO:BP - Enrichment map Down syndrome demented") 

dev.off()

# Enrichment Map: DOWNREGULATED
png(filename = "NPJ_dementia/3_Figures/7_Functional_pathway_analyses/DSD/GO_BP_DSD_DOWN_emapplot.png",
    width     = 40,
    height    = 30,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

emapplot(GO_BP_DSD_dotplot, showCategory = terms[31:50], size_category = 1, 
         layout = "kk", color = "p.adjust") +
  theme(strip.text = element_text(face = "bold", size = 10, color = "black"),
        strip.background = element_rect(fill = "lightgray", colour = "black", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text = element_blank(), 
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.5, units = "cm")) +
  ggtitle("Downregulated GO:BP - Enrichment map Down syndrome demented") 

dev.off()

## GSEA with Gene Ontology: Molecular Function (MF)
gseDSD_MF <- gseGO(geneList = gene_list_DSD, 
                     by = "fgsea",
                     ont = "MF", 
                     keyType = "ENSEMBL",
                     nPermSimple = 1000,
                     minGSSize = 3, 
                     maxGSSize = 10000, 
                     eps = 0,
                     pvalueCutoff = 0.05, 
                     verbose = TRUE, 
                     OrgDb = "org.Hs.eg.db", 
                     pAdjustMethod = "fdr")

# Dotplot: Gene Ontology:Molecular Function
png(filename = "NPJ_dementia/3_Figures/7_Functional_pathway_analyses/DSD/GO_MF_DSD_dotplot.png",
    width     = 30,
    height    = 50,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

dotplot(gseDSD_MF, showCategory = 13, split = ".sign") + 
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 20, color = "black"),
        strip.background = element_rect(fill = "lightgray", colour = "black", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.5, units = "cm")) +
  ggtitle("GO:MF - Down syndrome demented") +  
  facet_grid(.~.sign, labeller = labeller(.sign = 
                                            c(activated = "Upregulated",
                                              suppressed = "Downregulated")))

dev.off()

## GSEA with Kyoto Encyclopedia of Genes and Genomes (KEGG)
# Extract entrez gene IDs for KEGG
DSD2 <- DSD_HC[!is.na(DSD_HC$entrezgene_id),]
DSD2 <- DSD2[!duplicated(DSD2$entrezgene_id),]

# Create a vector with the metric
kegg_gene_list_DSD <- DSD2$metric

# Name vector with ENTREZ ids
names(kegg_gene_list_DSD) <- DSD2$entrezgene_id

# omit any NA values 
kegg_gene_list_DSD <- na.omit(kegg_gene_list_DSD)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list_DSD <- sort(kegg_gene_list_DSD, decreasing = TRUE)

gseKEGG_DSD <- gseKEGG(geneList = kegg_gene_list_DSD,
                       by = "fgsea",
                       organism = "hsa",
                       nPermSimple = 10000,
                       minGSSize = 40,
                       maxGSSize = 1000,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "fdr",
                       keyType = "kegg",
                       eps = 0,
                       verbose = TRUE)

# Dotplot: Kyoto Encyclopedia of Genes and Genomes (KEGG)
png(filename = "NPJ_dementia/3_Figures/7_Functional_pathway_analyses/DSD/KEGG_DSD_dotplot.png",
    width     = 30,
    height    = 50,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

dotplot(gseKEGG_DSD, showCategory = 10, split = ".sign") + 
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 20, color = "black"),
        strip.background = element_rect(fill = "lightgray", colour = "black", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.5, units = "cm")) +
  ggtitle("KEGG - Down syndrome demented") + 
  facet_grid(.~.sign, labeller = labeller(.sign = 
                                            c(activated = "Upregulated",
                                              suppressed = "Downregulated")))

dev.off()

## GSEA with Disease Ontology (DO)
gseDO_DSD <- gseDO(kegg_gene_list_DSD,
                     by = "fgsea",
                     nPermSimple = 10000,
                     minGSSize = 40,
                     maxGSSize = 1000,
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "fdr",
                     eps = 0,
                     verbose = TRUE)

# Dotplot: Disease Ontology (DO)
png(filename = "NPJ_dementia/3_Figures/7_Functional_pathway_analyses/DSD/DO_DSD_dotplot.png",
    width     = 30,
    height    = 50,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

dotplot(gseDO_DSD, showCategory = 10, split = ".sign") + 
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 20, color = "black"),
        strip.background = element_rect(fill = "lightgray", colour = "black", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.5, units = "cm")) +
  ggtitle("DO - Down syndrome demented") + 
  facet_grid(.~.sign, labeller = labeller(.sign = 
                                            c(activated = "Upregulated",
                                              suppressed = "Downregulated")))

dev.off()

# Compare Clusters PDD, AD and DSD ----------------------------------------
# ENSEMBL IDs LIST
PDD_list_ensembl <- gene_list_PDD[abs(gene_list_PDD) > 5]
AD_list_ensembl <- gene_list_AD[abs(gene_list_AD) > 5]
DSD_list_ensembl <- gene_list_DSD[abs(gene_list_DSD) > 5]

length(PDD_list_ensembl)
length(AD_list_ensembl)
length(DSD_list_ensembl)

list_ensembl <- list(PDD = gene_list_PDD,
             AD = gene_list_AD,
             DSD = gene_list_DSD)

# ENTREZ IDs LIST
PDD_list_entrez <- kegg_gene_list_PDD[abs(kegg_gene_list_PDD) > 5]
AD_list_entrez <- kegg_gene_list_AD[abs(kegg_gene_list_AD) > 5]
DSD_list_entrez <- kegg_gene_list_DSD[abs(kegg_gene_list_DSD) > 5]

length(PDD_list_entrez)
length(AD_list_entrez)
length(DSD_list_entrez)

list_entrez <- list(PDD = kegg_gene_list_PDD,
                    AD = kegg_gene_list_AD,
                    DSD = kegg_gene_list_DSD)

# DEMENTIA CLUSTER - GO:BP
ck1 <- compareCluster(geneCluster = list_ensembl, fun = gseGO,
                      ont = "BP",
                      OrgDb = "org.Hs.eg.db",
                      keyType = "ENSEMBL",
                      by = "fgsea",
                      nPermSimple = 10000,
                      eps = 0,
                      minGSSize = 50,
                      maxGSSize = 1000,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "fdr",
                      verbose = TRUE)

png(filename = "NPJ_dementia/3_Figures/7_Functional_pathway_analyses/Dementia_GO_BP_dotplot.png",
    width     = 30,
    height    = 50,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

dotplot(ck1, showCategory = 5, split = ".sign") +
  theme(strip.text = element_text(face = "bold", size = 22, color = "black"),
        strip.background = element_rect(fill = "lightgray", colour = "black", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18),
        legend.key.size = unit(0.5, units = "cm")) +
  ggtitle("GO:BP - Dementia") + 
  facet_grid(.~.sign, labeller = labeller(.sign = 
                                            c(activated = "Upregulated",
                                              suppressed = "Downregulated")))

dev.off()

# DEMENTIA CLUSTER - GO:MF
ck2 <- compareCluster(geneCluster = list_ensembl, fun = gseGO,
                      ont = "MF",
                      OrgDb = "org.Hs.eg.db",
                      keyType = "ENSEMBL",
                      by = "fgsea",
                      nPermSimple = 10000,
                      eps = 0,
                      minGSSize = 50,
                      maxGSSize = 1000,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "fdr",
                      verbose = TRUE)

png(filename = "NPJ_dementia/3_Figures/7_Functional_pathway_analyses/Dementia_GO_MF_dotplot.png",
    width     = 30,
    height    = 50,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

dotplot(ck2, showCategory = 5, split = ".sign") +
  theme(strip.text = element_text(face = "bold", size = 22, color = "black"),
        strip.background = element_rect(fill = "lightgray", colour = "black", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18),
        legend.key.size = unit(0.5, units = "cm")) +
  ggtitle("GO:MF - Dementia") + 
  facet_grid(.~.sign, labeller = labeller(.sign = 
                                            c(activated = "Upregulated",
                                              suppressed = "Downregulated")))

dev.off()

# DEMENTIA CLUSTER - GO:CC
ck3 <- compareCluster(geneCluster = list_ensembl, fun = gseGO,
                      ont = "CC",
                      OrgDb = "org.Hs.eg.db",
                      keyType = "ENSEMBL",
                      by = "fgsea",
                      nPermSimple = 10000,
                      eps = 0,
                      minGSSize = 50,
                      maxGSSize = 1000,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "fdr",
                      verbose = TRUE)

png(filename = "NPJ_dementia/3_Figures/7_Functional_pathway_analyses/Dementia_GO_CC_dotplot.png",
    width     = 30,
    height    = 50,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

dotplot(ck3, showCategory = 5, split = ".sign") +
  theme(strip.text = element_text(face = "bold", size = 22, color = "black"),
        strip.background = element_rect(fill = "lightgray", colour = "black", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18),
        legend.key.size = unit(0.5, units = "cm")) +
  ggtitle("GO:CC - Dementia") + 
  facet_grid(.~.sign, labeller = labeller(.sign = 
                                            c(activated = "Upregulated",
                                              suppressed = "Downregulated")))

dev.off()

# DEMENTIA CLUSTER - KEGG
ck4 <- compareCluster(geneCluster = list_entrez, fun = gseKEGG,
                      by = "fgsea",
                      nPermSimple = 100000,
                      eps = 0,
                      minGSSize = 40,
                      maxGSSize = 1000,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "fdr",
                      verbose = TRUE)

png(filename = "NPJ_dementia/3_Figures/7_Functional_pathway_analyses/Dementia_KEGG_dotplot.png",
    width     = 30,
    height    = 50,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

dotplot(ck4, showCategory = 5, split = ".sign") +
  theme(strip.text = element_text(face = "bold", size = 22, color = "black"),
        strip.background = element_rect(fill = "lightgray", colour = "black", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18),
        legend.key.size = unit(0.5, units = "cm")) +
  ggtitle("KEGG - Dementia") + 
  facet_grid(.~.sign, labeller = labeller(.sign = 
                                            c(activated = "Upregulated",
                                              suppressed = "Downregulated")))

dev.off()

# DEMENTIA CLUSTER - DISEASE ONTOLOGY
ck5 <- compareCluster(geneCluster = list_entrez, fun = gseDO,
                      by = "fgsea",
                      nPermSimple = 100000,
                      eps = 0,
                      minGSSize = 40,
                      maxGSSize = 1000,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "fdr",
                      verbose = TRUE)

png(filename = "NPJ_dementia/3_Figures/7_Functional_pathway_analyses/Dementia_DO_dotplot.png",
    width     = 30,
    height    = 50,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

dotplot(ck5, showCategory = 6, split = ".sign") +
  theme(strip.text = element_text(face = "bold", size = 22, color = "black"),
        strip.background = element_rect(fill = "lightgray", colour = "black", size = 0.5), 
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18),
        legend.key.size = unit(0.5, units = "cm")) +
  ggtitle("Disease Ontology - Dementia") + 
  facet_grid(.~.sign, labeller = labeller(.sign = 
                                            c(activated = "Upregulated",
                                              suppressed = "Downregulated")))

dev.off()

# DisGeNet Common DEGs ----------------------------------------------------
# Extract common differential expressed genes
DEGs <- read.csv("NPJ_dementia/3_Figures/6_Common expressed_genes/DEGs_Genes.csv", header = TRUE)
myListOfGenes <- as.character(DEGs$X)

results <- gene2disease(
  gene     = myListOfGenes,
  vocabulary = "HGNC",
  database = "CURATED",
  score =c(0.1, 1),
  verbose  = TRUE
)

# Heatmap: DisGeNet (DiseaseClass)
png(filename = "NPJ_dementia/3_Figures/7_Functional_pathway_analyses/HEATMAP_DEGs.png",
    width     = 35,
    height    = 20,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

gene_disease <- plot(results, 
     type = "Heatmap",
     class = "DiseaseClass", 
     nchars = 100, 
     interactive = FALSE)

gene_disease +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 26, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.key.size = unit(1, units = "cm"))

dev.off()

# Heatmap: DisGeNet 
png(filename = "NPJ_dementia/3_Figures/7_Functional_pathway_analyses/HEATMAP2_DEGs.png",
    width     = 30,
    height    = 40,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

plot(results, 
     type = "Heatmap",
     limit = 200, 
     nchars = 50, 
     interactive = FALSE)

dev.off()

# Network: DisGeNet 
png(filename = "NPJ_dementia/3_Figures/7_Functional_pathway_analyses/NETWORK_DEGs.png",
    width     = 100,
    height    = 70,
    units     = "cm",
    res       = 800,
    pointsize = 25)

disgenet_igraph <- plot(results,
     type = "Network",
     prop = 7,
     eprop = 10, 
     nchars = 3,
     limit = 140,
     verbose = TRUE,
     layout = "layout.kamada.kawai")


createNetworkFromIgraph(
  results,
  title = "DisGeNET Network",
  collection = "DisGeNET"
)
dev.off()
