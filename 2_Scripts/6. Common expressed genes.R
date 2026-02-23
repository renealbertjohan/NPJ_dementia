# 6. Common differential expressed genes (Venn diagrams)
# Load libraries
library(VennDiagram)

# Step 1: Load phenotype data
phenoHEROES <- read.delim("NPJ_dementia/1_Datasets/PhenoData/phenoHEROES_RNAAge_FPKM.csv", sep = ",", header = TRUE)
row.names(phenoHEROES) <- phenoHEROES[,1]
phenoHEROES$X <- NULL

# Step 2: Load in gene information ----------------------------------------
genedata <- read.delim("NPJ_dementia/1_Datasets/Counts_Matrices/counts_length.csv", sep = ",", header = TRUE)
genedata$X.1 <- NULL
colnames(genedata)[1] <- "Geneid" 
genedata <- genedata[,-c(2:21)]

# Step 3: Load previously filtered DEGs results (see script DEGs)
# Parkinson's disease demented
HC_PDD_vs_Control <- read.delim("NPJ_dementia/1_Datasets/DEG_results/PDD_vs_control_with_Genedata.csv", sep = ",", header = TRUE)
row.names(HC_PDD_vs_Control) <- HC_PDD_vs_Control[,1]
colnames(HC_PDD_vs_Control)[1] <- "ID"
dim(HC_PDD_vs_Control) # 2897 Differential expressed genes

# Alzheimer's disease
HC_AD_vs_Control <- read.delim("NPJ_dementia/1_Datasets/DEG_results/AD_vs_control_with_Genedata.csv", sep = ",", header = TRUE)
row.names(HC_AD_vs_Control) <- HC_AD_vs_Control[,1]
colnames(HC_AD_vs_Control)[1] <- "ID"
dim(HC_AD_vs_Control) # 657 Differential expressed genes

# Down syndrome demented
HC_DSD_vs_Control <- read.delim("NPJ_dementia/1_Datasets/DEG_results/DSD_vs_control_with_Genedata.csv", sep = ",", header = TRUE)
row.names(HC_DSD_vs_Control) <- HC_DSD_vs_Control[,1]
colnames(HC_DSD_vs_Control)[1] <- "ID"
dim(HC_DSD_vs_Control) # 196 Differential expressed genes

# Step 4: Venn diagrams ---------------------------------------------------
# Venn diagram of UPREGULATED genes
HC_PDD_vs_Control_adj_p_log <- HC_PDD_vs_Control[HC_PDD_vs_Control$adj.P.Val <= 0.05 & HC_PDD_vs_Control$logFC >= 0.3,]
UP_HC_PDD <- HC_PDD_vs_Control_adj_p_log[HC_PDD_vs_Control_adj_p_log$DifExpressed == "UP",]
length(UP_HC_PDD$ID)

HC_AD_vs_Control_adj_p_log <- HC_AD_vs_Control[HC_AD_vs_Control$adj.P.Val <= 0.05 & HC_AD_vs_Control$logFC >= 0.3,]
UP_HC_AD <- HC_AD_vs_Control_adj_p_log[HC_AD_vs_Control_adj_p_log$DifExpressed == "UP",]
length(UP_HC_AD$ID)

HC_DSD_vs_Control_adj_p_log <- HC_DSD_vs_Control[HC_DSD_vs_Control$adj.P.Val <= 0.05 & HC_DSD_vs_Control$logFC >= 0.3,]
UP_HC_DSD <- HC_DSD_vs_Control_adj_p_log[HC_DSD_vs_Control_adj_p_log$DifExpressed == "UP",]
length(UP_HC_DSD$ID)

# Select the common Upregulated DEGs 
my_list_up <- list(UP_HC_PDD$ID, UP_HC_AD$ID, UP_HC_DSD$ID)
names(my_list_up) <- c("PDD", "AD", "DSD")
subset(stack(my_list_up), duplicated(values))$values
names(which(table(unlist(my_list_up)) > 2))
length(names(which(table(unlist(my_list_up)) > 2)))

my_list_up12 <- list(UP_HC_PDD$ID, UP_HC_AD$ID)
length(names(which(table(unlist(my_list_up12)) > 1)))
my_list_up23 <- list(UP_HC_AD$ID, UP_HC_DSD$ID)
length(names(which(table(unlist(my_list_up23)) > 1)))
my_list_up13 <- list(UP_HC_PDD$ID, UP_HC_DSD$ID)
length(names(which(table(unlist(my_list_up13)) > 1)))

my_list_up123 <- list(UP_HC_PDD$ID, UP_HC_AD$ID, UP_HC_DSD$ID)
length(names(which(table(unlist(my_list_up123)) > 2)))

# Venn diagram: Upregulated
png(filename="NPJ_dementia/3_Figures/6_Common expressed_genes/Venn_diagram_Upregulated.png",
    width     = 20,
    height    = 10,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

draw.triple.venn(area1 = 1296, area2 = 357, area3 = 133, 
                 n12 = 270, n23 = 42, n13 = 51,
                 n123 = 29,
                 category=c("Parkinson's disease dementia",
                            "Alzheimer's disease", "Down syndrome dementia"),
                 fill = c("navyblue", "gold3", "darkred"),
                 cex = 2.5, fontfamily = "Helvetica", fontface = "bold",
                 cat.col = c("navyblue", "gold3", "darkred"), cat.fontface = "bold", 
                 cat.cex = 2.5, cat.dist = 0.09, cat.fontfamily = "Helvetica",
                 margin = 0.05, lty = 1, alpha = 0.5, lwd = 0.5)
dev.off()

# Venn diagram of DOWNREGULATED genes
HC_PDD_vs_Control_adj_p_log <- HC_PDD_vs_Control[HC_PDD_vs_Control$adj.P.Val <= 0.05 & HC_PDD_vs_Control$logFC <= 0.3,]
DOWN_HC_PDD <- HC_PDD_vs_Control_adj_p_log[HC_PDD_vs_Control_adj_p_log$DifExpressed == "DOWN",]
length(DOWN_HC_PDD$ID)

HC_AD_vs_Control_adj_p_log <- HC_AD_vs_Control[HC_AD_vs_Control$adj.P.Val <= 0.05 & HC_AD_vs_Control$logFC <= 0.3,]
DOWN_HC_AD <- HC_AD_vs_Control_adj_p_log[HC_AD_vs_Control_adj_p_log$DifExpressed == "DOWN",]
length(DOWN_HC_AD$ID)

HC_DSD_vs_Control_adj_p_log <- HC_DSD_vs_Control[HC_DSD_vs_Control$adj.P.Val <= 0.05 & HC_DSD_vs_Control$logFC <= 0.3,]
DOWN_HC_DSD <- HC_DSD_vs_Control_adj_p_log[HC_DSD_vs_Control_adj_p_log$DifExpressed == "DOWN",]
length(DOWN_HC_DSD$ID)

# Select the common Downregulated DEGs 
my_list_down <- list(DOWN_HC_PDD$ID, DOWN_HC_AD$ID, DOWN_HC_DSD$ID)
names(my_list_down) <- c("PDD", "AD", "DSD")
subset(stack(my_list_down), duplicated(values))$values
names(which(table(unlist(my_list_down)) > 2))
length(names(which(table(unlist(my_list_down)) > 2)))

my_list_down12 <- list(DOWN_HC_PDD$ID, DOWN_HC_AD$ID)
length(names(which(table(unlist(my_list_down12)) > 1)))
my_list_down23 <- list(DOWN_HC_AD$ID, DOWN_HC_DSD$ID)
length(names(which(table(unlist(my_list_down23)) > 1)))
my_list_down13 <- list(DOWN_HC_PDD$ID, DOWN_HC_DSD$ID)
length(names(which(table(unlist(my_list_down13)) > 1)))

my_list_down123 <- list(DOWN_HC_PDD$ID, DOWN_HC_AD$ID, DOWN_HC_DSD$ID)
length(names(which(table(unlist(my_list_down123)) > 2)))

# Venn diagram: Downregulated
png(filename="NPJ_dementia/3_Figures/6_Common expressed_genes/Venn_diagram_Downregulated.png",
    width     = 20,
    height    = 10,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

draw.triple.venn(area1 = 1601, area2 = 300, area3 = 63, 
                 n12 = 198, n23 = 21, n13 = 23,
                 n123 = 16,
                 category=c("Parkinson's disease dementia",
                            "Alzheimer's disease", "Down syndrome dementia"),
                 fill = c("navyblue", "gold3", "darkred"),
                 cex = 2.5, fontfamily = "Helvetica", fontface = "bold",
                 cat.col = c("navyblue", "gold3", "darkred"), cat.fontface = "bold", 
                 cat.cex = 2.5, cat.dist = 0.09, cat.fontfamily = "Helvetica",
                 margin = 0.05, lty = 1, alpha = 0.5, lwd = 0.5)
dev.off()

# Venn diagram of ALL genes
HC_PDD_vs_Control_adj_p_log <- HC_PDD_vs_Control[HC_PDD_vs_Control$adj.P.Val <= 0.05,]
length(HC_PDD_vs_Control_adj_p_log$ID)

HC_AD_vs_Control_adj_p_log <- HC_AD_vs_Control[HC_AD_vs_Control$adj.P.Val <= 0.05,]
length(HC_AD_vs_Control_adj_p_log$ID)

HC_DSD_vs_Control_adj_p_log <- HC_DSD_vs_Control[HC_DSD_vs_Control$adj.P.Val <= 0.05,]
length(HC_DSD_vs_Control_adj_p_log$ID)

# Select the ALL common DEGs 
my_list <- list(HC_PDD_vs_Control_adj_p_log$ID, HC_AD_vs_Control_adj_p_log$ID, HC_DSD_vs_Control_adj_p_log$ID)
names(my_list) <- c("PDD", "AD", "DS")
subset(stack(my_list), duplicated(values))$values
names(which(table(unlist(my_list)) > 2))
length(names(which(table(unlist(my_list)) > 2)))

my_list12 <- list(HC_PDD_vs_Control_adj_p_log$ID, HC_AD_vs_Control_adj_p_log$ID)
length(names(which(table(unlist(my_list12)) > 1)))
my_list23 <- list(HC_AD_vs_Control_adj_p_log$ID, HC_DSD_vs_Control_adj_p_log$ID)
length(names(which(table(unlist(my_list23)) > 1)))
my_list13 <- list(HC_PDD_vs_Control_adj_p_log$ID, HC_DSD_vs_Control_adj_p_log$ID)
length(names(which(table(unlist(my_list13)) > 1)))

my_list123 <- list(HC_PDD_vs_Control_adj_p_log$ID, HC_AD_vs_Control_adj_p_log$ID, HC_DSD_vs_Control_adj_p_log$ID)
length(names(which(table(unlist(my_list123)) > 2)))

# Venn diagram: ALL
png(filename="NPJ_dementia/3_Figures/6_Common expressed_genes/Venn_diagram_ALL.png",
    width     = 20,
    height    = 10,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

draw.triple.venn(area1 = 2897, area2 = 657, area3 = 196, 
                 n12 = 468, n23 = 63, n13 = 74,
                 n123 = 45,
                 category=c("Parkinson's disease dementia",
                            "Alzheimer's disease", "Down syndrome dementia"),
                 fill = c("navyblue", "gold3", "darkred"),
                 cex = 2.5, fontfamily = "Helvetica", fontface = "bold",
                 cat.col = c("navyblue", "gold3", "darkred"), cat.fontface = "bold", 
                 cat.cex = 2.5, cat.dist = 0.09, cat.fontfamily = "Helvetica",
                 margin = 0.05, lty = 1, alpha = 0.5, lwd = 0.5)

dev.off()

# Load the information about the genes
GeneInfo <- read.csv("NPJ_dementia/1_Datasets/Counts_Matrices/counts_length.csv", header = TRUE, row.names = 1)
colnames(GeneInfo)[1] <- "ID"
GeneInfo <- GeneInfo[,-c(2:21)]

# Write away the interesting common differential expressed genes 
info_up <- HC_PDD_vs_Control[HC_PDD_vs_Control$ID %in% names(which(table(unlist(my_list_up123)) > 2)),]
write.csv(info_up, "NPJ_dementia/3_Figures/6_Common expressed_genes/Upregulated_Genes.csv")

info_down <- HC_PDD_vs_Control[HC_PDD_vs_Control$ID %in% names(which(table(unlist(my_list_down123)) > 2)),]
write.csv(info_down, "NPJ_dementia/3_Figures/6_Common expressed_genes/Downregulated_Genes.csv")

info_all <- HC_PDD_vs_Control[HC_PDD_vs_Control$ID %in% names(which(table(unlist(my_list123)) > 2)),]
write.csv(info_all, "NPJ_dementia/3_Figures/6_Common expressed_genes/DEGs_Genes.csv")

# Write away differential expressed genes information for manuscript
gene_up <- GeneInfo[GeneInfo$ID %in% names(which(table(unlist(my_list_up123)) > 2)),]
gene_up <- gene_up[order(gene_up$chromosome_name, gene_up$start_position),]
write.csv(gene_up, "NPJ_dementia/3_Figures/6_Common expressed_genes/Upregulated_Table.csv")

gene_down <- GeneInfo[GeneInfo$ID %in% names(which(table(unlist(my_list_down123)) > 2)),]
gene_down <- gene_down[order(gene_down$chromosome_name, gene_down$start_position),]
write.csv(gene_down, "NPJ_dementia/3_Figures/6_Common expressed_genes/Downregulated_Table.csv")

gene_all <- GeneInfo[GeneInfo$ID %in% names(which(table(unlist(my_list123)) > 2)),]
gene_all <- gene_all[order(gene_all$chromosome_name, gene_all$start_position),]
write.csv(gene_all, "NPJ_dementia/3_Figures/6_Common expressed_genes/All_Table.csv")
