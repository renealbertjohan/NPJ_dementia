# 2. RNAAgeCalc: A multi-tissue transcriptional age calculator
# Load libraries
library(biomaRt)
library(RNAAgeCalc)
library(tidyverse)
library(ggpubr) # For function stat_compare_means
library(ggh4x) # For function facet_grid2
library(ggtext) # For adding text to figure, if wanted/needed
library(rstatix) # For tukeyhsd() function in figures
library(car) # qqPlot function for Two-way ANOVA

# LOAD DATA ---------------------------------------------------------------
# Load phenotype data
phenoHEROES <- read.delim("NPJ_dementia/1_Datasets/PhenoData/phenoHEROES.csv", sep = ",", dec = ",", header = TRUE)
row.names(phenoHEROES) <- phenoHEROES$Tube_code
phenoHEROES$X <- NULL

# Load counts matrix
countsHEROES <- read.delim("NPJ_dementia/1_Datasets/Counts.matrix.csv", sep = ",", header = TRUE)
row.names(countsHEROES) <- countsHEROES$X
countsHEROES$X <- NULL

# EXTRACT GENE DATA WITH BIOMART (NEW STYLE) ------------------------------
# Step 1: Selecting an Ensembl BioMart database and dataset ---------------
# I tried option = , but it seems that the annotation in org.Hs.eg.db is not as updated as biomaRt
biomartCacheClear() # Remove previous queries from your local computer

# Select mart
ensembl_114 <- useEnsembl("genes", 
                          dataset = "hsapiens_gene_ensembl", 
                          version = 114) # Select version 114


# Step 2: Query gene information from counts matrix -----------------------
# Query gene data for the Counts Matrix, such chromosome, start, end, description
genes <- rownames(countsHEROES) # Mixture of ENSEMBL IDs and HUGO symbols

## Select HUGO symbols from counts matrix
hugo <- genes[-grep("ENSG", genes)]

## Select ENSEMBL IDs from counts matrix
ens <- genes[grep("ENSG", genes)]

# HUGO + ENSEMBL
length(hugo) + length(ens) # 60645 genes

# Query the ENSEMBL IDs for the HUGO symbols
genes_hugo <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id", 
                                   "chromosome_name", "start_position", "end_position", 
                                   "strand", "percentage_gene_gc_content", "description"), 
                    filters = c("hgnc_symbol"), 
                    values = list(hugo), 
                    mart = ensembl_114,
                    uniqueRows = TRUE)

# Query the transcript/cds and ENSEMBL IDs for the HUGO symbols
genes_length_hugo <- getBM(attributes = c("ensembl_gene_id", "transcript_length", "cds_length"), 
                    filters = c("hgnc_symbol"), 
                    values = list(hugo), 
                    mart = ensembl_114,
                    uniqueRows = TRUE)

# Delete duplicated ENSEMBL IDs from the query (genes_length_hugo)
genes_length_hugo <- genes_length_hugo[!duplicated(genes_length_hugo$ensembl_gene_id),]

# Merge the two HUGO biomaRt queries (i.e., feature page and sequence page) with each other
HUGO_genes <- merge(genes_hugo, genes_length_hugo, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")

# Delete duplicated ENSEMBL IDs and HUGO symbols from the query
genes_hugo <- HUGO_genes[!duplicated(HUGO_genes$ensembl_gene_id) & !duplicated(HUGO_genes$hgnc_symbol),]

# Query the HUGO symbols for the ENSEMBL IDs (OTHERWAY AROUND)
genes_ens <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id",
                                "chromosome_name", "start_position", "end_position", 
                                "strand", "percentage_gene_gc_content", "description"), 
                   filters = c("ensembl_gene_id"),
                   values = list(ens), 
                   mart = ensembl_114,
                   uniqueRows = TRUE)

# Query the transcript/cds and ENSEMBL IDs for the HUGO symbols
genes_length_ens <- getBM(attributes = c("ensembl_gene_id", "transcript_length", "cds_length"), 
                           filters = c("ensembl_gene_id"), 
                           values = list(ens), 
                           mart = ensembl_114,
                           uniqueRows = TRUE)

# Delete duplicated ENSEMBL IDs from the query (genes_length_hugo)
genes_length_ens <- genes_length_ens[!duplicated(genes_length_ens$ensembl_gene_id),]

# Merge the two ENSEMBL biomaRt queries (i.e., feature page and sequence page) with each other
ENS_genes <- merge(genes_ens, genes_length_ens, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")

# Select ENSEMBL IDs that contain a HUGO symbol
hgn <- which(ENS_genes$hgnc_symbol != "")
ensemblWITH <- ENS_genes[hgn,]
ensemblWITH  <- ensemblWITH [!duplicated(ensemblWITH$ensembl_gene_id) & !duplicated(ensemblWITH$hgnc_symbol),]

# Select ENSEMBL IDs that does not contain a HUGO symbol
ensemblWITHOUT <- ENS_genes[-hgn,]
ensemblWITHOUT  <- ensemblWITHOUT[!duplicated(ensemblWITHOUT$ensembl_gene_id),]

ensemblWITHOUT <- cbind(ensemblWITHOUT[,1], ensemblWITHOUT[,1], ensemblWITHOUT[,3], 
                        ensemblWITHOUT[,4], ensemblWITHOUT[,5], ensemblWITHOUT[,6],
                        ensemblWITHOUT[,7], ensemblWITHOUT[,8], ensemblWITHOUT[,9],
                        ensemblWITHOUT[,10], ensemblWITHOUT[,11])

ensemblWITHOUT <- as.data.frame(ensemblWITHOUT)

colnames(ensemblWITHOUT)[1] <- "ensembl_gene_id"
colnames(ensemblWITHOUT)[2] <- "hgnc_symbol"
colnames(ensemblWITHOUT)[3] <- "entrezgene_id"
colnames(ensemblWITHOUT)[4] <- "chromosome_name"
colnames(ensemblWITHOUT)[5] <- "start_position"
colnames(ensemblWITHOUT)[6] <- "end_position"
colnames(ensemblWITHOUT)[7] <- "strand"
colnames(ensemblWITHOUT)[8] <- "percentage_gene_gc_content"
colnames(ensemblWITHOUT)[9] <- "description"
colnames(ensemblWITHOUT)[10] <- "transcript_length"
colnames(ensemblWITHOUT)[11] <- "cds_length"

# Check if colnames for all data frames are the same
colnames(ensemblWITH)
colnames(ensemblWITHOUT)
colnames(genes_hugo)

QUERY <- rbind(genes_hugo, ensemblWITH, ensemblWITHOUT)

countsHEROES$X <- rownames(countsHEROES)
counts_length <- merge(countsHEROES, QUERY, by.x = "X", by.y = "hgnc_symbol")
dim(counts_length[!duplicated(counts_length$X),]) # There are 12 duplications

# Take out the duplicates
counts_length <- counts_length[!duplicated(counts_length$X),]

# Change the strand code to + and - 
counts_length$strand <- gsub("-1", "-", counts_length$strand)
counts_length$strand <- gsub("1", "+", counts_length$strand)

# Trimming the description (i.e. remove the [Source: .....])
counts_length$description <- sub("\\[Source.*", "", counts_length$description)
counts_length$description <- gsub("^\\s+|\\s+$","",counts_length$description)

# Quality check with same genes
countsHEROES[countsHEROES$X == "A1CF",]
counts_length[counts_length$X == "A1CF",]

countsHEROES[countsHEROES$X == "DYRK1A",]
counts_length[counts_length$X == "DYRK1A",]

countsHEROES[countsHEROES$X == "APP",]
counts_length[counts_length$X == "APP",]

countsHEROES[countsHEROES$X == "SOD1",]
counts_length[counts_length$X == "SOD1",]

# Data object counts_length contains the counts and the gene information
write.csv(counts_length, "NPJ_dementia/1_Datasets/Counts_Matrices/counts_length.csv")

# Step 3: RNAAge Calculation ----------------------------------------------
# Option 1: RNAAgeCalc with FPKM ------------------------------------------
# To perform FPKM (for paired-end reads) or RPKM (for single-end reads), 
# We first divide by the library size and then by gene length (i.e., transcript_length). 
# Notice that the sum of each sample after FPKM normalization is different.
FPKM <- counts_length
row.names(FPKM) <- FPKM$X
FPKM <- FPKM[,c(2:21,30)]
FPKM$transcript_length <- as.numeric(FPKM$transcript_length) # Have to make the transcript_length numeric

# Step 1: normalize for read depth and multiply by million
FKPM2 <- apply(FPKM[,1:20], 2, function(x) x / sum(x) * 10^6) 

# Step 2. scale by gene length and multiply by thousand
countsFPKM <- FKPM2 / FPKM$transcript_length * 10^3

# Check if rownames(phenoHEROES) and colnames(countsFPKM) are identical, including order
identical(phenoHEROES$Tube_code, colnames(countsFPKM))

# Extract the chronological age
chronage <- data.frame(sampleid = colnames(countsFPKM), age = phenoHEROES$Age)

# Calculate RNAAge with predictage from the FPKM counts matrix
set.seed(123456)
resultAGE <- predict_age(exprdata = countsFPKM, tissue = "brain", exprtype = "FPKM", 
                  chronage = chronage, idtype = "SYMBOL", stype = "caucasian",
                  signature = "DESeq2", maxp = 30000) # 6.5539% genes in the gene signature are not included in the supplied gene expression (imputed).

# Option 2: RNAAgeCalc with counts ----------------------------------------
# If using 'exprtype = "counts"', the raw count will be converted to FPKM. 
# If 'genelength' is provided, the function will convert raw count to FPKM by the user-supplied gene length. 
# Otherwise, gene length is obtained from the internal database.

# Calculate RNAAge with predictage from the original counts matrix (it uses a internal data base to convert to FPKM)
COUNTS <- counts_length
row.names(COUNTS) <- COUNTS$X
COUNTS <- COUNTS[,c(2:21)]
head(COUNTS)

# Providing our own gene_length will give the same results as previously 
# (use argument genelength = gene_length in predict_age function) 
gene_length <- as.numeric(counts_length[,30])

# Check if rownames(phenoHEROES) and colnames(countsFPKM) are identical, including order
identical(phenoHEROES$Tube_code, colnames(COUNTS))

# Extract the chronological age
chronage_2 <- data.frame(Sample_id = colnames(COUNTS), Age = phenoHEROES$Age)

# Calculate RNAAge with predictage from the original counts matrix (without providing genelength)
set.seed(123456)
resultAGE_2 <- predict_age(exprdata = COUNTS, tissue = "brain", exprtype = "counts", 
                         chronage = chronage_2, idtype = "SYMBOL", stype = "caucasian",
                         signature = "DESeq2", maxp = 30000) # 6.5539% genes in the gene signature are not included in the supplied gene expression (imputed).
# Warning message:
# In count2FPKM(exprdata, genelength = genelength, idtype = idtype) :
# Can't find gene length for 44.2697% genes when converting raw count to FPKM.

# Results option 1: FPKM (biomaRt) or option 2: FPKM (internal data base)
resultAGE # Option 1
resultAGE_2 # Option 2
# Option 1 provides a better result, as more transcript lengths are obtained through biomaRt data base.
# The internal data base of the RNAAgeCalc package missed 45%!! of the gene length for the genes in our counts matrix to calculate the FKPM.
# Therefore, these lengths are computed and, thus, not reliable per se.

# Save results RNAAgeCalc
AgeDif <- resultAGE$RNAAge - resultAGE$ChronAge # Add also the calculated age difference (i.e., RNAAge - ChronAge)
phenoHEROES_RNAAge_FPKM <- as.data.frame(append(phenoHEROES, list(RNAAge = resultAGE$RNAAge, AgeDif = AgeDif), after = 3))
colnames(phenoHEROES_RNAAge_FPKM)[3] <- "ChronAge"
row.names(phenoHEROES_RNAAge_FPKM) <- phenoHEROES_RNAAge_FPKM$Tube_code

# Check if phenoHEROES_RNAAge_FPKM$Tube_code and phenoHEROES$Tube_code are identical, including order
identical(phenoHEROES_RNAAge_FPKM$Tube_code, phenoHEROES$Tube_code)
write.csv(phenoHEROES_RNAAge_FPKM, "NPJ_dementia/1_Datasets/PhenoData/phenoHEROES_RNAAge_FPKM.csv")

# Visualization of data and statistical testing
ChronAge <- as.data.frame(cbind(phenoHEROES_RNAAge_FPKM$Status, phenoHEROES_RNAAge_FPKM$Tube_code, phenoHEROES_RNAAge_FPKM$Condition, phenoHEROES_RNAAge_FPKM$ChronAge))
ChronAge$Timepoint <- rep("ChronAge", 20)
RNAAge <- as.data.frame(cbind(phenoHEROES_RNAAge_FPKM$Status, phenoHEROES_RNAAge_FPKM$Tube_code, phenoHEROES_RNAAge_FPKM$Condition, phenoHEROES_RNAAge_FPKM$RNAAge))
RNAAge$Timepoint <- rep("RNAAge", 20)
RNAAGECAL <- rbind(ChronAge, RNAAge)

colnames(RNAAGECAL)[1] <- "Status"
colnames(RNAAGECAL)[2] <- "Tube_code"
colnames(RNAAGECAL)[3] <- "Condition"
colnames(RNAAGECAL)[4] <- "Age"
RNAAGECAL$Age <- as.numeric(RNAAGECAL$Age)
RNAAGECAL$Status <- factor(RNAAGECAL$Status, levels = c("Control", "PDD", "AD", "DSD"))

# Calculations
# Control
mean(RNAAGECAL[RNAAGECAL$Status == "Control" & RNAAGECAL$Timepoint == "RNAAge", "Age"])
sd(RNAAGECAL[RNAAGECAL$Status == "Control" & RNAAGECAL$Timepoint == "RNAAge", "Age"])

# PDD
mean(RNAAGECAL[RNAAGECAL$Status == "PDD" & RNAAGECAL$Timepoint == "RNAAge", "Age"])
sd(RNAAGECAL[RNAAGECAL$Status == "PDD" & RNAAGECAL$Timepoint == "RNAAge", "Age"])

# AD
mean(RNAAGECAL[RNAAGECAL$Status == "AD" & RNAAGECAL$Timepoint == "RNAAge", "Age"])
sd(RNAAGECAL[RNAAGECAL$Status == "AD" & RNAAGECAL$Timepoint == "RNAAge", "Age"])

# DSD
mean(RNAAGECAL[RNAAGECAL$Status == "DSD" & RNAAGECAL$Timepoint == "RNAAge", "Age"])
sd(RNAAGECAL[RNAAGECAL$Status == "DSD" & RNAAGECAL$Timepoint == "RNAAge", "Age"])

## Comparison plot --> Chronological age versus Transcriptional Age
stat.test_RNAAGECAL <- aov(Age ~ Status * Timepoint, data = RNAAGECAL) %>%
  tukey_hsd()
stat.test_RNAAGECAL <- stat.test_RNAAGECAL[c(11,18,24,29),] 

png(filename="NPJ_dementia/3_Figures/2_RNAAgeCalc/Age_Calculation.png",
    width     = 50,
    height    = 15,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

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
  facet_grid2(. ~ Status, scales="free", space="free",
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
anno4 <- annotate("text", label = c("***"), size = 12, x = c(1.5), y = 120)
p + at_panel(anno2, PANEL == 2) +
  at_panel(anno3, PANEL == 3) + at_panel(anno4, PANEL == 4)

dev.off()

# Factorize
phenoHEROES_RNAAge_FPKM$Status <- factor(phenoHEROES_RNAAge_FPKM$Status, levels = c("Control", "PDD", "AD", "DSD"))

## Chronological Age
stat.test_ChronAge <- aov(ChronAge ~ Status, data = phenoHEROES_RNAAge_FPKM) %>%
  tukey_hsd()
stat.test_ChronAge

png(filename="NPJ_dementia/3_Figures/2_RNAAgeCalc/Chron_Age.png",
    width     = 20,
    height    = 15,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(data = phenoHEROES_RNAAge_FPKM, aes(x = Status, y = ChronAge, fill = factor(Status))) +
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

## Transcriptional Age
stat.test_RNAAge <- aov(RNAAge ~ Status, data = phenoHEROES_RNAAge_FPKM) %>%
  tukey_hsd()
stat.test_RNAAge

png(filename="NPJ_dementia/3_Figures/2_RNAAgeCalc/RNA_Age.png",
    width     = 20,
    height    = 15,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(data = phenoHEROES_RNAAge_FPKM, aes(x = Status, y = RNAAge, fill = factor(Status))) +
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
stat.test_AgeDif <- aov(AgeDif ~ Status, data = phenoHEROES_RNAAge_FPKM) %>%
  tukey_hsd()
stat.test_AgeDif

png(filename="NPJ_dementia/3_Figures/2_RNAAgeCalc/Age_Dif.png",
    width     = 20,
    height    = 15,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(data = phenoHEROES_RNAAge_FPKM, aes(x = Status, y = AgeDif, fill = factor(Status))) +
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
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) +
  ylab(expression(bold("Age acceleration (years)"))) +
  ggtitle(paste("Differential age")) + 
  stat_pvalue_manual(stat.test_AgeDif, hide.ns = TRUE, y.position = 50, size = 12)

dev.off()
