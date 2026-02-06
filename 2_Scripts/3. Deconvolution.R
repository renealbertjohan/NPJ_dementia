# 3. MuSiC2: Cell type deconvolution for multi-condition bulk RNA-seq data
# Load libraries
library(clusterProfiler)
library(MuSiC)
library(MuSiC2) # Deconvolution method
library(scRNAseq) # To load the snRNA-seq data from Darmanis et., 2015 for brain data
library(edgeR)
library(ggpubr)

# Load phenotype data
phenoHEROES <- read.delim("NPJ_dementia/1_Datasets/PhenoData/phenoHEROES_RNAAge_FPKM.csv", sep = ",", dec = ",", header = TRUE)
row.names(phenoHEROES) <- phenoHEROES$Tube_code
phenoHEROES$X <- NULL

# Check samples
table(phenoHEROES$Brain_region, factor(phenoHEROES$Status, levels = c("Control", "PDD", "AD", "DSD")))
table(phenoHEROES$Sex, factor(phenoHEROES$Status, levels = c("Control", "PDD", "AD", "DSD")))

# Load the original Counts Matrix
countsHERO <- read.delim("NPJ_dementia/1_Datasets/Counts.matrix.csv", sep = ",", header = TRUE)
countsHERO <- cbind(countsHERO$X, countsHERO[, phenoHEROES$Tube_code])
row.names(countsHERO) <- countsHERO[,1]
countsHEROES <- countsHERO[,2:21]

# Place bulk RNA-seq in Expression object and snRNA-seq data in SingleCellExperiment object for MuSiC2 deconvolution
# Step 1: Darmanis et al., 2015 -- > The counts of all genes for any given cell where converted 
# to counts per million (CPM) by diving with the total number of reads and multiplying 
# by 10^6 followed by conversion to a log base 10.
sn_RNA_seq <- DarmanisBrainData(ensembl = FALSE,
                                location = TRUE,
                                remove.htseq = TRUE,
                                legacy = FALSE)

# Step 2: Make my own data in an ExpressionSet object
# Phenotype data for Demented for the brain regions Hippocampus
phenoHEROES_DEMENTED_HC <- phenoHEROES

# Make assaydata objects
countsHEROES_DEMENTED_HC <- countsHEROES
countsHEROES_DEMENTED_HC <- edgeR::cpm(countsHEROES_DEMENTED_HC, log = FALSE)
summary(countsHEROES_DEMENTED_HC)

# Check for same names in counts matrices and phenotype data frame
all(rownames(phenoHEROES_DEMENTED_HC) == colnames(countsHEROES_DEMENTED_HC))
colnames(countsHEROES_DEMENTED_HC)
rownames(phenoHEROES_DEMENTED_HC)

# Create the metadata
metadata <- data.frame(labelDescription = c("Case/control status", "Case number", "Chronological age", "Biological age",
                                            "Difference between Biological age and Chronological age", "Sex of the subject",
                                            "APOE alleles", "APOE4 status", 
                                            "BioBank ID", "Biobank", "Batch number", "Postmortem interval in hours",
                                            "Braak and Braak NFT stage", "Thal Pahse for Aβ plaques", "CERAD neuritic plaque score", "Spread of Lewy pathology (α-synuclein protein)", 
                                            "Brain region", "Labelled code", "Sequencing library",
                                            "RNA intergrity number",
                                            "RNA concentration", "Condition of the patient",
                                            "PC1", "PC2", "PC3"),
                       row.names = c("Status", "Case", "ChronAge", "RNAAge",
                                     "AgeDif", "Sex",
                                     "APOE", "APOE4", "BrainID", "BrainBank",
                                     "Batch", "PMD", "AD_Braak", "Thal_phase", "CERAD_score", "PD_Braak",
                                     "Brain_region", "Tube_code", "Library_code",
                                     "RIN", "Concentration", "Condition",
                                     "The first principal component",
                                     "The second principal component",
                                     "The third principal component"))

# Create AnnotatedDataFrame
phenoData_DEMENTED_HC <- new("AnnotatedDataFrame", data = phenoHEROES_DEMENTED_HC, varMetadata = metadata)

# Create annotation
annotation <- "NovaSeq 6000 Sequencing System (Illumina)"

experimentData <- new("MIAME",
                      name = "René A.J. Crans",
                      lab="Cellular & Systems Neurobiology Lab",
                      contact = "rene.crans@crg.eu",
                      title = "Hippocampal brain samples of control and demented individuals",
                      abstract = "ExpressionSet",
                      url = "https://www.crg.eu/en/programmes-groups/dierssen-lab")

# ExpressionSet
DEMENTED_HC_Set <- ExpressionSet(assayData = countsHEROES_DEMENTED_HC,
                           phenoData = phenoData_DEMENTED_HC, 
                           experimentData = experimentData,
                           annotation = "NovaSeq 6000 Sequencing System")

exprs(DEMENTED_HC_Set)

# Step 3: MuSiC2 deconvolution 
# Use T-statistics, because it is more robust when dealing with smaller sample sizes compared to TOAST.
unique(colData(sn_RNA_seq)$cell.type) # Cell-types present in the snRNA-seq experiment
counts(sn_RNA_seq) <- edgeR::cpm(counts(sn_RNA_seq), log = FALSE)

# Demented subjects -----------------------------------------------
bulk_control_DEMENTED_HC <- exprs(DEMENTED_HC_Set)[, DEMENTED_HC_Set$Status == "Control"]
bulk_case_DEMENTED_HC <- exprs(DEMENTED_HC_Set)[, DEMENTED_HC_Set$Status != "Control"]

set.seed(1234)
est_DEMENTED_HC <- music2_prop_t_statistics(bulk.control.mtx = bulk_control_DEMENTED_HC, 
                                            bulk.case.mtx = bulk_case_DEMENTED_HC, 
                                            sc.sce = sn_RNA_seq, clusters = 'cell.type', samples = 'experiment_sample_name', 
                                            select.ct = c("astrocytes", "endothelial", "microglia", "neurons", "oligodendrocytes"),
                                            n_resample = 1000, sample_prop = 0.5, cutoff_c = 0.05, cutoff_r = 0.01)

est.prop_HC <- est_DEMENTED_HC$Est.prop

# Plot estimated cell type proportions
prop_HC <- cbind("proportion" = c(est.prop_HC), 
                  "sampleID" = rep(rownames(est.prop_HC), 
                                   times=ncol(est.prop_HC)), 
                  "celltype" = rep(colnames(est.prop_HC), 
                                   each = nrow(est.prop_HC)))

prop_HC <- as.data.frame(prop_HC)
prop_HC$proportion <- as.numeric(as.character(prop_HC$proportion))

CONTROL_DEMENTED_HC <- phenoHEROES_DEMENTED_HC$Status

prop_HC$group <- rep(CONTROL_DEMENTED_HC, 5) # As I will select 5 cell-types
prop_HC$group <- factor(prop_HC$group, levels = c("Control", "PDD", "AD", "DSD"))
cols <-c("astrocytes" = "cadetblue2", 
         "endothelial" = "lightsalmon1", 
         "microglia" = "palegreen2", 
         "neurons" = "goldenrod1",
         "oligodendrocytes" = "steelblue3")

# Plot samples cell deconvolution
stat.test_DEMENTED_HC <- prop_HC %>%
  group_by(celltype) %>%
  tukey_hsd(proportion ~ group) 
stat.test_DEMENTED_HC

png(filename="NPJ_dementia/3_Figures/3_Deconvolution/MuSiC2_deconvolution.png",
    width     = 30,
    height    = 15,
    units     = "cm",
    res       = 1200,
    pointsize = 4)

ggplot(prop_HC, aes(x = group, y = proportion, colour = celltype, fill = celltype)) + xlab('')+
  geom_violin(width = 0.5, alpha = 0.75, size = 5) +
  geom_jitter(position = position_dodge2(width = 0.5, preserve = "total"), size = 4 , alpha = 0.5, color = "black") +
  stat_summary(fun = mean,
               geom = "crossbar", width = 0.5,size=0.5,color='gray36') +
  theme_bw() +
  facet_grid(.~celltype, scales = "free_y", 
             labeller = labeller(celltype = c(astrocytes = "Astrocytes",
                                              endothelial = "Endothelia",
                                              microglia = "Microglia",
                                              neurons = "Neurons",
                                              oligodendrocytes = "Oligodendrocytes"))) + 
  theme(plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.line = element_line(colour = "black"),
        strip.text.x = element_text(size = 16, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')+
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  ylab('Cell-type proportions') +
  ggtitle(paste("MuSiC2 cell-type deconvolution")) +
  stat_pvalue_manual(stat.test_DEMENTED_HC, "p = {p.adj}", hide.ns = TRUE, y.position = 0.15, size = 5)

dev.off()
