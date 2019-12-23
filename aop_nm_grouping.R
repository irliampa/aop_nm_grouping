
# Import libraries --------------------------------------------------------

library(GEOquery)
library(readxl)

library(limma)
library(AnnotationDbi)
library(mgug4122a.db)

library(enrichR)
library(factoextra)

library(bmd)
library(drc)

library(tidyverse)
library(stringr)
library(reshape2)
library(RColorBrewer)

# Download data -----------------------------------------------------------

gse_rahman16 <- getGEO("GSE81570")
save(gse_rahman16, file = "gse_rahman16.RData") 

#load("gse_rahman16.RData")
expr_rahman16 <- exprs(gse_rahman16[[1]])
pheno_rahman16 <- pData(gse_rahman16[[1]])
sampleNames(gse_rahman16[[1]])
experimentData(gse_rahman16[[1]])

# QC ----------------------------------------------------------------------
# take 40 random samples to check dataset status
pdf("qc.plots.pdf")
boxplot(expr_rahman16[, sample(1:ncol(expr_rahman16), size = 50)], pch = 16, cex = .9)
title(main = "Rahman et al, 2016 GSE81570 Dataset")
dev.off()

# Data are normalized!!!

# Data disctributions
range(expr_rahman16)
pdf("Rahman16.distribution.plots.pdf")
hist(expr_rahman16, breaks = 50, freq = FALSE, main = "Rahman et al, 2016 GSE81570 Dataset")
lines(density(expr_rahman16, na.rm=TRUE), col="red", lwd=2)
lines(density(expr_rahman16, adjust=2, na.rm=TRUE), lty="dotted", col="blue", lwd=2) 
dev.off()

# Number of samples in the complete datasets
n_samples_rahman16 <- ncol(expr_rahman16)
n_samples_rahman16 #233

## clean phenotype data
table_rahman16 <- apply(pheno_rahman16, 2, table)
pheno_rahman16 <- pheno_rahman16[, c(1, 2, 8, 13:15, 17, 24, 27, 36, 48:55)]
write_csv(pheno_rahman16, path = "pheno_rahman16.csv")

table_rahman16 <- apply(pheno_rahman16, 2, table)

gpl <- featureData(gse_rahman16[[1]])

## Keep only data for 1day after exposure
p_rah <- pheno_rahman16[ pheno_rahman16$characteristics_ch1.5 == "time point: 1d post-exposure", ]
e_rah <- expr_rahman16[, p_rah$geo_accession]

# Inspect
table_p_rah <- apply(p_rah, 2, table)
table_p_rah

treatment_colors <- as.factor(p_rah$`exposed to:ch1`)
treatment_colors_palette <- brewer.pal(nlevels(treatment_colors), "RdYlBu")
treatment_colors <- treatment_colors_palette[treatment_colors]
dose_colors <- as.factor(p_rah$`dosage:ch1`)
dose_colors_palette <- brewer.pal(nlevels(dose_colors), "RdYlBu")
dose_colors <- dose_colors_palette[dose_colors]

# plotMDS
# "pairwise" to choose the top genes separately for each pairwise comparison between the 
# samples or "common" to select the same genes for all comparisons
pdf("rahman_gene_expression_plots.pdf")
boxplot(e_rah) 
title(main = "Rahman et al, 2016 GSE81570 Boxplot")
plotMDS(e_rah, top= 1000, labels = p_rah$`exposed to:ch1`, col = treatment_colors)
title(main = "top 1000 pairwise genes, coloured by treatment")
plotMDS(e_rah, top= 1000, labels = p_rah$`dosage:ch1`, col = dose_colors)
title(main = "top 1000 pairwise genes, coloured by dose")
plotMDS(e_rah, top= 1000, labels = p_rah$`exposed to:ch1`, gene.selection = "common", 
        col = treatment_colors)
title(main = "top 1000 common genes, coloured by treatment")
plotMDS(e_rah, top= 1000, labels = p_rah$`dosage:ch1`, gene.selection = "common", 
        col = dose_colors)
title(main = "top 1000 common genes, coloured by dose")
dev.off()

# Reform dataset 
p_rah$Array <- as.factor(as.character(1:nrow(p_rah)))

p_rah$Dose <- p_rah$`dosage:ch1`

# remove the "control" sample, since missing info
p_rah <- p_rah[!(p_rah$`dosage:ch1` == "control"), ]

p_rah$Dose <- ifelse(p_rah$Dose == "0ug", "negative", (ifelse(p_rah$Dose == "54ug", "middle", 
                                                              ifelse(p_rah$Dose == "162ug" , "high", "very_high"))))
p_rah$Dose <- factor(p_rah$Dose, levels = c("negative", "middle", "high", "very_high"))
p_rah$Dose_n <- ifelse(p_rah$`dosage:ch1` == "0ug", 0, (ifelse(p_rah$`dosage:ch1` == "54ug", 54, 
                                                               ifelse(p_rah$`dosage:ch1` == "162ug" , 162, 486))))
p_rah$Dose_n <-  as.factor(p_rah$Dose_n)

p_rah$NM <- p_rah$`exposed to:ch1`
p_rah$NM <- ifelse(p_rah$NM == "anatase TiO2NPs of 8 nm diameter", "anatase_8", 
                   ifelse(p_rah$NM == "anatase TiO2NPs of 20nm diameter", "anatase_20", 
                          ifelse(p_rah$NM == "anatase TiO2NPs of 300 nm diameter", "anatase_300", 
                                 ifelse(p_rah$NM == "mix rutile/anatase TiO2NPs of 20 nm diameter" , "mix_20", 
                                        ifelse(p_rah$NM == "rutile hydrophilic TiO2NPs of 20 nm diameter" , "rutile_hydrophilic_20", 
                                               "rutile_hydrophobic_20")))))
p_rah$NM <- as.factor(p_rah$NM)
p_rah <- p_rah[, (ncol(p_rah)-3):ncol(p_rah)]

# include interaction
p_rah$Interaction <- paste(p_rah$NM, p_rah$Dose_n, sep = "_")
p_rah$Interaction <- as.factor(p_rah$Interaction)

e_rah <- e_rah[, rownames(p_rah)]

write.csv(p_rah, file = "p_rah.csv")
write.csv(e_rah, file = "e_rah.csv")

treatment_colors <- as.factor(p_rah$NM)
treatment_colors_palette <- brewer.pal(nlevels(treatment_colors), "RdYlBu")
treatment_colors <- treatment_colors_palette[treatment_colors]
dose_colors <- as.factor(p_rah$Dose)
dose_colors_palette <- brewer.pal(nlevels(dose_colors), "RdYlBu")
dose_colors <- dose_colors_palette[dose_colors]
pdf("rahman_gene_expression_plots_filtered.pdf")
boxplot(e_rah) 
title(main = "Rahman et al, 2016 GSE81570 Boxplot")
plotMDS(e_rah, top= 1000, labels = p_rah$NM, col = treatment_colors)
title(main = "top 1000 pairwise genes, coloured by treatment")
plotMDS(e_rah, top= 1000, labels = p_rah$Dose, col = dose_colors)
title(main = "top 1000 pairwise genes, coloured by dose")
plotMDS(e_rah, top= 1000, labels = p_rah$NM, gene.selection = "common", 
        col = treatment_colors)
title(main = "top 1000 common genes, coloured by treatment")
plotMDS(e_rah, top= 1000, labels = p_rah$Dose, gene.selection = "common", 
        col = dose_colors)
title(main = "top 1000 common genes, coloured by dose")
dev.off()

rm(expr_rahman16, table_p_rah, dose_colors, dose_colors_palette, treatment_colors, treatment_colors_palette)
save.image("rah16_reformed.RData")


# A. Bioinformatics Analysis ----------------------------------------------

#rm(list = ls())
#load("rah16_reformed.RData")

levels(p_rah$NM) <- c(levels(p_rah$NM), "control") 
p_rah$NM[p_rah$Dose == "negative"] <- "control"

p_rah$NM <- factor(p_rah$NM, levels = c("control", "anatase_8", "anatase_20", "anatase_300", "mix_20", "rutile_hydrophilic_20", "rutile_hydrophobic_20"))

model.rah.all <- model.matrix(~ 0 + Dose_n + NM, data = p_rah)
model.rah.all.intercept <- model.matrix(~ Dose + NM, data = p_rah)

fit.rah.limma <- lmFit(e_rah, model.rah.all[, 3:10])
fit.rah.limma.intercept <- lmFit(e_rah, model.rah.all.intercept[, c(1, 3:10)])

# filtering unexpressed probes when background probes are not available
# unexpressed probes 

hist(fit.rah.limma$Amean, breaks = 100)

# Plot residual standard deviation versus average log expression for a fitted microarray linear model.
# sigma is the estimated residual standard deviation
plotSA(fit.rah.limma)  
keep <- abs(fit.rah.limma$Amean) > 0.2
plotSA(fit.rah.limma[keep, ])

fit.rah.limma2 <- eBayes(fit.rah.limma)
fit.rah.limma3 <- eBayes(fit.rah.limma[keep,], trend=TRUE) 

length(fit.rah.limma$Amean)
length(keep[keep == T])

plotSA(fit.rah.limma.intercept)  
keep <- abs(fit.rah.limma.intercept$Amean) > 0.2
plotSA(fit.rah.limma.intercept[keep, ])

fit.rah.limma.intercept2 <- eBayes(fit.rah.limma.intercept)
fit.rah.limma.intercept3 <- eBayes(fit.rah.limma.intercept[keep,], trend=TRUE) 

length(fit.rah.limma.intercept$Amean)
length(keep[keep == T])

# full model
# each category versus controls
number_genes <- 99999
pval <- 0.05
logfc <- 0.5

# FULL, FILTERED DATASET
result.rah.limma3.NManatase_8 <- topTable(fit.rah.limma.intercept3, coef = 4, adjust="fdr", p.value = pval, number = number_genes)
result.rah.limma3.NManatase_20 <- topTable(fit.rah.limma.intercept3, coef = 5, adjust="fdr", p.value = pval, number = number_genes)
result.rah.limma3.NManatase_300 <- topTable(fit.rah.limma.intercept3, coef = 6, adjust="fdr", p.value = pval, number = number_genes)
result.rah.limma3.NMmix_20 <- topTable(fit.rah.limma.intercept3, coef = 7, adjust="fdr", p.value = pval, number = number_genes)
result.rah.limma3.NMrutile_hydrophilic_20 <- topTable(fit.rah.limma.intercept3, coef = 8, adjust="fdr", p.value = pval, number = number_genes)
result.rah.limma3.NMrutile_hydrophobic_20 <- topTable(fit.rah.limma.intercept3, coef = 9, adjust="fdr", p.value = pval, number = number_genes)

# filter log2FC
result.rah.limma3.lfcNManatase_8 <- topTable(fit.rah.limma.intercept3, coef = 4, adjust="fdr", p.value = pval, lfc = logfc, number = number_genes)
result.rah.limma3.lfcNManatase_20 <- topTable(fit.rah.limma.intercept3, coef = 5, adjust="fdr", p.value = pval, lfc = logfc, number = number_genes)
result.rah.limma3.lfcNManatase_300 <- topTable(fit.rah.limma.intercept3, coef = 6, adjust="fdr", p.value = pval, lfc = logfc, number = number_genes)
result.rah.limma3.lfcNMmix_20 <- topTable(fit.rah.limma.intercept3, coef = 7, adjust="fdr", p.value = pval, lfc = logfc, number = number_genes)
result.rah.limma3.lfcNMrutile_hydrophilic_20 <- topTable(fit.rah.limma.intercept3, coef = 8, adjust="fdr", p.value = pval, lfc = logfc, number = number_genes)
result.rah.limma3.lfcNMrutile_hydrophobic_20 <- topTable(fit.rah.limma.intercept3, coef = 9, adjust="fdr", p.value = pval, lfc = logfc, number = number_genes)

# decideTests
result.rah.limma.intercept2<- decideTests(fit.rah.limma.intercept3, method = 'global', adjust.method = "fdr", p.value = 0.05)
result.rah.limma.intercept3<- decideTests(fit.rah.limma.intercept3, method = 'global', adjust.method = "fdr", p.value = 0.05)

# Annotation of limma results ---------------------------------------------

# keep only limma3 model & log2FC filtered results

keytypes(mgug4122a.db)
keys <- keytypes(mgug4122a.db)[c(3, 6, 10, 20, 23)]

keys3.lfcNManatase_8 <- rownames(result.rah.limma3.lfcNManatase_8)
keys3.lfcNManatase_20 <- rownames(result.rah.limma3.lfcNManatase_20)
keys3.lfcNManatase_300 <- rownames(result.rah.limma3.lfcNManatase_300)
keys3.lfcNMmix_20 <- rownames(result.rah.limma3.lfcNMmix_20)
keys3.lfcNMrutile_hydrophilic_20 <- rownames(result.rah.limma3.lfcNMrutile_hydrophilic_20)
keys3.lfcNMrutile_hydrophobic_20 <- rownames(result.rah.limma3.lfcNMrutile_hydrophobic_20)

annotate3.lfcNManatase_8 <- AnnotationDbi::select(mgug4122a.db, keys=keys3.lfcNManatase_8, columns=keys, keytype="PROBEID")
annotate3.lfcNManatase_20 <- AnnotationDbi::select(mgug4122a.db, keys=keys3.lfcNManatase_20, columns=keys, keytype="PROBEID")
annotate3.lfcNManatase_300 <- AnnotationDbi::select(mgug4122a.db, keys=keys3.lfcNManatase_300, columns=keys, keytype="PROBEID")
annotate3.lfcNMmix_20 <- AnnotationDbi::select(mgug4122a.db, keys=keys3.lfcNMmix_20, columns=keys, keytype="PROBEID")
annotate3.lfcNMrutile_hydrophilic_20 <- AnnotationDbi::select(mgug4122a.db, keys=keys3.lfcNMrutile_hydrophilic_20, columns=keys, keytype="PROBEID")
annotate3.lfcNMrutile_hydrophobic_20 <- AnnotationDbi::select(mgug4122a.db, keys=keys3.lfcNMrutile_hydrophobic_20, columns=keys, keytype="PROBEID")

symbols3.lfcNManatase_8 <- unique(annotate3.lfcNManatase_8$SYMBOL)
symbols3.lfcNManatase_20 <- unique(annotate3.lfcNManatase_20$SYMBOL)
symbols3.lfcNManatase_300 <- unique(annotate3.lfcNManatase_300$SYMBOL)
symbols3.lfcNMmix_20 <- unique(annotate3.lfcNMmix_20$SYMBOL)
symbols3.lfcNMrutile_hydrophilic_20 <- unique(annotate3.lfcNMrutile_hydrophilic_20$SYMBOL)
symbols3.lfcNMrutile_hydrophobic_20 <- unique(annotate3.lfcNMrutile_hydrophobic_20$SYMBOL)

# Merge information
resultTable3.lfcNManatase_8 <- merge(result.rah.limma3.lfcNManatase_8, annotate3.lfcNManatase_8, by.x=0, by.y="PROBEID")
resultTable3.lfcNManatase_20 <- merge(result.rah.limma3.lfcNManatase_20, annotate3.lfcNManatase_20, by.x=0, by.y="PROBEID")
resultTable3.lfcNManatase_300 <- merge(result.rah.limma3.lfcNManatase_300, annotate3.lfcNManatase_300, by.x=0, by.y="PROBEID")
resultTable3.lfcNMmix_20 <- merge(result.rah.limma3.lfcNMmix_20, annotate3.lfcNMmix_20, by.x=0, by.y="PROBEID")
resultTable3.lfcNMrutile_hydrophilic_20 <- merge(result.rah.limma3.lfcNMrutile_hydrophilic_20, annotate3.lfcNMrutile_hydrophilic_20, by.x=0, by.y="PROBEID")
resultTable3.lfcNMrutile_hydrophobic_20 <- merge(result.rah.limma3.lfcNMrutile_hydrophobic_20, annotate3.lfcNMrutile_hydrophobic_20, by.x=0, by.y="PROBEID")

write.table(resultTable3.lfcNManatase_8, file = "resultTable3.lfcNManatase_8.txt", sep = "\t", dec = ".")
write.table(resultTable3.lfcNManatase_20, file = "resultTable3.lfcNManatase_20.txt", sep = "\t", dec = ".")
write.table(resultTable3.lfcNManatase_300, file = "resultTable3.lfcNManatase_300.txt", sep = "\t", dec = ".")
write.table(resultTable3.lfcNMmix_20 , file = "resultTable3.lfcNMmix_20 .txt", sep = "\t", dec = ".")
write.table(resultTable3.lfcNMrutile_hydrophilic_20, file = "resultTable3.lfcNMrutile_hydrophilic_20.txt", sep = "\t", dec = ".")
write.table(resultTable3.lfcNMrutile_hydrophobic_20, file = "resultTable3.lfcNMrutile_hydrophobic_20.txt", sep = "\t", dec = ".")

save.image("rah16_dea.RData")

# Center Intensities ------------------------------------------------------

# get log2 intensities (ratios) for the DEGs
probes3.lfc <- unique(c(rownames(result.rah.limma3.lfcNManatase_8), rownames(result.rah.limma3.lfcNManatase_20),
                        rownames(result.rah.limma3.lfcNManatase_300), rownames(result.rah.limma3.lfcNMmix_20), 
                        rownames(result.rah.limma3.lfcNMrutile_hydrophilic_20), rownames(result.rah.limma3.lfcNMrutile_hydrophobic_20)))
length(probes3.lfc)
sign.intensities3.lfc <- e_rah[probes3.lfc, ]

# From now on work only with the filtered data **3** & **lfc**

# group average using MEDIAN INTENSITY
control.samples <- rownames(p_rah[which(p_rah$NM == "control"), ])
NManatase_8.samples <- rownames(p_rah[which(p_rah$NM =="anatase_8"), ])
NManatase_20.samples <- rownames(p_rah[which(p_rah$NM =="anatase_20"), ])
NManatase_300.samples <- rownames(p_rah[which(p_rah$NM =="anatase_300"), ])
NMmix_20.samples <- rownames(p_rah[which(p_rah$NM =="mix_20"), ])
NMrutile_hydrophilic_20.samples <- rownames(p_rah[which(p_rah$NM =="rutile_hydrophilic_20"), ])
NMrutile_hydrophobic_20.samples <- rownames(p_rah[which(p_rah$NM =="rutile_hydrophobic_20"), ])

control.intensities <- sign.intensities3.lfc[, control.samples]
NManatase_8.intensities <- sign.intensities3.lfc[, NManatase_8.samples]
NManatase_20.intensities <- sign.intensities3.lfc[, NManatase_20.samples]
NManatase_300.intensities <- sign.intensities3.lfc[, NManatase_300.samples]
NMmix_20.intensities <- sign.intensities3.lfc[, NMmix_20.samples]
NMrutile_hydrophilic_20.intensities <- sign.intensities3.lfc[, NMrutile_hydrophilic_20.samples]
NMrutile_hydrophobic_20.intensities <- sign.intensities3.lfc[, NMrutile_hydrophobic_20.samples]

control.means <- apply(control.intensities, MARGIN = 1, mean)
NManatase_8.means <- apply(NManatase_8.intensities, MARGIN = 1, mean)
NManatase_20.means <- apply(NManatase_20.intensities, MARGIN = 1, mean)
NManatase_300.means <- apply(NManatase_300.intensities, MARGIN = 1, mean)
NMmix_20.means <- apply(NMmix_20.intensities, MARGIN = 1, mean)
NMrutile_hydrophilic_20.means <- apply(NMrutile_hydrophilic_20.intensities, MARGIN = 1, mean)
NMrutile_hydrophobic_20.means <- apply(NMrutile_hydrophobic_20.intensities, MARGIN = 1, mean)

# Center log2 transformed intensities using median of controls (median subtraction)

NManatase_8.intensities.centered <- (NManatase_8.intensities - control.means)
NManatase_20.intensities.centered <- (NManatase_20.intensities - control.means)
NManatase_300.intensities.centered <- (NManatase_300.intensities - control.means)
NMmix_20.intensities.centered <- (NMmix_20.intensities - control.means)
NMrutile_hydrophilic_20.intensities.centered <- (NMrutile_hydrophilic_20.intensities - control.means)
NMrutile_hydrophobic_20.intensities.centered <- (NMrutile_hydrophobic_20.intensities - control.means)

save.image("rah16_centered.RData")

# GSEA --------------------------------------------------------------------

dbs <- listEnrichrDbs()
dbs.go <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018")
dbs.pathways<- c("KEGG_2019_Mouse", "KEGG_2019_Human", "WikiPathways_2019_Mouse", "WikiPathways_2019_Human")
dbs.diseases <- c("OMIM_Disease", "DisGeNET")
dbs.phenotypes <- "MGI_Mammalian_Phenotype_Level_4_2019"
dbs.tissues <- c("Tissue_Protein_Expression_from_ProteomicsDB", "Tissue_Protein_Expression_from_Human_Proteome_Map", 
                 "Jensen_TISSUES", "ARCHS4_Tissues")

enriched.go.3.lfcNManatase_8 <- enrichr(symbols3.lfcNManatase_8, dbs.go)
enriched.go.3.lfcNManatase_20 <- enrichr(symbols3.lfcNManatase_20, dbs.go)
enriched.go.3.lfcNManatase_300 <- enrichr(symbols3.lfcNManatase_300, dbs.go)
enriched.go.3.lfcNMmix_20 <- enrichr(symbols3.lfcNMmix_20, dbs.go)
enriched.go.3.lfcNMrutile_hydrophilic_20 <- enrichr(symbols3.lfcNMrutile_hydrophilic_20, dbs.go)
enriched.go.3.lfcNMrutile_hydrophobic_20 <- enrichr(symbols3.lfcNMrutile_hydrophobic_20, dbs.go)

enriched.pathways.3.lfcNManatase_8 <- enrichr(symbols3.lfcNManatase_8, dbs.pathways)
enriched.pathways.3.lfcNManatase_20 <- enrichr(symbols3.lfcNManatase_20, dbs.pathways)
enriched.pathways.3.lfcNManatase_300 <- enrichr(symbols3.lfcNManatase_300, dbs.pathways)
enriched.pathways.3.lfcNMmix_20 <- enrichr(symbols3.lfcNMmix_20, dbs.pathways)
enriched.pathways.3.lfcNMrutile_hydrophilic_20 <- enrichr(symbols3.lfcNMrutile_hydrophilic_20, dbs.pathways)
enriched.pathways.3.lfcNMrutile_hydrophobic_20 <- enrichr(symbols3.lfcNMrutile_hydrophobic_20, dbs.pathways)

enriched.diseases.3.lfcNManatase_8 <- enrichr(symbols3.lfcNManatase_8, dbs.diseases)
enriched.diseases.3.lfcNManatase_20 <- enrichr(symbols3.lfcNManatase_20, dbs.diseases)
enriched.diseases.3.lfcNManatase_300 <- enrichr(symbols3.lfcNManatase_300, dbs.diseases)
enriched.diseases.3.lfcNMmix_20 <- enrichr(symbols3.lfcNMmix_20, dbs.diseases)
enriched.diseases.3.lfcNMrutile_hydrophilic_20 <- enrichr(symbols3.lfcNMrutile_hydrophilic_20, dbs.diseases)
enriched.diseases.3.lfcNMrutile_hydrophobic_20 <- enrichr(symbols3.lfcNMrutile_hydrophobic_20, dbs.diseases)

enriched.phenotypes.3.lfcNManatase_8 <- enrichr(symbols3.lfcNManatase_8, dbs.phenotypes)
enriched.phenotypes.3.lfcNManatase_20 <- enrichr(symbols3.lfcNManatase_20, dbs.phenotypes)
enriched.phenotypes.3.lfcNManatase_300 <- enrichr(symbols3.lfcNManatase_300, dbs.phenotypes)
enriched.phenotypes.3.lfcNMmix_20 <- enrichr(symbols3.lfcNMmix_20, dbs.phenotypes)
enriched.phenotypes.3.lfcNMrutile_hydrophilic_20 <- enrichr(symbols3.lfcNMrutile_hydrophilic_20, dbs.phenotypes)
enriched.phenotypes.3.lfcNMrutile_hydrophobic_20 <- enrichr(symbols3.lfcNMrutile_hydrophobic_20, dbs.phenotypes)

enriched.tissues.3.lfcNManatase_8 <- enrichr(symbols3.lfcNManatase_8, dbs.tissues)
enriched.tissues.3.lfcNManatase_20 <- enrichr(symbols3.lfcNManatase_20, dbs.tissues)
enriched.tissues.3.lfcNManatase_300 <- enrichr(symbols3.lfcNManatase_300, dbs.tissues)
enriched.tissues.3.lfcNMmix_20 <- enrichr(symbols3.lfcNMmix_20, dbs.tissues)
enriched.tissues.3.lfcNMrutile_hydrophilic_20 <- enrichr(symbols3.lfcNMrutile_hydrophilic_20, dbs.tissues)
enriched.tissues.3.lfcNMrutile_hydrophobic_20 <- enrichr(symbols3.lfcNMrutile_hydrophobic_20, dbs.tissues)

save.image("rah16_gsea.RData")

# filter out terms with less than 5 genes
# GO Biological Process
GO.BP.lfcNManatase_8 <- enriched.go.3.lfcNManatase_8$GO_Biological_Process_2018
GO.BP.lfcNManatase_8$Genes_count <- sapply(strsplit(as.character(GO.BP.lfcNManatase_8$Overlap),'/'), "[", 1)
GO.BP.lfcNManatase_8.5 <- GO.BP.lfcNManatase_8[GO.BP.lfcNManatase_8$Genes_count >4, ]

GO.BP.lfcNManatase_20  <- enriched.go.3.lfcNManatase_20$GO_Biological_Process_2018
GO.BP.lfcNManatase_20$Genes_count <- sapply(strsplit(as.character(GO.BP.lfcNManatase_20$Overlap),'/'), "[", 1)
GO.BP.lfcNManatase_20.5 <- GO.BP.lfcNManatase_20[GO.BP.lfcNManatase_20$Genes_count >4, ]

GO.BP.lfcNManatase_300 <- enriched.go.3.lfcNManatase_300$GO_Biological_Process_2018
GO.BP.lfcNManatase_300$Genes_count <- sapply(strsplit(as.character(GO.BP.lfcNManatase_300$Overlap),'/'), "[", 1)
GO.BP.lfcNManatase_300.5 <- GO.BP.lfcNManatase_300[GO.BP.lfcNManatase_300$Genes_count >4, ]

GO.BP.lfcNMmix_20  <- enriched.go.3.lfcNMmix_20$GO_Biological_Process_2018
GO.BP.lfcNMmix_20 $Genes_count <- sapply(strsplit(as.character(GO.BP.lfcNMmix_20$Overlap),'/'), "[", 1)
GO.BP.lfcNMmix_20.5 <- GO.BP.lfcNMmix_20[GO.BP.lfcNMmix_20$Genes_count >4, ]

GO.BP.lfcNMrutile_hydrophilic_20   <- enriched.go.3.lfcNMrutile_hydrophilic_20$GO_Biological_Process_2018
GO.BP.lfcNMrutile_hydrophilic_20$Genes_count <- sapply(strsplit(as.character(GO.BP.lfcNMrutile_hydrophilic_20$Overlap),'/'), "[", 1)
GO.BP.lfcNMrutile_hydrophilic_20.5 <- GO.BP.lfcNMrutile_hydrophilic_20 [GO.BP.lfcNMrutile_hydrophilic_20$Genes_count >4, ]

GO.BP.lfcNMrutile_hydrophobic_20  <- enriched.go.3.lfcNMrutile_hydrophobic_20$GO_Biological_Process_2018
GO.BP.lfcNMrutile_hydrophobic_20$Genes_count <- sapply(strsplit(as.character(GO.BP.lfcNMrutile_hydrophobic_20$Overlap),'/'), "[", 1)
GO.BP.lfcNMrutile_hydrophobic_20.5 <- GO.BP.lfcNMrutile_hydrophobic_20[GO.BP.lfcNMrutile_hydrophobic_20$Genes_count >4, ]

# KEGG Mouse 
KEGG.Mm.lfcNManatase_8 <- enriched.go.3.lfcNManatase_8$GO_Biological_Process_2018
KEGG.Mm.lfcNManatase_8$Genes_count <- sapply(strsplit(as.character(KEGG.Mm.lfcNManatase_8$Overlap),'/'), "[", 1)
KEGG.Mm.lfcNManatase_8.5 <- KEGG.Mm.lfcNManatase_8[KEGG.Mm.lfcNManatase_8$Genes_count >4, ]

KEGG.Mm.lfcNManatase_20  <- enriched.go.3.lfcNManatase_20$GO_Biological_Process_2018
KEGG.Mm.lfcNManatase_20$Genes_count <- sapply(strsplit(as.character(KEGG.Mm.lfcNManatase_20$Overlap),'/'), "[", 1)
KEGG.Mm.lfcNManatase_20.5 <- KEGG.Mm.lfcNManatase_20[KEGG.Mm.lfcNManatase_20$Genes_count >4, ]

KEGG.Mm.lfcNManatase_300 <- enriched.go.3.lfcNManatase_300$GO_Biological_Process_2018
KEGG.Mm.lfcNManatase_300$Genes_count <- sapply(strsplit(as.character(KEGG.Mm.lfcNManatase_300$Overlap),'/'), "[", 1)
KEGG.Mm.lfcNManatase_300.5 <- KEGG.Mm.lfcNManatase_300[KEGG.Mm.lfcNManatase_300$Genes_count >4, ]

KEGG.Mm.lfcNMmix_20  <- enriched.go.3.lfcNMmix_20$GO_Biological_Process_2018
KEGG.Mm.lfcNMmix_20 $Genes_count <- sapply(strsplit(as.character(KEGG.Mm.lfcNMmix_20$Overlap),'/'), "[", 1)
KEGG.Mm.lfcNMmix_20.5 <- KEGG.Mm.lfcNMmix_20[KEGG.Mm.lfcNMmix_20$Genes_count >4, ]

KEGG.Mm.lfcNMrutile_hydrophilic_20   <- enriched.go.3.lfcNMrutile_hydrophilic_20$GO_Biological_Process_2018
KEGG.Mm.lfcNMrutile_hydrophilic_20$Genes_count <- sapply(strsplit(as.character(KEGG.Mm.lfcNMrutile_hydrophilic_20$Overlap),'/'), "[", 1)
KEGG.Mm.lfcNMrutile_hydrophilic_20.5 <- KEGG.Mm.lfcNMrutile_hydrophilic_20 [KEGG.Mm.lfcNMrutile_hydrophilic_20$Genes_count >4, ]

KEGG.Mm.lfcNMrutile_hydrophobic_20  <- enriched.go.3.lfcNMrutile_hydrophobic_20$GO_Biological_Process_2018
KEGG.Mm.lfcNMrutile_hydrophobic_20$Genes_count <- sapply(strsplit(as.character(KEGG.Mm.lfcNMrutile_hydrophobic_20$Overlap),'/'), "[", 1)
KEGG.Mm.lfcNMrutile_hydrophobic_20.5 <- KEGG.Mm.lfcNMrutile_hydrophobic_20[KEGG.Mm.lfcNMrutile_hydrophobic_20$Genes_count >4, ]

# WikiPathways Mouse 
WikiPathways.Mm.lfcNManatase_8 <- enriched.go.3.lfcNManatase_8$GO_Biological_Process_2018
WikiPathways.Mm.lfcNManatase_8$Genes_count <- sapply(strsplit(as.character(WikiPathways.Mm.lfcNManatase_8$Overlap),'/'), "[", 1)
WikiPathways.Mm.lfcNManatase_8.5 <- WikiPathways.Mm.lfcNManatase_8[WikiPathways.Mm.lfcNManatase_8$Genes_count >4, ]

WikiPathways.Mm.lfcNManatase_20  <- enriched.go.3.lfcNManatase_20$GO_Biological_Process_2018
WikiPathways.Mm.lfcNManatase_20$Genes_count <- sapply(strsplit(as.character(WikiPathways.Mm.lfcNManatase_20$Overlap),'/'), "[", 1)
WikiPathways.Mm.lfcNManatase_20.5 <- WikiPathways.Mm.lfcNManatase_20[WikiPathways.Mm.lfcNManatase_20$Genes_count >4, ]

WikiPathways.Mm.lfcNManatase_300 <- enriched.go.3.lfcNManatase_300$GO_Biological_Process_2018
WikiPathways.Mm.lfcNManatase_300$Genes_count <- sapply(strsplit(as.character(WikiPathways.Mm.lfcNManatase_300$Overlap),'/'), "[", 1)
WikiPathways.Mm.lfcNManatase_300.5 <- WikiPathways.Mm.lfcNManatase_300[WikiPathways.Mm.lfcNManatase_300$Genes_count >4, ]

WikiPathways.Mm.lfcNMmix_20  <- enriched.go.3.lfcNMmix_20$GO_Biological_Process_2018
WikiPathways.Mm.lfcNMmix_20 $Genes_count <- sapply(strsplit(as.character(WikiPathways.Mm.lfcNMmix_20$Overlap),'/'), "[", 1)
WikiPathways.Mm.lfcNMmix_20.5 <- WikiPathways.Mm.lfcNMmix_20[WikiPathways.Mm.lfcNMmix_20$Genes_count >4, ]

WikiPathways.Mm.lfcNMrutile_hydrophilic_20   <- enriched.go.3.lfcNMrutile_hydrophilic_20$GO_Biological_Process_2018
WikiPathways.Mm.lfcNMrutile_hydrophilic_20$Genes_count <- sapply(strsplit(as.character(WikiPathways.Mm.lfcNMrutile_hydrophilic_20$Overlap),'/'), "[", 1)
WikiPathways.Mm.lfcNMrutile_hydrophilic_20.5 <- WikiPathways.Mm.lfcNMrutile_hydrophilic_20 [WikiPathways.Mm.lfcNMrutile_hydrophilic_20$Genes_count >4, ]

WikiPathways.Mm.lfcNMrutile_hydrophobic_20  <- enriched.go.3.lfcNMrutile_hydrophobic_20$GO_Biological_Process_2018
WikiPathways.Mm.lfcNMrutile_hydrophobic_20$Genes_count <- sapply(strsplit(as.character(WikiPathways.Mm.lfcNMrutile_hydrophobic_20$Overlap),'/'), "[", 1)
WikiPathways.Mm.lfcNMrutile_hydrophobic_20.5 <- WikiPathways.Mm.lfcNMrutile_hydrophobic_20[WikiPathways.Mm.lfcNMrutile_hydrophobic_20$Genes_count >4, ]

# KEGG Human
KEGG.Hs.lfcNManatase_8 <- enriched.go.3.lfcNManatase_8$GO_Biological_Process_2018
KEGG.Hs.lfcNManatase_8$Genes_count <- sapply(strsplit(as.character(KEGG.Hs.lfcNManatase_8$Overlap),'/'), "[", 1)
KEGG.Hs.lfcNManatase_8.5 <- KEGG.Hs.lfcNManatase_8[KEGG.Hs.lfcNManatase_8$Genes_count >4, ]

KEGG.Hs.lfcNManatase_20  <- enriched.go.3.lfcNManatase_20$GO_Biological_Process_2018
KEGG.Hs.lfcNManatase_20$Genes_count <- sapply(strsplit(as.character(KEGG.Hs.lfcNManatase_20$Overlap),'/'), "[", 1)
KEGG.Hs.lfcNManatase_20.5 <- KEGG.Hs.lfcNManatase_20[KEGG.Hs.lfcNManatase_20$Genes_count >4, ]

KEGG.Hs.lfcNManatase_300 <- enriched.go.3.lfcNManatase_300$GO_Biological_Process_2018
KEGG.Hs.lfcNManatase_300$Genes_count <- sapply(strsplit(as.character(KEGG.Hs.lfcNManatase_300$Overlap),'/'), "[", 1)
KEGG.Hs.lfcNManatase_300.5 <- KEGG.Hs.lfcNManatase_300[KEGG.Hs.lfcNManatase_300$Genes_count >4, ]

KEGG.Hs.lfcNMmix_20  <- enriched.go.3.lfcNMmix_20$GO_Biological_Process_2018
KEGG.Hs.lfcNMmix_20$Genes_count <- sapply(strsplit(as.character(KEGG.Hs.lfcNMmix_20$Overlap),'/'), "[", 1)
KEGG.Hs.lfcNMmix_20.5 <- KEGG.Hs.lfcNMmix_20[KEGG.Hs.lfcNMmix_20$Genes_count >4, ]

KEGG.Hs.lfcNMrutile_hydrophilic_20   <- enriched.go.3.lfcNMrutile_hydrophilic_20$GO_Biological_Process_2018
KEGG.Hs.lfcNMrutile_hydrophilic_20$Genes_count <- sapply(strsplit(as.character(KEGG.Hs.lfcNMrutile_hydrophilic_20$Overlap),'/'), "[", 1)
KEGG.Hs.lfcNMrutile_hydrophilic_20.5 <- KEGG.Hs.lfcNMrutile_hydrophilic_20 [KEGG.Hs.lfcNMrutile_hydrophilic_20$Genes_count >4, ]

KEGG.Hs.lfcNMrutile_hydrophobic_20  <- enriched.go.3.lfcNMrutile_hydrophobic_20$GO_Biological_Process_2018
KEGG.Hs.lfcNMrutile_hydrophobic_20$Genes_count <- sapply(strsplit(as.character(KEGG.Hs.lfcNMrutile_hydrophobic_20$Overlap),'/'), "[", 1)
KEGG.Hs.lfcNMrutile_hydrophobic_20.5 <- KEGG.Hs.lfcNMrutile_hydrophobic_20[KEGG.Hs.lfcNMrutile_hydrophobic_20$Genes_count >4, ]

# WikiPathways Human
WikiPathways.Hs.lfcNManatase_8 <- enriched.go.3.lfcNManatase_8$GO_Biological_Process_2018
WikiPathways.Hs.lfcNManatase_8$Genes_count <- sapply(strsplit(as.character(WikiPathways.Hs.lfcNManatase_8$Overlap),'/'), "[", 1)
WikiPathways.Hs.lfcNManatase_8.5 <- WikiPathways.Hs.lfcNManatase_8[WikiPathways.Hs.lfcNManatase_8$Genes_count >4, ]

WikiPathways.Hs.lfcNManatase_20  <- enriched.go.3.lfcNManatase_20$GO_Biological_Process_2018
WikiPathways.Hs.lfcNManatase_20$Genes_count <- sapply(strsplit(as.character(WikiPathways.Hs.lfcNManatase_20$Overlap),'/'), "[", 1)
WikiPathways.Hs.lfcNManatase_20.5 <- WikiPathways.Hs.lfcNManatase_20[WikiPathways.Hs.lfcNManatase_20$Genes_count >4, ]

WikiPathways.Hs.lfcNManatase_300 <- enriched.go.3.lfcNManatase_300$GO_Biological_Process_2018
WikiPathways.Hs.lfcNManatase_300$Genes_count <- sapply(strsplit(as.character(WikiPathways.Hs.lfcNManatase_300$Overlap),'/'), "[", 1)
WikiPathways.Hs.lfcNManatase_300.5 <- WikiPathways.Hs.lfcNManatase_300[WikiPathways.Hs.lfcNManatase_300$Genes_count >4, ]

WikiPathways.Hs.lfcNMmix_20  <- enriched.go.3.lfcNMmix_20$GO_Biological_Process_2018
WikiPathways.Hs.lfcNMmix_20 $Genes_count <- sapply(strsplit(as.character(WikiPathways.Hs.lfcNMmix_20$Overlap),'/'), "[", 1)
WikiPathways.Hs.lfcNMmix_20.5 <- WikiPathways.Hs.lfcNMmix_20[WikiPathways.Hs.lfcNMmix_20$Genes_count >4, ]

WikiPathways.Hs.lfcNMrutile_hydrophilic_20   <- enriched.go.3.lfcNMrutile_hydrophilic_20$GO_Biological_Process_2018
WikiPathways.Hs.lfcNMrutile_hydrophilic_20$Genes_count <- sapply(strsplit(as.character(WikiPathways.Hs.lfcNMrutile_hydrophilic_20$Overlap),'/'), "[", 1)
WikiPathways.Hs.lfcNMrutile_hydrophilic_20.5 <- WikiPathways.Hs.lfcNMrutile_hydrophilic_20 [WikiPathways.Hs.lfcNMrutile_hydrophilic_20$Genes_count >4, ]

WikiPathways.Hs.lfcNMrutile_hydrophobic_20  <- enriched.go.3.lfcNMrutile_hydrophobic_20$GO_Biological_Process_2018
WikiPathways.Hs.lfcNMrutile_hydrophobic_20$Genes_count <- sapply(strsplit(as.character(WikiPathways.Hs.lfcNMrutile_hydrophobic_20$Overlap),'/'), "[", 1)
WikiPathways.Hs.lfcNMrutile_hydrophobic_20.5 <- WikiPathways.Hs.lfcNMrutile_hydrophobic_20[WikiPathways.Hs.lfcNMrutile_hydrophobic_20$Genes_count >4, ]

save.image("rah16_gsea.RData")

# BMD analysis ------------------------------------------------------------

# calculate BMD for all significantly affected probes

df.intensities <- data.frame(NManatase_8.intensities.centered, NManatase_20.intensities.centered, 
                             NManatase_300.intensities.centered, NMmix_20.intensities.centered,
                             NMrutile_hydrophilic_20.intensities.centered, NMrutile_hydrophobic_20.intensities.centered)
  
dose.treatment <- p_rah[colnames(df.intensities), c(3,4)]
df.dose.response <- data.frame(dose.treatment, t(df.intensities))

rm(df.intensities, dose.treatment)

# stratify data
DR.control <- df.dose.response[which(df.dose.response$NM == "control"), -2]
DR.anatase_8 <- df.dose.response[which(df.dose.response$NM == "anatase_8"), -2]
DR.anatase_20 <- df.dose.response[which(df.dose.response$NM == "anatase_20"), -2]
DR.anatase_300 <- df.dose.response[which(df.dose.response$NM == "anatase_300"), -2]
DR.mix_20 <- df.dose.response[which(df.dose.response$NM == "mix_20"), -2]
DR.rutile_hydrophilic_20  <- df.dose.response[which(df.dose.response$NM == "rutile_hydrophilic_20"), -2]
DR.rutile_hydrophobic_20 <- df.dose.response[which(df.dose.response$NM == "rutile_hydrophobic_20"), -2]

getMeanFunctions()
fct.list = list(LL.2(), LL.3(), LL.3u(), LL.4(), LL.5(), LL2.2(), LL2.3(), LL2.3u(), LL2.4(), LL2.5())

DR.anatase_8$Dose_n <- as.character(DR.anatase_8$Dose_n)
DR.anatase_8$Dose_n <- as.numeric(DR.anatase_8$Dose_n)

DR.anatase_20$Dose_n <- as.character(DR.anatase_20$Dose_n)
DR.anatase_20$Dose_n <- as.numeric(DR.anatase_20$Dose_n)

DR.anatase_300$Dose_n <- as.character(DR.anatase_300$Dose_n)
DR.anatase_300$Dose_n <- as.numeric(DR.anatase_300$Dose_n)

DR.mix_20$Dose_n <- as.character(DR.mix_20$Dose_n)
DR.mix_20$Dose_n <- as.numeric(DR.mix_20$Dose_n)

DR.rutile_hydrophilic_20$Dose_n <- as.character(DR.rutile_hydrophilic_20$Dose_n)
DR.rutile_hydrophilic_20$Dose_n <- as.numeric(DR.rutile_hydrophilic_20$Dose_n)

DR.rutile_hydrophobic_20$Dose_n <- as.character(DR.rutile_hydrophobic_20$Dose_n)
DR.rutile_hydrophobic_20$Dose_n <- as.numeric(DR.rutile_hydrophobic_20$Dose_n)

save.image("rah16_bmd.RData")

anatase_8.data <- list()
for (i in 2:ncol(DR.anatase_8)){
  anatase_8.data[[i]] <- DR.anatase_8[, c(1, i)] }

# remove first element
anatase_8.data[sapply(anatase_8.data, is.null)] <- NULL
# rename elements
element_names <- list()
for (i in 1:length(anatase_8.data)){
  element_names[i] <- colnames(anatase_8.data[[i]])[2]
  colnames(anatase_8.data[[i]]) <- c("Dose", "Response_expr")
  }
names(anatase_8.data) <- element_names

# Model Selection using AIC
# For Akaike's information criterion and the residual standard error: 
# the smaller the better and for lack-of-fit test (against a one-way ANOVA model): 
# the larger (the p-value) the better

BMDt.anatase_8.ll2 <- list()
anatase_8.model <- list()
for (i in 1:length(anatase_8.data)){
tryCatch(
    expr = {
      BMDt.anatase_8.ll2[[i]] <- drm(Response_expr ~ Dose, 
                               data = anatase_8.data[[i]],  type = "continuous" , fct = LL.2()) 
      anatase_8.model[[i]] <- mselect(BMDt.anatase_8.ll2[[i]], fctList = list(LL.3(), LL.3u(), LL.4(), LL.5(), 
                                                                        LL2.2(), LL2.3(), LL2.3u(), LL2.4(), LL2.5()))

    },
    error = function(e){ 
      BMDt.anatase_8.ll2[[i]] <- drm(Response_expr ~ Dose, 
                                     data = anatase_8.data[[i]],  type = "continuous" , fct = LL2.2()) 
      anatase_8.model[[i]] <- mselect(BMDt.anatase_8.ll2[[i]], fctList = list(LL.2(), LL.3(), LL.3u(), LL.4(), LL.5(), 
                                                                              LL2.3(), LL2.3u(), LL2.4(), LL2.5()))
      
    }
)
}    
    
anatase_8.model.best <- list()
for (i in 1:length(anatase_8.model)){
  anatase_8.model.best[i] <- rownames(anatase_8.model[[i]])[1]
}
for (i in 1:length(anatase_8.model.best)){
  anatase_8.model.best[[i]] <- ifelse(is.null(anatase_8.model.best[[i]]), paste("Unavailable"),
        anatase_8.model.best[[i]])
}

anatase_8.model.best <- unlist(anatase_8.model.best)
anatase_8.model.best <- lapply(anatase_8.model.best, as.factor)

anatase_8.drc <- list()
anatase_8.bmd <- list()
anatase_8.bmdl <- list()

  
for (i in 1:length(anatase_8.model.best)){  
  
  tryCatch(
    expr = {
      anatase_8.drc[[i]] <- drm(Response_expr ~ Dose, data =  anatase_8.data[[i]],  type = "continuous" , 
                                fct = LL2.3u())
      anatase_8.bmd[[i]] <- bmd(anatase_8.drc[[i]], bmr = 0.1, backg = 0.5831192, def = "hybrid")[1]
      anatase_8.bmdl[[i]] <- bmd(anatase_8.drc[[i]], bmr = 0.1, backg = 0.5831192, def = "hybrid")[2]
    },
    error = function(e){ 
      anatase_8.drc[[i]] <- drm(Response_expr ~ Dose, data =  anatase_8.data[[i]],  type = "continuous" , 
                                fct = LL2.3())
      anatase_8.bmd[[i]] <- bmd(anatase_8.drc[[i]], bmr = 0.1, backg = 0.5831192, def = "hybrid")[1]
      anatase_8.bmdl[[i]] <- bmd(anatase_8.drc[[i]], bmr = 0.1, backg = 0.5831192, def = "hybrid")[2]    
    }
  )
}    

names(anatase_8.drc) <- element_names
names(anatase_8.bmd) <- element_names
names(anatase_8.bmdl) <- element_names

##
anatase_20.data <- list()
for (i in 2:ncol(DR.anatase_20)){
  anatase_20.data[[i]] <- DR.anatase_20[, c(1, i)] }

# remove first element
anatase_20.data[sapply(anatase_20.data, is.null)] <- NULL
# rename elements
element_names <- list()
for (i in 1:length(anatase_20.data)){
  element_names[i] <- colnames(anatase_20.data[[i]])[2]
  colnames(anatase_20.data[[i]]) <- c("Dose", "Response_expr")
}
names(anatase_20.data) <- element_names

BMDt.anatase_20.ll2 <- list()
anatase_20.model <- list()
for (i in 1:length(anatase_20.data)){
  tryCatch(
    expr = {
      BMDt.anatase_20.ll2[[i]] <- drm(Response_expr ~ Dose, 
                                     data = anatase_20.data[[i]],  type = "continuous" , fct = LL.2()) 
      anatase_20.model[[i]] <- mselect(BMDt.anatase_20.ll2[[i]], fctList = list(LL.3(), LL.3u(), LL.4(), LL.5(), 
                                                                              LL2.2(), LL2.3(), LL2.3u(), LL2.4(), LL2.5()))
      
    },
    error = function(e){ 
      BMDt.anatase_20.ll2[[i]] <- drm(Response_expr ~ Dose, 
                                     data = anatase_20.data[[i]],  type = "continuous" , fct =  LL2.5()) 
      anatase_20.model[[i]] <- mselect(BMDt.anatase_20.ll2[[i]], fctList = list(LL.2(), LL.3(), LL.3u(), LL.4(), LL.5(), 
                                                                           LL2.2(), LL2.3(), LL2.3u(), LL2.4(), LL2.5()))
      
    }
  )
}    

anatase_20.model.best <- list()
for (i in 1:length(anatase_20.model)){
  anatase_20.model.best[i] <- rownames(anatase_20.model[[i]])[1]
}
for (i in 1:length(anatase_20.model.best)){
  anatase_20.model.best[[i]] <- ifelse(is.null(anatase_20.model.best[[i]]), paste("Unavailable"),
                                      anatase_20.model.best[[i]])
}

anatase_20.model.best <- unlist(anatase_20.model.best)
anatase_20.model.best <- lapply(anatase_20.model.best, as.factor)

anatase_20.drc <- list()
anatase_20.bmd <- list()
anatase_20.bmdl <- list()

for (i in 1:length(anatase_20.model.best)){  
  
  tryCatch(
    expr = {
      anatase_20.drc[[i]] <- drm(Response_expr ~ Dose, data =  anatase_20.data[[i]],  type = "continuous" , 
                                fct = LL.5())
      anatase_20.bmd[[i]] <- bmd(anatase_20.drc[[i]], bmr = 0.1, backg = 0.5831192, def = "hybrid")[1]
      anatase_20.bmdl[[i]] <- bmd(anatase_20.drc[[i]], bmr = 0.1, backg = 0.5831192, def = "hybrid")[2]
    },
    error = function(e){ 
      anatase_20.drc[[i]] <- drm(Response_expr ~ Dose, data =  anatase_20.data[[i]],  type = "continuous" , 
                                fct = LL2.3())
      anatase_20.bmd[[i]] <- bmd(anatase_20.drc[[i]], bmr = 0.1, backg = 0.5831192, def = "hybrid")[1]
      anatase_20.bmdl[[i]] <- bmd(anatase_20.drc[[i]], bmr = 0.1, backg = 0.5831192, def = "hybrid")[2]    
    }
  )
}    

names(anatase_20.drc) <- element_names
names(anatase_20.bmd) <- element_names
names(anatase_20.bmdl) <- element_names

##
anatase_300.data <- list()
for (i in 2:ncol(DR.anatase_300)){
  anatase_300.data[[i]] <- DR.anatase_300[, c(1, i)] }

# remove first element
anatase_300.data[sapply(anatase_300.data, is.null)] <- NULL
# rename elements
element_names <- list()
for (i in 1:length(anatase_300.data)){
  element_names[i] <- colnames(anatase_300.data[[i]])[2]
  colnames(anatase_300.data[[i]]) <- c("Dose", "Response_expr")
}
names(anatase_300.data) <- element_names

BMDt.anatase_300.ll2 <- list()
anatase_300.model <- list()
for (i in 1:length(anatase_300.data)){
  tryCatch(
    expr = {
      BMDt.anatase_300.ll2[[i]] <- drm(Response_expr ~ Dose, 
                                      data = anatase_300.data[[i]],  type = "continuous" , fct = LL2.3u()) 
      anatase_300.model[[i]] <- mselect(BMDt.anatase_300.ll2[[i]], fctList = list(LL.2(), LL.3(), LL.3u(), LL.4(), LL.5(), 
                                                                                LL2.2(), LL2.3(), LL2.3u(), LL2.4(), LL2.5()))
      
    },
    error = function(e){ 
      BMDt.anatase_300.ll2[[i]] <- drm(Response_expr ~ Dose, 
                                      data = anatase_300.data[[i]],  type = "continuous" , fct =  LL.2()) 
      anatase_300.model[[i]] <- mselect(BMDt.anatase_300.ll2[[i]], fctList = list(LL.2(), LL.3(), LL.3u(), LL.4(), LL.5(), 
                                                                                LL2.2(), LL2.3(), LL2.3u(), LL2.4(), LL2.5()))
      
    }
  )
}    

anatase_300.model.best <- list()
for (i in 1:length(anatase_300.model)){
  anatase_300.model.best[i] <- rownames(anatase_300.model[[i]])[1]
}
for (i in 1:length(anatase_300.model.best)){
  anatase_300.model.best[[i]] <- ifelse(is.null(anatase_300.model.best[[i]]), paste("Unavailable"),
                                       anatase_300.model.best[[i]])
}

anatase_300.model.best <- unlist(anatase_300.model.best)
anatase_300.model.best <- lapply(anatase_300.model.best, as.factor)

anatase_300.drc <- list()
anatase_300.bmd <- list()
anatase_300.bmdl <- list()

for (i in 1:length(anatase_300.model.best)){  
  
  tryCatch(
    expr = {
      anatase_300.drc[[i]] <- drm(Response_expr ~ Dose, data =  anatase_300.data[[i]],  type = "continuous" , 
                                 fct = LL.5())
      anatase_300.bmd[[i]] <- bmd(anatase_300.drc[[i]], bmr = 0.1, backg = 0.5831192, def = "hybrid")[1]
      anatase_300.bmdl[[i]] <- bmd(anatase_300.drc[[i]], bmr = 0.1, backg = 0.5831192, def = "hybrid")[2]
      names(anatase_300.drc)[[i]] <- element_names[[i]]
      names(anatase_300.bmd)[[i]] <- element_names[[i]]
      names(anatase_300.bmdl)[[i]] <- element_names[[i]]    
      },
    error = function(e){ 
      anatase_300.drc[[i]] <- drm(Response_expr ~ Dose, data =  anatase_300.data[[i]],  type = "continuous" , 
                                 fct = LL2.3())
      anatase_300.bmd[[i]] <- bmd(anatase_300.drc[[i]], bmr = 0.1, backg = 0.5831192, def = "hybrid")[1]
      anatase_300.bmdl[[i]] <- bmd(anatase_300.drc[[i]], bmr = 0.1, backg = 0.5831192, def = "hybrid")[2]    
      names(anatase_300.drc)[[i]] <- element_names[[i]]
      names(anatase_300.bmd)[[i]] <- element_names[[i]]
      names(anatase_300.bmdl)[[i]] <- element_names[[i]]
      }
  )
}    


##
mix_20.data <- list()
for (i in 2:ncol(DR.mix_20)){
  mix_20.data[[i]] <- DR.mix_20[, c(1, i)] }

# remove first element
mix_20.data[sapply(mix_20.data, is.null)] <- NULL
# rename elements
element_names <- list()
for (i in 1:length(mix_20.data)){
  element_names[i] <- colnames(mix_20.data[[i]])[2]
  colnames(mix_20.data[[i]]) <- c("Dose", "Response_expr")
}
names(mix_20.data) <- element_names

BMDt.mix_20.ll2 <- list()
mix_20.model <- list()
for (i in 1:length(mix_20.data)){
  tryCatch(
    expr = {
      BMDt.mix_20.ll2[[i]] <- drm(Response_expr ~ Dose, 
                                      data = mix_20.data[[i]],  type = "continuous" , fct = LL2.3()) 
      mix_20.model[[i]] <- mselect(BMDt.mix_20.ll2[[i]], fctList = list(LL.2(), LL.3(), LL.3u(), LL.4(), LL.5(), 
                                                                                LL2.2(), LL2.3(), LL2.3u(), LL2.4(), LL2.5()))
      
    },
    error = function(e){ 
      BMDt.mix_20.ll2[[i]] <- drm(Response_expr ~ Dose, 
                                      data = mix_20.data[[i]],  type = "continuous" , fct =  LL2.2()) 
      mix_20.model[[i]] <- mselect(BMDt.mix_20.ll2[[i]], fctList = list(LL.2(), LL.3(), LL.3u(), LL.4(), LL.5(), 
                                                                                LL2.2(), LL2.3(), LL2.3u(), LL2.4(), LL2.5()))
      
    }
  )
}    

mix_20.model.best <- list()
for (i in 1:length(mix_20.model)){
  mix_20.model.best[i] <- rownames(mix_20.model[[i]])[1]
}
for (i in 1:length(mix_20.model.best)){
  mix_20.model.best[[i]] <- ifelse(is.null(mix_20.model.best[[i]]), paste("Unavailable"),
                                       mix_20.model.best[[i]])
}

mix_20.model.best <- unlist(mix_20.model.best)
mix_20.model.best <- lapply(mix_20.model.best, as.factor)

mix_20.drc <- list()
mix_20.bmd <- list()
mix_20.bmdl <- list()

for (i in 1:length(mix_20.model.best)){  
  
  tryCatch(
    expr = {
      mix_20.drc[[i]] <- drm(Response_expr ~ Dose, data =  mix_20.data[[i]],  type = "continuous" , 
                                 fct = LL.5())
      mix_20.bmd[[i]] <- bmd(mix_20.drc[[i]], bmr = 0.1, backg = 0.5831192, def = "hybrid")[1]
      mix_20.bmdl[[i]] <- bmd(mix_20.drc[[i]], bmr = 0.1, backg = 0.5831192, def = "hybrid")[2]
      names(mix_20.drc)[[i]] <- element_names[[i]]
      names(mix_20.bmd)[[i]] <- element_names[[i]]
      names(mix_20.bmdl)[[i]] <- element_names[[i]]
          },
    error = function(e){ 
      mix_20.drc[[i]] <- drm(Response_expr ~ Dose, data =  mix_20.data[[i]],  type = "continuous" , 
                                 fct = LL2.3())
      mix_20.bmd[[i]] <- bmd(mix_20.drc[[i]], bmr = 0.1, backg = 0.5831192, def = "hybrid")[1]
      mix_20.bmdl[[i]] <- bmd(mix_20.drc[[i]], bmr = 0.1, backg = 0.5831192, def = "hybrid")[2]    
      names(mix_20.drc)[[i]] <- element_names[[i]]
      names(mix_20.bmd)[[i]] <- element_names[[i]]
      names(mix_20.bmdl)[[i]] <- element_names[[i]]
      }
  )
}    


##
rutile_hydrophilic_20.data <- list()
for (i in 2:ncol(DR.rutile_hydrophilic_20)){
  rutile_hydrophilic_20.data[[i]] <- DR.rutile_hydrophilic_20[, c(1, i)] }

# remove first element
rutile_hydrophilic_20.data[sapply(rutile_hydrophilic_20.data, is.null)] <- NULL
# rename elements
element_names <- list()
for (i in 1:length(rutile_hydrophilic_20.data)){
  element_names[i] <- colnames(rutile_hydrophilic_20.data[[i]])[2]
  colnames(rutile_hydrophilic_20.data[[i]]) <- c("Dose", "Response_expr")
}
names(rutile_hydrophilic_20.data) <- element_names

BMDt.rutile_hydrophilic_20.ll2 <- list()
rutile_hydrophilic_20.model <- list()
for (i in 1:length(rutile_hydrophilic_20.data)){
  tryCatch(
    expr = {
      BMDt.rutile_hydrophilic_20.ll2[[i]] <- drm(Response_expr ~ Dose, 
                                      data = rutile_hydrophilic_20.data[[i]],  type = "continuous" , fct = LL.2()) 
      rutile_hydrophilic_20.model[[i]] <- mselect(BMDt.rutile_hydrophilic_20.ll2[[i]], fctList = list(LL.3(), LL.3u(), LL.4(), LL.5(), 
                                                                                LL2.2(), LL2.3(), LL2.3u(), LL2.4(), LL2.5()))
      
    },
    error = function(e){ 
      BMDt.rutile_hydrophilic_20.ll2[[i]] <- drm(Response_expr ~ Dose, 
                                      data = rutile_hydrophilic_20.data[[i]],  type = "continuous" , fct =  LL2.5()) 
      rutile_hydrophilic_20.model[[i]] <- mselect(BMDt.rutile_hydrophilic_20.ll2[[i]], fctList = list(LL.2(), LL.3(), LL.3u(), LL.4(), LL.5(), 
                                                                                LL2.2(), LL2.3(), LL2.3u(), LL2.4(), LL2.5()))
      
    }
  )
}    

rutile_hydrophilic_20.model.best <- list()
for (i in 1:length(rutile_hydrophilic_20.model)){
  rutile_hydrophilic_20.model.best[i] <- rownames(rutile_hydrophilic_20.model[[i]])[1]
}
for (i in 1:length(rutile_hydrophilic_20.model.best)){
  rutile_hydrophilic_20.model.best[[i]] <- ifelse(is.null(rutile_hydrophilic_20.model.best[[i]]), paste("Unavailable"),
                                       rutile_hydrophilic_20.model.best[[i]])
}

rutile_hydrophilic_20.model.best <- unlist(rutile_hydrophilic_20.model.best)
rutile_hydrophilic_20.model.best <- lapply(rutile_hydrophilic_20.model.best, as.factor)

rutile_hydrophilic_20.drc <- list()
rutile_hydrophilic_20.bmd <- list()
rutile_hydrophilic_20.bmdl <- list()

for (i in 1:length(rutile_hydrophilic_20.model.best)){  
  
  tryCatch(
    expr = {
      rutile_hydrophilic_20.drc[[i]] <- drm(Response_expr ~ Dose, data =  rutile_hydrophilic_20.data[[i]],  type = "continuous" , 
                                 fct = LL.5())
      rutile_hydrophilic_20.bmd[[i]] <- bmd(rutile_hydrophilic_20.drc[[i]], bmr = 0.1, backg = 0.5831192, def = "hybrid")[1]
      rutile_hydrophilic_20.bmdl[[i]] <- bmd(rutile_hydrophilic_20.drc[[i]], bmr = 0.1, backg = 0.5831192, def = "hybrid")[2]
      
      names(rutile_hydrophilic_20.drc)[[i]] <- element_names[[i]]
      names(rutile_hydrophilic_20.bmd)[[i]] <- element_names[[i]]
      names(rutile_hydrophilic_20.bmdl)[[i]] <- element_names[[i]]
      },
    error = function(e){ 
      rutile_hydrophilic_20.drc[[i]] <- drm(Response_expr ~ Dose, data =  rutile_hydrophilic_20.data[[i]],  type = "continuous" , 
                                 fct = LL2.4())
      rutile_hydrophilic_20.bmd[[i]] <- bmd(rutile_hydrophilic_20.drc[[i]], bmr = 0.1, backg = 0.5831192, def = "hybrid")[1]
      rutile_hydrophilic_20.bmdl[[i]] <- bmd(rutile_hydrophilic_20.drc[[i]], bmr = 0.1, backg = 0.5831192, def = "hybrid")[2]    
      
      names(rutile_hydrophilic_20.drc)[[i]] <- element_names[[i]]
      names(rutile_hydrophilic_20.bmd)[[i]] <- element_names[[i]]
      names(rutile_hydrophilic_20.bmdl)[[i]] <- element_names[[i]]
      }
  )
}    


##
rutile_hydrophobic_20.data <- list()
for (i in 2:ncol(DR.rutile_hydrophobic_20)){
  rutile_hydrophobic_20.data[[i]] <- DR.rutile_hydrophobic_20[, c(1, i)] }

# remove first element
rutile_hydrophobic_20.data[sapply(rutile_hydrophobic_20.data, is.null)] <- NULL
# rename elements
element_names <- list()
for (i in 1:length(rutile_hydrophobic_20.data)){
  element_names[i] <- colnames(rutile_hydrophobic_20.data[[i]])[2]
  colnames(rutile_hydrophobic_20.data[[i]]) <- c("Dose", "Response_expr")
}
names(rutile_hydrophobic_20.data) <- element_names

BMDt.rutile_hydrophobic_20.ll2 <- list()
rutile_hydrophobic_20.model <- list()
for (i in 1:length(rutile_hydrophobic_20.data)){
  tryCatch(
    expr = {
      BMDt.rutile_hydrophobic_20.ll2[[i]] <- drm(Response_expr ~ Dose, 
                                      data = rutile_hydrophobic_20.data[[i]],  type = "continuous" , fct = LL.2()) 
      rutile_hydrophobic_20.model[[i]] <- mselect(BMDt.rutile_hydrophobic_20.ll2[[i]], fctList = list(LL.3(), LL.3u(), LL.4(), LL.5(), 
                                                                                LL2.2(), LL2.3(), LL2.3u(), LL2.4(), LL2.5()))
      
    },
    error = function(e){ 
      BMDt.rutile_hydrophobic_20.ll2[[i]] <- drm(Response_expr ~ Dose, 
                                      data = rutile_hydrophobic_20.data[[i]],  type = "continuous" , fct =  LL2.5()) 
      rutile_hydrophobic_20.model[[i]] <- mselect(BMDt.rutile_hydrophobic_20.ll2[[i]], fctList = list(LL.2(), LL.3(), LL.3u(), LL.4(), LL.5(), 
                                                                                LL2.2(), LL2.3(), LL2.3u(), LL2.4(), LL2.5()))
      
    }
  )
}    

rutile_hydrophobic_20.model.best <- list()
for (i in 1:length(rutile_hydrophobic_20.model)){
  rutile_hydrophobic_20.model.best[i] <- rownames(rutile_hydrophobic_20.model[[i]])[1]
}
for (i in 1:length(rutile_hydrophobic_20.model.best)){
  rutile_hydrophobic_20.model.best[[i]] <- ifelse(is.null(rutile_hydrophobic_20.model.best[[i]]), paste("Unavailable"),
                                       rutile_hydrophobic_20.model.best[[i]])
}

rutile_hydrophobic_20.model.best <- unlist(rutile_hydrophobic_20.model.best)
rutile_hydrophobic_20.model.best <- lapply(rutile_hydrophobic_20.model.best, as.factor)

rutile_hydrophobic_20.drc <- list()
rutile_hydrophobic_20.bmd <- list()
rutile_hydrophobic_20.bmdl <- list()

for (i in 1:length(rutile_hydrophobic_20.model.best)){  
  
  tryCatch(
    expr = {
      rutile_hydrophobic_20.drc[[i]] <- drm(Response_expr ~ Dose, data =  rutile_hydrophobic_20.data[[i]],  type = "continuous" , 
                                 fct = LL2.3u())
      rutile_hydrophobic_20.bmd[[i]] <- bmd(rutile_hydrophobic_20.drc[[i]], bmr = 0.1, backg = 0.5831192, def = "hybrid")[1]
      rutile_hydrophobic_20.bmdl[[i]] <- bmd(rutile_hydrophobic_20.drc[[i]], bmr = 0.1, backg = 0.5831192, def = "hybrid")[2]
      names(rutile_hydrophobic_20.drc)[[i]] <- element_names[[i]]
      names(rutile_hydrophobic_20.bmd)[[i]] <- element_names[[i]]
      names(rutile_hydrophobic_20.bmdl)[[i]] <- element_names[[i]]   
      },
    error = function(e){ 
      rutile_hydrophobic_20.drc[[i]] <- drm(Response_expr ~ Dose, data =  rutile_hydrophobic_20.data[[i]],  type = "continuous" , 
                                 fct = LL2.5())
      rutile_hydrophobic_20.bmd[[i]] <- bmd(rutile_hydrophobic_20.drc[[i]], bmr = 0.1, backg = 0.5831192, def = "hybrid")[1]
      rutile_hydrophobic_20.bmdl[[i]] <- bmd(rutile_hydrophobic_20.drc[[i]], bmr = 0.1, backg = 0.5831192, def = "hybrid")[2]    
      names(rutile_hydrophobic_20.drc)[[i]] <- element_names[[i]]
      names(rutile_hydrophobic_20.bmd)[[i]] <- element_names[[i]]
      names(rutile_hydrophobic_20.bmdl)[[i]] <- element_names[[i]] 
      }
  )
}    


save.image("rah16_bmd.RData")

# Collect BMDL ------------------------------------------------------------

bmdl <- list(anatase_8.bmdl, anatase_20.bmdl, anatase_300.bmdl, mix_20.bmdl, 
                   rutile_hydrophilic_20.bmdl, rutile_hydrophobic_20.bmdl)
bmd <- list(anatase_8.bmd, anatase_20.bmd, anatase_300.bmd, mix_20.bmd, 
             rutile_hydrophilic_20.bmd, rutile_hydrophobic_20.bmd)

common.probes <- Reduce(intersect, list(names(mix_20.bmdl), names(anatase_20.bmdl), names(anatase_300.bmdl),
                                        names(mix_20.bmdl), names(rutile_hydrophilic_20.bmdl), 
                                        names(rutile_hydrophobic_20.bmdl)))
anatase_8.bmdl.common <- anatase_8.bmdl[common.probes]
anatase_20.bmdl.common <- anatase_20.bmdl[common.probes]
anatase_300.bmdl.common <- anatase_300.bmdl[common.probes]
mix_20.bmdl.common <- mix_20.bmdl[common.probes]
rutile_hydrophilic_20.bmdl.common <- rutile_hydrophilic_20.bmdl[common.probes]
rutile_hydrophobic_20.bmdl.common <- rutile_hydrophobic_20.bmdl[common.probes]

bmdl.common <- list(anatase_8.bmdl.common, anatase_20.bmdl.common, anatase_300.bmdl.common, 
                    mix_20.bmdl.common, rutile_hydrophilic_20.bmdl.common, rutile_hydrophobic_20.bmdl.common)
names(bmdl.common) <- c("anatase_8", "anatase_20", "anatase_300", "mix_20", 
                        "rutile_hydrophilic_20", "rutile_hydrophobic_20")

# Annotate BMDL for common probes between NMs
keys.common <- common.probes
annotate <- AnnotationDbi::select(mgug4122a.db, keys=keys.common, columns=keys, keytype="PROBEID")
symbols <- unique(annotate$SYMBOL)
bmdl.common.annotated <- list()
for (i in 1:length(bmdl.common)){
  bmdl.common.annotated[[i]] <- unlist(bmdl.common[[i]]) 
  bmdl.common.annotated[[i]]<- merge(as.matrix(bmdl.common.annotated[[i]]), annotate, by.x=0, by.y="PROBEID")
}
names(bmdl.common.annotated) <- c("anatase_8", "anatase_20", "anatase_300", "mix_20", 
                        "rutile_hydrophilic_20", "rutile_hydrophobic_20")


# GSEA
bmdl.go <- "GO_Biological_Process_2018"
bmdl.pathways<- c("KEGG_2019_Mouse", "WikiPathways_2019_Mouse")

enriched.go.bmdl <- enrichr(symbols, bmdl.go)
enriched.pathways.bmdl <- enrichr(symbols, bmdl.pathways)

enriched.go.bmdl$GO_Biological_Process_2018$Genes_count <- sapply(strsplit(as.character(enriched.go.bmdl$GO_Biological_Process_2018$Overlap),'/'), "[", 1)
enriched.go.bmdl.5 <- enriched.go.bmdl$GO_Biological_Process_2018[enriched.go.bmdl$GO_Biological_Process_2018$Genes_count >4, ]

enriched.go.bmdl.5 <- enriched.go.bmdl.5[order(enriched.go.bmdl.5$Genes_count, decreasing = T), ]
enriched.go.bmdl.counts10 <- enriched.go.bmdl.5[1:10, ]
enriched.go.bmdl.5 <- enriched.go.bmdl.5[order(enriched.go.bmdl.5$Old.Adjusted.P.value, decreasing = F), ]
enriched.go.bmdl.pvalue10 <- enriched.go.bmdl.5[1:10, ]
enriched.go.bmdl.5 <- enriched.go.bmdl.5[order(enriched.go.bmdl.5$Combined.Score, decreasing = F), ]
enriched.go.bmdl.score10 <- enriched.go.bmdl.5[1:10, ]

enriched.kegg.bmdl <- enriched.pathways.bmdl$KEGG_2019_Mouse
enriched.kegg.bmdl$Genes_count <- sapply(strsplit(as.character(enriched.kegg.bmdl$Overlap),'/'), "[", 1)
enriched.kegg.bmdl.5 <- enriched.kegg.bmdl[enriched.kegg.bmdl$Genes_count >4, ]

enriched.kegg.bmdl.5 <- enriched.kegg.bmdl.5[order(enriched.kegg.bmdl.5$Genes_count, decreasing = T), ]
enriched.kegg.bmdl.counts10 <- enriched.kegg.bmdl.5[1:10, ]
enriched.kegg.bmdl.5 <- enriched.kegg.bmdl.5[order(enriched.kegg.bmdl.5$Old.Adjusted.P.value, decreasing = F), ]
enriched.kegg.bmdl.pvalue10 <- enriched.kegg.bmdl.5[1:10, ]
enriched.kegg.bmdl.5 <- enriched.kegg.bmdl.5[order(enriched.kegg.bmdl.5$Combined.Score, decreasing = F), ]
enriched.kegg.bmdl.score10 <- enriched.kegg.bmdl.5[1:10, ]

enriched.wikipathways.bmdl <- enriched.pathways.bmdl$WikiPathways_2019_Mouse
enriched.wikipathways.bmdl$Genes_count <- sapply(strsplit(as.character(enriched.wikipathways.bmdl$Overlap),'/'), "[", 1)
enriched.wikipathways.bmdl.5 <- enriched.wikipathways.bmdl[enriched.wikipathways.bmdl$Genes_count >3, ] 
# only 6 when taking only the pathways with more than 5 genes perturbed, so a lower threshold by 1 gene, namely 4 genes/pathway, is used

enriched.wikipathways.bmdl.5 <- enriched.wikipathways.bmdl.5[order(enriched.wikipathways.bmdl.5$Genes_count, decreasing = T), ]
enriched.wikipathways.bmdl.counts10 <- enriched.wikipathways.bmdl.5[1:10, ]
enriched.wikipathways.bmdl.5 <- enriched.wikipathways.bmdl.5[order(enriched.wikipathways.bmdl.5$Old.Adjusted.P.value, decreasing = F), ]
enriched.wikipathways.bmdl.pvalue10 <- enriched.wikipathways.bmdl.5[1:10, ]
enriched.wikipathways.bmdl.5 <- enriched.wikipathways.bmdl.5[order(enriched.wikipathways.bmdl.5$Combined.Score, decreasing = F), ]
enriched.wikipathways.bmdl.score10 <- enriched.wikipathways.bmdl.5[1:10, ]

# CTD ---------------------------------------------------------------------

# search for LUNG FIBROSIS term to get disease associated gene list 
ctd.fibrosis <- CTDquerier::query_ctd_dise(terms = "Pulmonary Fibrosis", verbose = TRUE )

#ctd.pathways <- read.csv( file = "CTD_D011658_pathways_20191218242937.csv")
#ctd.pathways <- as.matrix(ctd.pathways)
#ctd.kegg <- ctd.pathways[str_detect(ctd.pathways[ ,"Pathway.ID"], pattern = "KEGG"), ]

fibrosis.genes <- ctd.fibrosis@gene_interactions
fibrosis.chemicals <- ctd.fibrosis@chemicals_interactions
fibrosis.pathways <- ctd.fibrosis@kegg
fibrosis.kegg <- fibrosis.pathways[str_detect(fibrosis.pathways[ ,"Pathway.ID"], pattern = "KEGG"), ]
fibrosis.reactome <- fibrosis.pathways[str_detect(fibrosis.pathways[ ,"Pathway.ID"], pattern = "REACT"), ]

# BMD for pathways --------------------------------------------------------

go.genes <- strsplit(enriched.go.bmdl.score10$Genes, ";")
names(go.genes) <- enriched.go.bmdl.score10$Term
go.genes.all <- strsplit(enriched.go.bmdl.5$Genes, ";")
names(go.genes.all) <- enriched.go.bmdl.5$Term

kegg.genes <- strsplit(enriched.kegg.bmdl.score10$Genes, ";")
names(kegg.genes) <- enriched.kegg.bmdl.score10$Term
kegg.genes.all <- strsplit(enriched.kegg.bmdl.5$Genes, ";")
names(kegg.genes.all) <- enriched.kegg.bmdl.5$Term

wikipathways.genes <- strsplit(enriched.wikipathways.bmdl.score10$Genes, ";")
names(wikipathways.genes) <- enriched.wikipathways.bmdl.score10$Term
wikipathways.genes.all <- strsplit(enriched.wikipathways.bmdl.5$Genes, ";")
names(wikipathways.genes.all) <- enriched.wikipathways.bmdl.5$Term

bmdl.anatase_8 <- bmdl.common.annotated$anatase_8
bmdl.anatase_20 <- bmdl.common.annotated$anatase_20
bmdl.anatase_300 <- bmdl.common.annotated$anatase_300
bmdl.mix_20 <- bmdl.common.annotated$mix_20
bmdl.rutile_hydrophilic_20 <- bmdl.common.annotated$rutile_hydrophilic_20
bmdl.rutile_hydrophobic_20 <- bmdl.common.annotated$rutile_hydrophobic_20

bmdl.anatase_8$SYMBOL <- toupper(bmdl.anatase_8$SYMBOL)
bmdl.anatase_20$SYMBOL <- toupper(bmdl.anatase_20$SYMBOL)
bmdl.anatase_300$SYMBOL <- toupper(bmdl.anatase_300$SYMBOL)
bmdl.mix_20$SYMBOL <- toupper(bmdl.mix_20$SYMBOL)
bmdl.rutile_hydrophilic_20$SYMBOL <- toupper(bmdl.rutile_hydrophilic_20$SYMBOL)
bmdl.rutile_hydrophobic_20$SYMBOL <- toupper(bmdl.rutile_hydrophobic_20$SYMBOL)

colnames(bmdl.anatase_8)[2] <- "BMDL"
colnames(bmdl.anatase_20)[2] <- "BMDL"
colnames(bmdl.anatase_300)[2] <- "BMDL"
colnames(bmdl.mix_20)[2] <- "BMDL"
colnames(bmdl.rutile_hydrophilic_20)[2] <- "BMDL"
colnames(bmdl.rutile_hydrophobic_20)[2] <- "BMDL"

bmdl.anatase_8 <- na.omit(bmdl.anatase_8)
bmdl.anatase_20 <- na.omit(bmdl.anatase_20)
bmdl.anatase_300 <- na.omit(bmdl.anatase_300)
bmdl.mix_20 <- na.omit(bmdl.mix_20)
bmdl.rutile_hydrophilic_20 <- na.omit(bmdl.rutile_hydrophilic_20)
bmdl.rutile_hydrophobic_20 <- na.omit(bmdl.rutile_hydrophobic_20)

bmdl.anatase_8 <- bmdl.anatase_8[bmdl.anatase_8$BMDL>= 0 | bmdl.anatase_8 < 486, ]
bmdl.anatase_20 <- bmdl.anatase_20[bmdl.anatase_20$BMDL>= 0 | bmdl.anatase_20 < 486, ]
bmdl.anatase_300 <- bmdl.anatase_300[bmdl.anatase_300$BMDL>= 0 | bmdl.anatase_300 < 486, ]
bmdl.mix_20 <- bmdl.mix_20[bmdl.mix_20$BMDL>= 0 | bmdl.mix_20 < 486, ]
bmdl.rutile_hydrophilic_20 <- bmdl.rutile_hydrophilic_20[bmdl.rutile_hydrophilic_20 >= 0 |bmdl.rutile_hydrophilic_20 < 486, ]
bmdl.rutile_hydrophobic_20 <- bmdl.rutile_hydrophobic_20[bmdl.rutile_hydrophobic_20 >= 0 |bmdl.rutile_hydrophobic_20 < 486, ]

bmdl.anatase_8 <- bmdl.anatase_8[, -c(6, 8, 9)]
bmdl.anatase_8 <- bmdl.anatase_8[!duplicated(bmdl.anatase_8), ]

bmdl.anatase_20 <- bmdl.anatase_20[, -c(6, 8, 9)]
bmdl.anatase_20 <- bmdl.anatase_20[!duplicated(bmdl.anatase_20), ]

bmdl.anatase_300 <- bmdl.anatase_300[, -c(6, 8, 9)]
bmdl.anatase_300 <- bmdl.anatase_300[!duplicated(bmdl.anatase_300), ]

bmdl.mix_20 <- bmdl.mix_20[, -c(6, 8, 9)]
bmdl.mix_20 <- bmdl.mix_20[!duplicated(bmdl.mix_20), ]

bmdl.rutile_hydrophilic_20 <- bmdl.rutile_hydrophilic_20[, -c(6, 8, 9)]
bmdl.rutile_hydrophilic_20 <- bmdl.rutile_hydrophilic_20[!duplicated(bmdl.rutile_hydrophilic_20), ]

bmdl.rutile_hydrophobic_20 <- bmdl.rutile_hydrophobic_20[, -c(6, 8, 9)]
bmdl.rutile_hydrophobic_20 <- bmdl.rutile_hydrophobic_20[!duplicated(bmdl.rutile_hydrophobic_20), ]

go.bmdl.anatase_8 <- list()
go.bmdl.anatase_8.medians <- list()
for (i in 1:length(go.genes)){
  go.bmdl.anatase_8[[i]] <- bmdl.anatase_8[which(bmdl.anatase_8$SYMBOL %in% go.genes[[i]]), 2]
  go.bmdl.anatase_8.medians[i] <- median(go.bmdl.anatase_8[[i]])
  }
names(go.bmdl.anatase_8) <- names(go.genes)
names(go.bmdl.anatase_8.medians) <- names(go.genes)

go.all.bmdl.anatase_8 <- list()
go.all.bmdl.anatase_8.medians <- list()
for (i in 1:length(go.genes.all)){
  go.all.bmdl.anatase_8[[i]] <- bmdl.anatase_8[which(bmdl.anatase_8$SYMBOL %in% go.genes.all[[i]]), 2]
  go.all.bmdl.anatase_8.medians[i] <- median(go.all.bmdl.anatase_8[[i]])
}
names(go.all.bmdl.anatase_8) <- names(go.genes.all)
names(go.all.bmdl.anatase_8.medians) <- names(go.genes.all)

kegg.bmdl.anatase_8 <- list()
kegg.bmdl.anatase_8.medians <- list()
for (i in 1:length(kegg.genes)){
  kegg.bmdl.anatase_8[[i]] <- bmdl.anatase_8[which(bmdl.anatase_8$SYMBOL %in% kegg.genes[[i]]), 2]
  kegg.bmdl.anatase_8.medians[i] <- median(kegg.bmdl.anatase_8[[i]])
}
names(kegg.bmdl.anatase_8) <- names(kegg.genes)
names(kegg.bmdl.anatase_8.medians) <- names(kegg.genes)

kegg.all.bmdl.anatase_8 <- list()
kegg.all.bmdl.anatase_8.medians <- list()
for (i in 1:length(kegg.genes.all)){
  kegg.all.bmdl.anatase_8[[i]] <- bmdl.anatase_8[which(bmdl.anatase_8$SYMBOL %in% kegg.genes.all[[i]]), 2]
  kegg.all.bmdl.anatase_8.medians[i] <- median(kegg.all.bmdl.anatase_8[[i]])
}
names(kegg.all.bmdl.anatase_8) <- names(kegg.genes.all)
names(kegg.all.bmdl.anatase_8.medians) <- names(kegg.genes.all)

wikipathways.bmdl.anatase_8 <- list()
wikipathways.bmdl.anatase_8.medians <- list()
for (i in 1:length(wikipathways.genes)){
  wikipathways.bmdl.anatase_8[[i]] <- bmdl.anatase_8[which(bmdl.anatase_8$SYMBOL %in% wikipathways.genes[[i]]), 2]
  wikipathways.bmdl.anatase_8.medians[i] <- median(wikipathways.bmdl.anatase_8[[i]])
}
names(wikipathways.bmdl.anatase_8) <- names(wikipathways.genes)
names(wikipathways.bmdl.anatase_8.medians) <- names(wikipathways.genes)

wikipathways.all.bmdl.anatase_8 <- list()
wikipathways.all.bmdl.anatase_8.medians <- list()
for (i in 1:length(wikipathways.genes.all)){
  wikipathways.all.bmdl.anatase_8[[i]] <- bmdl.anatase_8[which(bmdl.anatase_8$SYMBOL %in% wikipathways.genes.all[[i]]), 2]
  wikipathways.all.bmdl.anatase_8.medians[i] <- median(wikipathways.all.bmdl.anatase_8[[i]])
}
names(wikipathways.all.bmdl.anatase_8) <- names(wikipathways.genes.all)
names(wikipathways.all.bmdl.anatase_8.medians) <- names(wikipathways.genes.all)

kegg.fibrosis.bmdl.anatase_8 <- kegg.all.bmdl.anatase_8[names(kegg.all.bmdl.anatase_8) %in% fibrosis.kegg$Pathway]
kegg.fibrosis.bmdl.anatase_8.medians <- kegg.all.bmdl.anatase_8.medians[names(kegg.all.bmdl.anatase_8.medians) %in% fibrosis.kegg$Pathway]

go.bmdl.anatase_20 <- list()
go.bmdl.anatase_20.medians <- list()
for (i in 1:length(go.genes)){
  go.bmdl.anatase_20[[i]] <- bmdl.anatase_20[which(bmdl.anatase_20$SYMBOL %in% go.genes[[i]]), 2]
  go.bmdl.anatase_20.medians[i] <- median(go.bmdl.anatase_20[[i]])
}
names(go.bmdl.anatase_20) <- names(go.genes)
names(go.bmdl.anatase_20.medians) <- names(go.genes)

go.all.bmdl.anatase_20 <- list()
go.all.bmdl.anatase_20.medians <- list()
for (i in 1:length(go.genes.all)){
  go.all.bmdl.anatase_20[[i]] <- bmdl.anatase_20[which(bmdl.anatase_20$SYMBOL %in% go.genes.all[[i]]), 2]
  go.all.bmdl.anatase_20.medians[i] <- median(go.all.bmdl.anatase_20[[i]])
}
names(go.all.bmdl.anatase_20) <- names(go.genes.all)
names(go.all.bmdl.anatase_20.medians) <- names(go.genes.all)

kegg.bmdl.anatase_20 <- list()
kegg.bmdl.anatase_20.medians <- list()
for (i in 1:length(kegg.genes)){
  kegg.bmdl.anatase_20[[i]] <- bmdl.anatase_20[which(bmdl.anatase_20$SYMBOL %in% kegg.genes[[i]]), 2]
  kegg.bmdl.anatase_20.medians[i] <- median(kegg.bmdl.anatase_20[[i]])
}
names(kegg.bmdl.anatase_20) <- names(kegg.genes)
names(kegg.bmdl.anatase_20.medians) <- names(kegg.genes)

kegg.all.bmdl.anatase_20 <- list()
kegg.all.bmdl.anatase_20.medians <- list()
for (i in 1:length(kegg.genes.all)){
  kegg.all.bmdl.anatase_20[[i]] <- bmdl.anatase_20[which(bmdl.anatase_20$SYMBOL %in% kegg.genes.all[[i]]), 2]
  kegg.all.bmdl.anatase_20.medians[i] <- median(kegg.all.bmdl.anatase_20[[i]])
}
names(kegg.all.bmdl.anatase_20) <- names(kegg.genes.all)
names(kegg.all.bmdl.anatase_20.medians) <- names(kegg.genes.all)

wikipathways.bmdl.anatase_20 <- list()
wikipathways.bmdl.anatase_20.medians <- list()
for (i in 1:length(wikipathways.genes)){
  wikipathways.bmdl.anatase_20[[i]] <- bmdl.anatase_20[which(bmdl.anatase_20$SYMBOL %in% wikipathways.genes[[i]]), 2]
  wikipathways.bmdl.anatase_20.medians[i] <- median(wikipathways.bmdl.anatase_20[[i]])
}
names(wikipathways.bmdl.anatase_20) <- names(wikipathways.genes)
names(wikipathways.bmdl.anatase_20.medians) <- names(wikipathways.genes)

wikipathways.all.bmdl.anatase_20 <- list()
wikipathways.all.bmdl.anatase_20.medians <- list()
for (i in 1:length(wikipathways.genes.all)){
  wikipathways.all.bmdl.anatase_20[[i]] <- bmdl.anatase_20[which(bmdl.anatase_20$SYMBOL %in% wikipathways.genes.all[[i]]), 2]
  wikipathways.all.bmdl.anatase_20.medians[i] <- median(wikipathways.all.bmdl.anatase_20[[i]])
}
names(wikipathways.all.bmdl.anatase_20) <- names(wikipathways.genes.all)
names(wikipathways.all.bmdl.anatase_20.medians) <- names(wikipathways.genes.all)

kegg.fibrosis.bmdl.anatase_20 <- kegg.all.bmdl.anatase_20[names(kegg.all.bmdl.anatase_20) %in% fibrosis.kegg$Pathway]
kegg.fibrosis.bmdl.anatase_20.medians <- kegg.all.bmdl.anatase_20.medians[names(kegg.all.bmdl.anatase_20.medians) %in% fibrosis.kegg$Pathway]

go.bmdl.anatase_300 <- list()
go.bmdl.anatase_300.medians <- list()
for (i in 1:length(go.genes)){
  go.bmdl.anatase_300[[i]] <- bmdl.anatase_300[which(bmdl.anatase_300$SYMBOL %in% go.genes[[i]]), 2]
  go.bmdl.anatase_300.medians[i] <- median(go.bmdl.anatase_300[[i]])
}
names(go.bmdl.anatase_300) <- names(go.genes)
names(go.bmdl.anatase_300.medians) <- names(go.genes)

go.all.bmdl.anatase_300 <- list()
go.all.bmdl.anatase_300.medians <- list()
for (i in 1:length(go.genes.all)){
  go.all.bmdl.anatase_300[[i]] <- bmdl.anatase_300[which(bmdl.anatase_300$SYMBOL %in% go.genes.all[[i]]), 2]
  go.all.bmdl.anatase_300.medians[i] <- median(go.all.bmdl.anatase_300[[i]])
}
names(go.all.bmdl.anatase_300) <- names(go.genes.all)
names(go.all.bmdl.anatase_300.medians) <- names(go.genes.all)

kegg.bmdl.anatase_300 <- list()
kegg.bmdl.anatase_300.medians <- list()
for (i in 1:length(kegg.genes)){
  kegg.bmdl.anatase_300[[i]] <- bmdl.anatase_300[which(bmdl.anatase_300$SYMBOL %in% kegg.genes[[i]]), 2]
  kegg.bmdl.anatase_300.medians[i] <- median(kegg.bmdl.anatase_300[[i]])
}
names(kegg.bmdl.anatase_300) <- names(kegg.genes)
names(kegg.bmdl.anatase_300.medians) <- names(kegg.genes)

kegg.all.bmdl.anatase_300 <- list()
kegg.all.bmdl.anatase_300.medians <- list()
for (i in 1:length(kegg.genes.all)){
  kegg.all.bmdl.anatase_300[[i]] <- bmdl.anatase_300[which(bmdl.anatase_300$SYMBOL %in% kegg.genes.all[[i]]), 2]
  kegg.all.bmdl.anatase_300.medians[i] <- median(kegg.all.bmdl.anatase_300[[i]])
}
names(kegg.all.bmdl.anatase_300) <- names(kegg.genes.all)
names(kegg.all.bmdl.anatase_300.medians) <- names(kegg.genes.all)

wikipathways.bmdl.anatase_300 <- list()
wikipathways.bmdl.anatase_300.medians <- list()
for (i in 1:length(wikipathways.genes)){
  wikipathways.bmdl.anatase_300[[i]] <- bmdl.anatase_300[which(bmdl.anatase_300$SYMBOL %in% wikipathways.genes[[i]]), 2]
  wikipathways.bmdl.anatase_300.medians[i] <- median(wikipathways.bmdl.anatase_300[[i]])
}
names(wikipathways.bmdl.anatase_300) <- names(wikipathways.genes)
names(wikipathways.bmdl.anatase_300.medians) <- names(wikipathways.genes)

wikipathways.all.bmdl.anatase_300 <- list()
wikipathways.all.bmdl.anatase_300.medians <- list()
for (i in 1:length(wikipathways.genes.all)){
  wikipathways.all.bmdl.anatase_300[[i]] <- bmdl.anatase_300[which(bmdl.anatase_300$SYMBOL %in% wikipathways.genes.all[[i]]), 2]
  wikipathways.all.bmdl.anatase_300.medians[i] <- median(wikipathways.all.bmdl.anatase_300[[i]])
}
names(wikipathways.all.bmdl.anatase_300) <- names(wikipathways.genes.all)
names(wikipathways.all.bmdl.anatase_300.medians) <- names(wikipathways.genes.all)

kegg.fibrosis.bmdl.anatase_300 <- kegg.all.bmdl.anatase_300[names(kegg.all.bmdl.anatase_300) %in% fibrosis.kegg$Pathway]
kegg.fibrosis.bmdl.anatase_300.medians <- kegg.all.bmdl.anatase_300.medians[names(kegg.all.bmdl.anatase_300.medians) %in% fibrosis.kegg$Pathway]

go.bmdl.mix_20 <- list()
go.bmdl.mix_20.medians <- list()
for (i in 1:length(go.genes)){
  go.bmdl.mix_20[[i]] <- bmdl.mix_20[which(bmdl.mix_20$SYMBOL %in% go.genes[[i]]), 2]
  go.bmdl.mix_20.medians[i] <- median(go.bmdl.mix_20[[i]])
}
names(go.bmdl.mix_20) <- names(go.genes)
names(go.bmdl.mix_20.medians) <- names(go.genes)

go.all.bmdl.mix_20 <- list()
go.all.bmdl.mix_20.medians <- list()
for (i in 1:length(go.genes.all)){
  go.all.bmdl.mix_20[[i]] <- bmdl.mix_20[which(bmdl.mix_20$SYMBOL %in% go.genes.all[[i]]), 2]
  go.all.bmdl.mix_20.medians[i] <- median(go.all.bmdl.mix_20[[i]])
}
names(go.all.bmdl.mix_20) <- names(go.genes.all)
names(go.all.bmdl.mix_20.medians) <- names(go.genes.all)

kegg.bmdl.mix_20 <- list()
kegg.bmdl.mix_20.medians <- list()
for (i in 1:length(kegg.genes)){
  kegg.bmdl.mix_20[[i]] <- bmdl.mix_20[which(bmdl.mix_20$SYMBOL %in% kegg.genes[[i]]), 2]
  kegg.bmdl.mix_20.medians[i] <- median(kegg.bmdl.mix_20[[i]])
}
names(kegg.bmdl.mix_20) <- names(kegg.genes)
names(kegg.bmdl.mix_20.medians) <- names(kegg.genes)

kegg.all.bmdl.mix_20 <- list()
kegg.all.bmdl.mix_20.medians <- list()
for (i in 1:length(kegg.genes.all)){
  kegg.all.bmdl.mix_20[[i]] <- bmdl.mix_20[which(bmdl.mix_20$SYMBOL %in% kegg.genes.all[[i]]), 2]
  kegg.all.bmdl.mix_20.medians[i] <- median(kegg.all.bmdl.mix_20[[i]])
}
names(kegg.all.bmdl.mix_20) <- names(kegg.genes.all)
names(kegg.all.bmdl.mix_20.medians) <- names(kegg.genes.all)

wikipathways.bmdl.mix_20 <- list()
wikipathways.bmdl.mix_20.medians <- list()
for (i in 1:length(wikipathways.genes)){
  wikipathways.bmdl.mix_20[[i]] <- bmdl.mix_20[which(bmdl.mix_20$SYMBOL %in% wikipathways.genes[[i]]), 2]
  wikipathways.bmdl.mix_20.medians[i] <- median(wikipathways.bmdl.mix_20[[i]])
}
names(wikipathways.bmdl.mix_20) <- names(wikipathways.genes)
names(wikipathways.bmdl.mix_20.medians) <- names(wikipathways.genes)

wikipathways.all.bmdl.mix_20 <- list()
wikipathways.all.bmdl.mix_20.medians <- list()
for (i in 1:length(wikipathways.genes.all)){
  wikipathways.all.bmdl.mix_20[[i]] <- bmdl.mix_20[which(bmdl.mix_20$SYMBOL %in% wikipathways.genes.all[[i]]), 2]
  wikipathways.all.bmdl.mix_20.medians[i] <- median(wikipathways.all.bmdl.mix_20[[i]])
}
names(wikipathways.all.bmdl.mix_20) <- names(wikipathways.genes.all)
names(wikipathways.all.bmdl.mix_20.medians) <- names(wikipathways.genes.all)

kegg.fibrosis.bmdl.mix_20 <- kegg.all.bmdl.mix_20[names(kegg.all.bmdl.mix_20) %in% fibrosis.kegg$Pathway]
kegg.fibrosis.bmdl.mix_20.medians <- kegg.all.bmdl.mix_20.medians[names(kegg.all.bmdl.mix_20.medians) %in% fibrosis.kegg$Pathway]

go.bmdl.rutile_hydrophilic_20 <- list()
go.bmdl.rutile_hydrophilic_20.medians <- list()
for (i in 1:length(go.genes)){
  go.bmdl.rutile_hydrophilic_20[[i]] <- bmdl.rutile_hydrophilic_20[which(bmdl.rutile_hydrophilic_20$SYMBOL %in% go.genes[[i]]), 2]
  go.bmdl.rutile_hydrophilic_20.medians[i] <- median(go.bmdl.rutile_hydrophilic_20[[i]])
}
names(go.bmdl.rutile_hydrophilic_20) <- names(go.genes)
names(go.bmdl.rutile_hydrophilic_20.medians) <- names(go.genes)

go.all.bmdl.rutile_hydrophilic_20 <- list()
go.all.bmdl.rutile_hydrophilic_20.medians <- list()
for (i in 1:length(go.genes.all)){
  go.all.bmdl.rutile_hydrophilic_20[[i]] <- bmdl.rutile_hydrophilic_20[which(bmdl.rutile_hydrophilic_20$SYMBOL %in% go.genes.all[[i]]), 2]
  go.all.bmdl.rutile_hydrophilic_20.medians[i] <- median(go.all.bmdl.rutile_hydrophilic_20[[i]])
}
names(go.all.bmdl.rutile_hydrophilic_20) <- names(go.genes.all)
names(go.all.bmdl.rutile_hydrophilic_20.medians) <- names(go.genes.all)

kegg.bmdl.rutile_hydrophilic_20 <- list()
kegg.bmdl.rutile_hydrophilic_20.medians <- list()
for (i in 1:length(kegg.genes)){
  kegg.bmdl.rutile_hydrophilic_20[[i]] <- bmdl.rutile_hydrophilic_20[which(bmdl.rutile_hydrophilic_20$SYMBOL %in% kegg.genes[[i]]), 2]
  kegg.bmdl.rutile_hydrophilic_20.medians[i] <- median(kegg.bmdl.rutile_hydrophilic_20[[i]])
}
names(kegg.bmdl.rutile_hydrophilic_20) <- names(kegg.genes)
names(kegg.bmdl.rutile_hydrophilic_20.medians) <- names(kegg.genes)

kegg.all.bmdl.rutile_hydrophilic_20 <- list()
kegg.all.bmdl.rutile_hydrophilic_20.medians <- list()
for (i in 1:length(kegg.genes.all)){
  kegg.all.bmdl.rutile_hydrophilic_20[[i]] <- bmdl.rutile_hydrophilic_20[which(bmdl.rutile_hydrophilic_20$SYMBOL %in% kegg.genes.all[[i]]), 2]
  kegg.all.bmdl.rutile_hydrophilic_20.medians[i] <- median(kegg.all.bmdl.rutile_hydrophilic_20[[i]])
}
names(kegg.all.bmdl.rutile_hydrophilic_20) <- names(kegg.genes.all)
names(kegg.all.bmdl.rutile_hydrophilic_20.medians) <- names(kegg.genes.all)

wikipathways.bmdl.rutile_hydrophilic_20 <- list()
wikipathways.bmdl.rutile_hydrophilic_20.medians <- list()
for (i in 1:length(wikipathways.genes)){
  wikipathways.bmdl.rutile_hydrophilic_20[[i]] <- bmdl.rutile_hydrophilic_20[which(bmdl.rutile_hydrophilic_20$SYMBOL %in% wikipathways.genes[[i]]), 2]
  wikipathways.bmdl.rutile_hydrophilic_20.medians[i] <- median(wikipathways.bmdl.rutile_hydrophilic_20[[i]])
}
names(wikipathways.bmdl.rutile_hydrophilic_20) <- names(wikipathways.genes)
names(wikipathways.bmdl.rutile_hydrophilic_20.medians) <- names(wikipathways.genes)

wikipathways.all.bmdl.rutile_hydrophilic_20 <- list()
wikipathways.all.bmdl.rutile_hydrophilic_20.medians <- list()
for (i in 1:length(wikipathways.genes.all)){
  wikipathways.all.bmdl.rutile_hydrophilic_20[[i]] <- bmdl.rutile_hydrophilic_20[which(bmdl.rutile_hydrophilic_20$SYMBOL %in% wikipathways.genes.all[[i]]), 2]
  wikipathways.all.bmdl.rutile_hydrophilic_20.medians[i] <- median(wikipathways.all.bmdl.rutile_hydrophilic_20[[i]])
}
names(wikipathways.all.bmdl.rutile_hydrophilic_20) <- names(wikipathways.genes.all)
names(wikipathways.all.bmdl.rutile_hydrophilic_20.medians) <- names(wikipathways.genes.all)

kegg.fibrosis.bmdl.rutile_hydrophilic_20 <- kegg.all.bmdl.rutile_hydrophilic_20[names(kegg.all.bmdl.rutile_hydrophilic_20) %in% fibrosis.kegg$Pathway]
kegg.fibrosis.bmdl.rutile_hydrophilic_20.medians <- kegg.all.bmdl.rutile_hydrophilic_20.medians[names(kegg.all.bmdl.rutile_hydrophilic_20.medians) %in% fibrosis.kegg$Pathway]

go.bmdl.rutile_hydrophobic_20 <- list()
go.bmdl.rutile_hydrophobic_20.medians <- list()
for (i in 1:length(go.genes)){
  go.bmdl.rutile_hydrophobic_20[[i]] <- bmdl.rutile_hydrophobic_20[which(bmdl.rutile_hydrophobic_20$SYMBOL %in% go.genes[[i]]), 2]
  go.bmdl.rutile_hydrophobic_20.medians[i] <- median(go.bmdl.rutile_hydrophobic_20[[i]])
}
names(go.bmdl.rutile_hydrophobic_20) <- names(go.genes)
names(go.bmdl.rutile_hydrophobic_20.medians) <- names(go.genes)

go.all.bmdl.rutile_hydrophobic_20 <- list()
go.all.bmdl.rutile_hydrophobic_20.medians <- list()
for (i in 1:length(go.genes.all)){
  go.all.bmdl.rutile_hydrophobic_20[[i]] <- bmdl.rutile_hydrophobic_20[which(bmdl.rutile_hydrophobic_20$SYMBOL %in% go.genes.all[[i]]), 2]
  go.all.bmdl.rutile_hydrophobic_20.medians[i] <- median(go.all.bmdl.rutile_hydrophobic_20[[i]])
}
names(go.all.bmdl.rutile_hydrophobic_20) <- names(go.genes.all)
names(go.all.bmdl.rutile_hydrophobic_20.medians) <- names(go.genes.all)

kegg.bmdl.rutile_hydrophobic_20 <- list()
kegg.bmdl.rutile_hydrophobic_20.medians <- list()
for (i in 1:length(kegg.genes)){
  kegg.bmdl.rutile_hydrophobic_20[[i]] <- bmdl.rutile_hydrophobic_20[which(bmdl.rutile_hydrophobic_20$SYMBOL %in% kegg.genes[[i]]), 2]
  kegg.bmdl.rutile_hydrophobic_20.medians[i] <- median(kegg.bmdl.rutile_hydrophobic_20[[i]])
}
names(kegg.bmdl.rutile_hydrophobic_20) <- names(kegg.genes)
names(kegg.bmdl.rutile_hydrophobic_20.medians) <- names(kegg.genes)

kegg.all.bmdl.rutile_hydrophobic_20 <- list()
kegg.all.bmdl.rutile_hydrophobic_20.medians <- list()
for (i in 1:length(kegg.genes.all)){
  kegg.all.bmdl.rutile_hydrophobic_20[[i]] <- bmdl.rutile_hydrophobic_20[which(bmdl.rutile_hydrophobic_20$SYMBOL %in% kegg.genes.all[[i]]), 2]
  kegg.all.bmdl.rutile_hydrophobic_20.medians[i] <- median(kegg.all.bmdl.rutile_hydrophobic_20[[i]])
}
names(kegg.all.bmdl.rutile_hydrophobic_20) <- names(kegg.genes.all)
names(kegg.all.bmdl.rutile_hydrophobic_20.medians) <- names(kegg.genes.all)

wikipathways.bmdl.rutile_hydrophobic_20 <- list()
wikipathways.bmdl.rutile_hydrophobic_20.medians <- list()
for (i in 1:length(wikipathways.genes)){
  wikipathways.bmdl.rutile_hydrophobic_20[[i]] <- bmdl.rutile_hydrophobic_20[which(bmdl.rutile_hydrophobic_20$SYMBOL %in% wikipathways.genes[[i]]), 2]
  wikipathways.bmdl.rutile_hydrophobic_20.medians[i] <- median(wikipathways.bmdl.rutile_hydrophobic_20[[i]])
}
names(wikipathways.bmdl.rutile_hydrophobic_20) <- names(wikipathways.genes)
names(wikipathways.bmdl.rutile_hydrophobic_20.medians) <- names(wikipathways.genes)

wikipathways.all.bmdl.rutile_hydrophobic_20 <- list()
wikipathways.all.bmdl.rutile_hydrophobic_20.medians <- list()
for (i in 1:length(wikipathways.genes.all)){
  wikipathways.all.bmdl.rutile_hydrophobic_20[[i]] <- bmdl.rutile_hydrophobic_20[which(bmdl.rutile_hydrophobic_20$SYMBOL %in% wikipathways.genes.all[[i]]), 2]
  wikipathways.all.bmdl.rutile_hydrophobic_20.medians[i] <- median(wikipathways.all.bmdl.rutile_hydrophobic_20[[i]])
}
names(wikipathways.all.bmdl.rutile_hydrophobic_20) <- names(wikipathways.genes.all)
names(wikipathways.all.bmdl.rutile_hydrophobic_20.medians) <- names(wikipathways.genes.all)

kegg.fibrosis.bmdl.rutile_hydrophobic_20 <- kegg.all.bmdl.rutile_hydrophobic_20[names(kegg.all.bmdl.rutile_hydrophobic_20) %in% fibrosis.kegg$Pathway]
kegg.fibrosis.bmdl.rutile_hydrophobic_20.medians <- kegg.all.bmdl.rutile_hydrophobic_20.medians[names(kegg.all.bmdl.rutile_hydrophobic_20.medians) %in% fibrosis.kegg$Pathway]

save.image("rah16_bmdl.RData")


# Clustering --------------------------------------------------------------

df.medians <- matrix(c(unlist(go.all.bmdl.anatase_8.medians),
                       unlist(go.all.bmdl.anatase_20.medians),
                       unlist(go.all.bmdl.anatase_300.medians),
                       unlist(go.all.bmdl.mix_20.medians),
                       unlist(go.all.bmdl.rutile_hydrophilic_20.medians),
                       unlist(go.all.bmdl.rutile_hydrophobic_20.medians)),
                     nrow = 6, ncol = 44)
rownames(df.medians) <- c("anatase_8", "anatase_20", "anatase_300", "mix_20", "rutile_hydrophilic_20", "rutile_hydrophobic_20")
colnames(df.medians) <- names(go.all.bmdl.anatase_8.medians)


# hierarchical 
distance.spearman <- get_dist(df.medians, method = "spearman")
distance.euclidean <- get_dist(df.medians, method = "euclidean")
distance.manhattan <- get_dist(df.medians, method = "manhattan")
distance.canberra <- get_dist(df.medians, method = "canberra")
distance.max <- get_dist(df.medians, method = "maximum")
distance.minkowski <- get_dist(df.medians, method = "minkowski")

clust.spearm.ward2 <- hclust(distance.spearman, method = "ward.D2")
clust.spearm.complete <- hclust(distance.spearman, method = "complete")
clust.spearm.single <- hclust(distance.spearman, method = "single")
clust.spearm.average <- hclust(distance.spearman, method = "average")
clust.spearm.mcquity <- hclust(distance.spearman, method = "mcquitty")
clust.spearm.median <- hclust(distance.spearman, method = "median")
clust.spearm.centroid <- hclust(distance.spearman, method = "centroid")

clust.euclidean.ward2 <- hclust(distance.euclidean, method = "ward.D2")
clust.euclidean.complete <- hclust(distance.euclidean, method = "complete")
clust.euclidean.single <- hclust(distance.euclidean, method = "single")
clust.euclidean.average <- hclust(distance.euclidean, method = "average")
clust.euclidean.mcquity <- hclust(distance.euclidean, method = "mcquitty")
clust.euclidean.median <- hclust(distance.euclidean, method = "median")
clust.euclidean.centroid <- hclust(distance.euclidean, method = "centroid")

clust.maximum.ward2 <- hclust(distance.max, method = "ward.D2")
clust.maximum.complete <- hclust(distance.max, method = "complete")
clust.maximum.single <- hclust(distance.max, method = "single")
clust.maximum.average <- hclust(distance.max, method = "average")
clust.maximum.mcquity <- hclust(distance.max, method = "mcquitty")
clust.maximum.median <- hclust(distance.max, method = "median")
clust.maximum.centroid <- hclust(distance.max, method = "centroid")

clust.manhattan.ward2 <- hclust(distance.manhattan, method = "ward.D2")
clust.manhattan.complete <- hclust(distance.manhattan, method = "complete")
clust.manhattan.single <- hclust(distance.manhattan, method = "single")
clust.manhattan.average <- hclust(distance.manhattan, method = "average")
clust.manhattan.mcquity <- hclust(distance.manhattan, method = "mcquitty")
clust.manhattan.median <- hclust(distance.manhattan, method = "median")
clust.manhattan.centroid <- hclust(distance.manhattan, method = "centroid")

clust.canberra.ward2 <- hclust(distance.canberra, method = "ward.D2")
clust.canberra.complete <- hclust(distance.canberra, method = "complete")
clust.canberra.single <- hclust(distance.canberra, method = "single")
clust.canberra.average <- hclust(distance.canberra, method = "average")
clust.canberra.mcquity <- hclust(distance.canberra, method = "mcquitty")
clust.canberra.median <- hclust(distance.canberra, method = "median")
clust.canberra.centroid <- hclust(distance.canberra, method = "centroid")

clust.minkowski.ward2 <- hclust(distance.minkowski, method = "ward.D2")
clust.minkowski.complete <- hclust(distance.minkowski, method = "complete")
clust.minkowski.single <- hclust(distance.minkowski, method = "single")
clust.minkowski.average <- hclust(distance.minkowski, method = "average")
clust.minkowski.mcquity <- hclust(distance.minkowski, method = "mcquitty")
clust.minkowski.median <- hclust(distance.minkowski, method = "median")
clust.minkowski.centroid <- hclust(distance.minkowski, method = "centroid")

pdf("go_trees.pdf")
plot(clust.spearm.ward2)
plot(clust.spearm.complete)
plot(clust.spearm.single)
plot(clust.spearm.average)
plot(clust.spearm.mcquity)
plot(clust.spearm.median)
plot(clust.spearm.centroid)

plot(clust.euclidean.ward2)
plot(clust.euclidean.complete)
plot(clust.euclidean.single)
plot(clust.euclidean.average)
plot(clust.euclidean.mcquity)
plot(clust.euclidean.median)
plot(clust.euclidean.centroid)

plot(clust.manhattan.ward2)
plot(clust.manhattan.complete)
plot(clust.manhattan.single)
plot(clust.manhattan.average)
plot(clust.manhattan.mcquity)
plot(clust.manhattan.median)
plot(clust.manhattan.centroid)

plot(clust.canberra.ward2)
plot(clust.canberra.complete)
plot(clust.canberra.single)
plot(clust.canberra.average)
plot(clust.canberra.mcquity)
plot(clust.canberra.median)
plot(clust.canberra.centroid)

plot(clust.maximum.ward2)
plot(clust.maximum.complete)
plot(clust.maximum.single)
plot(clust.maximum.average)
plot(clust.maximum.mcquity)
plot(clust.maximum.median)
plot(clust.maximum.centroid)

plot(clust.minkowski.ward2)
plot(clust.minkowski.complete)
plot(clust.minkowski.single)
plot(clust.minkowski.average)
plot(clust.minkowski.mcquity)
plot(clust.minkowski.median)
plot(clust.minkowski.centroid)
dev.off()

df.medians <- matrix(c(unlist(kegg.all.bmdl.anatase_8.medians),
                       unlist(kegg.all.bmdl.anatase_20.medians),
                       unlist(kegg.all.bmdl.anatase_300.medians),
                       unlist(kegg.all.bmdl.mix_20.medians),
                       unlist(kegg.all.bmdl.rutile_hydrophilic_20.medians),
                       unlist(kegg.all.bmdl.rutile_hydrophobic_20.medians)),
                     nrow = 6, ncol = 32)
rownames(df.medians) <- c("anatase_8", "anatase_20", "anatase_300", "mix_20", "rutile_hydrophilic_20", "rutile_hydrophobic_20")
colnames(df.medians) <- names(kegg.all.bmdl.anatase_8.medians)


# hierarchical 
distance.spearman <- get_dist(df.medians, method = "spearman")
distance.euclidean <- get_dist(df.medians, method = "euclidean")
distance.manhattan <- get_dist(df.medians, method = "manhattan")
distance.canberra <- get_dist(df.medians, method = "canberra")
distance.max <- get_dist(df.medians, method = "maximum")
distance.minkowski <- get_dist(df.medians, method = "minkowski")

clust.spearm.ward2 <- hclust(distance.spearman, method = "ward.D2")
clust.spearm.complete <- hclust(distance.spearman, method = "complete")
clust.spearm.single <- hclust(distance.spearman, method = "single")
clust.spearm.average <- hclust(distance.spearman, method = "average")
clust.spearm.mcquity <- hclust(distance.spearman, method = "mcquitty")
clust.spearm.median <- hclust(distance.spearman, method = "median")
clust.spearm.centroid <- hclust(distance.spearman, method = "centroid")

clust.euclidean.ward2 <- hclust(distance.euclidean, method = "ward.D2")
clust.euclidean.complete <- hclust(distance.euclidean, method = "complete")
clust.euclidean.single <- hclust(distance.euclidean, method = "single")
clust.euclidean.average <- hclust(distance.euclidean, method = "average")
clust.euclidean.mcquity <- hclust(distance.euclidean, method = "mcquitty")
clust.euclidean.median <- hclust(distance.euclidean, method = "median")
clust.euclidean.centroid <- hclust(distance.euclidean, method = "centroid")

clust.maximum.ward2 <- hclust(distance.max, method = "ward.D2")
clust.maximum.complete <- hclust(distance.max, method = "complete")
clust.maximum.single <- hclust(distance.max, method = "single")
clust.maximum.average <- hclust(distance.max, method = "average")
clust.maximum.mcquity <- hclust(distance.max, method = "mcquitty")
clust.maximum.median <- hclust(distance.max, method = "median")
clust.maximum.centroid <- hclust(distance.max, method = "centroid")

clust.manhattan.ward2 <- hclust(distance.manhattan, method = "ward.D2")
clust.manhattan.complete <- hclust(distance.manhattan, method = "complete")
clust.manhattan.single <- hclust(distance.manhattan, method = "single")
clust.manhattan.average <- hclust(distance.manhattan, method = "average")
clust.manhattan.mcquity <- hclust(distance.manhattan, method = "mcquitty")
clust.manhattan.median <- hclust(distance.manhattan, method = "median")
clust.manhattan.centroid <- hclust(distance.manhattan, method = "centroid")

clust.canberra.ward2 <- hclust(distance.canberra, method = "ward.D2")
clust.canberra.complete <- hclust(distance.canberra, method = "complete")
clust.canberra.single <- hclust(distance.canberra, method = "single")
clust.canberra.average <- hclust(distance.canberra, method = "average")
clust.canberra.mcquity <- hclust(distance.canberra, method = "mcquitty")
clust.canberra.median <- hclust(distance.canberra, method = "median")
clust.canberra.centroid <- hclust(distance.canberra, method = "centroid")

clust.minkowski.ward2 <- hclust(distance.minkowski, method = "ward.D2")
clust.minkowski.complete <- hclust(distance.minkowski, method = "complete")
clust.minkowski.single <- hclust(distance.minkowski, method = "single")
clust.minkowski.average <- hclust(distance.minkowski, method = "average")
clust.minkowski.mcquity <- hclust(distance.minkowski, method = "mcquitty")
clust.minkowski.median <- hclust(distance.minkowski, method = "median")
clust.minkowski.centroid <- hclust(distance.minkowski, method = "centroid")

pdf("kegg_trees.pdf")
plot(clust.spearm.ward2)
plot(clust.spearm.complete)
plot(clust.spearm.single)
plot(clust.spearm.average)
plot(clust.spearm.mcquity)
plot(clust.spearm.median)
plot(clust.spearm.centroid)

plot(clust.euclidean.ward2)
plot(clust.euclidean.complete)
plot(clust.euclidean.single)
plot(clust.euclidean.average)
plot(clust.euclidean.mcquity)
plot(clust.euclidean.median)
plot(clust.euclidean.centroid)

plot(clust.manhattan.ward2)
plot(clust.manhattan.complete)
plot(clust.manhattan.single)
plot(clust.manhattan.average)
plot(clust.manhattan.mcquity)
plot(clust.manhattan.median)
plot(clust.manhattan.centroid)

plot(clust.canberra.ward2)
plot(clust.canberra.complete)
plot(clust.canberra.single)
plot(clust.canberra.average)
plot(clust.canberra.mcquity)
plot(clust.canberra.median)
plot(clust.canberra.centroid)

plot(clust.maximum.ward2)
plot(clust.maximum.complete)
plot(clust.maximum.single)
plot(clust.maximum.average)
plot(clust.maximum.mcquity)
plot(clust.maximum.median)
plot(clust.maximum.centroid)

plot(clust.minkowski.ward2)
plot(clust.minkowski.complete)
plot(clust.minkowski.single)
plot(clust.minkowski.average)
plot(clust.minkowski.mcquity)
plot(clust.minkowski.median)
plot(clust.minkowski.centroid)
dev.off()


df.medians <- matrix(c(unlist(wikipathways.all.bmdl.anatase_8.medians),
                     unlist(wikipathways.all.bmdl.anatase_20.medians),
                     unlist(wikipathways.all.bmdl.anatase_300.medians),
                     unlist(wikipathways.all.bmdl.mix_20.medians),
                     unlist(wikipathways.all.bmdl.rutile_hydrophilic_20.medians),
                     unlist(wikipathways.all.bmdl.rutile_hydrophobic_20.medians)),
                     nrow = 6, ncol = 10)
rownames(df.medians) <- c("anatase_8", "anatase_20", "anatase_300", "mix_20", "rutile_hydrophilic_20", "rutile_hydrophobic_20")
colnames(df.medians) <- names(wikipathways.all.bmdl.anatase_8.medians)
  

# hierarchical 
distance.spearman <- get_dist(df.medians, method = "spearman")
distance.euclidean <- get_dist(df.medians, method = "euclidean")
distance.manhattan <- get_dist(df.medians, method = "manhattan")
distance.canberra <- get_dist(df.medians, method = "canberra")
distance.max <- get_dist(df.medians, method = "maximum")
distance.minkowski <- get_dist(df.medians, method = "minkowski")

clust.spearm.ward2 <- hclust(distance.spearman, method = "ward.D2")
clust.spearm.complete <- hclust(distance.spearman, method = "complete")
clust.spearm.single <- hclust(distance.spearman, method = "single")
clust.spearm.average <- hclust(distance.spearman, method = "average")
clust.spearm.mcquity <- hclust(distance.spearman, method = "mcquitty")
clust.spearm.median <- hclust(distance.spearman, method = "median")
clust.spearm.centroid <- hclust(distance.spearman, method = "centroid")

clust.euclidean.ward2 <- hclust(distance.euclidean, method = "ward.D2")
clust.euclidean.complete <- hclust(distance.euclidean, method = "complete")
clust.euclidean.single <- hclust(distance.euclidean, method = "single")
clust.euclidean.average <- hclust(distance.euclidean, method = "average")
clust.euclidean.mcquity <- hclust(distance.euclidean, method = "mcquitty")
clust.euclidean.median <- hclust(distance.euclidean, method = "median")
clust.euclidean.centroid <- hclust(distance.euclidean, method = "centroid")

clust.maximum.ward2 <- hclust(distance.max, method = "ward.D2")
clust.maximum.complete <- hclust(distance.max, method = "complete")
clust.maximum.single <- hclust(distance.max, method = "single")
clust.maximum.average <- hclust(distance.max, method = "average")
clust.maximum.mcquity <- hclust(distance.max, method = "mcquitty")
clust.maximum.median <- hclust(distance.max, method = "median")
clust.maximum.centroid <- hclust(distance.max, method = "centroid")

clust.manhattan.ward2 <- hclust(distance.manhattan, method = "ward.D2")
clust.manhattan.complete <- hclust(distance.manhattan, method = "complete")
clust.manhattan.single <- hclust(distance.manhattan, method = "single")
clust.manhattan.average <- hclust(distance.manhattan, method = "average")
clust.manhattan.mcquity <- hclust(distance.manhattan, method = "mcquitty")
clust.manhattan.median <- hclust(distance.manhattan, method = "median")
clust.manhattan.centroid <- hclust(distance.manhattan, method = "centroid")

clust.canberra.ward2 <- hclust(distance.canberra, method = "ward.D2")
clust.canberra.complete <- hclust(distance.canberra, method = "complete")
clust.canberra.single <- hclust(distance.canberra, method = "single")
clust.canberra.average <- hclust(distance.canberra, method = "average")
clust.canberra.mcquity <- hclust(distance.canberra, method = "mcquitty")
clust.canberra.median <- hclust(distance.canberra, method = "median")
clust.canberra.centroid <- hclust(distance.canberra, method = "centroid")

clust.minkowski.ward2 <- hclust(distance.minkowski, method = "ward.D2")
clust.minkowski.complete <- hclust(distance.minkowski, method = "complete")
clust.minkowski.single <- hclust(distance.minkowski, method = "single")
clust.minkowski.average <- hclust(distance.minkowski, method = "average")
clust.minkowski.mcquity <- hclust(distance.minkowski, method = "mcquitty")
clust.minkowski.median <- hclust(distance.minkowski, method = "median")
clust.minkowski.centroid <- hclust(distance.minkowski, method = "centroid")

pdf("wikipathways_trees.pdf")
plot(clust.spearm.ward2)
plot(clust.spearm.complete)
plot(clust.spearm.single)
plot(clust.spearm.average)
plot(clust.spearm.mcquity)
plot(clust.spearm.median)
plot(clust.spearm.centroid)

plot(clust.euclidean.ward2)
plot(clust.euclidean.complete)
plot(clust.euclidean.single)
plot(clust.euclidean.average)
plot(clust.euclidean.mcquity)
plot(clust.euclidean.median)
plot(clust.euclidean.centroid)

plot(clust.manhattan.ward2)
plot(clust.manhattan.complete)
plot(clust.manhattan.single)
plot(clust.manhattan.average)
plot(clust.manhattan.mcquity)
plot(clust.manhattan.median)
plot(clust.manhattan.centroid)

plot(clust.canberra.ward2)
plot(clust.canberra.complete)
plot(clust.canberra.single)
plot(clust.canberra.average)
plot(clust.canberra.mcquity)
plot(clust.canberra.median)
plot(clust.canberra.centroid)

plot(clust.maximum.ward2)
plot(clust.maximum.complete)
plot(clust.maximum.single)
plot(clust.maximum.average)
plot(clust.maximum.mcquity)
plot(clust.maximum.median)
plot(clust.maximum.centroid)

plot(clust.minkowski.ward2)
plot(clust.minkowski.complete)
plot(clust.minkowski.single)
plot(clust.minkowski.average)
plot(clust.minkowski.mcquity)
plot(clust.minkowski.median)
plot(clust.minkowski.centroid)
dev.off()

df.medians <- matrix(c(unlist(kegg.fibrosis.bmdl.anatase_8.medians),
                       unlist(kegg.fibrosis.bmdl.anatase_20.medians),
                       unlist(kegg.fibrosis.bmdl.anatase_300.medians),
                       unlist(kegg.fibrosis.bmdl.mix_20.medians),
                       unlist(kegg.fibrosis.bmdl.rutile_hydrophilic_20.medians),
                       unlist(kegg.fibrosis.bmdl.rutile_hydrophobic_20.medians)),
                     nrow = 6, ncol = 27)
rownames(df.medians) <- c("anatase_8", "anatase_20", "anatase_300", "mix_20", "rutile_hydrophilic_20", "rutile_hydrophobic_20")
colnames(df.medians) <- names(kegg.fibrosis.bmdl.anatase_8.medians)


# hierarchical 
distance.spearman <- get_dist(df.medians, method = "spearman")
distance.euclidean <- get_dist(df.medians, method = "euclidean")
distance.manhattan <- get_dist(df.medians, method = "manhattan")
distance.canberra <- get_dist(df.medians, method = "canberra")
distance.max <- get_dist(df.medians, method = "maximum")
distance.minkowski <- get_dist(df.medians, method = "minkowski")

clust.spearm.ward2 <- hclust(distance.spearman, method = "ward.D2")
clust.spearm.complete <- hclust(distance.spearman, method = "complete")
clust.spearm.single <- hclust(distance.spearman, method = "single")
clust.spearm.average <- hclust(distance.spearman, method = "average")
clust.spearm.mcquity <- hclust(distance.spearman, method = "mcquitty")
clust.spearm.median <- hclust(distance.spearman, method = "median")
clust.spearm.centroid <- hclust(distance.spearman, method = "centroid")

clust.euclidean.ward2 <- hclust(distance.euclidean, method = "ward.D2")
clust.euclidean.complete <- hclust(distance.euclidean, method = "complete")
clust.euclidean.single <- hclust(distance.euclidean, method = "single")
clust.euclidean.average <- hclust(distance.euclidean, method = "average")
clust.euclidean.mcquity <- hclust(distance.euclidean, method = "mcquitty")
clust.euclidean.median <- hclust(distance.euclidean, method = "median")
clust.euclidean.centroid <- hclust(distance.euclidean, method = "centroid")

clust.maximum.ward2 <- hclust(distance.max, method = "ward.D2")
clust.maximum.complete <- hclust(distance.max, method = "complete")
clust.maximum.single <- hclust(distance.max, method = "single")
clust.maximum.average <- hclust(distance.max, method = "average")
clust.maximum.mcquity <- hclust(distance.max, method = "mcquitty")
clust.maximum.median <- hclust(distance.max, method = "median")
clust.maximum.centroid <- hclust(distance.max, method = "centroid")

clust.manhattan.ward2 <- hclust(distance.manhattan, method = "ward.D2")
clust.manhattan.complete <- hclust(distance.manhattan, method = "complete")
clust.manhattan.single <- hclust(distance.manhattan, method = "single")
clust.manhattan.average <- hclust(distance.manhattan, method = "average")
clust.manhattan.mcquity <- hclust(distance.manhattan, method = "mcquitty")
clust.manhattan.median <- hclust(distance.manhattan, method = "median")
clust.manhattan.centroid <- hclust(distance.manhattan, method = "centroid")

clust.canberra.ward2 <- hclust(distance.canberra, method = "ward.D2")
clust.canberra.complete <- hclust(distance.canberra, method = "complete")
clust.canberra.single <- hclust(distance.canberra, method = "single")
clust.canberra.average <- hclust(distance.canberra, method = "average")
clust.canberra.mcquity <- hclust(distance.canberra, method = "mcquitty")
clust.canberra.median <- hclust(distance.canberra, method = "median")
clust.canberra.centroid <- hclust(distance.canberra, method = "centroid")

clust.minkowski.ward2 <- hclust(distance.minkowski, method = "ward.D2")
clust.minkowski.complete <- hclust(distance.minkowski, method = "complete")
clust.minkowski.single <- hclust(distance.minkowski, method = "single")
clust.minkowski.average <- hclust(distance.minkowski, method = "average")
clust.minkowski.mcquity <- hclust(distance.minkowski, method = "mcquitty")
clust.minkowski.median <- hclust(distance.minkowski, method = "median")
clust.minkowski.centroid <- hclust(distance.minkowski, method = "centroid")

pdf("kegg_fibrosis_trees.pdf")
plot(clust.spearm.ward2)
plot(clust.spearm.complete)
plot(clust.spearm.single)
plot(clust.spearm.average)
plot(clust.spearm.mcquity)
plot(clust.spearm.median)
plot(clust.spearm.centroid)

plot(clust.euclidean.ward2)
plot(clust.euclidean.complete)
plot(clust.euclidean.single)
plot(clust.euclidean.average)
plot(clust.euclidean.mcquity)
plot(clust.euclidean.median)
plot(clust.euclidean.centroid)

plot(clust.manhattan.ward2)
plot(clust.manhattan.complete)
plot(clust.manhattan.single)
plot(clust.manhattan.average)
plot(clust.manhattan.mcquity)
plot(clust.manhattan.median)
plot(clust.manhattan.centroid)

plot(clust.canberra.ward2)
plot(clust.canberra.complete)
plot(clust.canberra.single)
plot(clust.canberra.average)
plot(clust.canberra.mcquity)
plot(clust.canberra.median)
plot(clust.canberra.centroid)

plot(clust.maximum.ward2)
plot(clust.maximum.complete)
plot(clust.maximum.single)
plot(clust.maximum.average)
plot(clust.maximum.mcquity)
plot(clust.maximum.median)
plot(clust.maximum.centroid)

plot(clust.minkowski.ward2)
plot(clust.minkowski.complete)
plot(clust.minkowski.single)
plot(clust.minkowski.average)
plot(clust.minkowski.mcquity)
plot(clust.minkowski.median)
plot(clust.minkowski.centroid)
dev.off()

save.image("rah16_bmdl.RData")

rm(list = ls())
load("rah16_bmdl.RData")
