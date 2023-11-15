#Summary: Remove batch (scan run) effects, protect sample-genotype group effects.
#Conda: karaaym

#Design: 12 samples - 3 groups(4 WT, 4 BRCA1, 4 BRCA2)

library(data.table)
setDTthreads(1)

library(doFuture)
registerDoFuture()
plan("sequential")

library(ggplot2)
library(umap)
library(Rtsne)

library(GeomxTools)
library(sva)



options(width=150)



##
#Setup
##

#Directories
data.dir = "../Data"
graphs.dir = "../Graphs"
outputdata.dir = "../OutputData"
rdata.dir = "../RData"


##
#Load normalized, QC-filtered data
##
file.name = paste0(rdata.dir,"/","geomxset.genes.filtered.rds")
file.name
dat = readRDS(file.name)
dat




################################################################################################################
#Batch Removal
################################################################################################################

##
#Get normalized values.
##
dat.norm = assayDataElement(dat, elt = "log2_q3_norm")
dim(dat.norm)

##
#Get phenotype data.
##
dat.pheno = sData(dat)[,c("scan","genotype","patient")]
dat.pheno

all(rownames(dat.pheno) == colnames(dat.norm))

##
#Apply ComBat to log2-Q3 values, using the primary phenotype (genotype) as the grouping variable and scan as the batch variable.
##
dim(dat.norm)
dim(dat.pheno)

dat.norm.combat = ComBat(dat.norm, batch = dat.pheno$scan, mod = model.matrix(~genotype, data = dat.pheno))

dim(dat.norm)
dim(dat.norm.combat)
head(dat.norm.combat)




################################################################################################################
#Dimensionality Reduction - After Batch Correction
################################################################################################################

##
#Generate PCs on adjusted data.
##
dat.norm.combat.pcs = prcomp(t(dat.norm.combat))$x
dim(dat.norm.combat.pcs)

##
#Add first 4 PCs back to phenotype data
##
pData(dat)[, paste0("ComBat_",c("PC1", "PC2", "PC3", "PC4"))] = dat.norm.combat.pcs[, c(1:4)]
head(pData(dat))



##
#Set random number seed
##
seed = 42
set.seed(seed) 

##
#Update defaults for umap to contain a stable random_state (seed)
##
custom_umap = umap::umap.defaults
custom_umap$random_state = seed

##
#Perform UMAP on log2-Q3 values
##
umap.out = umap(t(dat.norm.combat), config = custom_umap)

##
#Add UMAP dimensions to phenotype data.
##
pData(dat)[, paste0("ComBat_UMAP",c(1:2))] = umap.out$layout[, c(1,2)]


##
#Perform tSNE on log2-Q3 values
##
tsne.out = Rtsne(t(dat.norm.combat), perplexity = ncol(dat.norm.combat) * .15)

##
#Add tSNE dimensions to phenotype data.
##
pData(dat)[, paste0("ComBat_tSNE",c(1:2))] = tsne.out$Y[, c(1,2)]
head(pData(dat))





##
#Plot 1 - PC1 vs. PC2, shape by segment
##
gg1 = ggplot(
        data = sData(dat),
        aes(x = ComBat_PC1, y = ComBat_PC2, color = genotype, shape = segment)
    ) + 
    geom_point(size = 3) + 
    theme_bw() + 
    ggtitle("Segment")

##
#Plot 2 - PC3 vs. PC4, shape by segment
##
gg2 = ggplot(
        data = sData(dat),
        aes(x = ComBat_PC3, y = ComBat_PC4, color = genotype, shape = segment)
    ) + 
    geom_point(size = 3) + 
    theme_bw() + 
    ggtitle("Segment")

##
#Plot 3 - PC1 vs. PC2, shape by scan
##
gg3 = ggplot(
        data = sData(dat),
        aes(x = ComBat_PC1, y = ComBat_PC2, color = genotype, shape = scan)
    ) + 
    geom_point(size = 3) + 
    theme_bw() + 
    ggtitle("Scan")

##
#Plot 4 - PC3 vs. PC4, shape by scan
##
gg4 = ggplot(
        data = sData(dat),
        aes(x = ComBat_PC3, y = ComBat_PC4, color = genotype, shape = scan)
    ) + 
    geom_point(size = 3) +
    theme_bw() + 
    ggtitle("Scan")

##
#Plot 5 - UMAP1 vs. UMAP2, shape by segment
##
gg5 = ggplot(
        data = sData(dat),
        aes(x = ComBat_UMAP1, y = ComBat_UMAP2, color = genotype, shape = segment)
    ) + 
    geom_point(size = 3) + 
    theme_bw() + 
    ggtitle("Segment")

##
#Plot 6 - UMAP1 vs. UMAP2, shape by scan
##
gg6 = ggplot(
        data = sData(dat),
        aes(x = ComBat_UMAP1, y = ComBat_UMAP2, color = genotype, shape = scan)
    ) + 
    geom_point(size = 3) + 
    theme_bw() + 
    ggtitle("Scan")

##
#Plot 7 - tSNE1 vs. tSNE2, shape by segment
##
gg7 = ggplot(
        data = sData(dat),
        aes(x = ComBat_tSNE1, y = ComBat_tSNE2, color = genotype, shape = segment)
    ) + 
    geom_point(size = 3) + 
    theme_bw() + 
    ggtitle("Segment")

##
#Plot 8 - tSNE1 vs. tSNE2, shape by scan
##
gg8 = ggplot(
        data = sData(dat),
        aes(x = ComBat_tSNE1, y = ComBat_tSNE2, color = genotype, shape = scan)
    ) + 
    geom_point(size = 3) + 
    theme_bw() + 
    ggtitle("Scan")

pdf(file = paste0(graphs.dir,"/ComBat_Scan,Genotype_Dimensionality_Reduction.pdf"))
print(gg1)
print(gg2)
print(gg3)
print(gg4)
print(gg5)
print(gg6)
print(gg7)
print(gg8)
dev.off()




##
#Save batch-corrected counts and PCs.
##
out.file.base = paste0(outputdata.dir,"/","Log2.Q3Normalized.ComBatBatchScanCorrected,GenotypeProtected")
write.csv(dat.norm.combat, file = paste0(out.file.base,".Counts.csv"))
write.csv(dat.norm.combat.pcs, file = paste0(out.file.base,".PCs.csv"))




##
#Save combat-adjusted values back to geomx object
##
assayDataElement(dat, elt = "combat_scan_genotype_log2_q3_norm") = dat.norm.combat
dat

out.file.name = paste0(rdata.dir,"/","geomxset.genes.filtered.combat.scan,genotype.rds")
out.file.name
saveRDS(dat, file = out.file.name)




