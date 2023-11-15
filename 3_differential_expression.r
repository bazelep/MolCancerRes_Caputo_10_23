#Summary: Differential expression testing for gene-aggregated, QC-filtered, Q3 normalized gene expression values. 
# Data have had batch (scan run) effect removed, while protecting for sample-genotype group variable genotype.
#Conda: karaaym

#Design: 12 samples - 3 groups(4 WT, 4 BRCA1, 4 BRCA2)

library(data.table)
setDTthreads(1)

library(doFuture)
registerDoFuture()
plan("sequential")

library(GeomxTools)

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
#Load normalized data
##
file.name = paste0(rdata.dir,"/","geomxset.genes.filtered.combat.scan,genotype.rds")
file.name
dat = readRDS(file.name)
dat

##
#Convert phenotype columns to factors
##
pData(dat)$segment = factor(sData(dat)$segment)
pData(dat)$genotype = factor(sData(dat)$tags)
pData(dat)$scan = factor(sData(dat)$scan)
pData(dat)$patient = factor(sData(dat)$patient)

##
#Remove negative control gene(s)
##
negative.control.gene.indicator = fData(dat)$CodeClass == "Negative"
print(table(negative.control.gene.indicator))

dat.noctrl = dat[!negative.control.gene.indicator,]
print(dim(dat))
print(dim(dat.noctrl)) 
dat = dat.noctrl


################################################################################################################
#Differential Expression
################################################################################################################

#genotype is the grouping variable - there are 11-12 segments per patient
xtabs(~patient+genotype,data=pData(dat))


##
#Analyze each ROI segment compartment separately.
##
segment.res = foreach(segment = c("stroma","panck")) %do% {
    print(segment)

    #select segment compartment
    segment.indicator = pData(dat)[,"segment"] %in% c(segment)
    
    dat.seg = dat[,segment.indicator]
    print(dim(dat))
    print(dim(dat.seg))    


    #fit each gene
    fits = mixedModelDE(
        dat.seg,
        elt = "combat_scan_genotype_log2_q3_norm",
        modelFormula = ~ genotype + (1|patient),
        groupVar = "genotype",
        nCores = 20,
        multiCore = T,
        pAdjust = "fdr"
    )

    #extract estimates
    genes.res = foreach(gene = colnames(fits)) %do% {

        gene.fits = fits[,gene]

        #pairwise contrast estimates
        gene.lsmeans = gene.fits$lsmeans

        gene.results = data.frame(gene.lsmeans)
        names(gene.results) = c("Estimate","P.Value")
        
        gene.results$Gene.Omnibus.P.Value = gene.fits$anova
        gene.results$Contrast = rownames(gene.lsmeans)
        gene.results$Gene = gene
        
        setDT(gene.results)    
        
        return(gene.results)
    }
    genes.results = rbindlist(genes.res)
    genes.results

    #multiple testing (FDR) correction, separate for each contrast
    genes.results[, Contrast.FDR.P.Value := p.adjust(P.Value,method="fdr"), by = c("Contrast")]

    genes.results[, Segment := segment]

    return(genes.results)
}
segment.results = rbindlist(segment.res)
segment.results


##
#Save results
##
out.file.name = paste0(outputdata.dir,"/",gsub("\\.rds","",basename(file.name)), ".DEList.RandomPatient,FixedGenotype.csv")
out.file.name

#re-order columns and sort
segment.results = segment.results[,mget(c("Segment","Contrast","Gene","Estimate","P.Value","Gene.Omnibus.P.Value","Contrast.FDR.P.Value"))][order(Contrast.FDR.P.Value)]
segment.results

write.csv(segment.results, file = out.file.name, row.names = F)














































































