#Summary: Correlation analysis for known receptor-ligand gene pairs.
# Data are gene-aggregated, QC-filtered, Q3 normalized gene expression values. 
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
#Remove negative control gene(s)
##
negative.control.gene.indicator = fData(dat)$CodeClass == "Negative"
print(table(negative.control.gene.indicator))

dat.noctrl = dat[!negative.control.gene.indicator,]
print(dim(dat))
print(dim(dat.noctrl)) 
dat = dat.noctrl









################################################################################################################
#Receptor-Ligand Correlation Analysis
################################################################################################################


##
#Ramilowski 2015 data (https://fantom.gsc.riken.jp/5/suppl/Ramilowski_et_al_2015/)
##

ligand.receptor.expression.file = paste0(data.dir,"/ExpressionLigRec.txt")
ligand.receptor.pairs.file = paste0(data.dir,"/PairsLigRec.txt")

##
#Load expression data
##
ligand.receptor.expression = fread(ligand.receptor.expression.file, header = T)
ligand.receptor.expression

##
#Keep genes with >0 expression in Mammary Epithelial tissue
##
mammary.epithelial.expressed = ligand.receptor.expression[`Mammary Epithelial` > 0, .(ApprovedSymbol,Type,`Mammary Epithelial`)]
dim(ligand.receptor.expression)
dim(mammary.epithelial.expressed)
setnames(mammary.epithelial.expressed, c("Gene","Type","Expression"))
mammary.epithelial.expressed


##
#Load pair data.
##
ligand.receptor.pairs = fread(ligand.receptor.pairs.file, header = T)
ligand.receptor.pairs = ligand.receptor.pairs[,.(Ligand.ApprovedSymbol,Receptor.ApprovedSymbol)]
setnames(ligand.receptor.pairs, c("ligand","receptor"))

ligand.receptor.pairs[, Pair.Name := paste0(ligand,".",receptor)]


##
#Keep only pairs with both ligand and receptor expressed in mammary epithelial data 
##
ligand.receptor.pairs.long = melt(ligand.receptor.pairs, id.vars = c("Pair.Name"), variable.name = "Type", value.name = "Gene")

pairs = merge(ligand.receptor.pairs.long, mammary.epithelial.expressed, by = c("Type","Gene"))
pairs


expressed.pairs = dcast(pairs, Pair.Name ~ Type, value.var = "Expression")[!is.na(ligand) & !is.na(receptor),]
expressed.pairs[, c("Ligand","Receptor") := tstrsplit(Pair.Name,"\\.")]
expressed.pairs


##
#Regress out scan (batch) effect via linear regression.
##
dat.pheno = sData(dat)[,c("segment","genotype","patient","scan")]
dim(dat.pheno)

dat.norm = t(assayDataElement(dat, elt = "log2_q3_norm"))
dim(dat.norm)

dat.norm.scan.resid = resid(lm(dat.norm ~ dat.pheno$scan))


pheno.expr = data.frame(merge(dat.pheno,dat.norm.scan.resid,by="row.names"),row.names=1)
pheno.expr[1:5,1:10]




##
#Correlation within tissue (segment compartment)
##
options(future.globals.maxSize = 10000 * 1024^2)
options(future.seed=TRUE)
plan(multisession, workers = 50)

pair.res = foreach(pair.i = c(1:nrow(expressed.pairs))) %dopar% {

    pair = expressed.pairs[pair.i,]

    segment.res = foreach(seg = c("panck","stroma")) %dopar% {

        genotype.res = foreach(geno = c("WT","BRCA1","BRCA2")) %dopar% {     
        
            segs = subset(pheno.expr, segment == seg & genotype == geno)
        
            if( (pair$Ligand %in% colnames(segs)) & (pair$Receptor %in% colnames(segs)) ){
        
                #create table of pair's expression values
                LRG = data.frame(
                    "L" = segs[,pair$Ligand],
                    "R" = segs[,pair$Receptor],
                    segs[,c("segment","scan","genotype","patient")]
                )
        
                #fit correlation
                fit = cor.test(LRG$L, LRG$R, alternative = c("two.sided"), method = c("pearson"),exact=F)
        
                results = data.frame(list(
                    "Ligand" = pair$Ligand,
                    "Receptor" = pair$Receptor,
                    "Nominal.P.Value" = fit$p.value,
                    "Estimate" = fit$estimate)
                )
                results$Segment = seg
                results$Genotype = geno 
                setDT(results)

                return(results)
            } 
        }
        return(rbindlist(genotype.res))
    }
    return(rbindlist(segment.res))
}
within.results = rbindlist(pair.res)
within.results

#multiple testing (FDR) correction, separate for each comparison
within.results[, FDR.P.Value := p.adjust(Nominal.P.Value,method="fdr"), by = c("Segment")]
within.results = within.results[order(FDR.P.Value),]





##
#Correlation across tissue (segment compartment)
##
pair.res = foreach(pair.i = c(1:nrow(expressed.pairs))) %dopar% {

    pair = expressed.pairs[pair.i,]

    genotype.res = foreach(geno = c("WT","BRCA1","BRCA2")) %dopar% {     
        
        segs = subset(pheno.expr, genotype == geno)
        
        if( (pair$Ligand %in% colnames(segs)) & (pair$Receptor %in% colnames(segs)) ){
        
            #create table of pair's expression values
            LRG = data.frame(
                "L" = segs[,pair$Ligand],
                "R" = segs[,pair$Receptor],
                segs[,c("segment","scan","genotype","patient")]
            )

            #order by segment        
            LRG = LRG[order(LRG$segment),]

            #average multiple segments per patient
            means = list()
            patient.res = foreach(pat = as.character(unique(LRG$patient))) %do% {
            
                segment.res = foreach(seg = as.character(unique(LRG$segment))) %do% {
            
                    type.res = foreach(type = c("L","R")) %do% {
            
                        means[[seg]][[type]][[pat]] = mean(subset(LRG, patient == pat & segment == seg)[,type])
                    }
                }
            }
            
            #test both ligand in epithelium with receptor in stroma, and vice versa
            Epi.Lig.Str.Rec.means = data.frame(unlist(means[["panck"]][["L"]]),unlist(means[["stroma"]][["R"]]))
            Epi.Rec.Str.Lig.means = data.frame(unlist(means[["panck"]][["R"]]),unlist(means[["stroma"]][["L"]]))

            #spearman correlation estimates
            Epi.Lig.Str.Rec.means.fit = cor.test(Epi.Lig.Str.Rec.means[,1],Epi.Lig.Str.Rec.means[,2],method="spearman",exact=F)
            Epi.Rec.Str.Lig.means.fit = cor.test(Epi.Rec.Str.Lig.means[,1],Epi.Rec.Str.Lig.means[,2],method="spearman",exact=F)
            

            Epi.Lig.Str.Rec.results = data.frame(list(
                "Genotype" = geno,
                "Comparison"="Epithelium.Ligand.vs.Stroma.Receptor",
                "Group1"=paste0("Epithelium.Ligand.",pair$Ligand),
                "Group2"=paste0("Stroma.Receptor.",pair$Receptor),
                "Nominal.P.Value"=Epi.Lig.Str.Rec.means.fit$p.value,
                "Estimate"=Epi.Lig.Str.Rec.means.fit$estimate)
            )

            Epi.Rec.Str.Lig.results = data.frame(list(
                "Genotype"=geno,
                "Comparison"="Epithelium.Receptor.vs.Stroma.Ligand",
                "Group1"=paste0("Epithelium.Receptor.",pair$Receptor),
                "Group2"=paste0("Stroma.Ligand.",pair$Ligand),
                "Nominal.P.Value"=Epi.Rec.Str.Lig.means.fit$p.value,
                "Estimate"=Epi.Rec.Str.Lig.means.fit$estimate)
            )

            results = rbind(Epi.Lig.Str.Rec.results,Epi.Rec.Str.Lig.results)
            setDT(results)
            return(results)
        }
    }
    return(rbindlist(genotype.res))
}
across.results = rbindlist(pair.res)
across.results

#multiple testing (FDR) correction, separate for each comparison
across.results[, FDR.P.Value := p.adjust(Nominal.P.Value,method="fdr"), by = c("Comparison")]
across.results = across.results[order(FDR.P.Value),]



##
#Save results
##
out.file.base = paste0(outputdata.dir,"/",gsub("\\.rds","",basename(file.name)), ".ReceptorLigandPairs.Correlation")

within.file.name = paste0(out.file.base,".WithinTissue.csv")
within.file.name
write.csv(within.results, file = within.file.name)

across.file.name = paste0(out.file.base,".AcrossTissue.csv")
across.file.name
write.csv(across.results, file = across.file.name)


































