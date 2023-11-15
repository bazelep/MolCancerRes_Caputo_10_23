#Summary: Quality control and preprocessing of DCC count data
#Conda: karaaym

#Design: 12 samples - 3 groups(4 WT, 4 BRCA1, 4 BRCA2)

library(data.table)
setDTthreads(1)

library(doFuture)
registerDoFuture()
plan("sequential")

library(rlist)
library(readxl)

library(ggplot2)
library(cowplot)
library(pheatmap)
library(umap)
library(Rtsne)

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


#PKC files
pkc.files = paste0(data.dir,"/",c("Hs_R_NGS_WTA_v1.0.pkc"))
pkc.files

#Phenotype sample sheet
sample.sheet.file = paste0(data.dir,"/","SampleAnnotation.xlsx")
sample.sheet.file
sample.sheet.sheet = "Sheet1"
sample.sheet.dcc.col = "Sample_ID"
sample.sheet.protocol.cols = c("aoi", "roi","segment","tags")
sample.sheet.experiment.cols = c("panel")

#DCC files
dcc.dir = paste0(data.dir,"/DCC")
dcc.files = list.files(path=dcc.dir,pattern="dcc",full.names=T)
dcc.files


##
#Initialize GeoMx dataset object
##
geomxset = readNanoStringGeoMxSet(
    dccFiles = dcc.files,
    pkcFiles = pkc.files,
    phenoDataFile = sample.sheet.file,
    phenoDataSheet = sample.sheet.sheet,
    phenoDataDccColName = sample.sheet.dcc.col,
    protocolDataColNames = sample.sheet.protocol.cols,
    experimentDataColNames = sample.sheet.experiment.cols
)
geomxset
head(sData(geomxset))

##
#Tags column represents genotype
##
pData(geomxset)$genotype = sData(geomxset)$tags

##
#Load phenotype groups from sample sheet
##
pheno = data.frame(read_excel(sample.sheet.file, sheet = "Sheet2"))
setDT(pheno)
pheno


##
#Add sample-genotype phenotype to geomx object
##
head(pData(geomxset))

pData(geomxset)$DCC = gsub(".dcc","",rownames(pData(geomxset)))

#make sure all segments listed in phenotype samplesheet are present in geomx object
pheno.keep = pheno[pheno$Sample_ID %in% pData(geomxset)$DCC,]
dim(pheno)
dim(pheno.keep)

#order by segment
pheno.keep = pheno.keep[order(pheno.keep$Sample_ID),]

#check that all segments are the same
all(pheno.keep$Sample_ID == pData(geomxset)$DCC)

#add sample-genotype grouping variable
pData(geomxset)$patient = factor(pheno.keep$patient)
head(pData(geomxset))






##
#Remove ROIs identified as ECM
##
ecm.roi = c("DSP-1001660006157-C-A04", "DSP-1001660006157-C-C03", "DSP-1001660006157-C-D04", "DSP-1001660006157-C-D05")
ecm.roi = paste0(ecm.roi,".dcc")
geomxset = geomxset[,!colnames(geomxset) %in% ecm.roi]

##
#Keep only Epithelial and Stromal segments (exclude "Full ROI" segments)
##
roi.flag = sData(geomxset)[,"segment"] %in% c("panck","stroma")
sum(roi.flag)

geomxset.roi = geomxset[,roi.flag]
dim(geomxset)
dim(geomxset.roi)

geomxset = geomxset.roi


##
#Add pseudo-count
##
geomxset = shiftCountsOne(geomxset, useDALogic = TRUE)










################################################################################################################
#Segment QC
################################################################################################################

##
#Segment QC parameters
##
segment.qc.thresholds =
    list(
        minSegmentReads = 1000, # Minimum number of reads (1000)
        percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
        percentStitched = 80,   # Minimum % of reads stitched (80%)
        percentAligned = 80,    # Minimum % of reads aligned (80%)
        percentSaturation = 50, # Minimum sequencing saturation (50%)
        minNegativeCount = 10,  # Minimum negative control counts (10)
        maxNTCCount = 1000,     # Maximum counts observed in NTC well (1000)
        minNuclei = 100,        # Minimum # of nuclei estimated (100)
        minArea = 5000          # Minimum segment area (5000)
    )

##
#Set Segment QC Flags
##
geomxset = setSegmentQCFlags(geomxset, qcCutoffs = segment.qc.thresholds)
geomxset

##
#Get Segment QC Metrics
##
segment.qc = list.cbind(sData(geomxset))
segment.qc$DCC.File = rownames(segment.qc)
setDT(segment.qc)
segment.qc
colnames(segment.qc) = gsub(" \\(%\\)",".Percent",colnames(segment.qc))

##
#Calculate the negative geometric means for each PKC module (only WTA used in this data)
##
#endogenousSubset - get only probes or targets of endogenous Code Class
#negativeControlSubset - get only the probes with Negative Code Class
#   You can see the Code Class information in the featureData slot.

segment.qc$negative.GeoMeans = esBy(
    negativeControlSubset(geomxset), 
    GROUP = "Module",
    FUN = function(x) { 
        assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
    })

segment.qc$endogenous.GeoMeans = esBy(
    endogenousSubset(geomxset), 
    GROUP = "Module",
    FUN = function(x) { 
        assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
    }) 

quantile(segment.qc$endogenous.GeoMeans, probs = seq(0,1,by=0.1))


##
#Function to plot QC statistics
##
qc_metric_hist = function(
    dat = NULL,
    metric = NULL,
    fill_by = NULL,
    thres = NULL,
    scale_trans = NULL
) {
    
    gg = ggplot(
            data = dat,
            aes_string(x = metric, fill = fill_by)
        ) +
        geom_histogram(bins = 50) +
        geom_vline(xintercept = thres, lty = "dashed", color = "black") +
        theme_bw() + 
        guides(fill = "none") +
        facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
        labs(x = metric, y = "Segments, #", title = metric)

    if(!is.null(scale_trans)) {
        gg = gg + scale_x_continuous(trans = scale_trans)
    }
    return(gg)
}

##
#QC metrics to plot
##
metrics = c(
"Trimmed.Percent",
"Stitched.Percent",
"Aligned.Percent"
)

##
#Metric groupings
##
facets = c("segment","genotype","patient")

##
#Plot histograms of metrics
##
pdf(file = paste0(graphs.dir,"/Segment_QC_metrics_histograms.pdf"))
facet.res = foreach(facet = facets) %do% {
    print(facet)

    metric.res = foreach(metric = metrics) %do% {
        print(metric)
        print(qc_metric_hist(dat = segment.qc, metric = metric, fill_by = facet, thres = 80))
    }
    print(qc_metric_hist(dat = segment.qc, metric = "Saturated.Percent", fill_by = facet, thres = 50))
    print(qc_metric_hist(dat = segment.qc, metric = "DeduplicatedReads", fill_by = facet, thres = 1000, scale_trans = "log10"))
    print(qc_metric_hist(dat = segment.qc, metric = "negative.GeoMeans", fill_by = facet, thres = 2, scale_trans = "log10"))
    print(qc_metric_hist(dat = segment.qc, metric = "endogenous.GeoMeans", fill_by = facet, thres = 2, scale_trans = "log10"))
}
dev.off()

##
#Save QC metrics
##
write.csv(segment.qc, file = paste0(outputdata.dir,"/Segment_QC_metrics.csv"), row.names = F)



##
#No Template Controls (NTCs)
##

##
#Show frequency of each NTC count
##
segment.qc[,.N,by=c("NTC")]


##
#Segment Filtering
##

#DSP-1001660006157-C-A03.dcc
# - low sequencing saturation
# - low negatives but also low endogenous
segment.remove = c("DSP-1001660007412-B-B10.dcc")

dim(geomxset)
geomxset = geomxset[,!colnames(geomxset) %in% segment.remove]
dim(geomxset)



















################################################################################################################
#Probe QC
################################################################################################################

#Expression values for probes that are found to be local Grubb's outliers will set as NA
sum(is.na(exprs(geomxset)))


##
#Probe QC Flags
##
#Generally keep the qcCutoffs parameters unchanged. 
#Set removeLocalOutliers to FALSE if you do not want to remove local outliers

probe.qc.thresholds =
    list(
        minProbeRatio = 0.1, 
        percentFailGrubbs = 20
    )

##
#Set Probe QC Flags
##
geomxset = setBioProbeQCFlags(
    geomxset, 
    qcCutoffs = probe.qc.thresholds,
    removeLocalOutliers = TRUE
)

#local Grubb's outliers are now set as NA in the expression table
sum(is.na(exprs(geomxset)))


##
#Get Probe QC Metric Booleans
##
probe.qc = fData(geomxset)[["QCFlags"]]
colnames(probe.qc)

##
#Summarize probe QC
##
probe.qc.summary = data.frame(
    Passed = sum(rowSums(probe.qc[, -1]) == 0),
    LowProbeRatio = sum(probe.qc[, 1]),
    Global = sum(probe.qc$GlobalGrubbsOutlier),
    Local = sum(rowSums(probe.qc[, -2:-1]) > 0 & !probe.qc$GlobalGrubbsOutlier),
    Local.or.Global = sum(rowSums(probe.qc[, -2:-1]) > 0)
)
head(probe.qc.summary)

##
#Exclude probes with low probe ratios or are global grubbs outliers, but keep local grubbs outliers
##
dim(geomxset)

geomxset = subset(
    geomxset, 
    fData(geomxset)[["QCFlags"]][,c("LowProbeRatio")] == FALSE 
    &
    fData(geomxset)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE
)

dim(geomxset)

#how many local grubbs outliers (i.e. TRUE) before filtering
sum(probe.qc[, -2:-1])

#local outliers remaining after filtering
sum(is.na(exprs(geomxset)))

#summarize count of segments with count of local outlier probes
table(rowSums(is.na(exprs(geomxset))))





################################################################################################################
#Aggregate Probe-level counts to Genes
################################################################################################################

#check how many unique gene targets the object has
length(unique(featureData(geomxset)[["TargetName"]]))

##
#Calculate gene-level counts
##
geomxset.genes = aggregateCounts(geomxset, FUN = ngeoMean)
dim(geomxset)
dim(geomxset.genes)

exprs(geomxset)[1:5, 1:2]
exprs(geomxset.genes)[1:5, 1:2]


################################################################################################################
#Limit of Quantification (LOQ) per Segment
################################################################################################################

#In addition to Segment and Probe QC, we also determine the limit of quantification (LOQ) per segment. 
#The LOQ is calculated based on the distribution of negative control probes and is intended to approximate the 
# quantifiable limit of gene expression per segment. Please note that this process is more stable in larger segments. 
#Likewise, the LOQ may not be as accurately reflective of true signal detection rates in segments with low negative 
# probe counts (ex: <2).

pkc.modules = gsub(".pkc","",annotation(geomxset.genes))
pkc.modules

##
#Define LOQ SD threshold and minimum value
##
loq.sd.cutoff = 2 #We typically use 2 geometric standard deviations (n=2) above the geometric mean as the LOQ
loq.min = 2 #We also recommend that a minimum LOQ of 2 be used if the LOQ calculated in a segment is below this threshold.


##
#Get geometric mean values for negative and endogenous probes for each module (only WTA in this case)
##
geomxset.genes.pheno = pData(geomxset.genes)
head(geomxset.genes.pheno)


##
#Calculate LOQ per module tested
##

#create empty df to store results
loq.df = data.frame(row.names = colnames(geomxset.genes))

module.res = foreach(module = pkc.modules) %do% {

    vars = paste0(c("NegGeoMean_", "NegGeoSD_"),module)
    
    if(all(vars[1:2] %in% colnames(geomxset.genes.pheno))) {
    
        segment.neg.geo.means = geomxset.genes.pheno[, vars[1]]
    
        segment.neg.geo.sds = geomxset.genes.pheno[, vars[2]]
    
        boundaries = segment.neg.geo.sds ^ loq.sd.cutoff
    
        loq = segment.neg.geo.means * boundaries
    
        max.values = pmax(loq.min, loq)
    
        loq.df[, module] = max.values
    }
}
loq.df

#add back to phenotype data
pData(geomxset.genes)$LOQ = loq.df


##
#Calculate LOQ flags
##
#After determining the limit of quantification (LOQ) per segment, we recommend filtering out either segments 
# and/or genes with abnormally low signal. 
#Filtering is an important step to focus on the true biological data of interest.
#We determine the number of genes detected in each segment across the dataset.

#empty matrix to hold LOQ flags
loq.mat = c()

module.res = foreach(module = pkc.modules) %do% {

    module.indicator = fData(geomxset.genes)$Module == module

    module.mat = t(
        esApply(
            geomxset.genes[module.indicator, ], 
            MARGIN = 1,
            FUN = function(x) {
                x > loq.df[, module]
            }
        )
    )  

    loq.mat = rbind(loq.mat, module.mat)
}
dim(loq.mat)

##
#Ensure ordering since this is stored outside of the geomxSet
##
loq.mat = loq.mat[fData(geomxset.genes)$TargetName, ]





################################################################################################################
#Segment-level detection rate (i.e. number of genes detected per segment)
################################################################################################################

##
#Save detection rate information to pheno data
##
pData(geomxset.genes)$Detected.Gene.Count = colSums(loq.mat, na.rm = TRUE)

pData(geomxset.genes)$Detected.Gene.Proportion = pData(geomxset.genes)$Detected.Gene.Count / nrow(geomxset.genes)

##
#Determine detection thresholds to plot
##
pData(geomxset.genes)$Detected.Gene.Threshold = cut(
    pData(geomxset.genes)$Detected.Gene.Proportion,
    breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
    labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%")
)

##
#Plot stacked bar plot of different cut points (1%, 5%, 10%, 15%)
##
gg = ggplot(
        data = sData(geomxset.genes),
        aes(x = Detected.Gene.Threshold)
    ) +
    geom_bar(aes(fill = `scan name`)) +
    geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
    theme_bw() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
        x = "Gene Detection Rate",
        y = "Segments, #",
        fill = "Scan"
    )
ggsave(file = paste0(graphs.dir,"/Segment_Gene_Detection_Rate.pdf"), plot = gg)


##
#Tabulate cut percent genes detected at 1, 5, 10, 15
##
xtabs(~Detected.Gene.Threshold+genotype,data = sData(geomxset.genes))

##
#Apply gene detection rate segment-level filtering
##

#We choose to remove segments with less than 1% of the genes detected. 
# Generally, 5-10% detection is a reasonable segment filtering threshold. 
# However, based on the experimental design (e.g. segment types, size, nuclei) and 
# tissue characteristics (e.g. type, age), these guidelines may require adjustment.

sum(pData(geomxset.genes)$Detected.Gene.Proportion >= .1)
sum(pData(geomxset.genes)$Detected.Gene.Proportion >= .05)
sum(pData(geomxset.genes)$Detected.Gene.Proportion >= .01)

#Here we are doing exploratory analysis, so use lenient criteria, and keep segments with >1% gene detection rates.

dim(geomxset.genes)
geomxset.genes = geomxset.genes[, pData(geomxset.genes)$Detected.Gene.Proportion >= .01]
dim(geomxset.genes)






################################################################################################################
#Gene-level detection rate (i.e. number of segments detected per gene)
################################################################################################################

##
#Update LOQ matrix in case segments filtered out above
##
loq.mat = loq.mat[, colnames(geomxset.genes)]

##
#For each gene, calculate count and percentage of detected segments
##
fData(geomxset.genes)$Detected.Segment.Count = rowSums(loq.mat, na.rm = TRUE)

table(fData(geomxset.genes)$Detected.Segment.Count)

fData(geomxset.genes)$Detected.Segment.Proportion = fData(geomxset.genes)$Detected.Segment.Count / nrow(pData(geomxset.genes))

t(t(table(fData(geomxset.genes)$Detected.Segment.Proportion)))

head(fData(geomxset.genes))

##
#Calculate detection rates to plot:
##
detection.rates = data.frame(Detected.Segment.Proportion = c(0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5))

detection.rates$Detected.Segment.Proportion.Gene.Count = unlist(
    lapply(
        detection.rates$Detected.Segment.Proportion,
        function(x) {
            sum(fData(geomxset.genes)$Detected.Segment.Proportion >= x)
        }
    )
)
detection.rates$Detected.Segment.Proportion.Gene.Proportion = detection.rates$Detected.Segment.Proportion.Gene.Count / nrow(fData(geomxset.genes))
detection.rates

rownames(detection.rates) = detection.rates$Detected.Segment.Proportion


##
#Plot selected segment detection rates against corresponding gene detection rates
##
gg = ggplot(
        data = detection.rates, 
        aes(
            x = as.factor(Detected.Segment.Proportion), 
            y = Detected.Segment.Proportion.Gene.Proportion, 
            fill = Detected.Segment.Proportion.Gene.Proportion
        )
    ) +
    geom_bar(stat = "identity") +
    geom_text(
        aes(label = formatC(Detected.Segment.Proportion.Gene.Count, format = "d", big.mark = ",")),
        vjust = 1.6, 
        color = "black", 
        size = 4
    ) +
    scale_fill_gradient2(
        low = "orange2", 
        mid = "lightblue",
        high = "dodgerblue3", 
        midpoint = 0.65,
        limits = c(0,1),
        labels = scales::percent
    ) +
    theme_bw() +
    scale_y_continuous(
        labels = scales::percent, 
        limits = c(0,1),
        expand = expansion(mult = c(0, 0))
    ) +
    labs(
        x = "% of Segments Detected",
        y = "% of Genes Detected for given % of Segments Detected (% of Panel > LOQ)"
    )
ggsave(file = paste0(graphs.dir,"/Genes_Detected_by_Segments_Detected.pdf"), plot = gg)






##
#Filter genes based on LOQ detection
##

#We typically set a cutoff ranging from 5-20% based on the biological diversity 
# of our dataset. But if we know that a key gene is represented in only a small number of segments 
# due to biological diversity, we may select a different cutoff or keep the target gene 
# by manually selecting it for inclusion in the data object.

#Here we are doing exploratory analysis, where we do not necessarily know all the target genes.
# Thus we use lenient criteria and subset to target genes detected in at least 1% of the segments.
# Also we manually include the negative control probe, for downstream use.

negative.probes = unique(subset(fData(geomxset.genes), CodeClass == "Negative")$TargetName)
negative.probes

unfiltered.gene.count = nrow(geomxset.genes)

geomxset.genes.filtered = geomxset.genes[
    fData(geomxset.genes)$Detected.Segment.Proportion >= 0.01 
    | fData(geomxset.genes)$TargetName %in% negative.probes
    , ]

filtered.gene.count = nrow(geomxset.genes.filtered)

unfiltered.gene.count - filtered.gene.count












################################################################################################################
#Q3 Normalization
################################################################################################################

##
#First compare normalized values to background signal (negative probe genes)
##

##
#Add simpler variable name for scan.
##
pData(geomxset.genes.filtered)[["scan"]] = sData(geomxset.genes.filtered)$`scan name`

##
#Get raw expression data (raw read counts if 1 probe per gene, geometric mean if >1 probe per gene)
##
expr = exprs(geomxset.genes.filtered)
segment.names = colnames(expr)

##
#Plot Q3 values of all genes vs. negative probe gene values, grouping by different phenotype groups
##
groupings = c("segment","genotype","scan")

pdf(paste0(graphs.dir,"/Q3_Normalization.pdf"))
grouping.res = foreach(grouping = groupings) %do% {
    print(grouping)

    Q3 = data.frame(
        Segment = segment.names,
        Annotation = sData(geomxset.genes.filtered)[, grouping],
        Q3 = unlist(apply(expr, 2, quantile, 0.75, na.rm = TRUE)),
        NegProbe = expr[negative.probes, ]
    )
    setDT(Q3)
    print(head(Q3))

    Q3.long = melt(Q3, measure.vars = c("Q3","NegProbe"), variable.name = "Statistic", value.name = "Value")
    print(head(Q3.long))

    #plot 1 - histogram of values, faceted by grouping
    gg1 = ggplot(
        data = Q3.long,
        aes(x = Value, fill = Statistic)
    ) +
    geom_histogram(bins = 40) + 
    theme_bw() +
    scale_x_continuous(trans = "log2") +
    facet_wrap(~Annotation, nrow = 1) + 
    scale_fill_brewer(palette = 3, type = "qual") +
    labs(x = "Counts", y = "Segments, #")

    #plot 2 - scatter plot of values, colored by grouping
    gg2 = ggplot(
        data = Q3,
        aes(x = NegProbe, y = Q3, color = Annotation)
    ) +
    geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
    geom_point() + 
    guides(color = "none") + 
    theme_bw() +
    scale_x_continuous(trans = "log2") + 
    scale_y_continuous(trans = "log2") +
    theme(aspect.ratio = 1) +
    labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")

    #plot 3 - scatter plot of negative control values vs. ratio of Q3 to negative control values
    gg3 = ggplot(
        data = Q3,
        aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)
    ) +
    geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
    geom_point() + 
    theme_bw() +
    scale_x_continuous(trans = "log2") + 
    scale_y_continuous(trans = "log2") +
    theme(aspect.ratio = 1) +
    labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")

    gg.btm_row = plot_grid(gg2, gg3, nrow = 1, labels = c("B", ""),rel_widths = c(0.43,0.57))

    print(plot_grid(gg1, gg.btm_row, ncol = 1, labels = c("A", "")))
}
dev.off()






##
#Perform Q3 normalization (75th percentile) - save to new expression slot "q_norm"
##

geomxset.genes.filtered = normalize(
    geomxset.genes.filtered,
    norm_method = "quant",
    desiredQuantile = .75,
    toElt = "q3_norm"
)


##
#Perform background normalization using negative control values, compare to Q3 normalization
##
geomxset.genes.filtered = normalize(
    geomxset.genes.filtered,
    norm_method = "neg", 
    fromElt = "exprs",
    toElt = "neg_norm"
)
    
##
#Visualize the first 10 segments with each normalization method
##

outlier.size = .1
segment.text.size = 3

##
#Plot 1 - raw values
##
raw = data.frame(exprs(geomxset.genes.filtered))
raw$Gene = rownames(raw)
setDT(raw)
raw.long = melt(raw, id.vars = c("Gene"), variable.name = "Segment", value.name = "Expression")
raw.long[, Segment := gsub("\\.","-",gsub("\\.dcc","",Segment))]

gg1 = ggplot(
        data = raw.long,
        aes(x = log(Expression), y = Segment)
    ) +
    geom_boxplot(outlier.size = outlier.size) +
    ggtitle("Raw Counts") + 
    ylab("Segment") + 
    xlab("Log (Raw Counts)") +
    theme(axis.text.y = element_text(size = segment.text.size))
ggsave(file = paste0(graphs.dir,"/Raw_Counts.pdf"), plot = gg1)

##
#Plot 2 - Q3 values
##
q3.norm = data.frame(assayDataElement(geomxset.genes.filtered, elt = "q3_norm"))
q3.norm$Gene = rownames(q3.norm)
setDT(q3.norm)
q3.norm.long = melt(q3.norm, id.vars = c("Gene"), variable.name = "Segment", value.name = "Expression")
q3.norm.long[, Segment := gsub("\\.","-",gsub("\\.dcc","",Segment))]

gg2 = ggplot(
        data = q3.norm.long,
        aes(x = log(Expression), y = Segment)
    ) +
    geom_boxplot(outlier.size = outlier.size) +
    ggtitle("Q3-Normalized Counts") + 
    ylab("Segment") + 
    xlab("Log (Q3-Normalized Counts)") +
    theme(axis.text.y = element_text(size = segment.text.size))
ggsave(file = paste0(graphs.dir,"/Q3_Normalized_Counts.pdf"), plot = gg2)

##
#Plot 3 - background normalized values
##
neg.norm = data.frame(assayDataElement(geomxset.genes.filtered, elt = "neg_norm"))
neg.norm$Gene = rownames(neg.norm)
setDT(neg.norm)
neg.norm.long = melt(neg.norm, id.vars = c("Gene"), variable.name = "Segment", value.name = "Expression")
neg.norm.long[, Segment := gsub("\\.","-",gsub("\\.dcc","",Segment))]

gg3 = ggplot(
        data = neg.norm.long,
        aes(x = log(Expression), y = Segment)
    ) +
    geom_boxplot(outlier.size = outlier.size) +
    ggtitle("Background-Normalized Counts") + 
    ylab("Segment") + 
    xlab("Log (Background-Normalized Counts)") +
    theme(axis.text.y = element_text(size = segment.text.size))
ggsave(file = paste0(graphs.dir,"/Background_Normalized_Counts.pdf"), plot = gg3)








################################################################################################################
#Dimensionality Reduction
################################################################################################################

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
#Create log2-transformed, Q3-normalized values
##
assayDataElement(object = geomxset.genes.filtered, elt = "log2_q3_norm") = assayDataApply(geomxset.genes.filtered, 2, FUN = log, base = 2, elt = "q3_norm")
norm = t(assayDataElement(geomxset.genes.filtered, elt = "log2_q3_norm"))
dim(norm)

##
#Perform PCA on log2-Q3 values
##
pca = prcomp(norm)$x

##
#Add first 4 PCs back to phenotype data
##
pData(geomxset.genes.filtered)[, c("PC1", "PC2", "PC3", "PC4")] = pca[, c(1:4)]


##
#Perform UMAP on log2-Q3 values
##
umap.out = umap(norm, config = custom_umap)

##
#Add UMAP dimensions to phenotype data.
##
pData(geomxset.genes.filtered)[, c("UMAP1", "UMAP2")] = umap.out$layout[, c(1,2)]


##
#Perform tSNE on log2-Q3 values
##
tsne.out = Rtsne(norm, perplexity = nrow(norm) * .15)

##
#Add tSNE dimensions to phenotype data.
##
pData(geomxset.genes.filtered)[, c("tSNE1", "tSNE2")] = tsne.out$Y[, c(1,2)]



##
#Plot 1 - PC1 vs. PC2, shape by segment
##
gg1 = ggplot(
        data = sData(geomxset.genes.filtered),
        aes(x = PC1, y = PC2, color = genotype, shape = segment)
    ) + 
    geom_point(size = 3) + 
    theme_bw() + 
    ggtitle("Segment")

##
#Plot 2 - PC3 vs. PC4, shape by segment
##
gg2 = ggplot(
        data = sData(geomxset.genes.filtered),
        aes(x = PC3, y = PC4, color = genotype, shape = segment)
    ) + 
    geom_point(size = 3) + 
    theme_bw() + 
    ggtitle("Segment")

##
#Plot 3 - PC1 vs. PC2, shape by scan
##
gg3 = ggplot(
        data = sData(geomxset.genes.filtered),
        aes(x = PC1, y = PC2, color = genotype, shape = scan)
    ) + 
    geom_point(size = 3) + 
    theme_bw() + 
    ggtitle("Scan")

##
#Plot 4 - PC3 vs. PC4, shape by scan
##
gg4 = ggplot(
        data = sData(geomxset.genes.filtered),
        aes(x = PC3, y = PC4, color = genotype, shape = scan)
    ) + 
    geom_point(size = 3) +
    theme_bw() + 
    ggtitle("Scan")

##
#Plot 5 - UMAP1 vs. UMAP2, shape by segment
##
gg5 = ggplot(
        data = sData(geomxset.genes.filtered),
        aes(x = UMAP1, y = UMAP2, color = genotype, shape = segment)
    ) + 
    geom_point(size = 3) + 
    theme_bw() + 
    ggtitle("Segment")

##
#Plot 6 - UMAP1 vs. UMAP2, shape by scan
##
gg6 = ggplot(
        data = sData(geomxset.genes.filtered),
        aes(x = UMAP1, y = UMAP2, color = genotype, shape = scan)
    ) + 
    geom_point(size = 3) + 
    theme_bw() + 
    ggtitle("Scan")

##
#Plot 7 - tSNE1 vs. tSNE2, shape by segment
##
gg7 = ggplot(
        data = sData(geomxset.genes.filtered),
        aes(x = tSNE1, y = tSNE2, color = genotype, shape = segment)
    ) + 
    geom_point(size = 3) + 
    theme_bw() + 
    ggtitle("Segment")

##
#Plot 8 - tSNE1 vs. tSNE2, shape by scan
##
gg8 = ggplot(
        data = sData(geomxset.genes.filtered),
        aes(x = tSNE1, y = tSNE2, color = genotype, shape = scan)
    ) + 
    geom_point(size = 3) + 
    theme_bw() + 
    ggtitle("Scan")

pdf(file = paste0(graphs.dir,"/Dimensionality_Reduction.pdf"))
print(gg1)
print(gg2)
print(gg3)
print(gg4)
print(gg5)
print(gg6)
print(gg7)
print(gg8)
dev.off()
















################################################################################################################
#Heatmap - top 1/3 genes by Coefficient of Variation
################################################################################################################

##
#Coefficient of Variantion (CV) function
##
calc_CV = function(x) {sd(x) / mean(x)}

##
#Calculate CV for log2-Q3 values
##
CV.dat = assayDataApply(geomxset.genes.filtered, elt = "log2_q3_norm", MARGIN = 1, calc_CV)

##
#Show the highest CV genes and their CV values
##
sort(CV.dat, decreasing = TRUE)[1:5]

##
#Identify genes in the top 3rd of the CV values
##
top.genes = names(CV.dat)[CV.dat > quantile(CV.dat, 0.8)]
length(top.genes)

anno = sData(geomxset.genes.filtered)[, c("genotype", "scan", "segment")]
head(anno)

top.genes.expr = assayDataElement(geomxset.genes.filtered[top.genes, ], elt = "log2_q3_norm")


##
#Plot heatmap
##
pdf(paste0(graphs.dir,"/Top Third of Genes by CV - Heatmap.pdf"))

pheatmap(
    mat = top.genes.expr,
    scale = "row", 
    show_rownames = FALSE, 
    show_colnames = FALSE,
    border_color = NA,
    clustering_method = "average",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    breaks = seq(-3, 3, 0.05),
    color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
    annotation_col = anno
)
dev.off()



################################################################################################################
#Save Log2-transformed, Q3-Normalized, filtered, aggregated, gene values
################################################################################################################

file.name = paste0(rdata.dir, "/","geomxset.genes.filtered.rds")
file.name

saveRDS(geomxset.genes.filtered, file = file.name)








