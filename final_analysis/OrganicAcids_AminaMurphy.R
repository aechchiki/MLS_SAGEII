### presets

## define working directory

setwd('../input_files/') 

# import packages for downstream analysis; or install them if necessary
# source("https://bioconductor.org/biocLite.R") # or: source("http://bioconductor.org/biocLite.R")
# biocLite()
# biocLite("edgeR")
library(edgeR)
# biocLite("rtracklayer")
library(rtracklayer)
# biocLite("pheatmap") 
library(pheatmap)

## read files

# read gtf
gtf <- import.gff('Pseudomonas_protegens_S5_genome.gtf')
gtf <- as.data.frame(gtf)
head(gtf) 
# note: idk what's the deal. I have an issue with the column format of the function "import.gff" that the assistants suggested to use. 

# I prepared a file extracted from the original gtf, called 'Pseudo.csv', containing three columns: $1 starting position on genome; $2 end position on the genome; $ locus tag
# then, calculate the length of the locus ($end-$start)
gtf_info=read.csv('Pseudo.csv', h=T)
head(gtf_info); dim(gtf_info) 

locus_length=gtf_info[["end"]]-gtf_info[["start"]]
gtf_info["length"]=locus_length

# define genes we're interested in
TCA=read.csv('TCA.csv', h=T) # TCA metabolism (preset color: blue)
DeG=read.csv('DeG.csv', h=T) # dicarboxylate and glyoxylate metabolism (preset color: green)
TRANS=read.csv('TRANS.csv', h=T) # other TRANS metabolism, grep in GenDB (preset color: orange)
UPTK=read.csv('UPTK.csv', h=T) # uptake (preset color: red)


### pre-analysis

# calculate gene length
#geneLength <- aggregate(gtf$width, list(gtf[["group"]]$gene_id), max) 
#colnames(geneLength) <- c("gene_name", "length")

# defne function for RPKM 
rpkm <- function(counts, lengths) {
 	 rate <- counts / lengths
	return(rate / sum(counts) * 1e9) 
}
# import abundances from kallisto
path_to_abundances <- '../input_files/abundances/' 
files <- dir(path_to_abundances, pattern=".tsv$")
files <- paste0(path_to_abundances, files)

# define sample tag
samples <- paste0(rep(c('LM','SA','WL','WR'), each = 4), rep(1:4,4))

# read kallisto files 
transcripts <- readDGE(files,
                       columns = c(1,4),
                       group = rep(c('LM','SA','WL','WR'), each = 4),
                       labels = samples)
# visualize 
head(transcripts[["counts"]])

# store counts 
tr_counts <- transcripts$counts
head(tr_counts)

# normalize and filter counts 
pseudoCount=log2(tr_counts)
cpms=cpm(tr_counts)
keep=rowSums(cpms>1) >=4
countsFilter=tr_counts[keep,]
dim(tr_counts); dim(countsFilter)
pseudocountsFilter=log2(countsFilter+1)

# compute RPKMs 
geneLengths <- gtf_info[gtf_info[["locus"]] %in% (rownames(tr_counts)),]
RPKMs <- rpkm(tr_counts, geneLengths$length); head(RPKMs)

## overall data visualisation
setwd("../output_files/")
# pseudocount histogram all conditions
par(mfrow=c(4,4))
hist(as.data.frame(pseudocountsFilter)$LM1, main=" Filtered pseudocounts LM1", xlim=(range(0,20)), ylim=(range(0,1500)), ylab="log2(count + 1)", xlab="WR1", col="white")
hist(as.data.frame(pseudocountsFilter)$LM2, main=" Filtered pseudocounts LM2", xlim=(range(0,20)), ylim=(range(0,1500)), ylab="log2(count + 1)", xlab="WR1", col="white")
hist(as.data.frame(pseudocountsFilter)$LM3, main=" Filtered pseudocounts LM3", xlim=(range(0,20)), ylim=(range(0,1500)), ylab="log2(count + 1)", xlab="WR1", col="white")
hist(as.data.frame(pseudocountsFilter)$LM4, main=" Filtered pseudocounts LM4", xlim=(range(0,20)), ylim=(range(0,1500)), ylab="log2(count + 1)", xlab="WR1", col="white")
hist(as.data.frame(pseudocountsFilter)$SA1, main=" Filtered pseudocounts SA1", xlim=(range(0,20)), ylim=(range(0,1500)), ylab="log2(count + 1)", xlab="WR1", col="khaki1")
hist(as.data.frame(pseudocountsFilter)$SA2, main=" Filtered pseudocounts SA2", xlim=(range(0,20)), ylim=(range(0,1500)), ylab="log2(count + 1)", xlab="WR1", col="khaki1")
hist(as.data.frame(pseudocountsFilter)$SA3, main=" Filtered pseudocounts SA3", xlim=(range(0,20)), ylim=(range(0,1500)), ylab="log2(count + 1)", xlab="WR1", col="khaki1")
hist(as.data.frame(pseudocountsFilter)$SA4, main=" Filtered pseudocounts SA4", xlim=(range(0,20)), ylim=(range(0,1500)), ylab="log2(count + 1)", xlab="WR1", col="khaki1")
hist(as.data.frame(pseudocountsFilter)$WR1, main=" Filtered pseudocounts WR1", xlim=(range(0,20)), ylim=(range(0,1500)), ylab="log2(count + 1)", xlab="WR1", col="khaki4")
hist(as.data.frame(pseudocountsFilter)$WR2, main=" Filtered pseudocounts WR2", xlim=(range(0,20)), ylim=(range(0,1500)), ylab="log2(count + 1)", xlab="WR1", col="khaki4")
hist(as.data.frame(pseudocountsFilter)$WR3, main=" Filtered pseudocounts WR3", xlim=(range(0,20)), ylim=(range(0,1500)), ylab="log2(count + 1)", xlab="WR1", col="khaki4")
hist(as.data.frame(pseudocountsFilter)$WR4, main=" Filtered pseudocounts WR4", xlim=(range(0,20)), ylim=(range(0,1500)), ylab="log2(count + 1)", xlab="WR1", col="khaki4")
hist(as.data.frame(pseudocountsFilter)$WL1, main=" Filtered pseudocounts WL1", xlim=(range(0,20)), ylim=(range(0,1500)), ylab="log2(count + 1)", xlab="WR1", col="darkgreen")
hist(as.data.frame(pseudocountsFilter)$WL2, main=" Filtered pseudocounts WL2", xlim=(range(0,20)), ylim=(range(0,1500)), ylab="log2(count + 1)", xlab="WR1", col="darkgreen")
hist(as.data.frame(pseudocountsFilter)$WL3, main=" Filtered pseudocounts WL3", xlim=(range(0,20)), ylim=(range(0,1500)), ylab="log2(count + 1)", xlab="WR1", col="darkgreen")
hist(as.data.frame(pseudocountsFilter)$WL4, main=" Filtered pseudocounts WL4", xlim=(range(0,20)), ylim=(range(0,1500)), ylab="log2(count + 1)", xlab="WR1", col="darkgreen")

# pseudocount boxpots all conditions
par(mfrow=c(1,1))
boxplot(pseudocountsFilter, col=c(rep(c("white", "khaki2", "darkgreen", "khaki4"), each=4)), ylab="Filtered log2(pseudocounts)", xlab="Samples", main="Filtered pseudocounts per sample")
summary(as.data.frame(pseudocountsFilter))

# plot MDS all conditions
conditions=c(rep("LM", 4), rep("SA", 4), rep("WL", 4), rep("WR", 4))
d=DGEList(counts=countsFilter, group=conditions)
d=calcNormFactors(d)
plotMDS(d, labels=colnames(cpms), col=c("black", "khaki2", "darkgreen", "khaki4")[factor(conditions)], main="MDS plot on samples")


## pairwise comparisons

## LM-WR

# select counts
tr_counts_LMWR=tr_counts[, grep("^LM|WR", colnames(tr_counts))]
head(tr_counts_LMWR)

# compute counts per million
tr_cpms_LMWR <- cpm(tr_counts_LMWR)
head(tr_cpms_LMWR); length(rownames(tr_cpms_LMWR))

# filter lowly expressed transcripts
keep_LMWR <- rowSums(tr_cpms_LMWR > 1) >= 4
tr_counts_clean_LMWR <- tr_counts_LMWR[keep_LMWR,]
head(tr_counts_clean_LMWR); length(rownames(tr_counts_clean_LMWR))

# save filtered counts  
dt_LMWR <- DGEList(counts = tr_counts_clean_LMWR, 
		   group = rep(c('LM','WR'), 
		   each=4))
# visualize 
head(dt_LMWR[["counts"]])

# compute normalization factors 
dt_LMWR <- calcNormFactors(dt_LMWR)
# normalize according to dispersion 
dt_LMWR <- estimateCommonDisp(dt_LMWR)
dt_LMWR <- estimateTagwiseDisp(dt_LMWR)

# test conditions
de.tr_LMWR <- exactTest(dt_LMWR, pair = c("LM","WR"))

# correction for multiple testing and sort de transcripts
tT.transcripts_LMWR <- topTags(de.tr_LMWR, n = nrow(dt_LMWR))
head(tT.transcripts_LMWR[["table"]])

# exract differentially expressed transcripts
table.tr_LMWR=tT.transcripts_LMWR$table
head(table.tr_LMWR); length(table.tr_LMWR[["logFC"]])

# select counts with FDR <0.1
top.det_LMWR <- table.tr_LMWR[table.tr_LMWR$FDR < 0.01 , ]
head(top.det_LMWR); length(top.det_LMWR[["logFC"]])

#compute RPKMs
RPKMs_LMWR <- rpkm(tr_counts_LMWR, geneLengths$length)

# prepare data to plot
TCA_LMWR=top.det_LMWR[row.names(top.det_LMWR) %in% TCA$transcript_id, ]
length(TCA_LMWR$logFC)
DeG_LMWR=top.det_LMWR[row.names(top.det_LMWR) %in% DeG$transcript_id, ]
length(DeG_LMWR$logFC)
TRANS_LMWR=top.det_LMWR[row.names(top.det_LMWR) %in% TRANS$transcript_id, ]
length(TRANS_LMWR$logFC)
UPTK_LMWR=top.det_LMWR[row.names(top.det_LMWR) %in% UPTK$transcript_id, ]
length(UPTK_LMWR$logFC)

# export logFC to csv
write.csv(TCA_LMWR, "TCA_LMWR.csv")
write.csv(DeG_LMWR, "DeG_LMWR.csv")
write.csv(TRANS_LMWR, "TRANS_LMWR.csv")
write.csv(UPTK_LMWR, "UPTK_LMWR.csv")

# create MA plot 
plot(x=table.tr_LMWR$logCPM, y=table.tr_LMWR$logFC, pch=19, main="MA-plot LM-WR", xlab="Average log(CPM)", ylab="log2(FC) WR-LM")
points(x=TCA_LMWR[["logCPM"]], y=TCA_LMWR[["logFC"]], pch=19, col="blue")
points(x=DeG_LMWR[["logCPM"]], y=DeG_LMWR[["logFC"]], pch=19, col="green")
points(x=TRANS_LMWR[["logCPM"]], y=TRANS_LMWR[["logFC"]], pch=19, col="orange")
points(x=UPTK_LMWR[["logCPM"]], y=UPTK_LMWR[["logFC"]], pch=19, col="red")

# create volcano plot
plot(y=-log10(table.tr_LMWR[["PValue"]]), x=table.tr_LMWR$logFC, pch=19, main=" Volcano plot LM-WR", xlab="log2(FC) WR-LM", ylab="log10(p-value)")
points(y=-log10(TCA_LMWR[["PValue"]]), x=TCA_LMWR[["logFC"]], pch=19, col="blue")
points(y=-log10(DeG_LMWR[["PValue"]]), x=DeG_LMWR[["logFC"]], pch=19, col="green")
points(y=-log10(TRANS_LMWR[["PValue"]]), x=TRANS_LMWR[["logFC"]], pch=19, col="orange")
points(y=-log10(UPTK_LMWR[["PValue"]]), x=UPTK_LMWR[["logFC"]], pch=19, col="red")

# heatmaps
# TCA metabolism
TCA_RPKMs_LMWR=RPKMs_LMWR[row.names(RPKMs_LMWR) %in% TCA$transcript_id, ]
pheatmap(log2(TCA_RPKMs_LMWR), main = 'Heatmap of TCA metabolism genes: LM vs WR', scale="column", cluster_rows=FALSE)
# DeG metaolism
DeG_RPKMs_LMWR=RPKMs_LMWR[row.names(RPKMs_LMWR) %in% DeG$transcript_id, ]
pheatmap(log2(DeG_RPKMs_LMWR), main = 'Heatmap of GeD metabolism genes: LM vs WR', scale="column", cluster_rows=FALSE)
# transport & uptake systems
UPTK_RPKMs_LMWR=RPKMs_LMWR[row.names(RPKMs_LMWR) %in% UPTK$transcript_id, ]
pheatmap(log2(UPTK_RPKMs_LMWR), main = 'Heatmap of uptake genes: LM vs WR', scale="column", cluster_rows=FALSE)
# other organic acids metabolism
TRANS_RPKMs_LMWR=RPKMs_LMWR[row.names(RPKMs_LMWR) %in% TRANS$transcript_id, ]
pheatmap(log2(TRANS_RPKMs_LMWR), main = 'Heatmap of transport genes: LM vs WR', scale="column", cluster_rows=FALSE)

## LM-WL

# select counts
tr_counts_LMWL=tr_counts[, grep("^LM|WL", colnames(tr_counts))]
head(tr_counts_LMWL)

# compute counts per million
tr_cpms_LMWL <- cpm(tr_counts_LMWL)
head(tr_cpms_LMWL); length(rownames(tr_cpms_LMWL))

# filter lowly expressed transcripts
keep_LMWL <- rowSums(tr_cpms_LMWL > 1) >= 4
tr_counts_clean_LMWL <- tr_counts_LMWL[keep_LMWL,]
head(tr_counts_clean_LMWL); length(rownames(tr_counts_clean_LMWL))

# save filtered counts  
dt_LMWL <- DGEList(counts = tr_counts_clean_LMWL, 
		   group = rep(c('LM','WL'), 
		   each=4))
# visualize 
head(dt_LMWL[["counts"]])

# compute normalization factors 
dt_LMWL <- calcNormFactors(dt_LMWL)
# normalize according to dispersion 
dt_LMWL <- estimateCommonDisp(dt_LMWL)
dt_LMWL <- estimateTagwiseDisp(dt_LMWL)

# test conditions
de.tr_LMWL <- exactTest(dt_LMWL, pair = c("LM","WL"))

# correction for multiple testing and sort de transcripts
tT.transcripts_LMWL <- topTags(de.tr_LMWL, n = nrow(dt_LMWL))
head(tT.transcripts_LMWL[["table"]])

# exract differentially expressed transcripts
table.tr_LMWL=tT.transcripts_LMWL$table
head(table.tr_LMWL); length(table.tr_LMWL[["logFC"]])

# select counts with FDR <0.1
top.det_LMWL <- table.tr_LMWL[table.tr_LMWL$FDR < 0.01 , ]
head(top.det_LMWL); length(top.det_LMWL[["logFC"]])

#compute RPKMs
RPKMs_LMWL <- rpkm(tr_counts_LMWL, geneLengths$length)

# prepare data to plot
TCA_LMWL=top.det_LMWL[row.names(top.det_LMWL) %in% TCA$transcript_id, ]
length(TCA_LMWL$logFC)
DeG_LMWL=top.det_LMWL[row.names(top.det_LMWL) %in% DeG$transcript_id, ]
length(DeG_LMWL$logFC)
TRANS_LMWL=top.det_LMWL[row.names(top.det_LMWL) %in% TRANS$transcript_id, ]
length(TRANS_LMWL$logFC)
UPTK_LMWL=top.det_LMWL[row.names(top.det_LMWL) %in% UPTK$transcript_id, ]
length(UPTK_LMWL$logFC)

# export logFC to csv
write.csv(TCA_LMWL, "TCA_LMWL.csv")
write.csv(DeG_LMWL, "DeG_LMWL.csv")
write.csv(TRANS_LMWL, "TRANS_LMWL.csv")
write.csv(UPTK_LMWL, "UPTK_LMWL.csv")

# create MA plot 
plot(x=table.tr_LMWL$logCPM, y=table.tr_LMWL$logFC, pch=19, main="MA-plot LM-WL", xlab="Average log(CPM)", ylab="log2(FC) WL-LM" )
points(x=TCA_LMWL[["logCPM"]], y=TCA_LMWL[["logFC"]], pch=19, col="blue")
points(x=DeG_LMWL[["logCPM"]], y=DeG_LMWL[["logFC"]], pch=19, col="green")
points(x=TRANS_LMWL[["logCPM"]], y=TRANS_LMWL[["logFC"]], pch=19, col="orange")
points(x=UPTK_LMWL[["logCPM"]], y=UPTK_LMWL[["logFC"]], pch=19, col="red")

# create volcano plot
plot(y=-log10(table.tr_LMWL[["PValue"]]), x=table.tr_LMWL$logFC, pch=19, main=" Volcano plot LM-WL", xlab="log2(FC) WL-LM", ylab="log10(p-value)")
points(y=-log10(TCA_LMWL[["PValue"]]), x=TCA_LMWL[["logFC"]], pch=19, col="blue")
points(y=-log10(DeG_LMWL[["PValue"]]), x=DeG_LMWL[["logFC"]], pch=19, col="green")
points(y=-log10(TRANS_LMWL[["PValue"]]), x=TRANS_LMWL[["logFC"]], pch=19, col="orange")
points(y=-log10(UPTK_LMWL[["PValue"]]), x=UPTK_LMWL[["logFC"]], pch=19, col="red")

# heatmaps
# TCA metabolism
TCA_RPKMs_LMWL=RPKMs_LMWL[row.names(RPKMs_LMWL) %in% TCA$transcript_id, ]
pheatmap(log2(TCA_RPKMs_LMWL), main = 'Heatmap of TCA metabolism genes: LM vs WL', scale="column", cluster_rows=FALSE)
# DeG metaolism
DeG_RPKMs_LMWL=RPKMs_LMWL[row.names(RPKMs_LMWL) %in% DeG$transcript_id, ]
pheatmap(log2(DeG_RPKMs_LMWL), main = 'Heatmap of G&D metabolism genes: LM vs WL', scale="column", cluster_rows=FALSE)
# transport & uptake systems
UPTK_RPKMs_LMWL=RPKMs_LMWL[row.names(RPKMs_LMWL) %in% UPTK$transcript_id, ]
pheatmap(log2(UPTK_RPKMs_LMWL), main = 'Heatmap of uptake genes: LM vs WL', scale="column", cluster_rows=FALSE)
# other organic acids metabolism
TRANS_RPKMs_LMWL=RPKMs_LMWL[row.names(RPKMs_LMWL) %in% TRANS$transcript_id, ]
pheatmap(log2(TRANS_RPKMs_LMWL), main = 'Heatmap of transport genes: LM vs WL', scale="column", cluster_rows=FALSE)

## WR-WL

# select counts
tr_counts_WRWL=tr_counts[, grep("^WR|WL", colnames(tr_counts))]
head(tr_counts_WRWL)

# compute counts per million
tr_cpms_WRWL <- cpm(tr_counts_WRWL)
head(tr_cpms_WRWL); length(rownames(tr_cpms_WRWL))

# filter lowly expressed transcripts
keep_WRWL <- rowSums(tr_cpms_WRWL > 1) >= 4
tr_counts_clean_WRWL <- tr_counts_WRWL[keep_WRWL,]
head(tr_counts_clean_WRWL); length(rownames(tr_counts_clean_WRWL))

# save filtered counts  
dt_WRWL <- DGEList(counts = tr_counts_clean_WRWL, 
		   group = rep(c('WR','WL'), 
		   each=4))
# visualize 
head(dt_WRWL[["counts"]])

# compute normalization factors 
dt_WRWL <- calcNormFactors(dt_WRWL)
# normalize according to dispersion 
dt_WRWL <- estimateCommonDisp(dt_WRWL)
dt_WRWL <- estimateTagwiseDisp(dt_WRWL)

# test conditions
de.tr_WRWL <- exactTest(dt_WRWL, pair = c("WR","WL"))

# correction for multiple testing and sort de transcripts
tT.transcripts_WRWL <- topTags(de.tr_WRWL, n = nrow(dt_WRWL))
head(tT.transcripts_WRWL[["table"]])

# exract differentially expressed transcripts
table.tr_WRWL=tT.transcripts_WRWL$table
head(table.tr_WRWL); length(table.tr_WRWL[["logFC"]])

# select counts with FDR <0.1
top.det_WRWL <- table.tr_WRWL[table.tr_WRWL$FDR < 0.01 , ]
head(top.det_WRWL); length(top.det_WRWL[["logFC"]])

#compute RPKMs
RPKMs_WRWL <- rpkm(tr_counts_WRWL, geneLengths$length)

# prepare data to plot
TCA_WRWL=top.det_WRWL[row.names(top.det_WRWL) %in% TCA$transcript_id, ]
length(TCA_WRWL$logFC)
DeG_WRWL=top.det_WRWL[row.names(top.det_WRWL) %in% DeG$transcript_id, ]
length(DeG_WRWL$logFC)
TRANS_WRWL=top.det_WRWL[row.names(top.det_WRWL) %in% TRANS$transcript_id, ]
length(TRANS_WRWL$logFC)
UPTK_WRWL=top.det_WRWL[row.names(top.det_WRWL) %in% UPTK$transcript_id, ]
length(UPTK_WRWL$logFC)

# export logFC to csv
write.csv(TCA_WRWL, "TCA_WRWL.csv")
write.csv(DeG_WRWL, "DeG_WRWL.csv")
write.csv(TRANS_WRWL, "TRANS_WRWL.csv")
write.csv(UPTK_WRWL, "UPTK_WRWL.csv")

# create MA plot 
plot(x=table.tr_WRWL$logCPM, y=table.tr_WRWL$logFC, pch=19, main="MA-plot WR-WL", xlab="Average log(CPM)", ylab="log2(FC) WL-WR" )
points(x=TCA_WRWL[["logCPM"]], y=TCA_WRWL[["logFC"]], pch=19, col="blue")
points(x=DeG_WRWL[["logCPM"]], y=DeG_WRWL[["logFC"]], pch=19, col="green")
points(x=TRANS_WRWL[["logCPM"]], y=TRANS_WRWL[["logFC"]], pch=19, col="orange")
points(x=UPTK_WRWL[["logCPM"]], y=UPTK_WRWL[["logFC"]], pch=19, col="red")

# create volcano plot
plot(y=-log10(table.tr_WRWL[["PValue"]]), x=table.tr_WRWL$logFC, pch=19, main=" Volcano plot WR-WL", xlab="log2(FC) WL-WR", ylab="log10(p-value)")
points(y=-log10(TCA_WRWL[["PValue"]]), x=TCA_WRWL[["logFC"]], pch=19, col="blue")
points(y=-log10(DeG_WRWL[["PValue"]]), x=DeG_WRWL[["logFC"]], pch=19, col="green")
points(y=-log10(TRANS_WRWL[["PValue"]]), x=TRANS_WRWL[["logFC"]], pch=19, col="orange")
points(y=-log10(UPTK_WRWL[["PValue"]]), x=UPTK_WRWL[["logFC"]], pch=19, col="red")

# heatmaps
# TCA metabolism
TCA_RPKMs_WRWL=RPKMs_WRWL[row.names(RPKMs_WRWL) %in% TCA$transcript_id, ]
pheatmap(log2(TCA_RPKMs_WRWL), main = 'Heatmap of TCA metabolism genes: WR vs WL', scale="column", cluster_rows=FALSE)
# DeG metaolism
DeG_RPKMs_WRWL=RPKMs_WRWL[row.names(RPKMs_WRWL) %in% DeG$transcript_id, ]
pheatmap(log2(DeG_RPKMs_WRWL), main = 'Heatmap of G&D metabolism genes: WR vs WL', scale="column", cluster_rows=FALSE)
# transport & uptake systems
UPTK_RPKMs_WRWL=RPKMs_WRWL[row.names(RPKMs_WRWL) %in% UPTK$transcript_id, ]
pheatmap(log2(UPTK_RPKMs_WRWL), main = 'Heatmap of uptake genes: WR vs WL', scale="column", cluster_rows=FALSE)
# other organic acids metabolism
TRANS_RPKMs_WRWL=RPKMs_WRWL[row.names(RPKMs_WRWL) %in% TRANS$transcript_id, ]
pheatmap(log2(TRANS_RPKMs_WRWL), main = 'Heatmap of transport genes: WR vs WL', scale="column", cluster_rows=FALSE)

## WR-SA

# select counts
tr_counts_WRSA=tr_counts[, grep("^WR|SA", colnames(tr_counts))]
head(tr_counts_WRSA)

# compute counts per million
tr_cpms_WRSA <- cpm(tr_counts_WRSA)
head(tr_cpms_WRSA); length(rownames(tr_cpms_WRSA))

# filter lowly expressed transcripts
keep_WRSA <- rowSums(tr_cpms_WRSA > 1) >= 4
tr_counts_clean_WRSA <- tr_counts_WRSA[keep_WRSA,]
head(tr_counts_clean_WRSA); length(rownames(tr_counts_clean_WRSA))

# save filtered counts  
dt_WRSA <- DGEList(counts = tr_counts_clean_WRSA, 
		   group = rep(c('WR','SA'), 
		   each=4))
# visualize 
head(dt_WRSA[["counts"]])

# compute normalization factors 
dt_WRSA <- calcNormFactors(dt_WRSA)
# normalize according to dispersion 
dt_WRSA <- estimateCommonDisp(dt_WRSA)
dt_WRSA <- estimateTagwiseDisp(dt_WRSA)

# test conditions
de.tr_WRSA <- exactTest(dt_WRSA, pair = c("WR","SA"))

# correction for multiple testing and sort de transcripts
tT.transcripts_WRSA <- topTags(de.tr_WRSA, n = nrow(dt_WRSA))
head(tT.transcripts_WRSA[["table"]])

# exract differentially expressed transcripts
table.tr_WRSA=tT.transcripts_WRSA$table
head(table.tr_WRSA); length(table.tr_WRSA[["logFC"]])

# select counts with FDR <0.1
top.det_WRSA <- table.tr_WRSA[table.tr_WRSA$FDR < 0.01 , ]
head(top.det_WRSA); length(top.det_WRSA[["logFC"]])

#compute RPKMs
RPKMs_WRSA <- rpkm(tr_counts_WRSA, geneLengths$length)

# prepare data to plot
TCA_WRSA=top.det_WRSA[row.names(top.det_WRSA) %in% TCA$transcript_id, ]
length(TCA_WRSA$logFC)
DeG_WRSA=top.det_WRSA[row.names(top.det_WRSA) %in% DeG$transcript_id, ]
length(DeG_WRSA$logFC)
TRANS_WRSA=top.det_WRSA[row.names(top.det_WRSA) %in% TRANS$transcript_id, ]
length(TRANS_WRSA$logFC)
UPTK_WRSA=top.det_WRSA[row.names(top.det_WRSA) %in% UPTK$transcript_id, ]
length(UPTK_WRSA$logFC)

# export logFC to csv
write.csv(TCA_WRSA, "TCA_WRSA.csv")
write.csv(DeG_WRSA, "DeG_WRSA.csv")
write.csv(TRANS_WRSA, "TRANS_WRSA.csv")
write.csv(UPTK_WRSA, "UPTK_WRSA.csv")

# create MA plot 
plot(x=table.tr_WRSA$logCPM, y=table.tr_WRSA$logFC, pch=19, main="MA-plot WR-SA", xlab="Average log(CPM)", ylab="log2(FC) SA-WR" )
points(x=TCA_WRSA[["logCPM"]], y=TCA_WRSA[["logFC"]], pch=19, col="blue")
points(x=DeG_WRSA[["logCPM"]], y=DeG_WRSA[["logFC"]], pch=19, col="green")
points(x=TRANS_WRSA[["logCPM"]], y=TRANS_WRSA[["logFC"]], pch=19, col="orange")
points(x=UPTK_WRSA[["logCPM"]], y=UPTK_WRSA[["logFC"]], pch=19, col="red")

# create volcano plot
plot(y=-log10(table.tr_WRSA[["PValue"]]), x=table.tr_WRSA$logFC, pch=19, main=" Volcano plot WR-SA", xlab="log2(FC) SA-WR", ylab="log10(p-value)")
points(y=-log10(TCA_WRSA[["PValue"]]), x=TCA_WRSA[["logFC"]], pch=19, col="blue")
points(y=-log10(DeG_WRSA[["PValue"]]), x=DeG_WRSA[["logFC"]], pch=19, col="green")
points(y=-log10(TRANS_WRSA[["PValue"]]), x=TRANS_WRSA[["logFC"]], pch=19, col="orange")
points(y=-log10(UPTK_WRSA[["PValue"]]), x=UPTK_WRSA[["logFC"]], pch=19, col="red")

# heatmaps
# TCA metabolism
TCA_RPKMs_WRSA=RPKMs_WRSA[row.names(RPKMs_WRSA) %in% TCA$transcript_id, ]
pheatmap(log2(TCA_RPKMs_WRSA), main = 'Heatmap of TCA metabolism genes: WR vs SA', scale="column", cluster_rows=FALSE)
# DeG metaolism
DeG_RPKMs_WRSA=RPKMs_WRSA[row.names(RPKMs_WRSA) %in% DeG$transcript_id, ]
pheatmap(log2(DeG_RPKMs_WRSA), main = 'Heatmap of G&D metabolism genes: WR vs SA', scale="column", cluster_rows=FALSE)
# transport & uptake systems
UPTK_RPKMs_WRSA=RPKMs_WRSA[row.names(RPKMs_WRSA) %in% UPTK$transcript_id, ]
pheatmap(log2(UPTK_RPKMs_WRSA), main = 'Heatmap of uptake genes: WR vs SA', scale="column", cluster_rows=FALSE)
# other organic acids metabolism
TRANS_RPKMs_WRSA=RPKMs_WRSA[row.names(RPKMs_WRSA) %in% TRANS$transcript_id, ]
pheatmap(log2(TRANS_RPKMs_WRSA), main = 'Heatmap of transport genes: WR vs SA', scale="column", cluster_rows=FALSE)

## WL-SA

# select counts
tr_counts_WLSA=tr_counts[, grep("^WL|SA", colnames(tr_counts))]
head(tr_counts_WLSA)

# compute counts per million
tr_cpms_WLSA <- cpm(tr_counts_WLSA)
head(tr_cpms_WLSA); length(rownames(tr_cpms_WLSA))

# filter lowly expressed transcripts
keep_WLSA <- rowSums(tr_cpms_WLSA > 1) >= 4
tr_counts_clean_WLSA <- tr_counts_WLSA[keep_WLSA,]
head(tr_counts_clean_WLSA); length(rownames(tr_counts_clean_WLSA))

# save filtered counts  
dt_WLSA <- DGEList(counts = tr_counts_clean_WLSA, 
		   group = rep(c('WL','SA'), 
		   each=4))
# visualize 
head(dt_WLSA[["counts"]])

# compute normalization factors 
dt_WLSA <- calcNormFactors(dt_WLSA)
# normalize according to dispersion 
dt_WLSA <- estimateCommonDisp(dt_WLSA)
dt_WLSA <- estimateTagwiseDisp(dt_WLSA)

# test conditions
de.tr_WLSA <- exactTest(dt_WLSA, pair = c("WL","SA"))

# correction for multiple testing and sort de transcripts
tT.transcripts_WLSA <- topTags(de.tr_WLSA, n = nrow(dt_WLSA))
head(tT.transcripts_WLSA[["table"]])

# exract differentially expressed transcripts
table.tr_WLSA=tT.transcripts_WLSA$table
head(table.tr_WLSA); length(table.tr_WLSA[["logFC"]])

# select counts with FDR <0.1
top.det_WLSA <- table.tr_WLSA[table.tr_WLSA$FDR < 0.01 , ]
head(top.det_WLSA); length(top.det_WLSA[["logFC"]])

#compute RPKMs
RPKMs_WLSA <- rpkm(tr_counts_WLSA, geneLengths$length)

# prepare data to plot
TCA_WLSA=top.det_WLSA[row.names(top.det_WLSA) %in% TCA$transcript_id, ]
length(TCA_WLSA$logFC)
DeG_WLSA=top.det_WLSA[row.names(top.det_WLSA) %in% DeG$transcript_id, ]
length(DeG_WLSA$logFC)
TRANS_WLSA=top.det_WLSA[row.names(top.det_WLSA) %in% TRANS$transcript_id, ]
length(TRANS_WLSA$logFC)
UPTK_WLSA=top.det_WLSA[row.names(top.det_WLSA) %in% UPTK$transcript_id, ]
length(UPTK_WLSA$logFC)

# export logFC to csv
write.csv(TCA_WLSA, "TCA_WLSA.csv")
write.csv(DeG_WLSA, "DeG_WLSA.csv")
write.csv(TRANS_WLSA, "TRANS_WLSA.csv")
write.csv(UPTK_WLSA, "UPTK_WLSA.csv")

# create MA plot 
plot(x=table.tr_WLSA$logCPM, y=table.tr_WLSA$logFC, pch=19, main="MA-plot WL-SA", xlab="Average log(CPM)", ylab="log2(FC) SA-WL" )
points(x=TCA_WLSA[["logCPM"]], y=TCA_WLSA[["logFC"]], pch=19, col="blue")
points(x=DeG_WLSA[["logCPM"]], y=DeG_WLSA[["logFC"]], pch=19, col="green")
points(x=TRANS_WLSA[["logCPM"]], y=TRANS_WLSA[["logFC"]], pch=19, col="orange")
points(x=UPTK_WLSA[["logCPM"]], y=UPTK_WLSA[["logFC"]], pch=19, col="red")

# create volcano plot
plot(y=-log10(table.tr_WLSA[["PValue"]]), x=table.tr_WLSA$logFC, pch=19, main=" Volcano plot WL-SA", xlab="log2(FC) SA-WL", ylab="log10(p-value)")
points(y=-log10(TCA_WLSA[["PValue"]]), x=TCA_WLSA[["logFC"]], pch=19, col="blue")
points(y=-log10(DeG_WLSA[["PValue"]]), x=DeG_WLSA[["logFC"]], pch=19, col="green")
points(y=-log10(TRANS_WLSA[["PValue"]]), x=TRANS_WLSA[["logFC"]], pch=19, col="orange")
points(y=-log10(UPTK_WLSA[["PValue"]]), x=UPTK_WLSA[["logFC"]], pch=19, col="red")

# heatmaps
# TCA metabolism
TCA_RPKMs_WLSA=RPKMs_WLSA[row.names(RPKMs_WLSA) %in% TCA$transcript_id, ]
pheatmap(log2(TCA_RPKMs_WLSA), main = 'Heatmap of TCA metabolism genes: WL vs SA', scale="column", cluster_rows=FALSE)
# DeG metaolism
DeG_RPKMs_WLSA=RPKMs_WLSA[row.names(RPKMs_WLSA) %in% DeG$transcript_id, ]
pheatmap(log2(DeG_RPKMs_WLSA), main = 'Heatmap of G&D metabolism genes: WL vs SA', scale="column", cluster_rows=FALSE)
# transport & uptake systems
UPTK_RPKMs_WLSA=RPKMs_WLSA[row.names(RPKMs_WLSA) %in% UPTK$transcript_id, ]
pheatmap(log2(UPTK_RPKMs_WLSA), main = 'Heatmap of uptake genes: WL vs SA', scale="column", cluster_rows=FALSE)
# other organic acids metabolism
TRANS_RPKMs_WLSA=RPKMs_WLSA[row.names(RPKMs_WLSA) %in% TRANS$transcript_id, ]
pheatmap(log2(TRANS_RPKMs_WLSA), main = 'Heatmap of transport genes: WL vs SA', scale="column", cluster_rows=FALSE)

## LM-SA

# select counts
tr_counts_LMSA=tr_counts[, grep("^LM|SA", colnames(tr_counts))]
head(tr_counts_LMSA)

# compute counts per million
tr_cpms_LMSA <- cpm(tr_counts_LMSA)
head(tr_cpms_LMSA); length(rownames(tr_cpms_LMSA))

# filter lowly expressed transcripts
keep_LMSA <- rowSums(tr_cpms_LMSA > 1) >= 4
tr_counts_clean_LMSA <- tr_counts_LMSA[keep_LMSA,]
head(tr_counts_clean_LMSA); length(rownames(tr_counts_clean_LMSA))

# save filtered counts  
dt_LMSA <- DGEList(counts = tr_counts_clean_LMSA, 
		   group = rep(c('LM','SA'), 
		   each=4))
# visualize 
head(dt_LMSA[["counts"]])

# compute normalization factors 
dt_LMSA <- calcNormFactors(dt_LMSA)
# normalize according to dispersion 
dt_LMSA <- estimateCommonDisp(dt_LMSA)
dt_LMSA <- estimateTagwiseDisp(dt_LMSA)

# test conditions
de.tr_LMSA <- exactTest(dt_LMSA, pair = c("LM","SA"))

# correction for multiple testing and sort de transcripts
tT.transcripts_LMSA <- topTags(de.tr_LMSA, n = nrow(dt_LMSA))
head(tT.transcripts_LMSA[["table"]])

# exract differentially expressed transcripts
table.tr_LMSA=tT.transcripts_LMSA$table
head(table.tr_LMSA); length(table.tr_LMSA[["logFC"]])

# select counts with FDR <0.1
top.det_LMSA <- table.tr_LMSA[table.tr_LMSA$FDR < 0.01 , ]
head(top.det_LMSA); length(top.det_LMSA[["logFC"]])

#compute RPKMs
RPKMs_LMSA <- rpkm(tr_counts_LMSA, geneLengths$length)

# prepare data to plot
TCA_LMSA=top.det_LMSA[row.names(top.det_LMSA) %in% TCA$transcript_id, ]
length(TCA_LMSA$logFC)
DeG_LMSA=top.det_LMSA[row.names(top.det_LMSA) %in% DeG$transcript_id, ]
length(DeG_LMSA$logFC)
TRANS_LMSA=top.det_LMSA[row.names(top.det_LMSA) %in% TRANS$transcript_id, ]
length(TRANS_LMSA$logFC)
UPTK_LMSA=top.det_LMSA[row.names(top.det_LMSA) %in% UPTK$transcript_id, ]
length(UPTK_LMSA$logFC)

# export logFC to csv
write.csv(TCA_LMSA, "TCA_LMSA.csv")
write.csv(DeG_LMSA, "DeG_LMSA.csv")
write.csv(TRANS_LMSA, "TRANS_LMSA.csv")
write.csv(UPTK_LMSA, "UPTK_LMSA.csv")

# create MA plot 
plot(x=table.tr_LMSA$logCPM, y=table.tr_LMSA$logFC, pch=19, main="MA-plot LM-SA", xlab="Average log(CPM)", ylab="log2(FC) SA-LM" )
points(x=TCA_LMSA[["logCPM"]], y=TCA_LMSA[["logFC"]], pch=19, col="blue")
points(x=DeG_LMSA[["logCPM"]], y=DeG_LMSA[["logFC"]], pch=19, col="green")
points(x=TRANS_LMSA[["logCPM"]], y=TRANS_LMSA[["logFC"]], pch=19, col="orange")
points(x=UPTK_LMSA[["logCPM"]], y=UPTK_LMSA[["logFC"]], pch=19, col="red")

# create volcano plot
plot(y=-log10(table.tr_LMSA[["PValue"]]), x=table.tr_LMSA$logFC, pch=19, main=" Volcano plot LM-SA", xlab="log2(FC) SA-LM", ylab="log10(p-value)")
points(y=-log10(TCA_LMSA[["PValue"]]), x=TCA_LMSA[["logFC"]], pch=19, col="blue")
points(y=-log10(DeG_LMSA[["PValue"]]), x=DeG_LMSA[["logFC"]], pch=19, col="green")
points(y=-log10(TRANS_LMSA[["PValue"]]), x=TRANS_LMSA[["logFC"]], pch=19, col="orange")
points(y=-log10(UPTK_LMSA[["PValue"]]), x=UPTK_LMSA[["logFC"]], pch=19, col="red")

# heatmaps
# TCA metabolism
TCA_RPKMs_LMSA=RPKMs_LMSA[row.names(RPKMs_LMSA) %in% TCA$transcript_id, ]
pheatmap(log2(TCA_RPKMs_LMSA), main = 'Heatmap of TCA metabolism genes: LM vs SA', scale="column", cluster_rows=FALSE)
# DeG metaolism
DeG_RPKMs_LMSA=RPKMs_LMSA[row.names(RPKMs_LMSA) %in% DeG$transcript_id, ]
pheatmap(log2(DeG_RPKMs_LMSA), main = 'Heatmap of G&D metabolism genes: LM vs SA', scale="column", cluster_rows=FALSE)
# transport & uptake systems
UPTK_RPKMs_LMSA=RPKMs_LMSA[row.names(RPKMs_LMSA) %in% UPTK$transcript_id, ]
pheatmap(log2(UPTK_RPKMs_LMSA), main = 'Heatmap of uptake genes: LM vs SA', scale="column", cluster_rows=FALSE)
# other organic acids metabolism
TRANS_RPKMs_LMSA=RPKMs_LMSA[row.names(RPKMs_LMSA) %in% TRANS$transcript_id, ]
pheatmap(log2(TRANS_RPKMs_LMSA), main = 'Heatmap of tranport genes: LM vs SA', scale="column", cluster_rows=FALSE)


