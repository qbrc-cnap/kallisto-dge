library(DESeq2)
library(tximport)

# This script assumes that the various abundance files have already been split
# out to various directories as required for sleuth 

args <- commandArgs(TRUE)
argsLen = length(args)
ANNOTATIONS <- args[1] # annotation file
TRANSCRIPT_TO_GENE_MAPPING_FILE <- args[2] # file mapping the transcript to gene names. Needs column headers of `enst` and `name`
GROUP_A <- args[3] # the 'base' condition
GROUP_B <- args[4] # the 'experimental' condition
OUTPUT_DESEQ_FILE <- args[5]
OUTPUT_NORMALIZED_COUNTS_FILE <- args[6]
PDX = FALSE
if (argsLen > 6) {
	PDX = TRUE
	HUMAN_TRANSCRIPT_PREFIX = args[7]
}

# read the annotations and set the row names to be the sample names"
ann_df = read.table(ANNOTATIONS, sep='\t', header=T, row.names='sample')

# ingest the transcript abundances and map to genes using tximport:
if (PDX){
    fname = paste('abundance', HUMAN_TRANSCRIPT_PREFIX, 'tsv', sep='.')
    files<-file.path(ann_df$path, fname)
} else {
    files<-file.path(ann_df$path, 'abundance.h5')
}
names(files) <- rownames(ann_df)

# prepare the annotations for DESeq2.  First keep only the condition column
ann_df <- ann_df['condition']
# relevel to get the comparison order correct:
ann_df$condition <- relevel(ann_df$condition, ref=GROUP_A)

mapping = read.table(TRANSCRIPT_TO_GENE_MAPPING_FILE, sep='\t', header=T)
txi.kallisto <- tximport(files, type='kallisto', tx2gene=mapping[c('enst','name')], ignoreTxVersion=T)


# now run the DESeq2 portion:
dds <- DESeqDataSetFromTximport(txi.kallisto, ann_df, ~condition)
dds <- DESeq(dds)
res <- results(dds, cooksCutoff=F, contrast=c("condition", GROUP_B, GROUP_A))
original_colnames = colnames(res)
n = length(original_colnames)
baseMeanA = rowMeans(counts(dds,normalized=TRUE)[,dds$condition == GROUP_A]) 
baseMeanB = rowMeans(counts(dds,normalized=TRUE)[,dds$condition == GROUP_B]) 
res = cbind(rownames(res), res[,1],baseMeanA, baseMeanB, as.data.frame(res[,2:n])) 
colnames(res) = c('Gene', 'overall_mean', GROUP_A, GROUP_B, original_colnames[2:n])
resOrdered <- res[order(res$padj),]
write.table(resOrdered, OUTPUT_DESEQ_FILE, sep='\t', quote=F, row.names=F)

#normalized counts:
dds <- estimateSizeFactors(dds)
nc <- counts(dds, normalized=TRUE)
nc <- cbind(gene=rownames(nc), nc)
write.table(nc, OUTPUT_NORMALIZED_COUNTS_FILE, sep='\t', quote=F, row.names=F)
