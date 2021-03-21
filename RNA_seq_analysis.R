outputdir = 'output'

########################################################################
# reading transcript counts and preprocessing
#
library(SummarizedExperiment)
seq_count_dirs = c("SRR2156848","SRR2156849","SRR2156850","SRR2156851")
m_seqcounts = readKallisto(paste('data',seq_count_dirs,"abundance.tsv",sep="/"), as="matrix") # matrix of shape (transcipts, samples)
print('shape of transcript matrix (transcripts, shapes):')
print(dim(m_seqcounts))
print('total transcript counts per sample')
print(colSums(m_seqcounts))
print('transcripts detected in at least one sample')
print(sum(rowSums(m_seqcounts) > 0))
print('transcript types per sample')
print(colSums(m_seqcounts>0))
print('filter out transcripts with no reads')
print('shape of filtered transcript matrix (transcripts present, shapes):')
is_present = rowSums(m_seqcounts) > 0
m_seqcounts_present = m_seqcounts[is_present,]
print(dim(m_seqcounts_present))
## scale by total transcript num by transcripts
#m_seqcounts_present = m_seqcounts_present / rowSums(m_seqcounts_present)

## scale by total transcript num by samples
#m_seqcounts_present = t(t(m_seqcounts_present) / colSums(m_seqcounts_present))

########################################################################
# PCA of samples
print('calculating PCA')
pca_result = prcomp(t(m_seqcounts_present), scale. = FALSE)
pca_scores = pca_result$x[,1:2]
PCA_x_scores = pca_scores[,1]
PCA_y_scores = pca_scores[,2]
png(file.path(outputdir, "PCA_samples.png"))
plot(PCA_x_scores, PCA_y_scores, type="p", main="PCA of samples (not scaled)", xlab="PC1", ylab="PC2", pch=16, col="blue")
text(PCA_x_scores, PCA_y_scores, labels=colnames(m_seqcounts_present))
dev.off()

########################################################################
# search differencially expressed transcripts between the treatment groups
## linear model fit
library(limma)
treatment_groups = c(0,0,1,1)
design_matrix = cbind(intercept=1, treatment=c(0,0,1,1))
voom_obj = voom(m_seqcounts_present, design_matrix)
modelfit = eBayes(lmFit(voom_obj, design=design_matrix))
## filter for differential expression
modelfit_output = topTable(modelfit, coef=2, adjust="fdr", number=nrow(m_seqcounts_present), p.value=0.05)
modelfit_output$ID = rownames(modelfit_output)
## assign gene symbols
library(biomaRt)
ensembl = useMart(host="feb2014.archive.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
lookup_table = getBM(attributes=c("ensembl_transcript_id", "hgnc_symbol"), mart=ensembl)
lookup_table = lookup_table[!duplicated(lookup_table[,1]),]
rownames(lookup_table) = as.character(lookup_table[,1])
modelfit_output$Symbol = lookup_table[rownames(modelfit_output), "hgnc_symbol"]
write.table(modelfit_output, file=file.path(outputdir, "modelfit_output.txt"), sep="\t", quote=FALSE, row.names=FALSE)

diffexp_transcript_IDs = as.character(modelfit_output$ID)
diffexp_transcript_symbols = as.character(modelfit_output$Symbol)
de_row_idxs = which(rownames(m_seqcounts_present) %in% diffexp_transcript_IDs)

library(gplots)
png(file.path(outputdir, "diff_exp_heatmap.png"))
heatmap(m_seqcounts_present[de_row_idxs,], labCol=colnames(m_seqcounts_present), labRow=diffexp_transcript_symbols)
dev.off()
