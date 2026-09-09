suppressPackageStartupMessages({library(ggplot2); library(pheatmap)})
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=3) stop('Usage: Rscript 20_transcriptome_descriptive_plots.R gene_TPM.tsv immune_panel.tsv outdir')
tpm <- read.delim(args[1], check.names=FALSE, row.names=1)
panel <- read.delim(args[2], stringsAsFactors=FALSE)
out <- args[3]; dir.create(out, showWarnings=FALSE, recursive=TRUE)
X <- log2(as.matrix(tpm)+1)
v <- apply(X,1,var); Xnz <- X[is.finite(v) & v>0,,drop=FALSE]
# PCA uses all retained genes with non-zero variance.
pc <- prcomp(t(Xnz), center=TRUE, scale.=FALSE)
var_exp <- pc$sdev^2/sum(pc$sdev^2)
pca_df <- data.frame(sample=rownames(pc$x), PC1=pc$x[,1], PC2=pc$x[,2])
write.table(pca_df,file.path(out,'PCA_scores.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(PC=paste0('PC',seq_along(var_exp)),variance_explained=var_exp),file.path(out,'PCA_variance.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
p <- ggplot(pca_df,aes(PC1,PC2,label=sample))+geom_point(size=3)+geom_text(vjust=-0.7)+theme_classic()+
  xlab(sprintf('PC1 (%.1f%%)',100*var_exp[1]))+ylab(sprintf('PC2 (%.1f%%)',100*var_exp[2]))
ggsave(file.path(out,'Figure5A_PCA.pdf'),p,width=6,height=5)
# Top 100 variable genes for Figure S2; row-scale only for visualization.
top <- names(sort(v[is.finite(v)],decreasing=TRUE))[seq_len(min(100,sum(is.finite(v))))]
top_raw <- X[top,,drop=FALSE]; top_scaled <- t(scale(t(top_raw))); top_scaled[!is.finite(top_scaled)] <- 0
write.table(top_raw,file.path(out,'top100_log2TPM.tsv'),sep='\t',quote=FALSE,col.names=NA)
write.table(top_scaled,file.path(out,'top100_row_scaled.tsv'),sep='\t',quote=FALSE,col.names=NA)
pdf(file.path(out,'FigureS2_top100_heatmap.pdf'),width=8,height=10); pheatmap(top_scaled,cluster_rows=TRUE,cluster_cols=TRUE,scale='none'); dev.off()
# Match the complete predefined panel to detected genes.
hit <- intersect(panel$gene, rownames(X)); detected <- panel[match(hit,panel$gene),,drop=FALSE]
write.table(detected,file.path(out,'immune_panel_detected.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
# Immune categories exclude Viral-entry.
immune_genes <- detected$gene[detected$category!='Viral-entry']
if(length(immune_genes)>0){
  z <- t(scale(t(X[immune_genes,,drop=FALSE]))); z[!is.finite(z)] <- 0
  write.table(z,file.path(out,'immune_panel_row_scaled.tsv'),sep='\t',quote=FALSE,col.names=NA)
  pdf(file.path(out,'Figure5B_immune_heatmap.pdf'),width=8,height=max(6,0.18*nrow(z))); pheatmap(z,cluster_rows=TRUE,cluster_cols=TRUE,scale='none'); dev.off()
}
entry <- detected$gene[detected$category=='Viral-entry']
if(length(entry)>0){
  z <- t(scale(t(X[entry,,drop=FALSE]))); z[!is.finite(z)] <- 0
  write.table(z,file.path(out,'viral_entry_row_scaled.tsv'),sep='\t',quote=FALSE,col.names=NA)
  pdf(file.path(out,'Figure5C_viral_entry_heatmap.pdf'),width=8,height=max(4,0.45*nrow(z))); pheatmap(z,cluster_rows=FALSE,cluster_cols=TRUE,scale='none'); dev.off()
}
