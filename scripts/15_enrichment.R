suppressPackageStartupMessages(library(clusterProfiler))
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 6) {
  stop('Usage: Rscript 15_enrichment.R candidate_genes.txt gene2go.tsv gene2ko.tsv ko2pathway.tsv out_prefix minGSSize')
}
candidate_file <- args[1]
gene2go_file <- args[2]
gene2ko_file <- args[3]
ko2path_file <- args[4]
out_prefix <- args[5]
minGS <- as.integer(args[6])

cand <- unique(scan(candidate_file, what='character', quiet=TRUE))
g2go <- read.delim(gene2go_file, stringsAsFactors=FALSE)
g2ko <- read.delim(gene2ko_file, stringsAsFactors=FALSE)
k2p <- read.delim(ko2path_file, stringsAsFactors=FALSE)
colnames(g2go)[1:2] <- c('gene','term')
colnames(g2ko)[1:2] <- c('gene','ko')
colnames(k2p)[1:2] <- c('ko','pathway')

# GO: all genes carrying GO annotation are the universe.
go_term2gene <- unique(g2go[,c('term','gene')])
go_bg <- unique(go_term2gene$gene)
go_cand <- intersect(cand, go_bg)
go <- enricher(go_cand, universe=go_bg, TERM2GENE=go_term2gene,
               pAdjustMethod='BH', minGSSize=minGS, pvalueCutoff=1, qvalueCutoff=1)
if (!is.null(go) && nrow(as.data.frame(go))>0) {
  write.table(as.data.frame(go), paste0(out_prefix,'.GO.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
} else {
  write.table(data.frame(), paste0(out_prefix,'.GO.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
}

# KEGG: gene -> KO -> pathway. Genes mappable to a pathway form the universe.
kg <- unique(merge(g2ko, k2p, by='ko')[,c('pathway','gene')])
kegg_bg <- unique(kg$gene)
kegg_cand <- intersect(cand, kegg_bg)
ke <- enricher(kegg_cand, universe=kegg_bg, TERM2GENE=kg,
               pAdjustMethod='BH', minGSSize=minGS, pvalueCutoff=1, qvalueCutoff=1)
if (!is.null(ke) && nrow(as.data.frame(ke))>0) {
  write.table(as.data.frame(ke), paste0(out_prefix,'.KEGG.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
} else {
  write.table(data.frame(), paste0(out_prefix,'.KEGG.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
}
