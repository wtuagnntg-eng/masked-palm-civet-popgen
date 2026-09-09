suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
a <- commandArgs(trailingOnly=TRUE); if(length(a)!=2) stop('Usage: Rscript 24_plot_figure4.R differentiation_dir outdir')
d <- read.delim(file.path(a[1],'matched_windows.tsv')); th <- read.delim(file.path(a[1],'thresholds.tsv')); out <- a[2]; dir.create(out,showWarnings=FALSE,recursive=TRUE)
d$CHROM <- factor(d$CHROM, levels=unique(d$CHROM)); lens <- aggregate(END~CHROM,d,max); off <- c(0,cumsum(head(lens$END,-1))); names(off)<-as.character(lens$CHROM); d$x <- d$START + off[as.character(d$CHROM)]
fst_thr <- th$threshold[grepl('weighted_FST',th$metric)][1]; ratio_thr <- th$threshold[grepl('piWild_piFarm',th$metric)][1]
p1 <- ggplot(d,aes(x,WEIGHTED_FST))+geom_point(size=.3)+geom_hline(yintercept=fst_thr,linetype=2)+theme_classic()+xlab(NULL)+ylab(expression(F[ST]))
p2 <- ggplot(d,aes(x,PI_WILD_DIV_FARM))+geom_point(size=.3)+geom_hline(yintercept=ratio_thr,linetype=2)+theme_classic()+xlab('Genomic position')+ylab(expression(pi[Wild]/pi[Farm]))
ggsave(file.path(out,'Figure4.pdf'),p1/p2,width=11,height=7)
