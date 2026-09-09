suppressPackageStartupMessages(library(ggplot2))
a <- commandArgs(trailingOnly=TRUE); if(length(a)!=2) stop('Usage: Rscript 21_plot_figure1.R mapping_qc.tsv outdir')
d <- read.delim(a[1]); out <- a[2]; dir.create(out,showWarnings=FALSE,recursive=TRUE)
metrics <- c(mean_depth='Mean depth',mapping_rate='Mapping rate (%)',properly_paired_rate='Properly paired (%)',cov_1x='Coverage >=1x',cov_5x='Coverage >=5x',cov_10x='Coverage >=10x')
for(m in names(metrics)){
  p <- ggplot(d,aes(x=ecotype,y=.data[[m]],fill=ecotype))+geom_boxplot(outlier.shape=NA)+geom_jitter(width=.08)+theme_classic()+labs(x=NULL,y=metrics[[m]])+theme(legend.position='none')
  ggsave(file.path(out,paste0('Figure1_',m,'.pdf')),p,width=4.5,height=4)
}
