args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop('Usage: Rscript wilcoxon_individual_metrics.R input.tsv value_column output.tsv')
d <- read.delim(args[1], check.names=FALSE)
v <- args[2]
x <- d[d$ecotype=='Farm', v]
y <- d[d$ecotype=='Wild', v]
tt <- wilcox.test(x, y, alternative='two.sided', exact=FALSE)
out <- data.frame(metric=v, farm_n=length(x), wild_n=length(y), farm_mean=mean(x,na.rm=TRUE), wild_mean=mean(y,na.rm=TRUE), p_value=tt$p.value)
write.table(out,args[3],sep='\t',quote=FALSE,row.names=FALSE)
