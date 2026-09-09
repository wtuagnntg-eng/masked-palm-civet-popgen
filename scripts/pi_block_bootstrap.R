args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 5) stop('Usage: Rscript pi_block_bootstrap.R wild.pi farm.pi nboot summary.tsv paired_blocks.tsv')
w <- read.delim(args[1], check.names=FALSE)
f <- read.delim(args[2], check.names=FALSE)
nboot <- as.integer(args[3])
m <- merge(w[,c('CHROM','BIN_START','BIN_END','PI')], f[,c('CHROM','BIN_START','BIN_END','PI')],
           by=c('CHROM','BIN_START','BIN_END'), suffixes=c('_Wild','_Farm'))
m <- m[is.finite(m$PI_Wild) & is.finite(m$PI_Farm),]
write.table(m, args[5], sep='\t', quote=FALSE, row.names=FALSE)
delta <- m$PI_Wild - m$PI_Farm
if (length(delta)==0) stop('No paired 5-Mb blocks were available.')
set.seed(20260909)
b <- replicate(nboot, mean(sample(delta, length(delta), replace=TRUE)))
out <- data.frame(n_blocks=nrow(m), wild_mean=mean(m$PI_Wild), farm_mean=mean(m$PI_Farm),
                  delta_wild_minus_farm=mean(delta), ci_2.5=unname(quantile(b,0.025)),
                  ci_97.5=unname(quantile(b,0.975)), n_boot=nboot)
write.table(out,args[4],sep='\t',quote=FALSE,row.names=FALSE)
cat('Paired 5-Mb blocks:', nrow(m), '\n')
print(out)
