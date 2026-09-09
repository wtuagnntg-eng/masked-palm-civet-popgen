args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=4) stop('Usage: Rscript summarize_ld.R wild.stat.gz farm.stat.gz summary.tsv curves.tsv')
read_ld <- function(path, group) {
  d <- read.table(gzfile(path), header=FALSE, stringsAsFactors=FALSE, comment.char='#')
  if (ncol(d)<2) stop(paste('LD stat file has fewer than two columns:', path))
  d <- d[,1:2]; names(d) <- c('Dist','r2'); d$Group <- group
  d[is.finite(d$Dist) & is.finite(d$r2),]
}
d <- rbind(read_ld(args[1],'Wild'), read_ld(args[2],'Farm'))
d30 <- d[d$Dist <= 30000,]
s <- aggregate(r2 ~ Group, d30, mean)
names(s)[2] <- 'Mean_r2_0_30kb'
write.table(s,args[3],sep='\t',quote=FALSE,row.names=FALSE)
write.table(d30,args[4],sep='\t',quote=FALSE,row.names=FALSE)
print(s)
