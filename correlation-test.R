library(usedist)
library(reshape2)
library(plyr)

args=commandArgs(trailingOnly=TRUE)

#file1: distance matrix of pairwise allelic mismatch from Genome Comparator plugin in PubMLST website
#file2: distance matrix of branch length converted using 'branchlength-to-matrix.py' script
#file3: user-defined output file (.txt)

allele_tab <- read.table(file=args[1], header=TRUE, sep="\t", row.names=1)
allele <- as.dist(as(allele_tab, "matrix"))
df_allele <- melt(as.matrix(allele), varnames = c("row", "col"))

nt_tab <- read.table(file=args[2], header=TRUE, sep="\t", row.names=1)
nt <- as.dist(as(nt_tab, "matrix"))
df_nt <- melt(as.matrix(nt), varnames = c("row", "col"))
df_nt <- df_nt[order(df_nt[,2], df_nt[,1]), ]
df_nt$log10value <- with(df_nt, -1*log10(value))

compare_data2 <- data.frame(df_allele$value, df_nt$log10value)

#correlation test
result2 = cor.test(compare_data2$df_allele.value, compare_data2$df_nt.log10value, method = "spearman")

write(result2, file=args[3])
