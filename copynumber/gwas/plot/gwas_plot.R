rm(list=ls())
library(ggplot2)
library(reshape2)
library(plyr)
library(data.table)
library(qqman)
options(java.parameters = "-Xmx8000m")
library('optparse')
theme_set(theme_bw((base_size=20)))

option_list <- list(
    make_option(c("--in_file"), type="character", help="file containing gwas result", default="/home/fs01/rz253/project/uk10k/analysis/copynumber/gwas/rvtest1/rvtest2.SingleWald.assoc"),
    make_option(c("--out_dir"), type="character", help="output dir", default="/home/fs01/rz253/project/uk10k/analysis/copynumber/gwas/plot/"),
    make_option(c("--file_pre"), type="character", help="output file prefix", default='test1'),
    make_option(c("--names"), type="character", help="useful column names", default="CHROM,POS,Pvalue"),
    make_option(c("--factors"), type="integer", help="testing variable numbers", default=5)
)

opt <- parse_args(OptionParser(option_list=option_list))
print(opt)

# input file
in_file=opt$in_file
# outdir
out_dir=opt$out_dir
# names
names=strsplit(opt$names,split = ',')[[1]]

#SNP CHR  BP      P zscore
data <- fread(in_file, header=TRUE,data.table = F,sep='\t')
data2=data[seq(1,nrow(data),opt$factors)]

df=data.frame(SNP=NA,CHR=as.character(data2[,names[1]]),BP=data2[,names[2]],P=data[,names[3]],zscore=NA)
df=na.omit(df)
df$CHR=as.character(df$CHR)
df$CHR[df$CHR=='X']=23
df$CHR=as.numeric(df$CHR)

jpeg(filename=paste0(opt$out_dir, '/',opt$file_pre, "manhattan_plot.jpeg"),width=15, height=9, units="in",res=300)
manhattan(df)
dev.off()

jpeg(filename=paste0(opt$out_dir,'/',opt$file_pre, "qq_plot.jpeg"),width=15, height=9, units="in",res=300)
qq(df$P)
dev.off()