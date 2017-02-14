rm(list=ls())
library(data.table)
library(qqman)
library('optparse')
library(ggplot2)
library(RColorBrewer)

theme_set(theme_bw((base_size=20)))

option_list <- list(
    make_option(c("--in_file"), type="character", help="file containing gwas result", default="/home/fs01/rz253/project/uk10k/analysis/copynumber/gwas/plink1/plink.all.pop.age.log.assoc.linear"),
    make_option(c("--out_dir"), type="character", help="output dir", default="/home/fs01/rz253/project/uk10k/analysis/copynumber/gwas/plot/"),
    make_option(c("--file_pre"), type="character", help="output file prefix", default='plink_age_log')
)

opt <- parse_args(OptionParser(option_list=option_list))
print(opt)

# input file
in_file=opt$in_file
# outdir
out_dir=opt$out_dir
# names

#SNP CHR  BP      P zscore
data <- fread(in_file, header=TRUE,data.table = F)
df = data[data$TEST=='GENO_2DF',]
df = df[!is.na(df$P),]


#df=data.frame(SNP=NA,CHR=as.character(data2[,names[1]]),BP=data2[,names[2]],P=data[,names[3]],zscore=NA)
#df=na.omit(df)
df$CHR=as.character(df$CHR)
df$CHR[df$CHR=='X']=23
df$CHR=as.numeric(df$CHR)

# top hit
top=df[order(df$P)[1:100],]
write.table(top,file=paste0(opt$out_dir, '/',opt$file_pre, ".tophit.txt"),row.names=F,quote=F)

#df.temp=df[1:1000000,]
tiff(filename=paste0(opt$out_dir, '/',opt$file_pre, ".manhattan_plot.jpeg"),width=174, height=100, units="mm",res=300)
manhattan(df,ylim=c(0,9),cex.axis = 1, col=brewer.pal(8, "Set2"))
dev.off()

tiff(filename=paste0(opt$out_dir,'/',opt$file_pre, ".qq_plot.tiff"),width=114, height=114, units="mm",res=300)
qq(df$P)
dev.off()
