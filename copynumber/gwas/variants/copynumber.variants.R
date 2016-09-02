rm(list=ls())
library(ggplot2)
library(reshape2)
library(plyr)
library(optparse)
options(java.parameters = "-Xmx8000m")

theme_set(theme_bw((base_size=20)))

option_list <- list(
    make_option(c("--multi_copy_file"), type="character", help="file containing copy number and phenotype for single samples", default="f:/Cornell/experiment/uk10k/uk10k/copynumber/calculate_copy_number/single.copy.txt"),
    make_option(c("--single_copy_file"), type="character", help="file containing copy number and phenotype for single samples", default="f:/Cornell/experiment/uk10k/uk10k/copynumber/calculate_copy_number/multi.copy.txt"),
    make_option(c("--indv_file"), type="character", help="dir to file containing individual id", default="f:/Cornell/experiment/uk10k/uk10k/copynumber/gwas/variants/22.min.012.indv"),
    make_option(c("--geno_file"), type="character", help="dir to file containing individual genotype", default="f:/Cornell/experiment/uk10k/uk10k/copynumber/gwas/variants/22.min.012"),
    make_option(c("--pos_file"), type="character", help="dir to file containing variants position", default="f:/Cornell/experiment/uk10k/uk10k/copynumber/gwas/variants/22.min.012.pos"),
    make_option(c("--out_dir"), type="character", help="output file", default="f:/Cornell/experiment/uk10k/uk10k/copynumber/gwas/plot/")
    
)

opt <- parse_args(OptionParser(option_list=option_list))
print(opt)

# single file
single_file=opt$single_copy_file
# multi file
multi_file=opt$multi_copy_file

# read single file
single=read.table(single_file,header=T,sep=";")
# read multi file
multi=read.table(multi_file,header=T,sep=";")

single$type="single"
multi$type="multi"

# cancat two file
all.sample=rbind(single,multi)

# genotypes
indv=read.table(opt$indv_file,header = F)
geno=read.table(opt$geno_file,header=F)
pos=read.table(opt$pos_file,header=F)

geno.df=data.frame(sampleID=indv[[1]],geno=geno$V2)

merge.df=merge(all.sample[,c('sampleID','copynumber')],geno.df)

name=paste0('chr',pos$V1,'_',pos$V2)

p=ggplot(merge.df,aes((geno),copynumber))+geom_point()+geom_smooth(method = 'lm')+labs(x="",y="mtDNA copy number",title=name)
ggsave(p,filename=paste0(opt$out_dir,"/",name,".jpeg"), width=8, height=5, dpi=300)
