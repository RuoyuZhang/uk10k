rm(list=ls())
options(java.parameters = "-Xmx8000m")
library('optparse')

option_list <- list(
    make_option(c("--ped_file"), type="character", help="dir to input ped file", default="/home/fs01/rz253/project/uk10k/analysis/copynumber/gwas/prepare/pheno.ped"),
    make_option(c("--pca_file"), type="character", help="dir to pca file", default="/home/fs01/rz253/project/uk10k/analysis/copynumber/gwas/prepare/pca.pca"),
    make_option(c("--haplogroup_file"), type="character", help="dir to haplogroup file", default="/home/fs01/rz253/project/uk10k/analysis/copynumber/gwas/prepare/haplogroups.txt"),
    make_option(c("--out_dir"), type="character", help="output file", default="/home/fs01/rz253/project/uk10k/analysis/copynumber/gwas/prepare/")
    
)

opt <- parse_args(OptionParser(option_list=option_list))
#print(opt)

# ped file
ped_file=opt$ped_file
# pca file
pca_file=opt$pca_file
# hapgroup file
hap_file=opt$haplogroup_file
# out dir
out_dir=opt$out_dir

# read files
ped=read.table(ped_file,header=T)
pca=read.table(pca_file,header=T)
hap=read.table(hap_file,header = T,fill=T)

# process 
hap$hap=substr(hap$Haplogroup,1,1)

# merge files
temp=merge(ped,hap[,c(1,3)],by.x="fid",by.y="sampleID")
temp2=merge(temp,pca[,c(1,4,5)],by.x="fid",by.y="FID")

write.table(temp2,file=paste0(out_dir,"/pheno.hap.ped"),row.names = F,quote = F)

