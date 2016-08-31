rm(list=ls())
options(java.parameters = "-Xmx8000m")
library('optparse')

option_list <- list(
    make_option(c("--multi_file"), type="character", help="dir to input ped file", default="/home/fs01/rz253/project/uk10k/analysis/copynumber/calculate_copy_number/single.copy.txt"),
    make_option(c("--single_file"), type="character", help="dir to input ped file", default="/home/fs01/rz253/project/uk10k/analysis/copynumber/calculate_copy_number/multi.copy.txt"),
    make_option(c("--pca_file"), type="character", help="dir to pca file", default="/home/fs01/rz253/project/uk10k/analysis/copynumber/gwas/prepare/pca.pca"),
    make_option(c("--haplogroup_file"), type="character", help="dir to haplogroup file", default="/home/fs01/rz253/project/uk10k/analysis/copynumber/gwas/prepare/haplogroups.txt"),
    make_option(c("--pheno_list_file"), type="character", help="file containing phenotype list", default="/home/fs01/rz253/project/uk10k/uk10k/copynumber/copy_phe/pheno_list.csv"),
    make_option(c("--out_dir"), type="character", help="output file", default="/home/fs01/rz253/project/uk10k/analysis/copynumber/gwas/prepare/")
    
)

opt <- parse_args(OptionParser(option_list=option_list))
#print(opt)

# single file
single_file=opt$single_file
# multi file
multi_file=opt$multi_file
# pca file
pca_file=opt$pca_file
# hapgroup file
hap_file=opt$haplogroup_file
# phe_list
pheno_list_file=opt$pheno_list_file
# out dir
out_dir=opt$out_dir

# read files
single=read.table(single_file,header = T,sep=';')
multi=read.table(multi_file,header=T,sep=';')
pca=read.table(pca_file,header=T)
hap=read.table(hap_file,header = T,fill=T)
pheno.list=read.table(pheno_list_file,header = T,sep=',')

# all samples
all.sample=rbind(single,multi)

age=all.sample$Hematology_Tests_Age

for (i in which(is.na(age))){
    age[i]=(all.sample[i,as.character(unique(pheno.list$age))][!is.na(all.sample[i,as.character(unique(pheno.list$age))])][1])
}

ped=cbind(as.character(all.sample$sampleID),as.character(all.sample$sampleID),rep("0",length(all.sample$sampleID)),rep("0",length(all.sample$sampleID)),
          rep("2",length(all.sample$sampleID)),all.sample$copynumber,age,age^2)

colnames(ped)=c("fid","iid","fatid","matid","sex","copynumber","age","age2")

# process 
hap$hap=substr(hap$Haplogroup,1,1)

# merge files
temp=merge(ped,hap[,c(1,3)],by.x="fid",by.y="SampleID")
temp2=merge(temp,pca[,c(1,4,5)],by.x="fid",by.y="FID")

write.table(temp2,file=paste0(out_dir,"/pheno.hap.ped"),row.names = F,quote = F)

