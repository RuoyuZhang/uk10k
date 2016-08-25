rm(list=ls())
library(ggplot2)
library(reshape2)
library(plyr)
options(java.parameters = "-Xmx8000m")
library('optparse')
theme_set(theme_bw((base_size=20)))

option_list <- list(
    make_option(c("--multi_copy_file"), type="character", help="file containing copy number and phenotype for single samples", default="f:/Cornell/experiment/uk10k/uk10k/copynumber/calculate_copy_number/single.copy.txt"),
    make_option(c("--single_copy_file"), type="character", help="file containing copy number and phenotype for single samples", default="f:/Cornell/experiment/uk10k/uk10k/copynumber/calculate_copy_number/multi.copy.txt"),
    make_option(c("--pheno_list_file"), type="character", help="file containing phenotype list", default="f:/Cornell/experiment/uk10k/uk10k/copynumber/copy_phe/pheno_list.csv"),
    make_option(c("--out_dir"), type="character", help="output file", default="f:/Cornell/experiment/uk10k/uk10k/copynumber/copy_phe/")
    
)

opt <- parse_args(OptionParser(option_list=option_list))
print(opt)

# single file
single_file=opt$single_copy_file
# multi file
multi_file=opt$multi_copy_file
# phenotype list
pheno_list=opt$pheno_list_file
# outdir
out_dir=opt$out_dir

# read single file
single=read.table(single_file,header=T,sep=";")
# read multi file
multi=read.table(multi_file,header=T,sep=";")

single$type="single"
multi$type="multi"

# cancat two file
all.sample=rbind(single,multi)

ggplot(all.sample,aes(type,copynumber))+geom_boxplot()

# phenotype to investigate
pheno.list=read.csv(pheno_list,header=T)

# ages
ages=all.sample[,as.character(unique(pheno.list$age))]

ages$id=1:nrow(ages)
ages.melt=melt(ages,id.vars = "id")

ggplot(ages.melt,aes(id,value,col=variable))+geom_point()

ages.sd=apply(ages[,-ncol(ages)],1,sd,na.rm=T)
ages.mean=apply(ages[,-ncol(ages)],1,mean,na.rm=T)

all.sample$age_mean=ages.mean

select.in=which(ages.sd<5)

select.sample=all.sample[select.in,]

# copy number and age
lm.age=summary(lm(copynumber~age_mean+nuclear_cover,data=select.sample))
p=ggplot(select.sample,aes(age_mean,copynumber))+geom_point()+geom_smooth(method = "lm")+
    ylab("mtDNA copy number")+xlab("mean Age")
ggsave(p,file=paste0(out_dir,"Alkaline_copynumber_rmout.jpeg"),width=9, height=5, dpi=300)

# do asscoation for each phenotype
sig.level=0.05/nrow(pheno.list)
# all phenotypes
all.pheno=NULL
# significant phenotypes
sig.pheno=NULL
for (i in 1:nrow(pheno.list)){
    lm.res=summary(lm(as.numeric(as.character(select.sample[,as.character(pheno.list[i,1])]))~select.sample$copynumber+
                          as.numeric(as.character(select.sample[,as.character(pheno.list[i,2])]))+select.sample$nuclear_cover))
    res=paste(pheno.list[i,1],lm.res$coefficients[2,4])
    all.pheno=c(all.pheno,res)
    if (lm.res$coefficients[2,4] < sig.level){
        sig.pheno=c(sig.pheno,res)
    }
}

sig.pheno

# Alkaline and Bicarbonate are significant; take a deep look
lm.alk=summary(lm(Alkaline~copynumber+Alkaline_Age+nuclear_cover,data=select.sample))
lm.bic=summary(lm(Bicarbonate~copynumber+Bicarbonate_Age+nuclear_cover,data=select.sample))

p=ggplot(select.sample,aes(copynumber,Alkaline))+geom_point()+geom_smooth(method = "lm")+
    ylab("Serum alkaline phosphatase level")+xlab("mtDNA copy number")
ggsave(p,file=paste0(out_dir,"Alkaline_copynumber.jpeg"),width=9, height=5, dpi=300)
p=ggplot(select.sample,aes(copynumber,Bicarbonate))+geom_point()+geom_smooth(method = "lm")+
    ylab("Serum bicarbonate level")+xlab("mtDNA copy number")
ggsave(p,file=paste0(out_dir,"Bicarbonate_copynumber.jpeg"),width=9, height=5, dpi=300)

# looks there is an outlier in Alkaline; remove that one 
lm.alk2=summary(lm(Alkaline~copynumber+Alkaline_Age+nuclear_cover,data=select.sample[-72,]))

p=ggplot(select.sample[-72,],aes(copynumber,Alkaline))+geom_point()+geom_smooth(method = "lm")+
    ylab("Serum alkaline phosphatase level")+xlab("mtDNA copy number")
ggsave(p,file=paste0(out_dir,"Alkaline_copynumber_rmout.jpeg"),width=9, height=5, dpi=300)


