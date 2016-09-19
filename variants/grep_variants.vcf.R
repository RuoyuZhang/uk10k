rm(list=ls())
options(java.parameters = "-Xmx8000m")
library('optparse')
library('reshape2')

option_list <- list(
    make_option(c("--in_dir"), type="character", help="dir to input data", default="f:/Cornell/experiment/uk10k/uk10k/variants/"),
    make_option(c("--variants_file"), type="character", help="variants file", default="f:/Cornell/experiment/uk10k/uk10k/copynumber/copy_phe/result/mtDNA_variants_list.txt"),
    make_option(c("--out_file"), type="character", help="output file", default="f:/Cornell/experiment/uk10k/uk10k/variants/variants.genotype.txt")
    
)

opt <- parse_args(OptionParser(option_list=option_list))
print(opt)

# rawfile
in_dir=opt$in_dir
# variants file
var_file=opt$variants_file
# outdir
out_file=opt$out_file

files=list.files(path = in_dir, pattern = ".vcf$")

var.ped=NULL
for (file in files){
    # get file name
    sample=basename(file)
    sample=gsub('.snp.vcf','',sample)
    # read sample
    absdir=paste0(in_dir,'/',file)
    f=try(read.table(absdir,header=F,fill=T))
    if (class(f)!='try-error'){
        data=read.table(absdir,header=F,fill=T)
    
        data.file=paste(sample,paste0(data$V4,data$V2,data$V5))
    
        var.ped=c(var.ped,data.file)
    }
    
}

var.all=strsplit(var.ped,split=" ")
var.all=do.call(rbind,var.all)
colnames(var.all)=c("sampleID","snp")
var.all=as.data.frame(var.all)

snps=unique(var.all$snp)
samples=unique((var.all$sampleID))

spread=NULL
for (sample in samples){
    temp=subset(var.all,sampleID==sample)
    line=sample
    for (snp in snps){
        if (snp %in% temp$snp){
            line=c(line,1)
        }else{
            line=c(line,0)
        }
    }
    spread = rbind(spread,line)
}

colnames(spread)=c("sampleID",as.character(snps))

write.table(spread,file=out_file,row.names = F,col.names = T,quote = F)

