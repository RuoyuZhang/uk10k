rm(list=ls())
options(java.parameters = "-Xmx8000m")
library('optparse')
library('reshape2')

option_list <- list(
    make_option(c("--in_dir"), type="character", help="dir to input data", default="f:/Cornell/experiment/uk10k/uk10k/variants/"),
    make_option(c("--out_file"), type="character", help="output file", default="f:/Cornell/experiment/uk10k/uk10k/variants/variants.genotype.txt")
    
)

opt <- parse_args(OptionParser(option_list=option_list))
print(opt)

# rawfile
in_dir=opt$in_dir
# outdir
out_file=opt$out_file

files=list.files(path = in_dir, pattern = ".vcf$")

var.ped=NULL
for (file in files){
    # get file name
    sample=basename(file)
    sample=gsub('.vcf','',sample)
    # read sample
    absdir=paste0(in_dir,'/',file)
    f=try(read.table(absdir,header=F,fill=T))
    if (class(f)!='try-error'){
        data=read.table(absdir,header=F,fill=T)
        data=data[data$V7=="PASS",]
        data=data[nchar(as.character(data$V5))<2,]
        mutation=paste0(data$V4,data$V2,data$V5)
        detail=do.call(rbind,strsplit(as.character(data$V10),split=':'))[,3]
        
        data.file=paste(sample,mutation,detail)
    
        var.ped=c(var.ped,data.file)
    }
    
}

var.all=strsplit(var.ped,split=" ")
var.all=do.call(rbind,var.all)
colnames(var.all)=c("sampleID","snp","frequency")
var.all=as.data.frame(var.all)

spread=dcast(var.all,sampleID ~ snp,value.var = 'frequency',fill = NA)

spread[which(is.na(spread),arr.ind = T)]=0

#colnames(spread)=c("sampleID",as.character(snps))

write.table(spread,file=out_file,row.names = F,col.names = T,quote = F)

