rm(list=ls())
options(java.parameters = "-Xmx8000m")
library('optparse')

option_list <- list(
    make_option(c("--in_dir"), type="character", help="dir to input data", default="f:/Cornell/experiment/uk10k/uk10k/copynumber/calculate_copy_number/"),
    make_option(c("--chr_length"), type="character", help="output dir", default="f:/Cornell/experiment/uk10k/uk10k/copynumber/calculate_copy_number/chr.length"),
    make_option(c("--phe_file"), type="character", help="phenotype file", default="f:/Cornell/experiment/uk10k/uk10k/copynumber/calculate_copy_number/_EGAZ00001016607_UK10K_TWINSUK_Phenotype_Data_October2013_1754samples.txt"),
    make_option(c("--out_file"), type="character", help="output file", default="f:/Cornell/experiment/uk10k/uk10k/copynumber/calculate_copy_number/copy.phe.txt")
    
)

opt <- parse_args(OptionParser(option_list=option_list))
print(opt)

# rawfile
in_dir=opt$in_dir
# chr length file
chr_l=opt$chr_length
# phenotype file
phe_file=opt$phe_file
# outdir
out_file=opt$out_file

# read chromsome length
chr=read.table(chr_l,header=F)
rownames(chr)=chr$V1


files=list.files(path = in_dir, pattern = ".depth")

copys=NULL
for (file in files){
    # get file name
    sample=basename(file)
    sample=gsub('.depth','',sample)
    # read sample
    absdir=paste0(in_dir,'/',file)
    data=read.table(absdir,header=F)
    rownames(data)=data$V1
    # autosome
    auto=1:22
    auto.all=data[auto,4]
    auto.cover=auto.all/chr[,2]

    # mt cover
    mt.cover=data['MT',4]/16569

    copy = 2*mt.cover/mean(auto.cover,na.rm = T)
    
    res=c(sample,copy)
    copys=rbind(copys,res)
}

colnames(copys)=c("sampleID","copynumber")

phe=read.table(phe_file,header=T,sep=';')

phe.copy=merge(copys,phe,by.x="sampleID",by.y="Public_ID")

write.table(res,file=out_file,row.names = F,col.names = F,quote = F,sep=;)

