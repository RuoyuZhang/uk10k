rm(list=ls())
options(java.parameters = "-Xmx8000m")
library('optparse')

option_list <- list(
    make_option(c("--in_file"), type="character", help="file of input data", default="f:/Cornell/experiment/uk10k/uk10k/variants/UK10K_TW5093554.raw"),
    make_option(c("--variants_file"), type="character", help="variants file", default="f:/Cornell/experiment/uk10k/uk10k/copynumber/copy_phe/result/mtDNA_het_list.txt"),
    make_option(c("--out_file"), type="character", help="output file", default="f:/Cornell/experiment/uk10k/uk10k/variants/UK10K_TW5093554.het.txt")
    
)

opt <- parse_args(OptionParser(option_list=option_list))
print(opt)

# rawfile
in_file=opt$in_file
# variants file
var_file=opt$variants_file
# outdir
out_file=opt$out_file

# read chromsome length
variants=as.character(read.table(var_file,header=F)[[1]])

    # get file name
    sample=basename(in_file)
    sample=gsub('.raw','',sample)
    # read sample
    data=read.table(in_file,header=F,fill=T)
    
    line = c(sample)
    for (variant in variants){
        pos = gsub('[ATGC]','',variant)
        splitv=strsplit(variant,split = '')[[1]]
        n1=splitv[1]
        n2=splitv[length(splitv)]
        if (pos %in% data$V2){
            use = data[which(data$V2==pos),c(4,5,6,7,8)]
            names(use)=c("dep","A","T","G","C")
 #           print(use)
            if (use['dep']<50){
                line=c(line,NA)
            }else{
                f1=as.numeric(strsplit(as.character(unlist(use[n1])),split = ':')[[1]][5])
                f2=as.numeric(strsplit(as.character(unlist(use[n2])),split = ':')[[1]][5])
                line=c(line,f2)
            }
        }else{
            line=c(line,NA)
        }
    }
    
line1=paste(line,collapse = " ")[1]
write.table(line1,file=out_file,row.names = F,col.names = F,quote = F)

