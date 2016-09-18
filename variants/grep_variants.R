rm(list=ls())
options(java.parameters = "-Xmx8000m")
library('optparse')

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

# read chromsome length
variants=as.character(read.table(var_file,header=F)[[1]])

files=list.files(path = in_dir, pattern = ".raw")

var.ped=NULL
for (file in files){
    # get file name
    sample=basename(file)
    sample=gsub('.raw','',sample)
    # read sample
    absdir=paste0(in_dir,'/',file)
    data=read.table(absdir,header=F,fill=T)
    
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
            if (use['dep']<200){
                line=c(line,NA)
            }else{
                f1=as.numeric(strsplit(as.character(unlist(use[n1])),split = ':')[[1]][5])
                f2=as.numeric(strsplit(as.character(unlist(use[n2])),split = ':')[[1]][5])
                if (f1 > 0.9){
                    line=c(line,0)
                }else if (f2 > 0.9){
                    line=c(line,1)
                }else{
                    line=c(line,NA)
                }
            }
        }else{
            line=c(line,NA)
        }
    }
    var.ped=rbind(var.ped,line)
    
}

colnames(var.ped)=c("sampleID",variants)

write.table(var.ped,file=out_file,row.names = F,col.names = T,quote = F)

