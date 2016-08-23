rm(list=ls())
options(java.parameters = "-Xmx8000m")
library('optparse')

option_list <- list(
    make_option(c("--raw_file"), type="character", help="File containing original depth information", default="f:/Cornell/experiment/uk10k/uk10k/copynumber/calculate_copy_number/QTL211893.depth"),
    make_option(c("--out_dir"), type="character", help="File containing original depth information", default="f:/Cornell/experiment/uk10k/uk10k/copynumber/calculate_copy_number/")
    
)

opt <- parse_args(OptionParser(option_list=option_list))
print(opt)

# rawfile
file=opt$raw_file
# outdir
out_dir=opt$out_dir

# get file name
sample=basename(file)
sample=gsub('.depth','',sample)

# read sample
data=read.table(file,header=F)
rownames(data)=data$V1

# autosome
auto=1:22
auto.cover=data[auto,2]

# mt cover
mt.cover=data['MT',2]

copy = 2*mt.cover/mean(auto.cover,na.rm = T)

res=c(sample,copy)
write.table(res,file=paste0(out_dir,'/',sample),row.names = F,col.names = F,quote = F)

