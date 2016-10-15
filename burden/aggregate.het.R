rm(list=ls())
library('optparse')
library('reshape2')

option_list <- list(
    make_option(c("--in_dir"), type="character", help="dir containing base count table", default="f:/Cornell/experiment/uk10k/uk10k/burden/mutation.table/"),
    make_option(c("--out_file"), type="character", help="output file", default="f:/Cornell/experiment/uk10k/uk10k/burden/mutation.table/test.out.txt"),
    make_option(c("--min"), type="integer", help="min minor allele cover", default=4)
    
)

opt <- parse_args(OptionParser(option_list=option_list))

files=list.files(path=opt$in_dir,pattern='.table')

Cal_H = function(x){
    x=as.numeric(as.character(x))
    i = x[2]; n=x[1]
    H=2*i*(n-i)/n/(n-1)
    return(H)
}

out.table=NULL
for (file in files){
    sample=basename(file)
    sample=sub('.mutation.table','',sample)
    filepath=paste0(opt$in_dir,file)
    data=read.table(filepath,header=T)
    data=data[data$coverage>100 & data$mutation>opt$min,]
    if (nrow(data)>0){
        Diversity=apply(data[,c(2,9)],1,Cal_H)
        out.sample=cbind(as.character(sample),paste0(data$consensus,data$pos,data$second),Diversity)
        out.table=rbind(out.table,out.sample)
    }else{
        out.table=rbind(out.table,c(sample,'A10306C',0))
    }
    print(sample)
}

colnames(out.table)=c('sampleID','variant','het')
out.table=as.data.frame(out.table)
data=dcast(out.table,sampleID~variant,value.var = "het",fill=0)
data[is.na(data)]=0

write.table(data,file=opt$out_file,row.names = F,col.names = T,quote = F)

