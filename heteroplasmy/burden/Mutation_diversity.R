rm(list=ls())
library('optparse')

option_list <- list(
    make_option(c("--in_file"), type="character", help="file containing mutation table", default="f:/Cornell/experiment/uk10k/uk10k/heteroplasmy/burden/test.mutation.table"),
    make_option(c("--in_vcf_file"), type="character", help="file containing vcf file", default="f:/Cornell/experiment/uk10k/uk10k/heteroplasmy/burden/QTL190044.vcf"),
    make_option(c("--common_file"), type="character", help="file containing common variants", default="f:/Cornell/experiment/uk10k/uk10k/heteroplasmy/burden/common.variants.txt"),
    make_option(c("--out_file"), type="character", help="output file", default="f:/Cornell/experiment/uk10k/uk10k/heteroplasmy/burden/test.diversity")
    
)

opt <- parse_args(OptionParser(option_list=option_list))

# common variants
common=read.table(opt$common_file,header=F)

data=read.table(opt$in_file,header=T,fill = T)

# diversity
Cal_H = function(x){
    x=as.numeric(as.character(x))
    i = x[2]; n=x[1]
    H=2*i*(n-i)/n/(n-1)
    return(H)
}

#16,024¨C576
select.pos=577:16023
data.select=data[data$position %in% select.pos,]
l=nrow(data.select)
l2=sum(data.select$coverage>=4)


data.select=data.select[data.select$coverage>100 & data.select$mutation>1 & !(paste0(data.select$position,data.select$second) %in% common),]
Diversity=apply(data.select[,c(2,9)],1,Cal_H)
D=(sum(Diversity,na.rm = T)/l)
sample_name=basename(opt$in_file)
sample_name=sub('.mutation.table','',sample_name)

# homoplasmy load
variant.persite=0
f=try(read.table(opt$in_vcf_file,header=F,fill=T))
if (class(f)!='try-error'){
    homo=read.table(opt$in_vcf_file,header=F,fill=T)
    homo=homo[homo$V7=="PASS",]       
    homo=homo[nchar(as.character(homo$V5))<2,]
    mutation=paste0(homo$V2,homo$V5)
    detail=do.call(rbind,strsplit(as.character(homo$V10),split=':'))[,3]
    homo.data=data.frame(pos=homo$V2,mutation,detail,stringsAsFactors = F)
    homo.data=homo.data[!homo.data$mutation %in% common & homo.data$pos %in% select.pos & homo.data$detail > 0.7,]
    variant.persite=nrow(homo.data)/l2
}



write.table(t(c(sample_name,D,variant.persite)),file=opt$out_file,row.names = F,col.names = F,quote = F)



