rm(list=ls())
library('optparse')

option_list <- list(
    make_option(c("--in_file"), type="character", help="file containing mutation table", default="f:/Cornell/experiment/uk10k/uk10k/heteroplasmy/burden/test.mutation.table"),
    make_option(c("--in_vcf_file"), type="character", help="file containing vcf file", default="f:/Cornell/experiment/uk10k/uk10k/heteroplasmy/burden/QTL190044.vcf"),
    make_option(c("--region_file"), type="character", help="file containing mtDNA region", default="f:/Cornell/experiment/uk10k/uk10k/burden/gene.region.csv"),
    make_option(c("--min"), type="integer", help="min minor allele count", default=4),
    make_option(c("--out_file"), type="character", help="output file", default="f:/Cornell/experiment/uk10k/uk10k/heteroplasmy/burden/test.diversity")
    
)

opt <- parse_args(OptionParser(option_list=option_list))

data=read.table(opt$in_file,header=T,fill = T)

# diversity
Cal_H = function(x){
    x=as.numeric(as.character(x))
    i = x[2]; n=x[1]
    H=2*i*(n-i)/n/(n-1)
    return(H)
}

region=read.csv(opt$region_file)

coding.list=list()
for (i in 1:nrow(region)){
    gene.name=as.character(region[i,1])
    st=region[i,2]
    end=region[i,3]
    pos.list=st:end
    coding.list[[gene.name]]=pos.list
}

# control region
coding.list[['MT-CR']]=c(1:576,16024:16569)

length(coding.list)

# overall heterogenicity
data.select=data[data$coverage>100 & data$mutation>opt$min,]
l = sum(data$coverage>100)
Diversity=apply(data.select[,c(2,9)],1,Cal_H)
D=(sum(Diversity,na.rm = T)/l)

# region diversity
region.d=NULL
for (i in 1:length(coding.list)){
    data.select=data[data$coverage>100 & data$mutation>opt$min & data$position %in% coding.list[[i]],]
    if (nrow(data.select)>0){
        l = sum(data$coverage>100 & data$position %in% coding.list[[i]])
        Diversity=apply(data.select[,c(2,9)],1,Cal_H)
        D=(sum(Diversity,na.rm = T)/l)
        region.d=c(region.d,D)
    }else{
        region.d=c(region.d,0)
    }
    
}



sample_name=basename(opt$in_file)
sample_name=sub('.mutation.table','',sample_name)

write.table(t(c(sample_name,D,variant.persite)),file=opt$out_file,row.names = F,col.names = F,quote = F)



