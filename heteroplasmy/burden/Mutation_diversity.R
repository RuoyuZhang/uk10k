rm(list=ls())
library('optparse')

option_list <- list(
    make_option(c("--in_file"), type="character", help="file containing mutation table", default="f:/Cornell/experiment/uk10k/uk10k/burden/QTL_RP333139.mutation.table"),
    make_option(c("--ns_file"), type="character", help="file containing vcf file", default="f:/Cornell/experiment/uk10k/uk10k/burden/pathogenic_score.csv"),
    make_option(c("--region_file"), type="character", help="file containing mtDNA region", default="f:/Cornell/experiment/uk10k/uk10k/burden/gene.region.csv"),
    make_option(c("--min"), type="integer", help="min minor allele count", default=4),
    make_option(c("--out_file"), type="character", help="output file", default="f:/Cornell/experiment/uk10k/uk10k/heteroplasmy/burden/test.diversity"),
    make_option(c("--out_file2"), type="character", help="output file2", default="f:/Cornell/experiment/uk10k/uk10k/heteroplasmy/burden/test.ns.diversity")
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
data.select=data[data$coverage>100 & data$mutation>=opt$min,]
l = sum(data$coverage>100)
Diversity=apply(data.select[,c(2,9)],1,Cal_H)
D1=(sum(Diversity,na.rm = T)/l)

# region diversity
region.d=NULL
for (i in 1:length(coding.list)){
    data.select=data[data$coverage>100 & data$mutation>=opt$min & data$position %in% coding.list[[i]],]
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

write.table(t(c(sample_name,D1,region.d)),file=opt$out_file,row.names = F,col.names = F,quote = F)


# only consider ns
ns.file=read.csv(opt$ns_file,header=F)
n.change=do.call(rbind,(lapply(ns.file$V2,function(x){unlist(strsplit(as.character(x),split = '>'))})))
ns.site=paste0(n.change[,1],ns.file$V1,n.change[,2])

het.name1=paste0(data$consensus,data$position,data$second)
het.name2=paste0(data$second,data$position,data$consensus)
ns.index=((het.name1 %in% ns.site) | (het.name2 %in% ns.site))

# overall heterogenicity
data.select=data[ns.index,]
data.select=data.select[data.select$coverage>100 & data.select$mutation>=opt$min,]

l = sum(data$coverage>100)
Diversity=apply(data.select[,c(2,9)],1,Cal_H)
D1=(sum(Diversity,na.rm = T)/l)

# region diversity
region.d=NULL
for (i in 1:length(coding.list)){
    data.select2=data.select[data.select$coverage>100 & data.select$mutation>opt$min & data.select$position %in% coding.list[[i]],]
    if (nrow(data.select)>0){
        l = sum(data$coverage>100 & data$position %in% coding.list[[i]])
        Diversity=apply(data.select2[,c(2,9)],1,Cal_H)
        D=(sum(Diversity,na.rm = T)/l)
        region.d=c(region.d,D)
    }else{
        region.d=c(region.d,0)
    }
    
}

sample_name=basename(opt$in_file)
sample_name=sub('.mutation.table','',sample_name)

write.table(t(c(sample_name,D1,region.d)),file=opt$out_file2,row.names = F,col.names = F,quote = F)
