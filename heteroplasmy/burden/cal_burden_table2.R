rm(list=ls())
library('optparse')

option_list <- list(
    make_option(c("--in_file"), type="character", help="file containing GATK base count", default="f:/Cornell/experiment/uk10k/uk10k/heteroplasmy/burden/UK10K_ALS5086237.raw"),
    make_option(c("--out_file"), type="character", help="output file", default="f:/Cornell/experiment/uk10k/uk10k/heteroplasmy/burden/test.mutation.table")
    
)

opt <- parse_args(OptionParser(option_list=option_list))

data=read.table(opt$in_file,header=F,fill = T)
data=data[data$V1=='chrM',]



# consense
count_base=function(x){
    counts=as.numeric(unlist(lapply(x,function(y) {unlist(strsplit(as.character(y),split=':'))[2]})))
    return(counts)
}

basecount=t(apply(data[,5:8],1,count_base))

count_mutation=function(x){
    bases=c('A','T','G','C')
    mutation=sort(x,decreasing = T)[2]
    consence=bases[which.max(x)]
    ts=0;tv=0
    second=NA
    if (mutation > 0){
        second = bases[order(x,decreasing = T)][2]
        if (paste0(consence,second) %in% c('AC','CA','TC','CT')){
            ts=ts+x[order(x,decreasing = T)][2]
        }else{
            tv=tv+x[order(x,decreasing = T)][2]
        }
    }
    return(c(consence,second,mutation,ts,tv))
}

mutations=t(apply(basecount,1,count_mutation))

data.new=cbind(data[,c(2,4)],basecount,mutations)

colnames(data.new)=c('position','coverage','A','T','G','C','consensus','second','mutation','ts','tv')
#data.new$position=sub('chrM:','',data.new$position)
write.table(data.new,opt$out_file,col.names = T,row.names = F,quote = F)


