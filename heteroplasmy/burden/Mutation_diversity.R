rm(list=ls())
library('optparse')

option_list <- list(
    make_option(c("--in_file"), type="character", help="file containing mutation table", default="f:/Cornell/experiment/uk10k/uk10k/heteroplasmy/burden/test.mutation.table"),
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

#16,024¨C576
select.pos=577:16023
data.select=data[data$position %in% select.pos,]
l=nrow(data.select)
data.select=data.select[data.select$coverage>100 & data.select$mutation>1,]
Diversity=apply(data.select[,c(2,8)],1,Cal_H)
D=(sum(Diversity,na.rm = T)/l)
write.table(D,file=opt$out_file,row.names = F,col.names = F,quote = F)

