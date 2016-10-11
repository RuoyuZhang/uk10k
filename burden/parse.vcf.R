
var_file=read.table('f:/Cornell/experiment/uk10k/uk10k/burden/all.sample.vcf',header=T,stringsAsFactors = F)
header.names=colnames(var_file)

x=var_file[5,]
x=var_file[1,]

parse.single=function(x){
    x.split=strsplit(as.character(x),split=':')[[1]]
    if (x.split[1]=='.'){
        return(NA)
    }else{
        genodep=as.numeric(strsplit(x.split[2],',')[[1]])
        genotype=rep(0,length(genodep))
        genotype[which.max(genodep)]=1
        return(genotype[-1])
    }
    
}

parse.line=function(x){
    # variants
    varinfo=strsplit(as.character(x[5]),split = ',')[[1]]
    var_head=unlist(lapply(varinfo,function(y){paste0(x[4],x[2],y)}))
    
    # variants info for each individual
    xsample=x[-(1:9)]
    x.matrix=do.call(rbind,lapply(xsample,parse.single))
    colnames(x.matrix)=var_head
    return(x.matrix)
}

m=do.call(cbind,apply(var_file,1,parse.line))

# only snp
variant_header=colnames(m)
variant_header=unlist(lapply(variant_header,function(x){gsub(' ','',x)}))
colnames(m)=variant_header

header.mat=do.call(rbind,lapply(variant_header,function(x){unlist(strsplit(x,split = '[0-9]+',fixed=F))}))

# select only snp
snp.in=(nchar(header.mat[,2]) ==1 & nchar(header.mat[,1])==1)
m.snp=m[,snp.in]

out_dir='f:/Cornell/experiment/uk10k/uk10k/burden/'

write.table(m,file=paste0(out_dir,"all.variants.table"),quote = F)
write.table(m.snp,file=paste0(out_dir,"snp.variants.table"),quote = F)

# for haplogrep
parse.hap=function(x,sampleID){
    x=x[-grep('\\*',names(x),perl=T)]
    vars=sub('[ATGC]','',names(x[which(x>0)]))
    if (length(vars)==0){
        vars=c('16A')
    }
    line=c(sampleID,'1-16569','?',vars)
    return(line)
}

haplogrep=NULL
for(i in 1:nrow(m.snp)){
    line = parse.hap(m.snp[i,],rownames(m.snp)[i])
    line=paste(line,sep="",collapse = "\t")
    
    haplogrep=c(haplogrep,line)
}

write.table(haplogrep,file=paste0(out_dir,'haplogrep.vcf.hsd'),row.names = F,col.names = F,quote = F)
