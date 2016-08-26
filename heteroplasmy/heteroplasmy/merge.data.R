rm(list=ls())
options(java.parameters = "-Xmx8000m")
library('optparse')

option_list <- list(
    make_option(c("--in_dir"), type="character", help="dir to input data", default="/home/fs01/rz253/project/uk10k/twins/heteroplasmy/het/"),
    make_option(c("--out_file"), type="character", help="output file", default="/home/fs01/rz253/project/uk10k/analysis/heteroplasmy/het.homo.txt")
    
)

opt <- parse_args(OptionParser(option_list=option_list))
#print(opt)

# file dir
in_dir=opt$in_dir
# out_dir
out_file=opt$out_file

# list files
multi.het.files=list.files(path = paste0(in_dir,"multi"), pattern = ".heteroplasmy",full.names = T)
multi.homo.files=list.files(path = paste0(in_dir,"multi"), pattern = ".homoplasmy",full.names = T)
single.het.files=list.files(path = paste0(in_dir,"single"), pattern = ".heteroplasmy",full.names = T)
single.homo.files=list.files(path = paste0(in_dir,"single"), pattern = ".homoplasmy",full.names = T)

single.het=NULL
for (file in single.het.files){
    if (file.size(file)>0){
        data=read.table(file,header=F,stringsAsFactors=F,fill=T,colClasses = c("character"))
        sampleID=basename(file)
        sampleID=gsub(".heteroplasmy","",sampleID)
        data$sampleID=sampleID
        data$type="heteroplasmy"
        data$lane="single"
        single.het=rbind(single.het,data)
    }
}

single.het=NULL
for (file in single.het.files){
    if (file.size(file)>0){
        data=read.table(file,header=F,stringsAsFactors=F,fill=T,colClasses = c("character"))
        sampleID=basename(file)
        sampleID=gsub(".heteroplasmy","",sampleID)
        data$sampleID=sampleID
        data$type="heteroplasmy"
        data$lane="single"
        single.het=rbind(single.het,data)
    }
}

multi.het=NULL
for (file in multi.het.files){
    if (file.size(file)>0){
        data=read.table(file,header=F,stringsAsFactors=F,fill=T,colClasses = c("character"))
        sampleID=basename(file)
        sampleID=gsub(".heteroplasmy","",sampleID)
        data$sampleID=sampleID
        data$type="heteroplasmy"
        data$lane="multi"
        multi.het=rbind(multi.het,data)
    }
}

single.homo=NULL
for (file in single.homo.files){
    if (file.size(file)>0){
        data=read.table(file,header=F,stringsAsFactors=F,fill=T,colClasses = c("character"))
        sampleID=basename(file)
        sampleID=gsub(".homoplasmy","",sampleID)
        data$V14="no"
        data$V15=0
        data$sampleID=sampleID
        data$type="homoplasmy"
        data$lane="single"
        single.homo=rbind(single.homo,data)
    }
}

multi.homo=NULL
for (file in multi.homo.files){
    if (file.size(file)>0){
        data=read.table(file,header=F,stringsAsFactors=F,fill=T,colClasses = c("character"))
        sampleID=basename(file)
        sampleID=gsub(".homoplasmy","",sampleID)
        data$V14="no"
        data$V15=0
        data$sampleID=sampleID
        data$type="homoplasmy"
        data$lane="multi"
        multi.homo=rbind(multi.homo,data)
    }
}

all.data=rbind(single.het,multi.het,single.homo,multi.homo)

colnames(all.data)=c("chr","pos","ref","cover","A","T","G","C","minorF","mle","tag","A1","A1f","A2","A2f","sampleID","type","lane")

# drop NA 
na.in=which((all.data$A1f)=="")
all.data=all.data[-na.in,]
all.data$A1f=as.numeric(all.data$A1f)
all.data$A2f=as.numeric(all.data$A2f)

for (i in 1:nrow(all.data)){
    if (all.data[i,"A1f"]>=all.data[i,"A2f"]){
        all.data[i,"MajorA"]=all.data[i,"A1"]
        all.data[i,"Majorf"]=all.data[i,"A1f"]
        all.data[i,"MinorA"]=all.data[i,"A2"]
        all.data[i,"Minorf"]=all.data[i,"A2f"]
        
    }else{
        all.data[i,"MajorA"]=all.data[i,"A2"]
        all.data[i,"Majorf"]=all.data[i,"A2f"]
        all.data[i,"MinorA"]=all.data[i,"A1"]
        all.data[i,"Minorf"]=all.data[i,"A1f"]
    }
    
    if (all.data[i,"A1"]==all.data[i,"ref"]){
        all.data[i,"reff"]=all.data[i,"A1f"]
        all.data[i,"altA"]=all.data[i,"A2"]
        all.data[i,"altf"]=all.data[i,"A2f"]
    }else{
        all.data[i,"reff"]=all.data[i,"A2f"]
        all.data[i,"altA"]=all.data[i,"A1"]
        all.data[i,"altf"]=all.data[i,"A1f"]
    }
}


write.table(all.data,file=out_file,row.names = F,col.names = T,quote = F,sep=';')

