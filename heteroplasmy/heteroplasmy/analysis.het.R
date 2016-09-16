rm(list=ls())
options(java.parameters = "-Xmx8000m")
library('optparse')
library(outliers)
library(ggplot2)
library(reshape)
library(plyr)
library(cowplot)

theme_set(theme_bw(base_size = 20))

option_list <- list(
    make_option(c("--in_file"), type="character", help="dir to input data", default="f:/Cornell/experiment/uk10k/uk10k/heteroplasmy/heteroplasmy/het.homo.annotate.txt"),
    make_option(c("--sample_name_file"), type="character", help="dir to sample name file", default="f:/Cornell/experiment/uk10k/uk10k/heteroplasmy/heteroplasmy/sampleID"),
    make_option(c("--hap_file"), type="character", help="dir to haplogroup file", default="f:/Cornell/experiment/uk10k/uk10k/heteroplasmy/heteroplasmy/result/heteroplasmy.haplogroups.hsd"),
    make_option(c("--disease_file"), type="character", help="dir to sample name file", default="f:/Cornell/experiment/uk10k/uk10k/heteroplasmy/heteroplasmy/disease.score"),
    make_option(c("--out_dir"), type="character", help="output file", default="f:/Cornell/experiment/uk10k/uk10k/heteroplasmy/heteroplasmy/result/")
    
)

opt <- parse_args(OptionParser(option_list=option_list))
#print(opt)

# in file
in_file=opt$in_file
# sampleID
sampleID_file=opt$sample_name_file
# haplogroup file
hap_file = opt$hap_file
# disease file
disease_file=opt$disease_file
# out_dir
out_dir=opt$out_dir

single.num=980
multi.num=733

sampleID=read.table(sampleID_file,header=F)[[1]]

data=read.table(in_file,header=T,sep='|',stringsAsFactors=F,quote = "")

# Basic statistics
# heteroplasmy distribution
ind.het=ddply(subset(data,tag=='H'),.(sampleID),nrow)
ind.homo=ddply(subset(data,tag=='M'),.(sampleID),nrow)

rownames(ind.het)=ind.het[,1]

ind.num=data.frame(sampleID=sampleID,het.num=rep(0,length(sampleID)),homo.num=rep(0,length(sampleID)))
rownames(ind.num)=ind.num$sampleID

ind.num[ind.het$sampleID,2]=ind.het$V1
ind.num[ind.homo$sampleID,3]=ind.homo$V1

ind.num=ind.num[-nrow(ind.num),]

# plot distribution
# heteroplasmy 
# individuals with many heteroplasmies may be containminated
# do outlier test by Q3 + 3IQR
out.cutoff=quantile(ind.num$het.num,probs = 0.75)+1.5*IQR(ind.num$het.num)
p=ggplot(ind.num,aes(het.num))+geom_histogram(binwidth = 1)+xlab("Number of Heteroplasmy")+ylab("Count of Individuals")+
    geom_vline(xintercept = out.cutoff,col="red")
ggsave(p,filename=paste0(out_dir,"/","heteroplasmy_dist_all_sample.jpeg"), width=8, height=5, dpi=300)

# prepare files for haplogrep
data.homo=data[data$tag=='M',]
data.het=data[data$tag=='H',]

haplogrep.het=NULL
haplogrep.homo=NULL
for (id in sampleID){
    homo.in=which(data.homo$sampleID==id)
    het.in=which(data.het$sampleID==id)
    
    if (length(homo.in)>0){
        variants=data.homo[homo.in,c('pos','altA')]
        variants=paste0(variants$pos,variants$altA)
        line=c(id,'1-16569','?',variants)
        line=paste(line,sep="",collapse = "\t")
        haplogrep.homo=c(haplogrep.homo,line)
    }else{
        line=c(id,'1-16569','?',"709G")
        line=paste(line,sep="",collapse = "\t")
        haplogrep.homo=c(haplogrep.homo,line)
    }
    
    if (length(het.in)>0){
        variants1=data.het[het.in,c('pos','MajorA')]
        variants1=paste0(variants1$pos,variants1$MajorA)
        line=c(paste0(id,"Ma"),'1-16569','?',variants,variants1)
        line=paste(line,sep="",collapse = "\t")
        haplogrep.het=c(haplogrep.het,line)
        
        variants2=data.het[het.in,c('pos','MinorA')]
        variants2=paste0(variants2$pos,variants2$MinorA)
        line=c(paste0(id,"Mi"),'1-16569','?',variants,variants2)
        line=paste(line,sep="",collapse = "\t")
        haplogrep.het=c(haplogrep.het,line)
    }
    

}

#write.table(haplogrep.homo,file=paste0(out_dir,'haplogrep.homoplasmy.txt'),row.names = F,col.names = F,quote = F)
#write.table(haplogrep.het,file=paste0(out_dir,'haplogrep.heteroplasmy.txt'),row.names = F,col.names = F,quote = F)

# minor allele has same haplogroup as major allele?
hap_data = read.table(hap_file,header = F, fill = T)
hap_data = hap_data[-1,c(1,3)]

diff_hap = NULL
for (i in seq(1,nrow(hap_data),2)){
    if (hap_data[i,2]!=hap_data[i+1,2]){
        sample_name = sub('Ma','',hap_data[i,1])
        diff_hap = c(diff_hap,sample_name)
    }
}

# samples with different haplogroup
length(diff_hap)

# heteroplasmy greater than cutoff 
outlier_sample = ind.num$sampleID[ind.num$het.num>out.cutoff]
length(outlier_sample)

#write.table(ind.num,file=paste0(out_dir,"het.number.txt"),row.names = F,quote = F)

# remove containmated samples
con_sample = intersect(diff_hap,outlier_sample)
keep_sample = setdiff(rownames(ind.num),con_sample)

write.table(con_sample,file=paste0(out_dir,'contaminated_sample.txt'),row.names = F,col.names = F,quote = F)

ind.num.rmcon=ind.num[keep_sample,]

# plot again
# heteroplasmy
p=ggplot(ind.num.rmcon,aes(het.num,fill=I("lightblue"),col=I("blue")))+geom_histogram(binwidth = 1)+xlab("Number of Heteroplasmy")+ylab("Count of Individuals")
ggsave(p,filename=paste0(out_dir,"/","heteroplasmy_dist_rm_containmination.jpeg"), width=8, height=5, dpi=300)
# homoplasmy
ggplot(ind.num.rmcon,aes(homo.num))+geom_histogram()

# remove containminated samples
rm.in.data=which(data$sampleID %in% con_sample)
data.rmcon=data[-rm.in.data,]

N = 1595

# heteroplasmy prevelence
het.data = subset(data.rmcon,tag=='H')
# total number of heteroplasmy
nrow(het.data)
length(unique(het.data$sampleID)); length(unique(het.data$sampleID))/N

# # of individuals haboring heteroplasmy
maf = seq(0.02,0.2,0.02)
individual.per = NULL
for (f in maf){
    het.temp = subset(het.data,Minorf>=f)
    nrow(het.temp)
    length(unique(het.temp$sampleID)); 
    ind.per=length(unique(het.temp$sampleID))/N
    individual.per = c(individual.per,ind.per)
}

p=qplot(maf,individual.per) + geom_line() + xlab("MAF cutoff") + ylab("Percentage of individuals with heteroplasmy")  
ggsave(p,filename=paste0(out_dir,"/","Percentage of individual with different cutoff.jpeg"), width=8, height=5, dpi=300)

# heteroplasmy > 5% 
het.0.05 = subset(het.data,Minorf>=0.05)
nrow(het.0.05)
length(unique(het.0.05$sampleID)); length(unique(het.0.05$sampleID))/N
# heteroplasmy > 10% 
het.0.1 = subset(het.data,Minorf>=0.1)
nrow(het.0.1)
length(unique(het.0.1$sampleID)); length(unique(het.0.1$sampleID))/N

# homoplasmy prevelence
homo.data = subset(data.rmcon,tag=='M')
nrow(homo.data)
length(unique(homo.data$sampleID)); length(unique(homo.data$sampleID))/N

#sharing of heteroplasmy and homoplasmy
#
het.share=ddply(subset(data.rmcon,tag=='H'),.(redetail),nrow)
dt=as.data.frame(table(cut(het.share$V1,breaks=c(0,1,2,3,10,max(het.share$V1)))))
colnames(dt)=c("B","A")
dt = dt[order(dt$A, decreasing = TRUE),]
myLabel = c("1","2","3","4-10",">10")  
myLabel = paste(myLabel, " (", round(dt$A / sum(dt$A) * 100, 2), "%)", sep = "")   

p1 = ggplot(dt, aes(x = "", y = A, fill = B)) +
    geom_bar(stat = "identity", width = 1) +    
    coord_polar(theta = "y") + 
    labs(x = "", y = "", title = "") + 
    theme(axis.ticks = element_blank()) + 
    theme(legend.title = element_blank(), legend.position = "right") + 
    scale_fill_discrete(breaks = dt$B, labels = myLabel) + 
    theme(axis.text.x = element_blank())+
    theme(panel.grid=element_blank()) +
    theme(panel.border=element_blank()) 

p1
ggsave(p1,filename=paste0(out_dir,"/","heteroplasmy_sharing.jpeg"), width=8, height=5, dpi=300)


#
het.share=ddply(subset(data.rmcon,tag=='M'),.(redetail),nrow)
dt=as.data.frame(table(cut(het.share$V1,breaks=c(0,1,2,3,10,max(het.share$V1)))))
colnames(dt)=c("B","A")
dt = dt[order(dt$A, decreasing = TRUE),]
myLabel = c("1","2","3","4-10",">10")  
myLabel = paste(myLabel, " (", round(dt$A / sum(dt$A) * 100, 2), "%)", sep = "")   

p1 = ggplot(dt, aes(x = "", y = A, fill = B)) +
    geom_bar(stat = "identity", width = 1) +    
    coord_polar(theta = "y") + 
    labs(x = "", y = "", title = "") + 
    theme(axis.ticks = element_blank()) + 
    theme(legend.title = element_blank(), legend.position = "right") + 
    scale_fill_discrete(breaks = dt$B, labels = myLabel) + 
    theme(axis.text.x = element_blank())+
    theme(panel.grid=element_blank()) +
    theme(panel.border=element_blank()) 

p1
ggsave(p1,filename=paste0(out_dir,"/","homoplasmy_sharing.jpeg"), width=8, height=5, dpi=300)

# How many heteroplasmy alleles are observed in homoplasmy?
het.unique=unique(het.data$redetail)
homo.unique=unique(homo.data$redetail)
length(intersect(het.unique,homo.unique))/length(het.unique)

#regions
data.rmcon$retype[is.na(data.rmcon$retype)]='Intergenic'
data.rmcon$mutype[data.rmcon$retype=="Intergenic"]='Intergenic'
het1=table(data.rmcon$mutype[data.rmcon$tag=='H'])[c('control_region','Intergenic','rRNA','tRNA')]
het2=table(data.rmcon$tRNA[data.rmcon$tag=='H'])[c('NS','SY')]

homo1=table(data.rmcon$mutype[data.rmcon$tag=='M'])[c('control_region','Intergenic','rRNA','tRNA')]
homo2=table(data.rmcon$tRNA[data.rmcon$tag=='M'])[c('NS','SY')]

heteroplasmy=c(het1,het2)/sum(c(het1,het2))
homoplasmy=c(homo1,homo2)/sum(c(homo1,homo2))

#region plot
region=c("Control Region","Intergenic","rRNA","tRNA","Nonsynonymous","Synonymous")

dat=data.frame(Percentage=c(heteroplasmy,homoplasmy),
               Group=c(rep("Heteroplasmy",6),rep("Homoplasmy",6)),
               region=c(region,region))
p=ggplot(dat,aes(x=Group,y=Percentage,fill=region))+
    geom_bar(stat="identity")+coord_flip()+ylab("")+xlab("")

ggsave(p,filename=paste0(out_dir,"/","region_distribution.jpeg"), width=8, height=5, dpi=300)

# gene location
het.gene=data.rmcon$retype[data.rmcon$tag=='H']
homo.gene=data.rmcon$retype[data.rmcon$tag=='M']

het.gene = (table(het.gene))/sum(table(het.gene))
homo.gene = (table(homo.gene))/sum(table(homo.gene))

het.gene=data.frame(freq=het.gene,tag=rep("Heteroplasmy",length(het.gene)));colnames(het.gene)=c("gene","frequency",'tag')
homo.gene=data.frame(freq=homo.gene,tag=rep("Homoplasmy",length(homo.gene)));colnames(homo.gene)=c("gene","frequency",'tag')
df = rbind(het.gene,homo.gene)

#gene.seq = homo.gene$gene
#gene.seq = c(as.character(gene.seq[c(1,2,8,10)]),as.character(gene.seq[-c(1,2,8,10)]))

#df$gene=factor(df$gene,levels = gene.seq)

p=ggplot(df, aes(x=gene,y=frequency, fill=tag)) + geom_bar(position="dodge",stat='identity') + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("")+ylab("Frequency in all variants")
ggsave(p,filename=paste0(out_dir,"/","variants.gene.distribution.jpeg"), width=12, height=6, dpi=300)

# Pathogenicity of heteroplasmy
# how many heteroplasmies are associatied with disease?
sum(!is.na(het.data$disease))
1-sum(is.na(het.data$disease))/nrow(het.data)

# how many individual harbor disease associated heteroplamy?
length(unique(het.data[!is.na(het.data$disease),]$sampleID))
length(unique(het.data[!is.na(het.data$disease),]$sampleID))/1595

# compare to homoplasmy
sum(!is.na(homo.data$disease))
1-sum(is.na(homo.data$disease))/nrow(homo.data)

# how many individual harbor disease associated heteroplamy?
length(unique(homo.data[!is.na(homo.data$disease),]$sampleID))
length(unique(homo.data[!is.na(homo.data$disease),]$sampleID))/1595


# Pathogenic score of heteroplasmy and homoplasmy sites
# All sites
# read disease associated variants
disease=read.table(disease_file,header=T,sep='|',quote="")
disease$Minorf=0

het.score=data.rmcon[data.rmcon$tag=='H',"phred"]
homo.score=data.rmcon[data.rmcon$tag=='M',"phred"]

# boxplot
df=data.frame(score=c(disease$phred,het.score,homo.score),
              tag=c(rep("Disease-assoicated sites",length(disease$phred)),rep("Heteroplasmy sites",length(het.score)),rep("Homoplasmy sites",length(homo.score))))
p=ggplot(data=df,aes(tag,score,color=tag))+geom_boxplot()+scale_x_discrete(breaks=NULL)+scale_color_discrete(name="")+
    ylab("Pathogenic Score")+theme(axis.text=element_text(size=20))+xlab("")
ggsave(p,filename=paste0(out_dir,"/","all.sites.pathogenic.boxplot.jpeg"), width=8, height=5, dpi=300)

# t test
t.test(het.score,homo.score)
t.test(het.score,disease$phred)

HM.all=data.rmcon[,c("pos","tag","phred","Minorf")]
HM.h=HM.all[HM.all$tag=='H',]
h.maxf=tapply(HM.h$Minorf,HM.h$pos,max)
h.ps=tapply(HM.h$phred,HM.h$pos,mean)
h.tag=rep('l',length(h.ps))
h.tag[h.maxf>0.1]='h'

# t test
t.test(h.ps[h.tag=='l'],h.ps[h.tag=='h'])
t.test(h.ps[h.tag=='l'],disease$phred)

HM.homo=HM.all[HM.all$tag=='M',]
HM.homo=unique(HM.homo)

plot.df=data.frame(Minorf=h.maxf,phred=h.ps,tag=h.tag)
ggplot(plot.df,aes(h.tag,h.ps))+geom_boxplot(notch = T)

ggplot(plot.df, aes(h.ps, colour = h.tag)) + stat_ecdf(size=2) +
    theme(text = element_text(size=30),legend.title=element_blank(),legend.position="top") + 
    xlab("Pathogenic score") + ylab("")+guides(col=guide_legend(ncol=2))

ALL=NULL
ALL=rbind(plot.df,HM.homo[,c('Minorf','phred','tag')],disease[,c('Minorf','phred','tag')])

ALL$tag2=NULL
ALL$tag2[ALL$tag=='M']='Homoplasmy'
ALL$tag2[ALL$tag=='D']='Disease associated mutation'
ALL$tag2[ALL$tag=="H"&ALL$Minorf<0.05]='Heteroplsamy 2%-10%'
ALL$tag2[ALL$Minorf>=0.1]='Heteroplasmy higher than 10%'

ALL$tag=factor(ALL$tag,levels=c("M",'h','l','D'))
p=ggplot(ALL, aes(phred, colour = tag)) + stat_ecdf(size=2) +
    theme(text = element_text(size=20),legend.title=element_blank(),legend.position="top") + 
    xlab("Pathogenic score") + ylab("")+guides(col=guide_legend(ncol=2))+
    scale_color_discrete(labels=c("Homoplasmy",'Heteroplasmy higher than 10%',"Heteroplsamy 2%-10%","Disease associated mutation"))
ggsave(p,filename=paste0(out_dir,"/","heteroplasmy_pathogenicity_cumulative.all.sites.jpeg"), width=8, height=5, dpi=300)

p=ggplot(ALL, aes(x=tag,y=phred,color=tag))+geom_boxplot(notch = T)+guides(fill=F)+ 
    scale_color_discrete(name="",labels=c("Homoplasmy",'Heteroplasmy higher than 10%',"Heteroplsamy 2%-10%","Disease associated mutation"))+
    ylab("Pathogenic Score")+theme(axis.text=element_text(size=20))+scale_x_discrete(breaks=NULL)+xlab("")
ggsave(p,filename=paste0(out_dir,"/","l_h_heteroplasmy.boxplot.all.sites.jpeg"), width=8, height=5, dpi=300)


# NS sites
disease.NS=disease[disease$tRNA=='NS',]

het.score=data.rmcon[which(data.rmcon$tag=='H' & data.rmcon$tRNA=='NS'),"phred"]
homo.score=data.rmcon[which(data.rmcon$tag=='M' & data.rmcon$tRNA=='NS'),"phred"]

# boxplot
df=data.frame(score=c(disease.NS$phred,het.score,homo.score),
              tag=c(rep("Disease-assoicated sites",length(disease.NS$phred)),rep("Heteroplasmy sites",length(het.score)),rep("Homoplasmy sites",length(homo.score))))
p=ggplot(data=df,aes(tag,score,color=tag))+geom_boxplot()+scale_x_discrete(breaks=NULL)+scale_color_discrete(name="")+
    ylab("Pathogenic Score")+theme(axis.text=element_text(size=20))+xlab("")
ggsave(p,filename=paste0(out_dir,"/","NS.sites.pathogenic.boxplot.jpeg"), width=8, height=5, dpi=300)

t.test(het.score,homo.score)
t.test(het.score,disease.NS$phred)

HM.all=data.rmcon[which(data.rmcon$tRNA=='NS'),c("pos","tag","phred","Minorf")]
HM.h=HM.all[HM.all$tag=='H',]
h.maxf=tapply(HM.h$Minorf,HM.h$pos,max)
h.ps=tapply(HM.h$phred,HM.h$pos,mean)
h.tag=rep('l',length(h.ps))
h.tag[h.maxf>0.1]='h'

t.test(h.ps[h.tag=='l'],h.ps[h.tag=='h'])

HM.homo=HM.all[HM.all$tag=='M',]
HM.homo=unique(HM.homo)

plot.df=data.frame(Minorf=h.maxf,phred=h.ps,tag=h.tag)
ggplot(plot.df,aes(h.tag,h.ps))+geom_boxplot(notch = T)

ggplot(plot.df, aes(h.ps, colour = h.tag)) + stat_ecdf(size=2) +
    theme(text = element_text(size=30),legend.title=element_blank(),legend.position="top") + 
    xlab("Pathogenic score") + ylab("")+guides(col=guide_legend(ncol=2))

ALL=NULL
ALL=rbind(plot.df,HM.homo[,c('Minorf','phred','tag')],disease.NS[,c('Minorf','phred','tag')])

ALL$tag2=NULL
ALL$tag2[ALL$tag=='M']='Homoplasmy'
ALL$tag2[ALL$tag=='D']='Disease associated mutation'
ALL$tag2[ALL$tag=="H"&ALL$Minorf<0.1]='Heteroplsamy 2%-10%'
ALL$tag2[ALL$Minorf>=0.1]='Heteroplasmy higher than 10%'

ALL$tag=factor(ALL$tag,levels=c("M",'h','l','D'))
p=ggplot(ALL, aes(phred, colour = tag)) + stat_ecdf(size=2) +
    theme(text = element_text(size=20),legend.title=element_blank(),legend.position="top") + 
    xlab("Pathogenic score") + ylab("")+guides(col=guide_legend(ncol=2))+
    scale_color_discrete(labels=c("Homoplasmy",'Heteroplasmy higher than 10%',"Heteroplsamy 2%-10%","Disease associated mutation"))
ggsave(p,filename=paste0(out_dir,"/","heteroplasmy_pathogenicity_cumulative.NS.jpeg"), width=8, height=5, dpi=300)


p=ggplot(ALL, aes(x=tag,y=phred,color=tag))+geom_boxplot(notch = T)+guides(fill=F)+ 
    scale_color_discrete(name="",labels=c("Homoplasmy",'Heteroplasmy higher than 10%',"Heteroplsamy 2%-10%","Disease associated mutation"))+
    ylab("Pathogenic Score")+theme(axis.text=element_text(size=20))+scale_x_discrete(breaks=NULL)+xlab("")
ggsave(p,filename=paste0(out_dir,"/","l_h_heteroplasmy.boxplot.NS.sites.jpeg"), width=8, height=5, dpi=300)

