rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
file.names <- args[1]
###Provide full path to the results from SSB###
ttype<-basename(file.names)


###Provide full path to the list of genes to test
gene2test <- args[2]

#Read main file and gene file
df<-read.table(file.names, header=T, sep="\t", quote="", stringsAsFactors=F, blank.lines.skip=T, dec=".")
genes<-read.csv(gene2test,header=F,sep="\t",stringsAsFactors = TRUE)
df <- df[df$Hugo_symbol %in% genes$V1,]

##Using real mutations 
m<-sum(df$nonsilent_variant)
s<-sum(df$synonymous_variant)
M<-sum(df$nonsynsites_SSB)
S<-sum(df$synsites_SSB)

##Calculate conf interval using Katz method
p1<-(m/(M+1))
p2<-(s/(S+1))
globaldnds<-p1/p2
N1 <- M
N2 <- S

SE = sqrt( (1-p1)/(N1*p1) + (1-p2)/(N2*p2) )

finalLowCI = globaldnds * exp(-1.96*SE)
finalHighCI = globaldnds * exp(1.96*SE)

N<-m+s
df.global<-as.data.frame(cbind(globaldnds,finalLowCI,finalHighCI,N))
rownames(df.global)<-NULL

a<-df.global
b<-cbind(ttype,a)

df$pval_SSB.adj<-p.adjust(df$pval_SSB, method = "BH", n=length(df$Hugo_symbol))

write.table(df,file=paste(file.names,"geneset.txt",sep="."), quote=F,sep='\t',row.names=F)
write.table(b,file="",quote=F,row.names=F)
