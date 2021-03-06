rm(list=ls())
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
file.names <- args[1]
###Provide full path to the results from SSB###
ttype<-basename(file.names)


###Provide full path to the list of genes to test
gene2test <- args[2]

#Read main file and gene file
df<-read.table(file.names, header=T, sep="\t", quote="", stringsAsFactors=F, blank.lines.skip=T, dec=".")
df<-read.table("~/Dropbox (Personal)/PI_position/Projects/NegativeSelection/NegativeSelection2021TCGA/TCGA_COAD.mutect.finalSSB.txt", header=T, sep="\t", quote="", stringsAsFactors=F, blank.lines.skip=T, dec=".")
genes<-read.csv(gene2test,header=F,sep="\t",stringsAsFactors = TRUE)
genes<-read.csv("~/Dropbox (Personal)/Postdoc/Projects/EPICC/results_TCGA/IntOGen-DriverGenes_COREAD.txt",header=F,sep="\t",stringsAsFactors = TRUE)
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

LowCI = globaldnds * exp(-1.96*SE)
HighCI = globaldnds * exp(1.96*SE)

N<-m+s
df.global<-as.data.frame(cbind(globaldnds,LowCI,HighCI,N))
rownames(df.global)<-NULL

a<-df.global
b<-cbind(ttype,a)

df$pval_SSB.adj<-p.adjust(df$pval_SSB, method = "BH", n=length(df$Hugo_symbol))

get_lowCI=function(m,s,M,S){
  p1<-(m/(M+1))
  p2<-(s/(S+1))
  globaldnds<-p1/p2
  N1 <- M
  N2 <- S
  
  SE = sqrt( (1-p1)/(N1*p1) + (1-p2)/(N2*p2) )
  
  LowCI = globaldnds * exp(-1.96*SE)
  HighCI = globaldnds * exp(1.96*SE)
  return(LowCI)
}
get_highCI=function(m,s,M,S){
  p1<-(m/(M+1))
  p2<-(s/(S+1))
  globaldnds<-p1/p2
  N1 <- M
  N2 <- S
  
  SE = sqrt( (1-p1)/(N1*p1) + (1-p2)/(N2*p2) )
  
  LowCI = globaldnds * exp(-1.96*SE)
  HighCI = globaldnds * exp(1.96*SE)
  return(HighCI)
}

df2<-df %>% rowwise() %>% dplyr::mutate(lowCI=get_lowCI(nonsilent_variant,synonymous_variant,nonsynsites_SSB,synsites_SSB))

df3<-df2 %>% rowwise() %>% dplyr::mutate(highCI=get_highCI(nonsilent_variant,synonymous_variant,nonsynsites_SSB,synsites_SSB))

write.table(df3,file=paste(file.names,"geneset.txt",sep="."), quote=F,sep='\t',row.names=F)
write.table(b,file=paste(file.names,"geneset.globaldNdS.txt",sep="."),quote=F,sep='\t',row.names=F)
