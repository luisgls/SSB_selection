rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
file.names <- args[1]
#file.names <- "~/Dropbox (Personal)/Projects/NegativeSelection/negDriver_tool/SSB_NatSel/results/pancan33.mutect.ssb192.all.finalSSB.txt"

gene2test <- args[2]
#gene2test <- "~/Dropbox (Personal)/Projects/NegativeSelection/Martincorena_dnds/Test_callsets/hg19/list_of_NOT_cancer_genes.txt"

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
colnames(df.global)<-c("globaldnds","low_CI","high_CI","Total_Muts")

a<-df.global
print(a)
#cat("\n")

#write.table(df.global,file=paste(file.names,"geneset_globaldnds.txt",sep="."), quote=F,sep='\t',row.names=F)

