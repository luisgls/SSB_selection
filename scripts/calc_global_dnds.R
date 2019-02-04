rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
file.names <- args[1]

#Read main file
df<-read.table(file.names, header=T, sep="\t", quote="", stringsAsFactors=F, blank.lines.skip=T, dec=".")

##Using real mutations SSB192
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

write.table(df.global,file=paste(file.names,"globaldnds.txt",sep="."), quote=F,sep='\t',row.names=F)

