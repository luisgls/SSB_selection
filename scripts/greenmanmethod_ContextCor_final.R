##greenman chi-square testing
rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
file.names <- args[1]

#Read main file
df.genes<-read.table(file.names, header=T, sep="\t", quote="", stringsAsFactors=F, blank.lines.skip=T, dec=".") 

#pseudocount to add to keep neutrality =1
pseudna<-df.genes$nonsynsites_SSB/df.genes$synsites_SSB

###Real data
pn <- (df.genes$nonsilent_variant+pseudna)/df.genes$nonsynsites_SSB
ps <- (df.genes$synonymous_variant+1)/df.genes$synsites_SSB
df.genes$dnds_SSB <- pn/ps

##Using real mutations SSB192
m<-df.genes$nonsilent_variant
s<-df.genes$synonymous_variant
M<-df.genes$nonsynsites_SSB
S<-df.genes$synsites_SSB

#Define values of expected and observed and transform into a chi-square distribution
t=(m+s)
Ta=(M+S)
U=m-(t*(M/Ta))
V=t*(M*((Ta-M)/Ta^2))
testscore=U^2/V

df.genes$pval_SSB=pchisq(testscore,df=1,lower.tail=FALSE)
df.genes<-df.genes[complete.cases(df.genes$pval_SSB),]
df.genes$pval_SSB.adj<-p.adjust(df.genes$pval_SSB, method = "BH", n=length(rownames(df.genes)))

df.genes.sort<-df.genes[order(df.genes$pval_SSB),]



write.table(df.genes.sort,file=paste(file.names,"res.txt",sep="_"), quote=F,sep='\t',row.names=F)
