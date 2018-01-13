##greenman chi-square testing
rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
file.names <- args[1]

#Read main file
df.genes<-read.table(file.names, header=T, sep="\t", quote="", stringsAsFactors=F, blank.lines.skip=T, dec=".") 

#pseudocount to add to keep neutrality =1
pseudna_SSB192<-df.genes$nonsynsites_SSB192/df.genes$synsites_SSB192

###Real data
pn_SSB192 <- (df.genes$nonsilent_variant+pseudna_SSB192)/df.genes$nonsynsites_SSB192
ps_SSB192 <- (df.genes$synonymous_variant+1)/df.genes$synsites_SSB192
df.genes$dnds_SSB192 <- pn_SSB192/ps_SSB192

##Using real mutations SSB192
m<-df.genes$nonsilent_variant
s<-df.genes$synonymous_variant
M<-df.genes$nonsynsites_SSB192
S<-df.genes$synsites_SSB192

#Define values of expected and observed and transform into a chi-square distribution
t=(m+s)
Ta=(M+S)
U=m-(t*(M/Ta))
V=t*(M*((Ta-M)/Ta^2))
testscore=U^2/V

df.genes$pval_SSB192=pchisq(testscore,df=1,lower.tail=FALSE)
df.genes<-df.genes[complete.cases(df.genes$pval_SSB192),]
df.genes$pval_SSB192.adj<-p.adjust(df.genes$pval_SSB192, method = "BH", n=length(rownames(df.genes)))

df.genes.sort<-df.genes[order(df.genes$pval_SSB192),]

write.table(df.genes.sort,file=paste(file.names,"res.txt",sep="_"), quote=F,sep='\t',row.names=F)
