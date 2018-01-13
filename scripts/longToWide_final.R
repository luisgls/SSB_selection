rm(list=ls())
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)

#read filename and column number to calculate
file.name <-  args[1] 

file.name2 <- args[2]

df<-read.table(file.name,header=F,sep="\t",quote="",stringsAsFactors=F, blank.lines.skip=T)

data_wide<- spread(df,V2,V3)

colnames(data_wide)<-c("EnsemblID","nonsilent_variant","synonymous_variant")

data_wide[is.na(data_wide)] <- 0

df2<-read.table(file.name2,header=F,sep="\t",quote="",stringsAsFactors=F, blank.lines.skip=T)

colnames(df2) <- c("ensemblID","nonsyn_sites","syn_sites","nonsynsites_SSB","synsites_SSB")

df.tmp<-merge(data_wide,df2,by.x="EnsemblID",by.y="ensemblID")

write.table(df.tmp,file=paste(as.character(file.name),".mgd",sep=""),sep="\t",quote=F,row.names=F)
