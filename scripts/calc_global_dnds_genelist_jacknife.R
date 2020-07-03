rm(list=ls())
#library(tidyverse)
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

get_glob_dnds<-function(df){
##Using real mutations 
m<-sum(df$nonsilent_variant)
s<-sum(df$synonymous_variant)
M<-sum(df$nonsynsites_SSB)
S<-sum(df$synsites_SSB)

 ##If synonymous are zero add pseudo count that keeps the neutrality balance
 if(s==0){
   N<-m+s
   s = 1
   m = m + (M/S)
 }else{
   N<-m+s
 }
# 
# ##Calculate conf interval using Katz method
 p1<-(m/(M+1))
 p2<-(s/(S+1))
 globaldnds<-p1/p2
# N1 <- M
# N2 <- S
# 
# SE = sqrt( (1-p1)/(N1*p1) + (1-p2)/(N2*p2) )
# 
# finalLowCI = globaldnds * exp(-1.96*SE)
# finalHighCI = globaldnds * exp(1.96*SE)
# 
# #Print results
# high = finalHighCI
# low = finalLowCI
# val = globaldnds
# SE = (log(high)-log(low) )/ (2 * 1.96)
# EST = log(val)
# z = EST/SE
# PVAL1=exp(-0.717*z - 0.416*z^(2))
# PVAL2=exp(-0.717*-z - 0.416*-z^(2))
# 
# if(PVAL2 > 0 & PVAL2 <= 1) {
#   PVAL = PVAL2
# } else {
#   PVAL = PVAL1
# }

return(globaldnds)
}

##Get all gene list dnds
all_dnds<-get_glob_dnds(df)

##Get gene list
genelist<-genes$V1

for(i in 1:length(genelist)) {
  df <- df[df$Hugo_symbol %in% genes$V1,]
  tmpgene<-as.character(genelist[i])
  # remove ith element from list
  genelist_tmp<-genelist[-i] 
  tmpdf <- df[df$Hugo_symbol %in% genelist_tmp,]

  tmp_dnds<-get_glob_dnds(tmpdf)
  
  caca<-as.data.frame(cbind(tmpgene,tmp_dnds,all_dnds,all_dnds-tmp_dnds))
  
  colnames(caca) <- NULL
  print(caca,quote=F,sep='\t',row.names=F,col.names=F,max=4)
  
  #cat(caca)
  }




#df.global<-as.data.frame(cbind(globaldnds,finalLowCI,finalHighCI,N,PVAL))
#rownames(df.global)<-NULL

#a<-df.global
#b<-cbind(ttype,a)

#df$pval_SSB.adj<-p.adjust(df$pval_SSB, method = "BH")

#write.table(df,file=paste(file.names,"geneset.txt",sep="."), quote=F,sep='\t',row.names=F)
#write.table(b,file="",quote=F,row.names=F)
