#This is the work I am doing to try to streamline the use of JASPAR for determination of regulatory elements.  If someone wants to use this
#keep in mind that TAIR has bulk downloads of promoter sequences available and biomart can streamline that even further (depends on your tech saavy)
#This pipeline was made function by Colin Kruse on 6/28/16.  Next step if I don't do it, would be to replace the manual entry of seq1, seq1, etc.
#to incorporate biomart, so that a list of locus IDs could be used in its place.  See my pipeline for HumanSNP analysis for a basis for this.  

library(TFBSTools)
library(JASPAR)
library(JASPAR2018)
library(Biostrings)
library(biomaRt)
library(seqinr)
#library(doBy)

#setwd("WHERE_YOU_WANT_FILES_TO_BE_LOADED/SAVE_TO")
setwd("~/Documents/GeneLab Hack-a-thon/")

marts=listMarts(host="plants.ensembl.org")
datasets=listDatasets(useMart(biomart="plants_mart",host="plants.ensembl.org"))
ensemblArab=useMart(biomart="plants_mart",dataset = "athaliana_eg_gene",host="plants.ensembl.org")
Attributes<-listAttributes(ensemblArab)
#Import the gene description, gene start position and the chromosome of each gene in the TAIR10 genome.  
#We wont use the end positon in the default mode of this script, but I've included it in case anyone wants to evaluate downstream regions as well.  
#As an example and thought, this could be interesting if different polyadenlyation sequences or other 5' features impact degradation or translation
Annotations<-getBM(attributes=c('ensembl_gene_id','description','start_position','chromosome_name','end_position'),mart=ensemblArab)
row.names(Annotations)<-Annotations$ensembl_gene_id
#import the chromosome sequences for TAIR10
TAIR10<-read.fasta(file = "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa")

genome<-as.data.frame(getSequence(TAIR10["1"]))

#Now with both the chromosomal sequences and each genes start position, we can obtain the upstream region for promoter element analysis

#The genes examined could be any list of genes (even row names or a column from a data frame).
LISTOFGENES<-c("AT1G08810","AT1G74310")

#Bulk data retrieval for upstream sequences
#https://www.arabidopsis.org/tools/bulk/sequences/index.jsp
#Downloading the TAIR10 Genome
#ftp://ftp.ensemblgenomes.org/pub/plants/release-47/fasta/arabidopsis_thaliana/dna/

Annotations["AT1G74310","start_position"]

for (genes in unique(LISTOFGENES)){
  endpos<-as.numeric(Annotations[genes,"start_position"])
  chr<-Annotations[genes,"chromosome_name"]
  distance<-5000
  startpos<-endpos-5000
  test<-as.data.frame(paste(unlist(getSequence(TAIR10[chr]))[startpos:endpos],collapse = ""))
  row.names(test)<-genes
  colnames(test)<-"Promot"
  if(exists("UpstreamSeqs")){
    UpstreamSeqs<-rbind(test,UpstreamSeqs)
  }else{
    UpstreamSeqs<-test
  }
}
DNAStrings<-DNAStringSet(UpstreamSeqs$PromoterSeq)
names(DNAStrings)<-row.names(UpstreamSeqs)


optshs <- list()
optshs[["tax_group"]] <- "plants"
#optshs[["species"]] <- "Arabidopsis thaliana"
#optshs[["type"]] <- c("SELEX","ChIP-seq","ChIP-chip")
#optshs[["all_versions"]] <- TRUE
AthalSiteList <- getMatrixSet(JASPAR2018, optshs)
AthalPWM<-toPWM(AthalSiteList)

AthalsitesetList <- searchSeq(AthalPWM, DNAStrings, strand="*", min.score="90%")

Athalupstreambindingfeatures<-as(AthalsitesetList, "data.frame")

TFBSCounts=as.data.frame(table(Athalupstreambindingfeatures$TF))
