library(plyr)

setwd("PATH/TO/Hack-a-thon_2020/DIRECTORY")
Counts<-read.csv(file = "AWG_MAtrix_rnaseq_microarray_cropped_V5.csv",stringsAsFactors = FALSE,header = TRUE)
Metadata<-read.csv(file = "AWG_MAtrix_rnaseq_microarray_cropped_V5_MetaData_V3.csv",stringsAsFactors = FALSE,header = TRUE)
tMetadata<-as.data.frame(t(Metadata))
colnames(tMetadata)<-Metadata$Sample.Name
#Remove first row (Redundant)
tMetadata<-tMetadata[-1,]

CountsAndMeta<-join(Counts,Metadata,type="full")



sort(as.numeric(unique(tMetadata$GLDS)))
unique(tMetadata$Genotype)
unique(tMetadata$Ecotype)
