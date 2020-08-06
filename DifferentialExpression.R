suppressMessages(require(edgeR))
suppressMessages(require(ggplot2))
suppressMessages(require(grid))
suppressMessages(require(limma))
suppressMessages(require(NMF))
suppressMessages(require(plyr))
suppressMessages(require(RColorBrewer))
suppressMessages(require(reshape2))
suppressMessages(require(scales))
suppressMessages(require(gridExtra))
suppressMessages(require(dplyr))
suppressMessages(require(pheatmap))
suppressMessages(require(tidyr))
suppressMessages(require(cowplot))
suppressMessages(require(gplots))
suppressMessages(require(WGCNA))
suppressMessages(require(ParamHelpers))
suppressMessages(require(S4Vectors))
suppressMessages(require(IRanges))
suppressMessages(require(AnnotationDbi))
source("~/Desktop/AndrewRNAseq_Pipeline/GO_BP_Functions.R")
options(stringsAsFactors = F)
allowWGCNAThreads()
suppressMessages(require(genefilter))
#***WHERE WE USE THE ARGUMENTS TAKEN IN THE MAIN SCRIPT FROM COMMAND LINE***

plant_species <- "arabidopsis"
#Printing plant species
print(paste("Species: ",plant_species,sep = ""))
### Create CPM Matrix:


#Check the rownames of the counts matrix to make sure you've got rownames as gene IDs and everything that isn't a row or column name is a number
rownames(dataset) <- dataset[,1]
dataset<-dataset[,-1]

#makes final_data from data in dataset, excluding rows 1-5



#It doesn't seem like the alignment log stats are in the files, this means we probably don't need to run the line below.
#final_data<-dataset[-(1:5),]
final_data<-dataset

if(plant_species=="mouse"){
  rownames(final_data)<- sapply(strsplit(row.names(final_data),split = "\\."), "[",1)
}
if(plant_species=="arabidopsis"){
  rownames(final_data)<- sapply(strsplit(row.names(final_data),split = "\\."), "[",1)
}
count_data <- final_data





# #creating model matrix, MAKE SURE YOU CHOOSE THE RIGHT VALUE TO REPLACE "YOURGROUPCOLUMN" WITH (in both places)!!!
expDesign <- model.matrix(~0 + MetaData$YOURGROUPCOLUMN)

colnames(expDesign)<-sapply(strsplit(colnames(expDesign),split = "YOURGROUPCOLUMN"),'[', 2) 

#Change the "3" here to be a calculated value from the metadata file for the smallest number in any group/condition
y <- DGEList(counts= count_data[rowSums(cpm(as.matrix(count_data)) > 1) >= 3,])

#Counts are already normalized
#y<- calcNormFactors(y, method = "TMM")
# 
#v<-voomWithQualityWeights(y,design = expDesign, method = "genebygene",plot=TRUE)




# 
mds <- plotMDS(y, ndim=3)
# 
mds.out <- as.data.frame(mds$cmdscale.out)
# 
mds.out$Tissue <- strtrim(colnames(count_data),9)
mds.out$Sample <- factor(rownames(lib_sizes), levels = rownames(lib_sizes))
# 
cols <- colorRampPalette(brewer.pal(7, "Paired"))
# 
ggplot(mds.out, aes(x=V1, y=V2, color=Tissue)) +
   geom_point(size=5) +
   ylim(-6,6) +
   xlim(-6,6) +
   ylab(bquote('Leading' ~Log[2]~ 'Fold Change Dim 2')) +
   xlab(bquote('Leading' ~Log[2]~ 'Fold Change Dim 1')) +
   theme_bw() +
   theme(
     text = element_text(size = 16)
   )
# 
 fit<-lmFit(y,expDesign)
# #fit<-lmFit(v,expDesign)

 
# #tweak, create contrast matrix
design.pairs<-function(levels){
  n <- length(levels)
  design <- matrix(0,n,choose(n,2))
  rownames(design) <- levels
  colnames(design) <- 1:choose(n,2)
  k <- 0
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      k <- k+1
      design[i,k] <- 1
      design[j,k] <- -1
      colnames(design)[k] <- paste(levels[i],"_v_",levels[j],sep="")
    }
  }
  design
}

#we can start again from here
levels=colnames(expDesign)
contMatrix<-design.pairs(levels = levels)
DE<-contrasts.fit(fit,contMatrix)
DE<-eBayes(DE)

#This part will need to be done separately for each species


if(plant_species == "arabidopsis"){
suppressMessages(require(biomaRt))
suppressMessages(require(xml2))
dbthale<-useEnsembl(biomart="plants_mart",dataset = "athaliana_eg_gene",host="plants.ensembl.org")
AthalianaAnnotations<-getBM(attributes = c('ensembl_gene_id','description','external_gene_name','entrezgene_id'),mart=dbthale)
AthalianaAnnotations[AthalianaAnnotations$external_gene_name=="","external_gene_name"]<-AthalianaAnnotations[AthalianaAnnotations$external_gene_name=="","ensembl_gene_id"]
Annotations<-AthalianaAnnotations
}

if(plant_species == "mouse"){
suppressMessages(require(biomaRt))
suppressMessages(require(xml2))
dbmus<-useMart("ensembl", dataset="mmusculus_gene_ensembl")
Atributes<-listAttributes(dbmus)
MouseAnnotations<-getBM(attributes=c('ensembl_gene_id','description','mgi_symbol','entrezgene_id'),mart=dbmus)
Annotations<-MouseAnnotations
}

#Add if exists comments before these rm commands to prevent warnings
if(exists("deframe")){
rm(deframe)
}
if(exists("dframe")){
  rm(dframe)
}


print("Let's get it started in ha")
for(contrasts in colnames(contMatrix)){
  print(paste("Beginning DE, GO, and Pathway Analysis for ", contrasts,".",sep = ""))
  assign("temp",topTable(DE,coef =paste(contrasts),adjust.method = "BH",p.value = 0.05,number=length(DE$coefficients)))
  #print(paste(contrasts,"_DE",sep = ""))
  #print(c("The number of significantly DE genes for ",paste(contrasts),sum(topTable(DE,coef =paste(contrasts),adjust.method = "BH",p.value = 1,number=length(DE$coefficients))[5]<0.05)))
  #Add if exists command here to check if Annotations exists, and only run if that is true
  if(exists("Annotations")){
    temp<-merge(Annotations,temp,by.x="ensembl_gene_id",by.y="row.names")
    temp$Row.names <- NULL
    if(length(temp$logFC)>0){
      assign(paste(contrasts,"_DE",sep = ""),temp)
      assign(paste(contrasts,"_GO",sep = ""),GO_Wrapper(temp,out.prefix=paste(contrasts),SpeciesInput=plant_species))
      assign(paste(contrasts,"_KEGG",sep = ""),PathwayD(temp,out.suffix = paste(contrasts),SpeciesInput = plant_species))
      print("Welp, at least we got this far....")
    }else{
      print("No significantly differentially expressed genes for current comparison")
    }
  }
  else{
    print("No annotation file available. Pathway and GO analysis will not be performed.")
    assign(paste(contrasts,"_DE",sep = ""),temp)
  }
  print(paste("Finished DE, GO, and Pathway Analysis for ", contrasts,".",sep = ""))
  write.csv(temp,file=paste(dir,contrasts,".csv",sep = ""))
  if(exists("dframe")){
    print("and again")
    deframe=data.frame(contrasts,sum(topTable(DE,coef =paste(contrasts),adjust.method = "BH",p.value = 1,number=length(DE$coefficients))[5]<0.05),sum(topTable(DE,coef =paste(contrasts),adjust.method = "BH",p.value = 1,number=length(DE$coefficients))[4]<0.05))
    names(deframe)=c("contrast","SigDEGenes_adj","SigDEGenes")
    dframe<-rbind(dframe,deframe)
  }
  if(!exists("dframe")){
    print("First time through")
    dframe=data.frame(contrasts,sum(topTable(DE,coef =paste(contrasts),adjust.method = "BH",p.value = 1,number=length(DE$coefficients))[5]<0.05),sum(topTable(DE,coef =paste(contrasts),adjust.method = "BH",p.value = 1,number=length(DE$coefficients))[4]<0.05))
    names(dframe)<-c("contrast","SigDEGenes_adj","SigDEGenes")
  }
}

save.image(paste(dir,"/","DE_workspace.Rdata",sep = ""))
