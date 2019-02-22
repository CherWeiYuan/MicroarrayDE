#####Installation of packages#####
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("GEOquery", version = "3.8")
BiocManager::install("affy")
BiocManager::install("limma")
BiocManager::install("hgu133plus2.db")
BiocManager::install("org.Hs.eg.db")

#####Loading packages#####
library(GEOquery)
library(affy)
library(hgu133plus2.db)
library(org.Hs.eg.db)



#####Importing data from GEO#####
gse <- getGEO('GSE50697', GSEMatrix = FALSE)
class(gse) #confirms if downloaded file has class "Geoquery"
GSE50697cel <- getGEOSuppFiles('GSE50697')
#Shift the downloaded "GSE50697_RAW.tar" into the working directory before the following step:
untar("GSE50697_RAW.tar", exdir="data") #extract gz files to folder named "data"
cels <- list.files("data/", pattern = "[gz]") 
cels
sapply(paste("data", cels, sep="/"), gunzip) #converts the gz files into cel files



#####Attaching phenoData to our assayData#####
gse <- getGEO('GSE50697', GSEMatrix = FALSE) #run if not done yet
gsm <- GSMList(gse)[[1]]
#names(GSMList(gse)) gets the names of all GSM objects in the GSE whereas GSMList(gse)[[1]] gets the first GSM object on the list
names(Meta(gsm))
names(Meta(gse))

names(Meta(gsm))[!(names(Meta(gsm)) %in% names(Meta(gse)))] #tells us what metadata names present in GSM and GSE so we know what phenoData from GSE applies to all our GSMs
Meta(gsm)[!(names(Meta(gsm)) %in% names(Meta(gse)))] #we can see in $characteristics_ch1 that the GSMs contain data for three samples: cell line SUM159, control and breast cancer tissue

for (gsm in GSMList(gse)) { 
  print(Meta(gsm)[['characteristics_ch1']])
} # to check the characteristics of each GSM's metadata

#Since treatment informs us whether a GSM was a control or treated sample, we append the 'treatment' metadata to our assayData
treatment_type <- function(gsm) { 
  Meta(gsm)[['characteristics_ch1']][2]
} #creates a function that refers to treatment in the characteristics metadata

sapply(GSMList(gse),treatment_type) #uses the new function made in the previous step to append treatment metadata in the list of GSMs acquired from the GSE

#Making a dataframe of the GSMs and the treatment type
pd <- data.frame(treatment=as.factor(sapply(GSMList(gse),treatment_type)))
pd

pd$treatment <- as.factor(pd$treatment)
levels(pd$treatment) <- c("control","treatment") #simplifies the treatment column 
pd

celfiles <- paste0('data/',rownames(pd),'.CEL') #since the row names of the dataframe (pd) is e.g. GSM1226581, can use it to load the cel. files with the same name 
#we need to change CEL files in the data/ folder to e.g. from GSM1226581_S1100958.MDA.01_HG-U133_Plus_2_.CEL to GSM1226581.CEL for this code to work
affydata <- read.affybatch(celfiles,phenoData = new("AnnotatedDataFrame",pd))
affydata



#RMA processing of microarray data
plotDensity.AffyBatch(affydata) #plotting unprocessed data
eset <- rma(affydata, background = TRUE, normalize = TRUE)
plotDensity(exprs(eset),xlab='log intensity',main="feature level densities after RMA",lwd=2)
class(affydata)
class(eset) #after RMA, AffyBatch object becomes an ExpressionSet so we can use the following functions to access each data
annotation(eset) #chip information
pData(eset) #pData
exprs(eset) #normalized gene expression
experimentData(eset) #Experimental information
featureData(eset) #feature data

#Linear models to identify differentially expressed genes
library(limma)
model <- model.matrix(~0 + eset$treatment)
model
colnames(model) <- levels(eset$treatment)
model
contrasts <- makeContrasts(treatment - control, levels=model)
contrasts


fit <- lmFit(eset, model)
fitted.contrast <- contrasts.fit(fit,contrasts)
fitted.ebayes <- eBayes(fitted.contrast)
topTable(fitted.ebayes)



#Annotating genes
ps <- rownames(topTable(fitted.ebayes))
ps
library(hgu133plus2.db)

#Method 1: reaching into database
ls('package:hgu133plus2.db') 
unlist(mget(ps,hgu133plus2SYMBOL))

#Method 2: AnnotationDbi interface
columns(hgu133plus2.db)
keytypes(hgu133plus2.db)
head(keys(hgu133plus2.db,keytype="PROBEID"))

AnnotationDbi::select(hgu133plus2.db,ps,c(
  "SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")

#Retrieving genes all genes that are differentially expressed with a 
#adjusted p-value of less than 0.05, with at fold change of 
#at least 1.8 (log fold change at least one).
ps2 <- topTable(fitted.ebayes,number=Inf,p.value = 0.05,lfc=1.8)
ps2_up <- rownames(ps2[ps2$logFC > 0,])
AnnotationDbi::select(hgu133plus2.db,ps2_up,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")



#Downstream analysis of microarray data
#Volcano plot
interesting_genes <- topTable(fitted.ebayes,number=Inf,p.value = 0.05,lfc=1.8)
interesting_genes

interesting_genes_annotated <- AnnotationDbi::select(hgu133plus2.db,rownames(interesting_genes),c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
interesting_genes_annotated

print(interesting_genes_annotated[3], row.names = TRUE) #counts of genes differentially expressed
print(interesting_genes_annotated[3], row.names = FALSE) #copy and paste the results here into David webpage for downstream analysis

volcanoplot(fitted.ebayes, main=sprintf("%d features pass our cutoffs",nrow(interesting_genes)))
points(interesting_genes[['logFC']],-log10(interesting_genes[['P.Value']]),col='red')

#Heatmap
eset_of_interest <- eset[rownames(interesting_genes),]
heatmap(exprs(eset_of_interest))

library(RColorBrewer)
eset_of_interest <- eset[rownames(interesting_genes),]
heatmap(exprs(eset_of_interest),labCol=eset$treatment, col = rev(brewer.pal(10, "RdBu")))
heatmap(exprs(eset_of_interest),labCol=eset$treatment, labRow=NA, col = rev(brewer.pal(10, "RdBu")), distfun = function(x) as.dist(1-cor(t(x))))