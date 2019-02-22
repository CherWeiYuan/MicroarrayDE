#GSE studied here: GSE45121

#####Installation of packages#####
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("GEOquery", version = "3.8")
BiocManager::install("affy")
BiocManager::install("limma")
BiocManager::install("mouse4302.db", version = "3.8")

#####Loading packages#####
library(GEOquery)
library(affy)
library(hgu133plus2.db)
library(org.Hs.eg.db)
library(mouse4302.db)


#####Importing data from GEO#####
gse <- getGEO('GSE45121', GSEMatrix = FALSE)
class(gse) #confirms if downloaded file has class "Geoquery"
GSE45121cel <- getGEOSuppFiles('GSE45121')
#Shift the downloaded "GSE45121_RAW.tar" into the working directory before the following step:
untar("GSE45121_RAW.tar", exdir="data") #extract gz files to folder named "data"
cels <- list.files("data/", pattern = "[gz]") 
cels
sapply(paste("data", cels, sep="/"), gunzip) #converts the gz files into cel files

#Pseudo-image of the chip to check quality of microarray
image(data[,-1])

#####Attaching phenoData to our assayData#####
gse <- getGEO('GSE45121', GSEMatrix = FALSE) #run if not done yet
gsm <- GSMList(gse)[[1]]
#names(GSMList(gse)) gets the names of all GSM objects in the GSE whereas GSMList(gse)[[1]] gets the first GSM object on the list
names(Meta(gsm))
names(Meta(gse))

names(Meta(gsm))[!(names(Meta(gsm)) %in% names(Meta(gse)))] #tells us what metadata names present in GSM and GSE so we know what phenoData from GSE applies to all our GSMs
Meta(gsm)[!(names(Meta(gsm)) %in% names(Meta(gse)))] #we can see in $characteristics_ch1 that the GSMs contain data for three samples: cell line SUM159, control and breast cancer tissue

for (gsm in GSMList(gse)) { 
  print(Meta(gsm)[['characteristics_ch1']])
} # to check the characteristics of each GSM's metadata

#we append the 'genotype' metadata to our assayData
genotype <- function(gsm) { 
  Meta(gsm)[['characteristics_ch1']][2]
} #creates a function that refers to treatment in the characteristics metadata

sapply(GSMList(gse),genotype) #uses the new function made in the previous step to append genotype metadata in the list of GSMs acquired from the GSE

#Making a dataframe of the GSMs and the treatment type
pd <- data.frame(genotype=as.factor(sapply(GSMList(gse),genotype)))
pd

pd$genotype <- as.factor(pd$genotype)
levels(pd$genotype) <- c("miR203","WT") #simplifies the treatment column 
pd

celfiles <- paste0('data/',rownames(pd),'.CEL') #since the row names of the dataframe (pd) is e.g. GSM1226581, can use it to load the cel. files with the same name 
#we need to change CEL files in the data/ folder to e.g. from GSM1097476_DF-P3_9.25.09.CEL to GSM1097476.CEL for this code to work
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
model <- model.matrix(~0 + eset$genotype)
model
colnames(model) <- levels(eset$genotype)
model
contrasts <- makeContrasts(miR203 - WT, levels=model)
contrasts


fit <- lmFit(eset, model)
fitted.contrast <- contrasts.fit(fit,contrasts)
fitted.ebayes <- eBayes(fitted.contrast)
topTable(fitted.ebayes)



#Annotating genes
ps <- rownames(topTable(fitted.ebayes))
ps
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("mouse4302.db", version = "3.8")
library(mouse4302.db)

#AnnotationDbi interface
columns(mouse4302.db)
keytypes(mouse4302.db)
head(keys(mouse4302.db,keytype="PROBEID"))

AnnotationDbi::select(mouse4302.db,ps,c(
  "SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")

#Retrieving genes all genes that are differentially expressed with a 
#adjusted p-value of less than 0.05, with at fold change of 
#at least 1.8 (log fold change at least one).
ps2 <- topTable(fitted.ebayes,number=Inf,p.value = 0.3,lfc=1.2)
ps2_up <- rownames(ps2[ps2$logFC > 0,])
AnnotationDbi::select(mouse4302.db,ps2_up,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")


## Due to the high adjusted p-value and the insufficient log-fold change,
## we abandoned this dataset for further analysis.