options(warn=-1)
##############################  RUVSEQ  Normalization ##############################

##### Installing Package RUVSeq ############
if (!require("devtools")) install.packages("devtools",repos="http://cran.us.r-project.org")
print ("Loading required dependencies for devtools")
library(withr)
library(httr)
library(curl)
library(devtools)
if (!require("RUVSeq")) biocLite("RUVSeq")
library(RUVSeq)
############################################


####### Reading Command Line Arguments ########
args = commandArgs(trailingOnly=TRUE)
path <- args[1]
filenameR <- args[2]
NumberofControl <- as.numeric(args[3])
NumberofTreatment <- as.numeric(args[4])
Counts <- as.numeric(args[5])
NumberofSamples <- as.numeric(args[6])
##############################################


setwd(path)
#file <- filenameR
countData <- as.matrix(read.table(filenameR,sep="\t",header=TRUE,row.names=1,check.names=F))
print ("This is how your input data looks")
head(countData)

print ("Filtering and exploratory data analysis")
## Filtering and exploratory data analysis
##########countData <- as.integer(countData)
filter <- apply(countData, 1, function(x) length(x[x>Counts])>=NumberofSamples)
filtered <- countData[filter,]
genes <- rownames(filtered)
aa <- length(genes)
Controlrep <- rep("Ctl",NumberofControl)
Treatmentrep <- rep("Trt",NumberofTreatment)
x <- as.factor(c(Controlrep,Treatmentrep))
x
set <- newSeqExpressionSet(as.matrix(filtered), phenoData = data.frame(x, row.names=colnames(filtered)))

print ("Performing Upper Quantile Normalization")
##upper quartlie normalization
set <- betweenLaneNormalization(set, which="upper")

print ("Performing DE using EDGER package")
##Differential expression analysis using Empirical RUVg
design <- model.matrix(~x, data=pData(set))
y <- DGEList(counts=counts(set), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, method="deviance", robust=TRUE, subset=NULL)
y <- estimateGLMTagwiseDisp(y)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
print ("Showing top DE genes")
topTags(lrt)
result1 <- topTags(lrt, n=dim(y)[1]+1, adjust.method="BH", sort.by="logFC")
write.table(result1,file = "temporaryfile.txt",sep = "\t",quote=FALSE)
########################### END RUVSEq #############################