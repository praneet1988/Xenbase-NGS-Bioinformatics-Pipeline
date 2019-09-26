###START
###############################Hierarchical Clustering on Columns only##############################
##Color options are YellowGreenOrange , BlueWhiteRed , YellowBlackBlue, GreenWhitePurple
####Taking command line options
args = commandArgs(trailingOnly=TRUE)
path <- args[1]
filenameR <- args[2]
user_color <- args[3]
####supress warnings
options(warn=-1) 
###check and install packages
print ("Checking and installing required packages")
if (!require("devtools")) install.packages("devtools",repos="http://cran.us.r-project.org")
print ("Loading required dependencies for devtools")
library(withr)
library(httr)
library(curl)
library(devtools)
if (!require("lme4")) devtools::install_github("lme4/lme4")
if (!require("d3heatmap")) devtools::install_github("rstudio/d3heatmap")
if (!require("flexmix")) install.packages("flexmix",repos="http://cran.us.r-project.org")
print ("Loading required dependencies for d3heatmap ")
library(lme4)
library(flexmix)
if (!require("gplots")) install.packages("gplots",repos="http://cran.us.r-project.org")
if (!require("R2HTML")) install.packages("R2HTML",repos="http://cran.us.r-project.org")
if (!require("RColorBrewer")) install.packages("RColorBrewer",repos="http://cran.us.r-project.org")
if (!require("htmlwidgets")) install.packages("htmlwidgets",repos="http://cran.us.r-project.org")
#####load Libraries
print ("Loading required libraries")
library(d3heatmap)
library(R2HTML)
library(htmlwidgets)
library(RColorBrewer)
library(gplots)
####setting path
setwd(path)
file <- filenameR
#### reading data
data <- as.matrix(read.table(file,sep="\t",header=TRUE,row.names=1,check.names=F))
print ("This is how input data looks like")
head(data)
####Based on user color preferences generating heatmap
print ("Based on User preference generating interactive and static heatmap")


###check1 user color
if(identical(user_color,"YellowGreenOrange")){
file <- strsplit(file,".txt")
outfile1 <- paste0(file, "_static", "_heatmap", ".png")
hmcols<-colorRampPalette(c("yellow","green","orange"))(256)
plot <- d3heatmap(data, scale = "row",color = scales::col_numeric(hmcols,c(0,1)), Rowv = FALSE, Colv=TRUE, show_grid=FALSE, theme="dark", anim_duration=0, brush_color="blue", na.rm=TRUE, cellnote_scale=TRUE, cexRow = 0.5)
png(outfile1, width = 8, height = 15, units = 'in', res = 600)
heatmap.2(data, scale = "row", col = hmcols, Rowv = FALSE, Colv = TRUE, show_grid=FALSE, na.rm = TRUE, cellnote_scale = TRUE, cexRow = 0.5, cexCol = 0.5,  trace = "none", density.info = c("histogram","density","none"), labRow=FALSE)
outfile <- paste0(file, "_interactive", "_heatmap", ".html")
saveWidget(plot, outfile)
dev.off()}



####check2
if(identical(user_color,"BlueWhiteRed")){
file <- strsplit(file,".txt")
outfile1 <- paste0(file, "_static", "_heatmap", ".png")
hmcols<-colorRampPalette(c("blue","white","red"))(256)
plot <- d3heatmap(data, scale = "row",color = scales::col_numeric(hmcols,c(0,1)), Rowv = FALSE, Colv=TRUE, show_grid=FALSE, theme="dark", anim_duration=0, brush_color="yellow", na.rm=TRUE, cellnote_scale=TRUE, cexRow = 0.5)
png(outfile1, width = 8, height = 15, units = 'in', res = 600)
heatmap.2(data, scale = "row", col = hmcols, Rowv = FALSE, Colv = TRUE, show_grid=FALSE, na.rm = TRUE, cellnote_scale = TRUE, cexRow = 0.5, cexCol = 0.5,  trace = "none", density.info = c("histogram","density","none"), labRow=FALSE)
outfile <- paste0(file, "_interactive", "_heatmap", ".html")
saveWidget(plot, outfile)
dev.off()}



###check3
if(identical(user_color,"YellowBlackBlue")){
file <- strsplit(file,".txt")
outfile1 <- paste0(file, "_static", "_heatmap", ".png")
hmcols<-colorRampPalette(c("yellow","black","blue"))(256)
plot <- d3heatmap(data, scale = "row",color = scales::col_numeric(hmcols,c(0,1)), Rowv = FALSE, Colv=TRUE, show_grid=FALSE, theme="dark", anim_duration=0, brush_color="green", na.rm=TRUE, cellnote_scale=TRUE, cexRow = 0.5)
png(outfile1, width = 8, height = 15, units = 'in', res = 600)
heatmap.2(data, scale = "row", col = hmcols, Rowv = FALSE, Colv = TRUE, show_grid=FALSE, na.rm = TRUE, cellnote_scale = TRUE, cexRow = 0.5, cexCol = 0.5,  trace = "none", density.info = c("histogram","density","none"), labRow=FALSE)
outfile <- paste0(file, "_interactive", "_heatmap", ".html")
saveWidget(plot, outfile)
dev.off()}



####check4
if(identical(user_color,"GreenWhitePurple")){
file <- strsplit(file,".txt")
outfile1 <- paste0(file, "_static", "_heatmap", ".png")
hmcols<-colorRampPalette(c("green","white","purple"))(256)
plot <- d3heatmap(data, scale = "row",color = scales::col_numeric(hmcols,c(0,1)), Rowv = FALSE, Colv=TRUE, show_grid=FALSE, theme="dark", anim_duration=0, brush_color="yellow", na.rm=TRUE, cellnote_scale=TRUE, cexRow = 0.5)
png(outfile1, width = 8, height = 15, units = 'in', res = 600)
heatmap.2(data, scale = "row", col = hmcols, Rowv = FALSE, Colv = TRUE, show_grid=FALSE, na.rm = TRUE, cellnote_scale = TRUE, cexRow = 0.5, cexCol = 0.5,  trace = "none", density.info = c("histogram","density","none"), labRow=FALSE)
outfile <- paste0(file, "_interactive", "_heatmap", ".html")
saveWidget(plot, outfile)
dev.off()}
####END