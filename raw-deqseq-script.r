# source("http://bioconductor.org/biocLite.R")
# biocLite("DESeq2")

#set your working directory- will need to add paths for each of your computers
setwd("C:/Users/Maddy/Documents/BI586/seaurchinupwelling")
getwd()

###conduct array quality metrics to detect and remove outliers
library(DESeq2) 
library(affycoretools)
library(arrayQualityMetrics)
library(genefilter)
library(Biobase)

#read in counts 
countData <- read.table("geneCounts_02122019.txt")
View(countData)
head(countData)
length(countData[,1])
#our length is 30284
names(countData)=c("NN1", "NN2", "NN3", "UU1", "UU2", "UU3")
row.names(countData)=sub("", "isogroup", rownames(countData))
head(countData)


#will have to change these paths for each one of your computers
setwd("C:/Users/Maddy/Documents/BI586/seaurchinupwelling/outlier")
getwd()
v=setwd("C:/Users/Maddy/Documents/BI586/seaurchinupwelling/outlier")

setwd("/usr4/bi594/vfrench3/assignment2/seaurchinupwelling")
v= setwd("/usr4/bi594/vfrench3/assignment2/seaurchinupwelling/outlier")

# # # #look for outliers
treat=c("NN1", "NN2", "NN3", "UU1", "UU2", "UU3")
#Why not NU and UN treatments in count data?

#Creating coldata; data frame associating sample with treatment
g=data.frame(treat)
g
colData=g

#Calling on DESeq2 to create a model; storing results of analysis of differential expression; asking how our design varies by treatment 
dds=DESeqDataSetFromMatrix(countData=countData,
                           colData = g,
                           design = ~treat)

#Normalizing Data 
vsd.ge=assay(vst(dds))
rl=vst(dds)


#Identifying outliers in data
e=ExpressionSet(assay(rl), AnnotatedDataFrame(as.data.frame(colData(rl))))
arrayQualityMetrics(e,outdir=v,intgroup=c("treat"),force=T)

##dev.off() 
#double-click index.html

#####you only ever need to run the above code once. Outliers are decided at the beginning. 
## I like to close R and restart with packages etc
##So, please save your script and restart R
#no outliers detected! 
  #Looks like outlier is UU3 detected under "Distances between arrays"; remove that whole replication? 
#low or high sequencing depth can create outliers 
#remove outliers from counts file and matrix of call data; HOW? 


#will need to change wd for each computer
setwd("C:/Users/Maddy/Documents/BI586/seaurchinupwelling")
getwd()

setwd("/usr4/bi594/vfrench3/assignment2/seaurchinupwelling")

library("DESeq2")
library("ggplot2")

#read in counts 
countData <- read.table("geneCounts_02122019.txt")
head(countData)
length(countData[,1])
#our length is 30284

names(countData)=c( "NN1", "NN2", "NN3", "UU1", "UU2", "UU3")
row.names(countData)=sub("", "isogroup", rownames(countData))
head(countData)

totalCounts=colSums(countData)

totalCounts
#NN1 #NN2 #NN3 #UU1 #UU2 #UU3 
#8555156 #8577700 #8948115 #8570174 #7455015 #6376636
barplot(totalCounts, col=c("coral", "coral", "coral", "red", "red", "red"), ylab="raw counts")
#not sure if we should have a different color for each group?
#I think a diff color for each treatment not necessarily each replicate 
#raw counts generally uniformly distributed! Good signal for normalization 



min(totalCounts) #our number is 6376636
max(totalCounts)  # our number is 8948115


treat=c( "NN1", "NN2", "NN3", "UU1", "UU2", "UU3")
#Treatment = upwelling vs. nonupwelling conditions 
g=data.frame(treat)
g
colData<- g

#creating DESeq object
dds<-DESeqDataSetFromMatrix(countData=countData, colData=colData, design= ~ treat) #can only test for the main effects of site, pco2, temp

#one step DESeq
dds<-DESeq(dds)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing

#Getting Error; think it has to do with design formula 

#Error in checkForExperimentalReplicates(object, modelMatrix) : 
  
  #The design matrix has the same number of samples and coefficients to fit,
#so estimation of dispersion is not possible. Treating samples
#as replicates was deprecated in v1.20 and no longer supported since v1.22.

head(dds)
res<- results(dds)

#***********this is where I stopped and everything after this still has old data from the deseq lab*************************


#Look at dispersions plot
plotDispEsts(dds, main="Dispersion plot Snails")

####################pH8 vs pH7.6 pairwise comparisons
colData$pH76<-factor(colData$treat, levels=c("pH7.6","pH8"))
respH76 <- results(dds, contrast=c("treat","pH7.6","pH8"))
#how many FDR < 10%?
table(respH76$padj<0.01) 
# 0.1=486
# 0.05=381
# 0.01=242
summary(respH76)

nrow(respH76[respH76$padj<0.05 & !is.na(respH76$padj),])   # Num significantly differentially expressed genes excluding the no/low count genes   #228

#dev.off()

plotMA(respH76, main="pH8 vs pH7.6")
plotMA(respH76, main="pH8 vs pH7.6", ylim=c(-2,2))

results <- as.data.frame(respH76)
head(results)

nrow(respH76[respH76$padj<0.1 & respH76$log2FoldChange > 0 & !is.na(respH76$padj),])
nrow(respH76[respH76$padj<0.1 & respH76$log2FoldChange < 0 & !is.na(respH76$padj),])
#UP in 7.6 66
#DOWN in 7.6 420


write.table(respH76, file="7.6_2016.txt", quote=F, sep="\t") #include in paper

cd <- read.table("7.6_2016.txt")
head(cd)

##make the GO table for MWU
head(cd)

library(dplyr)
cd
go_input_7.6 = cd %>%
  tibble::rownames_to_column(var = "iso") %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(iso, mutated_p_updown)
head(go_input_7.6)
colnames(go_input_7.6) <- c("gene", "pval")
head(go_input_7.6)
write.csv(go_input_7.6, file="7.6_GO.csv", quote=F, row.names=FALSE)

#########################################################################################################
###Pco2 pH8 vs pH7.5
summary(res)
resph75 <- results(dds, contrast=c("treat","pH7.5", "pH8"))
table(resph75$padj<0.05)
# 0.1=118
# 0.05=72
# 0.01=43
summary(resph75)

plotMA(resph75, main="pH8 vs pH7.5")
plotMA(resph75, main="pH8 vs pH7.5", ylim=c(-2,2))

results <- as.data.frame(resph75)
head(results)

nrow(resph75[resph75$padj<0.1 & resph75$log2FoldChange > 0 & !is.na(resph75$padj),])
nrow(resph75[resph75$padj<0.1 & resph75$log2FoldChange < 0 & !is.na(resph75$padj),])
#UP in 7.5 38
#DOWN in 7.5 80

write.table(resph75, file="7.5_2016.txt", quote=F, sep="\t")

cd2 <- read.table("7.5_2016.txt")
head(cd2)

##make the GO table for MWU for 7.5
head(cd2)

library(dplyr)
go_input_7.5 = cd2 %>%
  tibble::rownames_to_column(var = "iso") %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(iso, mutated_p_updown)
head(go_input_7.5)
colnames(go_input_7.5) <- c("gene", "pval")
head(go_input_7.5)
write.csv(go_input_7.5, file="7.5_GO.csv", quote=F, row.names=FALSE)

write.table(resph75, file="7.5_2016.txt", quote=F, sep="\t")
###############################################################################################
##############################################################################
#--------------get pvals
val76=cbind(respH76$pvalue, respH76$padj)
head(val76)
colnames(val76)=c("pval.76", "padj.76")
length(val76[,1])
table(complete.cases(val76))

val75=cbind(resph75$pvalue, resph75$padj)
head(val75)
colnames(val75)=c("pval.75", "padj.75")
length(val75[,1])
table(complete.cases(val75))

######-------------make rlogdata and pvals table
rlog=rlogTransformation(dds, blind=TRUE) 
rld=assay(rlog)
head(rld)
colnames(rld)=paste(colData$treat)
head(rld)
length(rld[,1])

rldpvals=cbind(rld,val75, val76)
head(rldpvals)
dim(rldpvals)
# [1] 19717    13
table(complete.cases(rldpvals))
#FALSE  TRUE 
#17202  2515 

write.csv(rldpvals, "Crep2016_RLDandPVALS.csv", quote=F)

colnames(rld)=paste(colData$treat)
head(rld)

library(RColorBrewer)
# Sample distance heatmap
sampleDists <- as.matrix(dist(t(rld)))
library(gplots)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          margin=c(10, 10), main="Sample Distance Matrix")



#################################################################################
# VENN Diagram to include both up and down regulated genes in common for PC02
library(VennDiagram)

p76_up=row.names(respH76[respH76$padj<0.1 & !is.na(respH76$padj) & respH76$log2FoldChange>0,])
length(p76_up) #66
p76_down=row.names(respH76[respH76$padj<0.1 & !is.na(respH76$padj) & respH76$log2FoldChange<0,])
length(p76_down) #420
p75_up=row.names(resph75[resph75$padj<0.1 & !is.na(resph75$padj) & resph75$log2FoldChange>0,])
length(p75_up) #38
p75_down=row.names(resph75[resph75$padj<0.1 & !is.na(resph75$padj) & resph75$log2FoldChange<0,])
length(p75_down) #80

p76=row.names(respH76[respH76$padj<0.1 & !is.na(respH76$padj),])
p75=row.names(resph75[resph75$padj<0.1 & !is.na(resph75$padj),])

#UP
pdegs05_up=union(p76_up,p75_up)
length(pdegs05_up)
#93

#DOWN
pdegs05_down=union(p76_down,p75_down)
length(pdegs05_down)
#432

#ALL
pdegs05=union(p76,p75)
length(pdegs05)
#524

###do UP, DOWN, ALL
candidates=list("7.6"=p76, "7.5"=p75)
quartz()
prettyvenn=venn.diagram(
  x = candidates,
  filename=NULL,
  col = "transparent",
  fill = c("coral2", "forestgreen"),
  alpha = 0.5,
  # label.col = c("darkred", "white", "darkgreen", "white", "white", "white", "blue4"),
  cex = 2.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col = c("darkred", "darkgreen"),
  cat.cex = 2.5,
  cat.fontfamily = "sans",
  cat.dist = c(0.08, 0.08),
  cat.pos = 1
);
grid.draw(prettyvenn)

###########################heat map of sample distances for pco2
rldpvals <- read.csv(file="Crep2016_RLDandPVALS.csv", row.names=1)
head(rldpvals)
rld=rldpvals[,1:9]
head(rld)

sampleDists <- dist(t(rld))
sampleDistMatrix <- as.matrix( sampleDists )
treat=c( "pH7.5", "pH7.5", "pH7.5", "pH7.6", "pH7.6", "pH7.6", "pH8", "pH8",  "pH8")
colnames(sampleDistMatrix)=paste(treat)
rownames(sampleDistMatrix)=paste(treat)

library("pheatmap")
heat.colors = colorRampPalette(rev(c("blue","yellow","red")),bias=0.3)(100)
pheatmap(sampleDistMatrix,color = heat.colors,cex=0.9,border_color=NA,cluster_rows=T,cluster_cols=T)

library(vegan)
library(ggplot2)
library(ggrepel)
library(tidyverse)

rld_t=t(rld)
pca <- prcomp(rld_t,center = TRUE, scale. = TRUE)
head(pca)
li <- pca$sdev^2 / sum(pca$sdev^2)
pc1v <- round(li[1] * 100, 1)
pc2v <- round(li[2] * 100, 1)
pca_s <- as.data.frame(pca$x)
head(pca_s)
pca_s <- pca_s[,c(1,2)]
pca_s$Samples = row.names(pca_s)
pca_s$treat=colData$treat
head(pca_s)

cbPalette <- c("darkgoldenrod2",  "darkolivegreen3", "dodgerblue3")
ggplot(pca_s, aes(PC1, PC2, color = treat, pch = treat)) +
  geom_point(size=3) +
  #  geom_text_repel(aes(label=Samples)) +
  scale_colour_manual(values=cbPalette)+
  theme_bw() +
  # geom_density2d(alpha=.5)+
  geom_polygon(alpha=.2)+
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance")) 
head(pca)
adonis(pca$x ~ treat, data = pca_s, method='eu', na.rm = TRUE)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# treat      2     40140 20069.9  13.048 0.81306  0.003 **
# Residuals  6      9229  1538.1         0.18694          
# Total      8     49369                 1.00000 


###################################heatmaps for genes NS vs FR
rldpvals <- read.csv(file="Crep2016_RLDandPVALS.csv", row.names=1)
head(rldpvals)
rld_site= rldpvals[,1:9]
head(rld_site)
gg=read.table("Crep454_iso2gene.tab",sep="\t", row.names=1)
head(gg)

nrow(rldpvals[rldpvals$padj.76<0.01& !is.na(rldpvals$padj.76),])
#242

topnum= 100 # number of DEGS
head(rldpvals)
top100=head(rldpvals[order(rldpvals$padj.76), ],topnum)
head(top100)
length(top100[,1])
summary(top100)
###
library(pheatmap)
head(top100)
p.val=0.1 # FDR cutoff
conds=top100[top100$padj.76<=p.val & !is.na(top100$padj.76),]
length(conds[,1])

exp=conds[,1:9] # change numbers to be your vsd data columns
means=apply(exp,1,mean) # means of rows
explc=exp-means # subtracting them
head(explc)

ccol=colorRampPalette(rev(c("red","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)

pheatmap(explc,cluster_cols=T,scale="row",color=col0, show_rownames = F)

###################################Heatmap for the genes in common
rldpvals <- read.csv(file="Crep2016_RLDandPVALS.csv", row.names=1)
head(rldpvals)
p.val=0.1 # FDR cutoff
conds=rldpvals[rldpvals$padj.76<=p.val & !is.na(rldpvals$padj.76) & rldpvals$padj.75<=p.val & !is.na(rldpvals$padj.75),]
rld_data= conds[,c(1:9)]
head(rld_data)
nrow(rld_data)
gg=read.table("Crep454_iso2gene.tab",sep="\t", row.names=1)
library(pheatmap)
means=apply(rld_data,1,mean) # means of rows
explc=rld_data-means # subtracting them

ccol=colorRampPalette(rev(c("red","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)

pheatmap(explc,cluster_cols=T,scale="row",color=col0, show_rownames = F)

# Make annotation table for pheatmap
ann = data.frame(cond = c('7.5', '7.5', '7.5', '7.6', '7.6', '7.6', '8', '8', '8'))
rownames(ann) <- names(explc)

# Set colors
Var1        <- c("darkgoldenrod2",  "darkolivegreen3", "dodgerblue3")
names(Var1) <- c("7.5", "7.6", "8")
anno_colors <- list(cond = Var1)

pheatmap(as.matrix(explc),annotation_col=ann,annotation_colors=anno_colors,cex=1.2,color=col0,border_color=NA,clustering_distance_rows="correlation",clustering_distance_cols="correlation", show_rownames=T)