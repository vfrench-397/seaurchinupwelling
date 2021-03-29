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
#our length is 30284 - this is the number of genes we have represented

#here we are renaming our data by treatment: N=nonupwelling and U=upwelling conditions
names(countData)=c("NN", "NN", "NN", "UU", "UU", "UU")
#row.names(countData)=sub("", "isogroup", rownames(countData)) this is not necessary for our data 

head(countData)


#will have to change these paths for each one of your computers
setwd("C:/Users/Maddy/Documents/BI586/seaurchinupwelling/outlier")
getwd()
v=setwd("C:/Users/Maddy/Documents/BI586/seaurchinupwelling/outlier")

setwd("/usr4/bi594/vfrench3/assignment2/seaurchinupwelling")
v= setwd("/usr4/bi594/vfrench3/assignment2/seaurchinupwelling/outlier")

# # # #look for outliers
treat=c("NN", "NN", "NN", "UU", "UU", "UU")


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


#####you only ever need to run the above code once. Outliers are decided at the beginning. 

#(can see outlier in index html file in outlier folder)
#We have an outlier in UU3 detected under "Distances between arrays." We will be keeping this outlier in our data because it was only offensive in distance between arrays, but not in boxplots or MA plots.
#Based on the barplot of our total counts data below, we see even distribution between our treatments, which is another reason we are keeping outlier.
#low or high sequencing depth can create outliers 


#will need to change wd for each computer
setwd("C:/Users/Maddy/Documents/BI586/seaurchinupwelling")
setwd("/usr4/bi594/vfrench3/assignment2/seaurchinupwelling")
getwd()

library("DESeq2")
library("ggplot2")


#read in counts 
countData <- as.matrix(read.table("geneCounts_02122019.txt"))
class(countData)
head(countData)
length(countData[,1])
#our length is 30284

names(countData)=c( "NN", "NN", "NN", "UU", "UU", "UU")
#row.names(countData)=sub("", "isogroup", rownames(countData))
head(countData)

totalCounts=colSums(countData)

totalCounts #How many reads associated with an isogroup for each treatment 
#NN1 #NN2 #NN3 #UU1 #UU2 #UU3 
#8555156 #8577700 #8948115 #8570174 #7455015 #6376636
#our raw counts range from 6mil to 8mil
barplot(totalCounts, col=c("slateblue", "slateblue", "slateblue", "royalblue4", "royalblue4", "royalblue4"), ylab="raw counts")
#raw counts generally uniformly distributed for each treatment! Good signal for normalization 


min(totalCounts) #our number is 6376636
max(totalCounts)  # our number is 8948115



treat=c( "NN", "NN", "NN", "UU", "UU", "UU")

#Treatment = upwelling vs. nonupwelling conditions 
g=data.frame(treat)
g
colData<- g

str(colData)
#Creating colData again like above 

#creating DESeq object, design is treatment group
dds<-DESeqDataSetFromMatrix(countData=countData, colData=colData, design= ~ treat) 

#one step DESeq
dds<-DESeq(dds)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing



head(dds) #this is deseq object
res<- results(dds)
res



#Look at dispersions plot
plotDispEsts(dds, main="Dispersion plot")
#this should look like hockey stick, this is visual representation of deseq
#Our data follows the hockey stick shape/dispersion fits well to curve

####################upwelling vs nonupwelling pairwise comparisons
#here we are doing analyses to look at differentially expressed genes among our two treatments; upwelling vs. nonupwelling
colData$UU<-factor(colData$treat, levels=c("UU","NN"))
#order of levels matters: the second term (nonupwelling) is our control; upwelling will be compared to nonupwelling
resUU <- results(dds, contrast=c("treat","UU","NN"))
#how many FDR < 10%?
table(resUU$padj<0.001) 
# 0.1=3039 - not valuable
# 0.05=2289
# 0.01=1344
# 0.001 = 715
summary(resUU)
#doing summary allows us to see that there are genes removed in pairwise analysis due to low counts (23%)
#also see 5.4% of genes were downregulated which means downregulated in UU relative to NN
#7.2% of genes were upregulated in UU relative to NN

nrow(resUU[resUU$padj<0.05 & !is.na(resUU$padj),])   # Num significantly differentially expressed genes excluding the no/low count genes
#2289, this is the same result we get above for 0.05. This is another way to look at differentially expressed genes


plotMA(resUU, main="NN vs UU") #leave at ylim 4,-4

results <- as.data.frame(resUU)
head(results)
#shows us summary of results

nrow(resUU[resUU$padj<0.1 & resUU$log2FoldChange > 0 & !is.na(resUU$padj),])
#this is looking at upregulated bc logfold > 0
nrow(resUU[resUU$padj<0.1 & resUU$log2FoldChange < 0 & !is.na(resUU$padj),])
#UP in UU is 1735
#DOWN in UU is 1304


write.table(resUU, file="UU_DEG.txt", quote=F, sep="\t") #include in paper

#this is results summary from above in .txt format and is saved in the outlier folder 

cd <- read.table("UU_DEG.txt")
head(cd)

##make the GO table for MWU
head(cd)
library(dplyr)
cd
go_input_UU = cd %>%
  tibble::rownames_to_column(var = "iso") %>%
  mutate(mutated_p = -log(pvalue)) %>%
  #making a ranked p-value with directionality here
  #if pos log fold change means upregulated, if neg means down regulated
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(iso, mutated_p_updown)
head(go_input_UU)
colnames(go_input_UU) <- c("gene", "pval")
head(go_input_UU)
write.csv(go_input_UU, file="UU_GO.csv", quote=F, row.names=FALSE)



###############################################################################################
##############################################################################
#--------------get pvals

#here we are binding two columns, the first column is the pval from res comparison of UU to NN and then the p adjusted values
valUU=cbind(resUU$pvalue, resUU$padj)
head(valUU)
#this is pvalues and padjusted values
colnames(valUU)=c("pval.UU", "padj.UU")
length(valUU[,1]) #this is the number of genes we are looking at = 30284 (same as above)
table(complete.cases(valUU))
#False = NAs

######-------------make rlogdata and pvals table

#now doing r log transformation, this is important for making heat maps. it is normalization method and is also how we are going to make PCAs 

rlog=rlogTransformation(dds, blind=TRUE) 
rld=assay(rlog)
head(rld)
#this shows us for each isogroup, the r log normalized values for each of our samples
colnames(rld)=paste(colData$treat)
head(rld)
length(rld[,1]) #length should be still same, indicating retention of all data

#binding rld data and pvalues; ranking pvalues due to significance for heat map
rldpvals=cbind(rld,valUU)
head(rldpvals)
#bound r log normalized data with p-values
dim(rldpvals) #looking at dimensions [1] 30284     8 , there are more columns here because we added columns of pvalues
table(complete.cases(rldpvals))
#FALSE  TRUE 
#11745  18539 , we still have the same number of NAs (false) here

write.csv(rldpvals, "RLDandPVALS.csv", quote=F)


colnames(rld)=paste(colData$treat)
head(rld)

library(RColorBrewer)
# making sample distance heatmap
sampleDists <- as.matrix(dist(t(rld)))
library(gplots)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          margin=c(10, 10), main="Sample Distance Matrix")


#heatmap shows clustering - how similar samples are 
#these are a good way to visualize overall expression between our samples
#i think this heatmap looks good, it is showing the clustering of our treatments, which we expect 


#################################################################################
# VENN Diagram to include both up and down regulated genes in common for PC02
#shows number of differentially expressed genes
library(VennDiagram)

#making series of up reg and down reg genes (p adjusted values of 0.1 fairly standard in literature)
UU_up=row.names(resUU[resUU$padj<0.1 & !is.na(resUU$padj) & resUU$log2FoldChange>0,])
#using res function and asking what are the row names that meet these requirements: 
#has to have padj of 0.1, cant be a NA, and because we want upreg gene in this case we need logfold change to be >0.. these are the 3 requirements to be in the up
length(UU_up) #1735 genes upregulated in UU compared to NN (control)
UU_down=row.names(resUU[resUU$padj<0.1 & !is.na(resUU$padj) & resUU$log2FoldChange<0,])
length(UU_down) #1304


UU=row.names(resUU[resUU$padj<0.1 & !is.na(resUU$padj),])
length(UU) #sanity check, should be the sum of up and down gene numbers from above (1735+1304), which it is


################################
#not sure if we need these parts next because they are combining all upregulated groups in 7.6 and 7.5 treatments and then getting rid of repetitive isogroups 
#since we only have one treatment compared to control (UU) we probably dont need to combine anything????
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
##################################



###do UP, DOWN, ALL
candidates=list("UP"=UU_up, "DOWN"=UU_down) #I am not sure what to compare here, should we compare UU and NN?
#I agree, I will reach out to Sarah on slack tomorrow! I can't see any instance where we would need this? 
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
rldpvals <- read.csv(file="RLDandPVALS.csv", row.names=1)
head(rldpvals)
rld=rldpvals[,1:6] 
#making rld which is just columns 1-6 and will cut off pvalues to just leave us with r log normalized isogroups
head(rld)
#leaving NN.1, NN.2 instead of just NN etc. not sure if this is ok??????? 

sampleDists <- dist(t(rld))
sampleDistMatrix <- as.matrix( sampleDists )
treat=c( "upwelling", "upwelling", "upwelling", "non-upwelling", "non-upwelling", "non-upwelling")
colnames(sampleDistMatrix)=paste(treat)
rownames(sampleDistMatrix)=paste(treat)

#creating a heat map from the sampledistMatrix 
library("pheatmap")
heat.colors = colorRampPalette(rev(c("navy","beige")),bias=0.3)(100)
pheatmap(sampleDistMatrix,color = heat.colors,cex=0.9,border_color=NA,cluster_rows=T,cluster_cols=T)

library(vegan)
library(ggplot2)
library(ggrepel)
library(tidyverse)


#now we are doing principle components analaysis PCA
#PCA is looking at distance between 2 dots


rld_t=t(rld) #transposing data frame 

#Removing columns with zero variance (all log fold change 0 meaning not differentially expressed)
which(apply(rld_t, 2, var)==0)
rld_t <- rld_t[ , which(apply(rld_t, 2, var) != 0)]

pca <- prcomp(rld_t,center = TRUE, scale. = TRUE) #up to 6 principle components because n=6 samples 
head(pca)

#we are interested in PC1 and PC2 bc these are the two principle components that explain the most variance
#defining amount of variance described by PC1 and PC2 

li <- pca$sdev^2 / sum(pca$sdev^2)
pc1v <- round(li[1] * 100, 1)
pc2v <- round(li[2] * 100, 1)
pca_s <- as.data.frame(pca$x)
head(pca_s)
pca_s <- pca_s[,c(1,2)] #pulling out first two PCs
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
#ANOVA analysis testing if distances between samples on PCA are significantly different 
#treatment (at pvalue 0.003) has a specific impact on gene differential expression 
adonis(pca$x ~ treat, data = pca_s, method='eu', na.rm = TRUE)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# treat      2     40140 20069.9  13.048 0.81306  0.003 **
# Residuals  6      9229  1538.1         0.18694          
# Total      8     49369                 1.00000 


###################################heatmaps for genes NS vs FR
rldpvals <- read.csv(file="RLDandPVALS.csv", row.names=1)
head(rldpvals)
rld_site= rldpvals[,1:6]
head(rld_site)
gg=read.table("goAnnot_spu.tab",sep="\t", row.names=1)
head(gg)

nrow(rldpvals[rldpvals$padj.UU<0.01& !is.na(rldpvals$padj.UU),])
#1344; isogroups extremely differently expressed (p-value <0.01) 

#Ranking p-values to find top 100 most differentially expressed genes 
topnum= 100 # number of DEGS in heatmap 
head(rldpvals)
top100=head(rldpvals[order(rldpvals$padj.UU), ],topnum)
head(top100)
length(top100[,1])
summary(top100)
###
library(pheatmap)
head(top100)
p.val=0.1 # FDR cutoff
conds=top100[top100$padj.UU<=p.val & !is.na(top100$padj.UU),] #all top 100 where the p- adjusted value is less than or equal to .1 (significant)
length(conds[,1])

exp=conds[,1:6] #removing associated p-values from log fold change data
means=apply(exp,1,mean) # means of rows
explc=exp-means # subtracting them
head(explc)

ccol=colorRampPalette(rev(c("red","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)

pheatmap(explc,cluster_cols=T,scale="row",color=col0, show_rownames = F)

###################################Heatmap for the genes in common ##We don't need this, genes in common for more than 1 treatment 
rldpvals <- read.csv(file="RLDandPVALS.csv", row.names=1)
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
ann = data.frame(cond = c('UU', 'UU', 'UU', 'NN', 'NN', 'NN'))
rownames(ann) <- names(explc)

# Set colors
Var1        <- c("darkgoldenrod2",  "darkolivegreen3", "black")
names(Var1) <- c("UU", "NN")
anno_colors <- list(cond = Var1)

pheatmap(as.matrix(explc),annotation_col=ann,annotation_colors=anno_colors,cex=.85,color=col0,border_color=NA,clustering_distance_rows="correlation",clustering_distance_cols="correlation", show_rownames=T)
