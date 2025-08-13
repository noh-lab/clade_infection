############ ############# clade infection HOST (202506 version) ############ #############

library(BiocManager)
libraries <- c("DESeq2", "ggplot2", "RColorBrewer", "GenomicFeatures", "topGO", "GSEABase", "gridExtra", "Rgraphviz")
lapply(libraries, require, character.only=T) 

#"pheatmap", "Gviz", "gage", "ShortRead", "reshape2", "plyr", "gplots", 
manuscript_theme = ggplot2::theme_bw() + ggplot2::theme(axis.line.x = element_line(color="black", linewidth = .7),axis.line.y = element_line(color="black", size = .7),text=element_text(size=16))
mypalette <- colorRampPalette(c("#ff5858","#53f1f1","#c5fd57")) # http://paletton.com
hist(discoveries, col=mypalette(3))

getwd()



############ ############# DE ANALYSIS by each symbiont ############# #############
directory.d <- "./dicty_5_star_count/"
sampleFiles.d <- grep("count",list.files(directory.d),value=TRUE)
sampleName.d <- sub("(*).count","\\1",sampleFiles.d)
sampleTable.d <- data.frame(sampleName = sampleName.d,
                            fileName = sampleFiles.d,
                            condition = c("control","infected","infected","control","infected","infected","infected","infected","infected","infected","control","control","control","control","infected","infected","infected","infected"),
                            dicty = c(rep("qs18",6), rep("qs864",6), rep("qs9",6)),
                            burk = c("none","b11","b159","none","b70","b859",
                                     "b11","b159","b70","b859","none","none",
                                     "none","none","b11","b159","b70","b859"),
                            batch = c("1966_1","1966_1","1966_1","1966_2","1966_2","1966_2",
                                      "1966_3","1966_3","1966_4","1966_4","1966_3","1966_4",
                                      "1952_7","1952_8","1952_7","1952_7","1952_8","1952_8"),
                            clade = c("none","reduced","nonreduced","none","nonreduced","reduced",
                                      "reduced","nonreduced","nonreduced","reduced","none","none",
                                      "none","none","reduced","nonreduced","nonreduced","reduced")
)
sampleTable.d$condition <- factor(sampleTable.d$condition)
sampleTable.d$burk <- factor(sampleTable.d$burk)
sampleTable.d$burk <- relevel(sampleTable.d$burk, ref="none")
sampleTable.d$dicty <- factor(sampleTable.d$dicty)


# DE host by each symbiont, outlier removed
dds.strain <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable.d[-6,],
                                         directory = directory.d,
                                         design= ~burk+dicty)

# checking samples for proportion of genes with expression
colSums(counts(dds.strain) >= 10) / nrow(counts(dds.strain)) # STAR 48-56 %; GMAP 64-73%

# checking reps for correlated counts before outlier removal
#cor(counts(dds.strain)[,c(2,7,15)])
#cor(counts(dds.strain)[,c(3,8,16)])
#cor(counts(dds.strain)[,c(5,9,17)]) 
#cor(counts(dds.strain)[,c(6,10,18)])

# checking reps for correlated counts after outlier removal
cor(counts(dds.strain)[,c(1,4)])
cor(counts(dds.strain)[,c(10,11)])
cor(counts(dds.strain)[,c(12,13)])

cor(counts(dds.strain)[,c(2,6,14)])
cor(counts(dds.strain)[,c(3,7,15)])
cor(counts(dds.strain)[,c(5,8,16)]) 
cor(counts(dds.strain)[,c(9,17)])

# at least three samples must have expression
keep <- rowSums(counts(dds.strain) >= 10) >= 3
table(keep) # 9885 genes

dds.strain <- dds.strain[keep,]

dds.strain <- DESeq(dds.strain)
resultsNames(dds.strain)


# host responses to each burk strain
res.strain159 <- results(dds.strain, contrast=c("burk","b159","none"), alpha = 0.05)
summary(res.strain159) # 15% up, 11% down
mcols(res.strain159)$description

res.strain70 <- results(dds.strain, contrast=c("burk","b70","none"), alpha = 0.05)
summary(res.strain70) # 16% up, 13% down
mcols(res.strain70)$description

res.strain859 <- results(dds.strain, contrast=c("burk","b859","none"), alpha = 0.05)
summary(res.strain859) # 4.5% up, 2.9% down
mcols(res.strain859)$description

res.strain11 <- results(dds.strain, contrast=c("burk","b11","none"), alpha = 0.05)
summary(res.strain11) # 1.2% up, 1.6% down
mcols(res.strain11)$description


# MA plot difference between most DE and least DE
lfc.strain11 <- lfcShrink(dds.strain, coef="burk_b11_vs_none")
plotMA(lfc.strain11)

lfc.strain70 <- lfcShrink(dds.strain, coef="burk_b70_vs_none")
plotMA(lfc.strain70)


# overall visualization with rlog PCA 
rld.strain <- rlogTransformation(dds.strain, blind=F)
plotPCA(rld.strain, intgroup=c("burk"), ntop=5000) # this plot does not change due to glm design

# replot with ggplot2
pca.strain <- plotPCA(rld.strain, intgroup=c("burk"), ntop=5000, returnData=T)
pca.strain$host <- c(rep("QS18",5), rep("QS864",6), rep("QS9",6))
pca.strain$burk <- factor(pca.strain$burk, levels = c("none", "b11", "b859", "b159", "b70"))

ggplot(pca.strain, aes(x=PC1, y=PC2, shape=host, fill=burk)) + geom_point(size =3) + scale_shape_manual(values=c(21,22,24)) + scale_fill_manual(values=c("black", mypalette(3)[1], mypalette(3)[2], "white", mypalette(3)[3])) + manuscript_theme

# how much overlap in DE genes
sig.d159 <- row.names(subset(res.strain159, padj<=0.05))
sig.d70 <- row.names(subset(res.strain70, padj<=0.05))
table(sig.d159 %in% sig.d70)
table(sig.d70 %in% sig.d159) # 1394 genes in common


sig.d859 <- row.names(subset(res.strain859, padj<=0.05))
sig.d11 <- row.names(subset(res.strain11, padj<=0.05))
table(sig.d859 %in% sig.d11)
table(sig.d11 %in% sig.d859) # 115 genes in common


# write to flat tables
write.table(res.strain159, file="dicty159.DEresults.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(res.strain70, file="dicty70.DEresults.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(res.strain859, file="dicty859.DEresults.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(res.strain11, file="dicty11.DEresults.txt", quote=F, sep="\t", row.names=T, col.names=T)

res.strain159 <- read.table(file="dicty159.DEresults.txt", sep="\t", h=T)
res.strain70 <- read.table(file="dicty70.DEresults.txt", sep="\t", h=T)
res.strain859 <- read.table(file="dicty859.DEresults.txt", sep="\t", h=T)
res.strain11 <- read.table(file="dicty11.DEresults.txt", sep="\t", h=T)


############ ############# DE COMPARISON of clades ############# #############
# e.g. are nonreduced clade sig genes also DE in reduced genome clade, just not significantly so

# combine expression data (DE true/false)
de.strain <- data.frame(matrix(NA, ncol=1, nrow=9885))
row.names(de.strain) <- row.names(res.strain859)
de.strain$sig.d859 <- row.names(de.strain) %in% sig.d859
de.strain$sig.d11 <- row.names(de.strain) %in% sig.d11
de.strain$sig.d159 <- row.names(de.strain) %in% sig.d159
de.strain$sig.d70 <- row.names(de.strain) %in% sig.d70
temp <- de.strain[,2:5]
de.strain <- temp
limma::vennDiagram(de.strain) # 80 genes shared by all

all.cand <- rownames(subset(de.strain, sig.d859==TRUE & sig.d11==TRUE & sig.d159==TRUE & sig.d70==TRUE))
plotCounts(dds.strain, gene=all.cand[16], intgroup="burk")

#o859.cand <- rownames(subset(de.strain, sig.859==TRUE & sig.11==FALSE & sig.159==FALSE & sig.70==FALSE))
#plotCounts(dds.strain, gene=o859.cand[9], intgroup="burk")#, returnData=T) # outlier qs18_b859 was driving significance each time

# record mean expression for nonreduced vs. reduced
de.strain$lfc.1 <- 1/2*(res.strain159$log2FoldChange + res.strain70$log2FoldChange) # additive fold change in nonreduced
de.strain$lfc.2 <- 1/2*(res.strain11$log2FoldChange + res.strain859$log2FoldChange) # additive fold change in reduced

require(stats)


# stratify DE genes by whether significant in at least one nonreduced and/or reduced clade
de.strain$sig1.1 <- ifelse(de.strain$sig.d159=="TRUE" | de.strain$sig.d70=="TRUE", "TRUE", "FALSE") # significant in at least one 
de.strain$sig1.2 <- ifelse(de.strain$sig.d11=="TRUE" | de.strain$sig.d859=="TRUE", "TRUE", "FALSE") # significant in at least one

de.strain$comb <- with(de.strain, ifelse(sig1.1=="TRUE" & sig1.2=="TRUE", "both", ifelse(sig1.1=="TRUE" & sig1.2=="FALSE", "nonreduced only", ifelse(sig1.1=="FALSE" & sig1.2=="TRUE", "reduced only", "neither"))))

table(de.strain$comb)


de.strain$comb <- factor(de.strain$comb, levels=c("neither","reduced only","nonreduced only","both"))
table(de.strain$comb)

#reg.on1 <- lm(lfc.2 ~ lfc.1 + 0, data = subset(de.strain, comb=="nonreduced only"))
#coef(reg.on1)[1]

#reg.on2 <- lm(lfc.2 ~ lfc.1 + 0, data = subset(de.strain, comb=="reduced only"))
#coef(reg.on2)[1]

reg.bo1 <- lm(lfc.2 ~ lfc.1 + 0, data = subset(de.strain, comb=="both"))
coef(reg.bo1)[1]
# coef is < 1, indicating that for the same significantly DE genes, expression is higher in nonreduced than reduced 


c1 <- ggplot(subset(de.strain, comb=="neither"), aes(x=lfc.1, y=lfc.2)) + geom_point() + geom_abline(intercept = 0, slope = 1, lty="dashed") + ggtitle("(c) not DE") + xlim(-8,12) + ylim(-10,8) + xlab("response to nonreduced") + ylab("response to reduced") + theme_bw()

#c2 <- ggplot(subset(de.strain, comb=="reduced only"), aes(x=lfc.1, y=lfc.2)) + geom_point() + geom_abline(intercept = 0, slope = 1, lty="dashed") + geom_abline(intercept = 0, slope = coef(reg.on2)[1], color=mypalette(3)[1]) + ggtitle("(c) reduced only DE") + xlim(-8,12) + ylim(-10,8) + xlab("response to nonreduced-genome symbiont") + ylab("response to reduced-genome symbiont") + theme_bw()

#c3 <- ggplot(subset(de.strain, comb=="nonreduced only"), aes(x=lfc.1, y=lfc.2)) + geom_point() + geom_abline(intercept = 0, slope = 1, lty="dashed") + geom_abline(intercept = 0, slope = coef(reg.on1)[1], color=mypalette(3)[3]) + ggtitle("(b) nonreduced only DE") + xlim(-8,12) + ylim(-10,8) + xlab("response to nonreduced-genome symbiont") + ylab("response to reduced-genome symbiont") + theme_bw()

c4 <- ggplot(subset(de.strain, comb=="both"), aes(x=lfc.1, y=lfc.2)) + geom_point() + geom_abline(intercept = 0, slope = 1, lty="dashed") + geom_abline(intercept = 0, slope = coef(reg.bo1)[1], color=mypalette(3)[2]) + ggtitle("(b) both DE (lenient)") + xlim(-8,12) + ylim(-10,8) + xlab("response to nonreduced") + ylab("response to reduced") + theme_bw()



# strict stratification
de.strain$sig2.1 <- ifelse(de.strain$sig.d159=="TRUE" & de.strain$sig.d70=="TRUE", "TRUE", "FALSE") # significant in at least one 
de.strain$sig2.2 <- ifelse(de.strain$sig.d11=="TRUE" & de.strain$sig.d859=="TRUE", "TRUE", "FALSE") # significant in at least one

de.strain$comb2 <- with(de.strain, ifelse(sig2.1=="TRUE" & sig2.2=="TRUE", "both", ifelse(sig2.1=="TRUE" & sig2.2=="FALSE", "nonreduced only", ifelse(sig2.1=="FALSE" & sig2.2=="TRUE", "reduced only", "neither"))))

table(de.strain$comb2)

de.strain$comb2 <- factor(de.strain$comb2, levels=c("neither","reduced only","nonreduced only","both"))
table(de.strain$comb2)


reg.bo2 <- lm(lfc.2 ~ lfc.1 + 0, data = subset(de.strain, comb2=="both"))
coef(reg.bo2)[1]
# coef is < 1, indicating that for the same significantly DE genes, expression is higher in nonreduced than reduced 


c5 <- ggplot(subset(de.strain, comb2=="both"), aes(x=lfc.1, y=lfc.2)) + geom_point() + geom_abline(intercept = 0, slope = 1, lty="dashed") + geom_abline(intercept = 0, slope = coef(reg.bo2)[1], color=mypalette(3)[2]) + ggtitle("(a) both DE (strict)") + xlim(-8,12) + ylim(-10,8) + xlab("response to nonreduced") + ylab("response to reduced") + theme_bw()

grid.arrange(c5, c4, c1, ncol=2)
# the two clades have fairly distinct gene expression patterns: better way to show than venn diagram


############ ############# GO TERM ANALYSIS ############# #############
# define which genes to consider in these tests
diffGene <- function(allScore) {
  return(allScore < 0.05)
}

upGene <- function(allScore) {
  return(names(geneList) %in% up)
}

downGene <- function(allScore) {
  return(names(geneList) %in% down)
}


# define the tests; in recent topGO citations, classic is not used, either default (mix) or elim
#fisher.stat <- new("classicCount", testStatistic=GOFisherTest, name="Fisher test")
elim.stat <- new("elimCount", testStatistic=GOFisherTest, name="Fisher (elim)")
#mix.stat <- new("weight01Count", testStatistic=GOFisherTest, name="Fisher (mix)")


# read in go mappings
dicty.map <- read.table("geneid2go.dictyBase.20230102.filter.topgo.txt", h=F)
colnames(dicty.map) <- c("gene_id","go_id")

go_map.dicty <- by(dicty.map$go_id, dicty.map$gene_id, function(x) as.character(x))
str(head(go_map.dicty))



####  BA159 ASSOCIATED ####
# define gene lists
up <- gsub("_gene","",row.names(subset(res.strain159, padj<=0.05 & log2FoldChange>0)))
down <- gsub("_gene","",row.names(subset(res.strain159, padj<=0.05 & log2FoldChange<0)))
universe <- gsub("_gene","",row.names(subset(res.strain159, is.na(padj)==FALSE)))

# define base data structure
geneList <- subset(res.strain159, is.na(padj)==FALSE)[,6]
names(geneList) <- universe
str(geneList)


# set up test data structure
BP.data159 <- new("topGOdata", description="GO BP analysis of  All Genes", ontology="BP", allGenes=geneList , geneSel=diffGene, annot = annFUN.gene2GO, gene2GO=go_map.dicty, nodeSize=10)
BP.up159 <- new("topGOdata", description="GO BP analysis of  All Genes", ontology="BP", allGenes=geneList , geneSel=upGene, annot = annFUN.gene2GO, gene2GO=go_map.dicty, nodeSize=10)
BP.down159 <- new("topGOdata", description="GO BP analysis of  All Genes", ontology="BP", allGenes=geneList , geneSel=downGene, annot = annFUN.gene2GO, gene2GO=go_map.dicty, nodeSize=10)


# apply test
all159.resultElim <- getSigGroups(BP.data159, elim.stat)
all159.resultElim # 57 @ 0.05 sig terms, 23 @ 0.01
GenTable(BP.data159, elim=all159.resultElim, topNodes=70) 

up159.resultElim <- getSigGroups(BP.up159, elim.stat)
up159.resultElim # 72 sig terms @ 0.05, 30 @ 0.01
GenTable(BP.up159, elim=up159.resultElim, topNodes=100)

down159.resultElim <- getSigGroups(BP.down159, elim.stat)
down159.resultElim # 48 sig terms @ 0.05, 25 @ 0.01
GenTable(BP.down159, elim=down159.resultElim, topNodes=70)


# create and save result
all159.Dicty <- GenTable(BP.data159, elim=all159.resultElim, orderBy="elim", ranksOf="elim", topNodes=57, numChar=99)
up159.Dicty <- GenTable(BP.up159, elim=up159.resultElim, orderBy="elim", ranksOf="elim", topNodes=72, numChar=99)
down159.Dicty <- GenTable(BP.down159, elim=down159.resultElim, orderBy="elim", ranksOf="elim", topNodes=48, numChar=99)

# save results files for go-figure
write.table(all159.Dicty, file="dicty159.allreg.topGO.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(up159.Dicty, file="dicty159.upreg.topGO.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(down159.Dicty, file="dicty159.downreg.topGO.txt", quote=F, sep="\t", row.names=T, col.names=T)


# to export gene lists, do the following and manually remove "- " that yaml attaches
yaml::write_yaml(up, file = "up159.txt")
yaml::write_yaml(down, file = "down159.txt")
yaml::write_yaml(universe, file = "universe.txt")

# visualize 
old.par <- par(mar = c(0, 0, 0, 0))
par(old.par)
showSigOfNodes(BP.up159, score(up159.resultElim), firstSigNodes = 7, useInfo = 'def', .NO.CHAR=10)



####  BA70 ASSOCIATED ####
# define gene lists
up <- gsub("_gene","",row.names(subset(res.strain70, padj<=0.05 & log2FoldChange>0)))
down <- gsub("_gene","",row.names(subset(res.strain70, padj<=0.05 & log2FoldChange<0)))
universe <- gsub("_gene","",row.names(subset(res.strain70, is.na(padj)==FALSE)))

# define base data structure
geneList <- subset(res.strain70, is.na(padj)==FALSE)[,6]
names(geneList) <- universe
str(geneList)


# set up test data structure
BP.data70 <- new("topGOdata", description="GO BP analysis of  All Genes", ontology="BP", allGenes=geneList , geneSel=diffGene, annot = annFUN.gene2GO, gene2GO=go_map.dicty, nodeSize=10)
BP.up70 <- new("topGOdata", description="GO BP analysis of  All Genes", ontology="BP", allGenes=geneList , geneSel=upGene, annot = annFUN.gene2GO, gene2GO=go_map.dicty, nodeSize=10)
BP.down70 <- new("topGOdata", description="GO BP analysis of  All Genes", ontology="BP", allGenes=geneList , geneSel=downGene, annot = annFUN.gene2GO, gene2GO=go_map.dicty, nodeSize=10)


# apply test
all70.resultElim <- getSigGroups(BP.data70, elim.stat)
all70.resultElim # 96 @ 0.05 sig terms, 38 @ 0.01
GenTable(BP.data70, elim=all70.resultElim, topNodes=200) 

up70.resultElim <- getSigGroups(BP.up70, elim.stat)
up70.resultElim # 124 sig terms @ 0.05, 63 @ 0.01
GenTable(BP.up70, elim=up70.resultElim, topNodes=200)

down70.resultElim <- getSigGroups(BP.down70, elim.stat)
down70.resultElim # 65 sig terms @ 0.05, 13 @ 0.01
GenTable(BP.down70, elim=down70.resultElim, topNodes=100)


# create and save result
all70.Dicty <- GenTable(BP.data70, elim=all70.resultElim, orderBy="elim", ranksOf="elim", topNodes=96, numChar=99)
up70.Dicty <- GenTable(BP.up70, elim=up70.resultElim, orderBy="elim", ranksOf="elim", topNodes=124, numChar=99)
down70.Dicty <- GenTable(BP.down70, elim=down70.resultElim, orderBy="elim", ranksOf="elim", topNodes=65, numChar=99)

# save results files for go-figure
write.table(all70.Dicty, file="dicty70.allreg.topGO.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(up70.Dicty, file="dicty70.upreg.topGO.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(down70.Dicty, file="dicty70.downreg.topGO.txt", quote=F, sep="\t", row.names=T, col.names=T)


####  BB859 ASSOCIATED ####
# define gene lists
up <- gsub("_gene","",row.names(subset(res.strain859, padj<=0.05 & log2FoldChange>0)))
down <- gsub("_gene","",row.names(subset(res.strain859, padj<=0.05 & log2FoldChange<0)))
universe <- gsub("_gene","",row.names(subset(res.strain859, is.na(padj)==FALSE)))

# define base data structure
geneList <- subset(res.strain859, is.na(padj)==FALSE)[,6]
names(geneList) <- universe
str(geneList)


# set up test data structure
BP.data859 <- new("topGOdata", description="GO BP analysis of  All Genes", ontology="BP", allGenes=geneList , geneSel=diffGene, annot = annFUN.gene2GO, gene2GO=go_map.dicty, nodeSize=10)
BP.up859 <- new("topGOdata", description="GO BP analysis of  All Genes", ontology="BP", allGenes=geneList , geneSel=upGene, annot = annFUN.gene2GO, gene2GO=go_map.dicty, nodeSize=10)
BP.down859 <- new("topGOdata", description="GO BP analysis of  All Genes", ontology="BP", allGenes=geneList , geneSel=downGene, annot = annFUN.gene2GO, gene2GO=go_map.dicty, nodeSize=10)


# apply test
all859.resultElim <- getSigGroups(BP.data859, elim.stat)
all859.resultElim # 77 @ 0.05 sig terms, 19 @ 0.01
GenTable(BP.data859, elim=all859.resultElim, topNodes=100) 

up859.resultElim <- getSigGroups(BP.up859, elim.stat)
up859.resultElim # 57 sig terms @ 0.05, 18 @ 0.01
GenTable(BP.up859, elim=up859.resultElim, topNodes=100)

down859.resultElim <- getSigGroups(BP.down859, elim.stat)
down859.resultElim # 23 sig terms @ 0.05, 12 @ 0.01
GenTable(BP.down859, elim=down859.resultElim, topNodes=100)


# create and save result
all859.Dicty <- GenTable(BP.data859, elim=all859.resultElim, orderBy="elim", ranksOf="elim", topNodes=77, numChar=99)
up859.Dicty <- GenTable(BP.up859, elim=up859.resultElim, orderBy="elim", ranksOf="elim", topNodes=57, numChar=99)
down859.Dicty <- GenTable(BP.down859, elim=down859.resultElim, orderBy="elim", ranksOf="elim", topNodes=23, numChar=99)

# save results files for go-figure
write.table(all859.Dicty, file="dicty859.allreg.topGO.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(up859.Dicty, file="dicty859.upreg.topGO.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(down859.Dicty, file="dicty859.downreg.topGO.txt", quote=F, sep="\t", row.names=T, col.names=T)


####  BH11 ASSOCIATED ####
# define gene lists
up <- gsub("_gene","",row.names(subset(res.strain11, padj<=0.05 & log2FoldChange>0)))
down <- gsub("_gene","",row.names(subset(res.strain11, padj<=0.05 & log2FoldChange<0)))
universe <- gsub("_gene","",row.names(subset(res.strain11, is.na(padj)==FALSE)))

# define base data structure
geneList <- subset(res.strain11, is.na(padj)==FALSE)[,6]
names(geneList) <- universe
str(geneList)


# set up test data structure
BP.data11 <- new("topGOdata", description="GO BP analysis of  All Genes", ontology="BP", allGenes=geneList , geneSel=diffGene, annot = annFUN.gene2GO, gene2GO=go_map.dicty, nodeSize=10)
BP.up11 <- new("topGOdata", description="GO BP analysis of  All Genes", ontology="BP", allGenes=geneList , geneSel=upGene, annot = annFUN.gene2GO, gene2GO=go_map.dicty, nodeSize=10)
BP.down11 <- new("topGOdata", description="GO BP analysis of  All Genes", ontology="BP", allGenes=geneList , geneSel=downGene, annot = annFUN.gene2GO, gene2GO=go_map.dicty, nodeSize=10)


# apply test
all11.resultElim <- getSigGroups(BP.data11, elim.stat)
all11.resultElim # 47 @ 0.05 sig terms, 11 @ 0.01
GenTable(BP.data11, elim=all11.resultElim, topNodes=70) 

up11.resultElim <- getSigGroups(BP.up11, elim.stat)
up11.resultElim # 17 sig terms @ 0.05, 6 @ 0.01
GenTable(BP.up11, elim=up11.resultElim, topNodes=100)

down11.resultElim <- getSigGroups(BP.down11, elim.stat)
down11.resultElim # 45 sig terms @ 0.05, 11 @ 0.01
GenTable(BP.down11, elim=down11.resultElim, topNodes=70)


# create and save result
all11.Dicty <- GenTable(BP.data11, elim=all11.resultElim, orderBy="elim", ranksOf="elim", topNodes=47, numChar=99)
up11.Dicty <- GenTable(BP.up11, elim=up11.resultElim, orderBy="elim", ranksOf="elim", topNodes=17, numChar=99)
down11.Dicty <- GenTable(BP.down11, elim=down11.resultElim, orderBy="elim", ranksOf="elim", topNodes=45, numChar=99)

# save results files for go-figure
write.table(all11.Dicty, file="dicty11.allreg.topGO.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(up11.Dicty, file="dicty11.upreg.topGO.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(down11.Dicty, file="dicty11.downreg.topGO.txt", quote=F, sep="\t", row.names=T, col.names=T)

rm(all11.Dicty, up11.Dicty, down11.Dicty, all159.Dicty, up159.Dicty, down159.Dicty, all70.Dicty, up70.Dicty, down70.Dicty, all859.Dicty, up859.Dicty, down859.Dicty)


## can run go-figure on desktop
#./gofigure -i dicty159.allreg.topGO.txt -o dicty159.allreg.topGO -j topgo -n bpo -m 50
#./gofigure -i dicty159.upreg.topGO.txt -o dicty159.upreg.topGO -j topgo -n bpo -m 50
#./gofigure -i dicty159.downreg.topGO.txt -o dicty159.downreg.topGO -j topgo -n bpo -m 50

#./gofigure -i dicty70.allreg.topGO.txt -o dicty70.allreg.topGO -j topgo -n bpo -m 50
#./gofigure -i dicty70.upreg.topGO.txt -o dicty70.upreg.topGO -j topgo -n bpo -m 50
#./gofigure -i dicty70.downreg.topGO.txt -o dicty70.downreg.topGO -j topgo -n bpo -m 50

#./gofigure -i dicty859.allreg.topGO.txt -o dicty859.allreg.topGO -j topgo -n bpo -m 50
#./gofigure -i dicty859.upreg.topGO.txt -o dicty859.upreg.topGO -j topgo -n bpo -m 50
#./gofigure -i dicty859.downreg.topGO.txt -o dicty859.downreg.topGO -j topgo -n bpo -m 50

#./gofigure -i dicty11.allreg.topGO.txt -o dicty11.allreg.topGO -j topgo -n bpo -m 50
#./gofigure -i dicty11.upreg.topGO.txt -o dicty11.upreg.topGO -j topgo -n bpo -m 50
#./gofigure -i dicty11.downreg.topGO.txt -o dicty11.downreg.topGO -j topgo -n bpo -m 50


#### COMPARISON ACROSS STRAINS ####
dicty.b159u <- read.table("dicty159.upreg.topGO.txt", h=T, sep="\t")
dicty.b70u <- read.table("dicty70.upreg.topGO.txt", h=T, sep="\t")
dicty.b859u <- read.table("dicty859.upreg.topGO.txt", h=T, sep="\t")
dicty.b11u <- read.table("dicty11.upreg.topGO.txt", h=T, sep="\t")

dicty.b159d <- read.table("dicty159.downreg.topGO.txt", h=T, sep="\t")
dicty.b70d <- read.table("dicty70.downreg.topGO.txt", h=T, sep="\t")
dicty.b859d <- read.table("dicty859.downreg.topGO.txt", h=T, sep="\t")
dicty.b11d <- read.table("dicty11.downreg.topGO.txt", h=T, sep="\t")


# combine GO data (enrichment true/false)
go.strain <- data.frame(matrix(NA, ncol=1, nrow=316))
row.names(go.strain) <- unique(sort(c(dicty.b159u$GO.ID, dicty.b70u$GO.ID, dicty.b859u$GO.ID, dicty.b11u$GO.ID, dicty.b159d$GO.ID, dicty.b70d$GO.ID, dicty.b859d$GO.ID, dicty.b11d$GO.ID)))
go.strain$up.d159 <- row.names(go.strain) %in% dicty.b159u$GO.ID
go.strain$up.d70 <- row.names(go.strain) %in% dicty.b70u$GO.ID
go.strain$up.d859 <- row.names(go.strain) %in% dicty.b859u$GO.ID
go.strain$up.d11 <- row.names(go.strain) %in% dicty.b11u$GO.ID
go.strain$down.d159 <- row.names(go.strain) %in% dicty.b159d$GO.ID
go.strain$down.d70 <- row.names(go.strain) %in% dicty.b70d$GO.ID
go.strain$down.d859 <- row.names(go.strain) %in% dicty.b859d$GO.ID
go.strain$down.d11 <- row.names(go.strain) %in% dicty.b11d$GO.ID
temp <- go.strain[,2:9]
go.strain <- temp

limma::vennDiagram(go.strain[,1:4]) 
limma::vennDiagram(go.strain[,5:8]) 

# 3 GO terms shared by all up, # 3 GO terms shared by all down
write.table(subset(dicty.b159u, GO.ID %in% row.names(subset(go.strain, up.d159+up.d70+up.d859+up.d11 == 4)))[,1:2], file="dictyGO.all.up.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(subset(dicty.b159d, GO.ID %in% row.names(subset(go.strain, down.d159+down.d70+down.d859+down.d11 == 4)))[,1:2], file="dictyGO.all.down.txt", quote=F, sep="\t", row.names=T, col.names=T)

# not hayl shared GO terms
write.table(subset(dicty.b159u, GO.ID %in% row.names(subset(go.strain, up.d159==1 & up.d70==1 & up.d859==1 & up.d11==0)))[,1:2], file="dictyGO.nothayl.up.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(subset(dicty.b159d, GO.ID %in% row.names(subset(go.strain, down.d159==1 & down.d70==1 & down.d859==1 & down.d11==0)))[,1:2], file="dictyGO.nothayl.down.txt", quote=F, sep="\t", row.names=T, col.names=T)

# nonreduced only shared GO terms
write.table(subset(dicty.b159u, GO.ID %in% row.names(subset(go.strain, up.d159==1 & up.d70==1 & up.d859==0 & up.d11==0)))[,1:2], file="dictyGO.nonreduced.up.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(subset(dicty.b159d, GO.ID %in% row.names(subset(go.strain, down.d159==1 & down.d70==1 & down.d859==0 & down.d11==0)))[,1:2], file="dictyGO.nonreduced.down.txt", quote=F, sep="\t", row.names=T, col.names=T)


## get gene IDs for specific GO terms
# create a map for a topGO dataset, should work for all
ann.genes <- genesInTerm(BP.up159)
str(ann.genes)

# select GO terms of interest from the results, but these are not necessarily DE 
go.all1 <- unlist(genesInTerm(BP.up159, whichGO="GO:0019835"), use.names = F) # all up
go.all2 <- unlist(genesInTerm(BP.up159, whichGO="GO:0042742"), use.names = F) # all up
go.all3 <- unlist(genesInTerm(BP.up159, whichGO="GO:0043652"), use.names = F) # all down
go.ha1 <- unlist(genesInTerm(BP.up159, whichGO="GO:0090114"), use.names = F) # not hayl up
go.ha2 <- unlist(genesInTerm(BP.up159, whichGO="GO:0170035"), use.names = F) # not hayl up
go.ha3 <- unlist(genesInTerm(BP.up159, whichGO="GO:1902600"), use.names = F) # not hayl up
go.no1 <- unlist(genesInTerm(BP.up159, whichGO="GO:0010952"), use.names = F) # nonreduced up
go.no2 <- unlist(genesInTerm(BP.up159, whichGO="GO:0015986"), use.names = F) # nonreduced up
go.no3 <- unlist(genesInTerm(BP.up159, whichGO="GO:0042026"), use.names = F) # nonreduced up
go.no4 <- unlist(genesInTerm(BP.up159, whichGO="GO:2001235"), use.names = F) # nonreduced up

sort(c(rownames(subset(res.strain159[go.all1,], padj<0.05)), rownames(subset(res.strain70[go.all,], padj<0.05)), rownames(subset(res.strain859[go.all,], padj<0.05)), rownames(subset(res.strain11[go.all,], padj<0.05))))
sort(c(rownames(subset(res.strain159[go.all2,], padj<0.05)), rownames(subset(res.strain70[go.all,], padj<0.05)), rownames(subset(res.strain859[go.all,], padj<0.05)), rownames(subset(res.strain11[go.all,], padj<0.05))))
sort(c(rownames(subset(res.strain159[go.all3,], padj<0.05)), rownames(subset(res.strain70[go.all,], padj<0.05)), rownames(subset(res.strain859[go.all,], padj<0.05)), rownames(subset(res.strain11[go.all,], padj<0.05))))

sort(c(rownames(subset(res.strain159[go.ha1,], padj<0.05)), rownames(subset(res.strain70[go.ha1,], padj<0.05)), rownames(subset(res.strain859[go.ha1,], padj<0.05))))
sort(c(rownames(subset(res.strain159[go.ha2,], padj<0.05)), rownames(subset(res.strain70[go.ha2,], padj<0.05)), rownames(subset(res.strain859[go.ha2,], padj<0.05))))
sort(c(rownames(subset(res.strain159[go.ha3,], padj<0.05)), rownames(subset(res.strain70[go.ha3,], padj<0.05)), rownames(subset(res.strain859[go.ha3,], padj<0.05))))

sort(c(rownames(subset(res.strain159[go.no1,], padj<0.05)), rownames(subset(res.strain70[go.no1,], padj<0.05))))
sort(c(rownames(subset(res.strain159[go.no2,], padj<0.05)), rownames(subset(res.strain70[go.no2,], padj<0.05))))
sort(c(rownames(subset(res.strain159[go.no3,], padj<0.05)), rownames(subset(res.strain70[go.no3,], padj<0.05))))
sort(c(rownames(subset(res.strain159[go.no4,], padj<0.05)), rownames(subset(res.strain70[go.no4,], padj<0.05))))




############ ############# PHAGOCYTOSIS GENES ANALYSIS ############# #############
## manual hypergeometric test, make sure consistent with: https://systems.crump.ucla.edu/hypergeometric/index.php
# order: (gene list in target list, universe genes in target list, universe genes not in target list, gene list size) 
hg.test <- function(a,b,c,d) { min(1-cumsum(dhyper(0:(a-1),b,c,d))) }

# phagocytosis genes
endo <- read.table("phagocytosis_genes.tsv", sep="\t", h=T, as.is = T)
phago <- endo$ddb

# gene lists
sig.d159 <- row.names(subset(res.strain159, padj<=0.05))
sig.d70 <- row.names(subset(res.strain70, padj<=0.05))
sig.d859 <- row.names(subset(res.strain859, padj<=0.05))
sig.d11 <- row.names(subset(res.strain11, padj<=0.05))

# universe (same for all)
universe <- row.names(subset(res.strain159, is.na(padj)==FALSE))


mylist.1 <- sig.d159[sig.d159 %in% phago]
mylist.2 <- sig.d70[sig.d70 %in% phago]
mylist.3 <- sig.d859[sig.d859 %in% phago]
mylist.4 <- sig.d11[sig.d11 %in% phago]
checklist <- universe[universe %in% phago]

dicty.hg.table <- data.frame(Term=c("host_159","host_70", "host_859", "host_11"),
                                ExpCount=c(length(sig.d159)*length(checklist)/length(universe), length(sig.d70)*length(checklist)/length(universe), length(sig.d859)*length(checklist)/length(universe), length(sig.d11)*length(checklist)/length(universe)),
                                ObsCount=c(length(mylist.1), length(mylist.2), length(mylist.3), length(mylist.4)),
                                SizeTerm=c(length(phago), length(phago), length(phago), length(phago)))

dicty.hg.table$raw.p <- c(hg.test(length(mylist.1), length(checklist), length(universe)-length(checklist), length(sig.d159)), 
                             hg.test(length(mylist.2), length(checklist), length(universe)-length(checklist), length(sig.d70)), 
                             hg.test(length(mylist.3), length(checklist), length(universe)-length(checklist), length(sig.d859)),
                          hg.test(length(mylist.4), length(checklist), length(universe)-length(checklist), length(sig.d11)))

dicty.hg.table$adj.p <- p.adjust(dicty.hg.table$raw.p, method="BH", n=length(dicty.hg.table$raw.p))

dicty.hg.table # b11 is not significant while others are


#### COMPARISON ACROSS STRAINS ####
endo$sig.d159 <- endo$ddb %in% sig.d159
endo$sig.d70 <- endo$ddb %in% sig.d70
endo$sig.d859 <- endo$ddb %in% sig.d859
endo$sig.d11 <- endo$ddb %in% sig.d11

limma::vennDiagram(endo[,4:7]) # 3 genes shared by all, 7 more shared by at least 3

endo$count <- endo$sig.d159 + endo$sig.d70 + endo$sig.d859 + endo$sig.d11
endo[order(endo$count, decreasing=T), ]

phago_cand <- subset(endo, count>1) # these are all upreg except for tirA

de.strain$lfc.d159 <- res.strain159$log2FoldChange
de.strain$lfc.d70 <- res.strain70$log2FoldChange
de.strain$lfc.d859 <- res.strain859$log2FoldChange
de.strain$lfc.d11 <- res.strain11$log2FoldChange

phago_lfc <- de.strain[phago_cand[,2],c(13:16)]
phago_res <- cbind(phago_cand, phago_lfc)

# DDB_G0279307 vacC, DDB_G0279191 vacB, alyA DDB_G0275123 highest upreg

phago_res[order(phago_res$count, decreasing=T),]

############ ############# NETWORK ANALYSIS ############# #############
library(igraph)

# exporting expression data from deseq2
dicty_for_igraph <- varianceStabilizingTransformation(dds.strain, blind=F)
head(assay(dicty_for_igraph),1)

# grabbing significantly differentially expressed genes with larger fold changes
xsig.d159 <- row.names(subset(res.strain159, padj<=0.05 & abs(log2FoldChange)>=2))
xsig.d70 <- row.names(subset(res.strain70, padj<=0.05 & abs(log2FoldChange)>=2))
xsig.d859 <- row.names(subset(res.strain859, padj<=0.05 & abs(log2FoldChange)>=2))
xsig.d11 <- row.names(subset(res.strain11, padj<=0.05 & abs(log2FoldChange)>=2))

# subset assay for sig genes and appropriate columns - #1 controls were with b11 and b159, #2 controls with b70 and b859
head(assay(dicty_for_igraph)[xsig.d159,c(1,3,7,10,12,15)],1)
head(assay(dicty_for_igraph)[xsig.d70,c(4,5,8,11,13,16)],1)
head(assay(dicty_for_igraph)[xsig.d859,c(9,11,13,17)],1)
head(assay(dicty_for_igraph)[xsig.d11,c(1,2,6,10,12,14)],1)

d159_igraph <- assay(dicty_for_igraph)[xsig.d159,c(1,3,7,10,12,15)]
d70_igraph <- assay(dicty_for_igraph)[xsig.d70,c(4,5,8,11,13,16)]
d859_igraph <- assay(dicty_for_igraph)[xsig.d859,c(9,11,13,17)]
d11_igraph <- assay(dicty_for_igraph)[xsig.d11,c(1,2,6,10,12,14)]

rm(dicty_for_igraph)


####  BA159 ASSOCIATED ####
g.d159 <- graph_from_adjacency_matrix(as.matrix(as.dist(cor(t(d159_igraph), method="pearson"))), mode="undirected", weighted=TRUE, diag=FALSE)
g.d159 <- simplify(g.d159, remove.multiple=T, remove.loops=T)

# Color negative correlation edges as blue
E(g.d159)[which(E(g.d159)$weight<0)]$color <- "slategrey"

# Color positive correlation edges as red
E(g.d159)[which(E(g.d159)$weight>0)]$color <- "red"

# Convert edge weights to absolute values
E(g.d159)$weight <- abs(E(g.d159)$weight)

# Remove edges below absolute Pearson correlation 0.8
g.d159 <- delete_edges(g.d159, E(g.d159)[which(E(g.d159)$weight<0.8)])

# Remove any vertices remaining that have no edges
g.d159 <- delete_vertices(g.d159, degree(g.d159)==0)

# Assign names to the graph vertices (optional)
V(g.d159)$name <- V(g.d159)$name

# Change color of vertex frames
V(g.d159)$vertex.frame.color <- "white"

# Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst.d159 <- mst(g.d159, algorithm="prim")

# Community detection with Newman-Girvan algorithm
mst.d159.communities <- cluster_edge_betweenness(mst.d159, weights=NULL)
modularity(mst.d159.communities)
mst.d159.clustering <- make_clusters(mst.d159, membership=mst.d159.communities$membership)

#
d159.commSummary <- data.frame(mst.d159.communities$names, mst.d159.communities$membership, mst.d159.communities$modularity)
colnames(d159.commSummary) <- c("Gene", "Community", "Modularity")
d159.commSummary[order(d159.commSummary$Community),]
degree(mst.d159)[degree(mst.d159)>4]

d159.list <- ifelse(degree(mst.d159)>4, V(mst.d159)$name, NA)



####  BA70 ASSOCIATED ####
g.d70 <- graph_from_adjacency_matrix(as.matrix(as.dist(cor(t(d70_igraph), method="pearson"))), mode="undirected", weighted=TRUE, diag=FALSE)
g.d70 <- simplify(g.d70, remove.multiple=T, remove.loops=T)

E(g.d70)[which(E(g.d70)$weight<0)]$color <- "slategrey"
E(g.d70)[which(E(g.d70)$weight>0)]$color <- "red"

E(g.d70)$weight <- abs(E(g.d70)$weight)
g.d70 <- delete_edges(g.d70, E(g.d70)[which(E(g.d70)$weight<0.8)])
g.d70 <- delete_vertices(g.d70, degree(g.d70)==0)

V(g.d70)$name <- V(g.d70)$name
V(g.d70)$vertex.frame.color <- "white"

# Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst.d70 <- mst(g.d70, algorithm="prim")

# clustering
mst.d70.communities <- cluster_edge_betweenness(mst.d70, weights=NULL)
modularity(mst.d70.communities)
mst.d70.clustering <- make_clusters(mst.d70, membership=mst.d70.communities$membership)

#
d70.commSummary <- data.frame(mst.d70.communities$names, mst.d70.communities$membership, mst.d70.communities$modularity)
colnames(d70.commSummary) <- c("Gene", "Community", "Modularity")
d70.commSummary[order(d70.commSummary$Community),]
degree(mst.d70)[degree(mst.d70)>3]

d70.list <- ifelse(degree(mst.d70)>3, V(mst.d70)$name, NA)



####  BB859 ASSOCIATED ####
g.d859 <- graph_from_adjacency_matrix(as.matrix(as.dist(cor(t(d859_igraph), method="pearson"))), mode="undirected", weighted=TRUE, diag=FALSE)
g.d859 <- simplify(g.d859, remove.multiple=T, remove.loops=T)

E(g.d859)[which(E(g.d859)$weight<0)]$color <- "slategrey"
E(g.d859)[which(E(g.d859)$weight>0)]$color <- "red"

E(g.d859)$weight <- abs(E(g.d859)$weight)
g.d859 <- delete_edges(g.d859, E(g.d859)[which(E(g.d859)$weight<0.8)])
g.d859 <- delete_vertices(g.d859, degree(g.d859)==0)

V(g.d859)$name <- V(g.d859)$name
V(g.d859)$vertex.frame.color <- "white"

# Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst.d859 <- mst(g.d859, algorithm="prim")

# clustering
mst.d859.communities <- cluster_edge_betweenness(mst.d859, weights=NULL)
modularity(mst.d859.communities)
mst.d859.clustering <- make_clusters(mst.d859, membership=mst.d859.communities$membership)

#
d859.commSummary <- data.frame(mst.d859.communities$names, mst.d859.communities$membership, mst.d859.communities$modularity)
colnames(d859.commSummary) <- c("Gene", "Community", "Modularity")
d859.commSummary[order(d859.commSummary$Community),]
degree(mst.d859)[degree(mst.d859)>3]

d859.list <- ifelse(degree(mst.d859)>3, V(mst.d859)$name, NA)


####  BH11 ASSOCIATED ####
g.d11 <- graph_from_adjacency_matrix(as.matrix(as.dist(cor(t(d11_igraph), method="pearson"))), mode="undirected", weighted=TRUE, diag=FALSE)
g.d11 <- simplify(g.d11, remove.multiple=T, remove.loops=T)

E(g.d11)[which(E(g.d11)$weight<0)]$color <- "slategrey"
E(g.d11)[which(E(g.d11)$weight>0)]$color <- "red"

E(g.d11)$weight <- abs(E(g.d11)$weight)
g.d11 <- delete_edges(g.d11, E(g.d11)[which(E(g.d11)$weight<0.8)])
g.d11 <- delete_vertices(g.d11, degree(g.d11)==0)

V(g.d11)$name <- V(g.d11)$name
V(g.d11)$vertex.frame.color <- "white"

# Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst.d11 <- mst(g.d11, algorithm="prim")

# clustering
mst.d11.communities <- cluster_edge_betweenness(mst.d11, weights=NULL)
modularity(mst.d11.communities)
mst.d11.clustering <- make_clusters(mst.d11, membership=mst.d11.communities$membership)

#
d11.commSummary <- data.frame(mst.d11.communities$names, mst.d11.communities$membership, mst.d11.communities$modularity)
colnames(d11.commSummary) <- c("Gene", "Community", "Modularity")
d11.commSummary[order(d11.commSummary$Community),]
degree(mst.d11)[degree(mst.d11)>2]

d11.list <- ifelse(degree(mst.d11)>2, V(mst.d11)$name, NA)




## check different layouts - fr seems best
layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1] 
# Remove layouts that do not apply to our graph.
layouts <- layouts[grepl("fr|kk|lgl", layouts)]
# 3 rows, 1 column
par(mfrow=c(3,1), mar=c(1,1,1,1))

for (layout in layouts) {
  print(layout)
  l <- do.call(layout, list(mst.d159)) 
  plot(mst.d159.clustering, mst.d159, layout=l, edge.arrow.mode=0, main=layout, vertex.label=d159.list, vertex.label.cex=0.75) }

for (layout in layouts) {
  print(layout)
  l <- do.call(layout, list(mst.d70)) 
  plot(mst.d70.clustering, mst.d70, layout=l, edge.arrow.mode=0, main=layout, vertex.label=d70.list, vertex.label.cex=0.75) }

for (layout in layouts) {
  print(layout)
  l <- do.call(layout, list(mst.d859)) 
  plot(mst.d859.clustering, mst.d859, layout=l, edge.arrow.mode=0, main=layout, vertex.label=d859.list, vertex.label.cex=0.75) }

for (layout in layouts) {
  print(layout)
  l <- do.call(layout, list(mst.d11)) 
  plot(mst.d11.clustering, mst.d11, layout=l, edge.arrow.mode=0, main=layout, vertex.label=d11.list, vertex.label.cex=0.75) }


#### NETWORK GRAPHS ####
# Plot the tree object using Fruchterman-Reingold force-directed algorithm
plot(mst.d159, layout=layout_with_fr, vertex.label=d159.list,
     vertex.color="SkyBlue2", vertex.size=4, vertex.label.dist=-0.5, 
     vertex.label.color="slategrey", asp=FALSE, vertex.label.cex=0.75,
     edge.arrow.mode=0, main="dicty(b159)")

# internal vertices become default color despite previous coding
graph.d159 <- plot(mst.d159.clustering, mst.d159, layout=layout_with_fr, vertex.label=d159.list,
                   vertex.size=4, vertex.label.dist=-0.5, vertex.label.color="red", 
                   asp=FALSE, vertex.label.cex=0.75, edge.arrow.mode=0, main="dicty(b159)")

# Plot the tree object using Fruchterman-Reingold force-directed algorithm
plot(mst.d70, layout=layout_with_fr, vertex.label=d70.list,
     vertex.color="SkyBlue2", vertex.size=4, vertex.label.dist=-0.5, 
     vertex.label.color="slategrey", asp=FALSE, vertex.label.cex=0.75,
     edge.arrow.mode=0, main="dicty(b159)")

# internal vertices become default color despite previous coding
graph.d70 <- plot(mst.d70.clustering, mst.d70, layout=layout_with_fr, vertex.label=d70.list,
                  vertex.size=4, vertex.label.dist=-0.5, vertex.label.color="red", 
                  asp=FALSE, vertex.label.cex=0.75, edge.arrow.mode=0, main="dicty(b70)")

# Plot the tree object using Fruchterman-Reingold force-directed algorithm
plot(mst.d859, layout=layout_with_fr, vertex.label=d859.list,
     vertex.color="SkyBlue2", vertex.size=4, vertex.label.dist=-0.5, 
     vertex.label.color="slategrey", asp=FALSE, vertex.label.cex=0.75,
     edge.arrow.mode=0, main="dicty(b159)")

# internal vertices become default color despite previous coding
graph.d859 <- plot(mst.d859.clustering, mst.d859, layout=layout_with_fr, vertex.label=d859.list,
                   vertex.label.color="red", vertex.size=4, vertex.label.dist=-0.5, 
                   asp=FALSE, vertex.label.cex=0.75, edge.arrow.mode=0, main="dicty(b859)")

# Plot the tree object using Fruchterman-Reingold force-directed algorithm
plot(mst.d11, layout=layout_with_fr, vertex.label=d11.list,
     vertex.color="SkyBlue2", vertex.size=4, vertex.label.dist=-0.5, 
     vertex.label.color="slategrey", asp=FALSE, vertex.label.cex=0.75,
     edge.arrow.mode=0, main="dicty(b159)")

# internal vertices become default color despite previous coding
graph.d11 <- plot(mst.d11.clustering, mst.d11, layout=layout_with_fr, vertex.label=d11.list,
                  vertex.size=4, vertex.label.dist=-0.5, vertex.label.color="red", 
                  asp=FALSE, vertex.label.cex=0.75, edge.arrow.mode=0, main="dicty(b11)")


#### COMPARE GENES WITH HIGH VERTEX DEGREES ####
# are these descriptive or is there something to evaluate? 
# Hub score (higher score = more influence), only works for directed graphs
# Vertex degree (higher score = more influence)

# dicty
X <- unique(sort(c(names(degree(mst.d159)[degree(mst.d159)>3]),
                   names(degree(mst.d70)[degree(mst.d70)>3]),
                   names(degree(mst.d859)[degree(mst.d859)>3]),
                   names(degree(mst.d11)[degree(mst.d11)>3]))))
vennX <- data.frame(d159=rep(NA,length(X)), d70=NA, d859=NA, d11=NA)
rownames(vennX) <- X
vennX$d159 <- row.names(vennX) %in% names(degree(mst.d159)[degree(mst.d159)>3])
vennX$d70 <- row.names(vennX) %in% names(degree(mst.d70)[degree(mst.d70)>3])
vennX$d859 <- row.names(vennX) %in% names(degree(mst.d859)[degree(mst.d859)>3])
vennX$d11 <- row.names(vennX) %in% names(degree(mst.d11)[degree(mst.d11)>3])
limma::vennDiagram(vennX)
rownames(subset(vennX, d159==TRUE & d70==TRUE))
rownames(subset(vennX, d859==TRUE & d70==TRUE))

# 4 in common between 159 and 70, or 70 and 859 infected dicty; checked on stringDB to see what was interesting
# "DDB_G0269176" - racF1, RhoGTPase involved in regulation of SCAR complex actin polymerization, transiently associated with phagosomes (Rivero et al 1999); SCAR ("DDB_G0285253") is DE only in dicty(b159)
# "DDB_G0291678" pseudogene
# "DDB_G3961538" deleted
# "DDB_G0280389" unknown function

"DDB_G0269176" %in% sig.159
"DDB_G0269176" %in% sig.70
"DDB_G0269176" %in% sig.859
"DDB_G0269176" %in% sig.11

