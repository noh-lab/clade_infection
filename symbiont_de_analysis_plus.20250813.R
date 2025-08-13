############ ############# clade infection SYMBIONT (202506 version) ############ #############

library(BiocManager)
libraries <- c("DESeq2", "ggplot2", "RColorBrewer", "GenomicFeatures", "topGO", "GSEABase", "gridExtra", "Rgraphviz")
lapply(libraries, require, character.only=T) 

#"pheatmap", "Gviz", "gage", "ShortRead", "reshape2", "plyr", "gplots", 
manuscript_theme = ggplot2::theme_bw() + ggplot2::theme(axis.line.x = element_line(color="black", linewidth = .7),axis.line.y = element_line(color="black", size = .7),text=element_text(size=16))
mypalette <- colorRampPalette(c("#ff5858","#53f1f1","#c5fd57")) # http://paletton.com
hist(discoveries, col=mypalette(3))
getwd()



############ ############# DE ANALYSIS by each symbiont ############# #############
directory.b <- "./burk_4_bt2_count/"
sampleFiles.b <- grep("count",list.files(directory.b),value=TRUE)
sampleName.b <- sub("(*).bt2.count","\\1",sampleFiles.b)
sampleTable.b <- data.frame(sampleName = sampleName.b,
                            fileName = sampleFiles.b,
                            condition = c("control","control","control","intracellular","intracellular","intracellular",
                                          "control","control","control","intracellular","intracellular","intracellular",
                                          "control","control","control","intracellular","intracellular","intracellular", 
                                          "control","control","control","intracellular","intracellular","intracellular"),
                            dicty = rep(c(rep("none",3), "qs18", "qs864", "qs9"),4),
                            burk = c(rep("b11",6), rep("b159",6), rep("b70",6), rep("b859",6)),
                            batch = c("1952_7","1966_1","1966_3","1966_1","1966_3","1952_7",
                                      "1952_7","1966_1","1966_3","1966_1","1966_3","1952_7",
                                      "1952_8","1966_2","1966_4","1966_2","1966_4","1952_8",
                                      "1952_8","1966_2","1966_4","1966_2","1966_4","1952_8"),
                            clade = c("none","none","none","reduced","reduced","reduced",
                                      "none","none","none","nonreduced","nonreduced","nonreduced",
                                      "none","none","none","reduced","reduced","reduced",
                                      "none","none","none","nonreduced","nonreduced","nonreduced"))

sampleTable.b$condition <- factor(sampleTable.b$condition)


# individual strain DE
dds.b159 <- DESeqDataSetFromHTSeqCount(sampleTable = subset(sampleTable.b, burk=="b159"),
                                       directory = directory.b,
                                       design= ~condition+batch)
colSums(counts(dds.b159) >= 10) / nrow(counts(dds.b159))
keep <- rowSums(counts(dds.b159) >= 10) >= 3
table(keep) # false 561, true 7250
dds.b159 <- dds.b159[keep,]

dds.b159 <- estimateSizeFactors(dds.b159)
dds.b159 <- DESeq(dds.b159)
res.b159 <- results(dds.b159, contrast=c("condition","intracellular","control"), alpha=0.05)
summary(res.b159) # 25% up, 22% down
#res.b159 <- res.b159[order(res.b159$padj),]


dds.b70 <- DESeqDataSetFromHTSeqCount(sampleTable = subset(sampleTable.b, burk=="b70"),
                                      directory = directory.b,
                                      design= ~condition+batch)
colSums(counts(dds.b70) >= 10) / nrow(counts(dds.b70))
keep <- rowSums(counts(dds.b70) >= 10) >= 3
table(keep) # false 980, true 6971
dds.b70 <- dds.b70[keep,]

dds.b70 <- estimateSizeFactors(dds.b70)
dds.b70 <- DESeq(dds.b70)
res.b70 <- results(dds.b70, contrast=c("condition","intracellular","control"), alpha=0.05)
summary(res.b70) # 28% up, 29% down
#res.b70 <- res.b70[order(res.b70$padj),]


dds.b859 <- DESeqDataSetFromHTSeqCount(sampleTable = subset(sampleTable.b, burk=="b859"),
                                       directory = directory.b,
                                       design= ~condition+batch)
colSums(counts(dds.b859) >= 10) / nrow(counts(dds.b859))
keep <- rowSums(counts(dds.b859) >= 10) >= 3
table(keep) # false 218, true 3382
dds.b859 <- dds.b859[keep,]

dds.b859 <- estimateSizeFactors(dds.b859)
dds.b859 <- DESeq(dds.b859)
res.b859 <- results(dds.b859, contrast=c("condition","intracellular","control"), alpha=0.05)
summary(res.b859) # 22% up, 20% down
#res.b859 <- res.b859[order(res.b859$padj),]


dds.b11 <- DESeqDataSetFromHTSeqCount(sampleTable = subset(sampleTable.b, burk=="b11"),
                                      directory = directory.b,
                                      design= ~condition+batch)
colSums(counts(dds.b11) >= 10) / nrow(counts(dds.b11))
keep <- rowSums(counts(dds.b11) >= 10) >= 3
table(keep) # false 133, true 3553
dds.b11 <- dds.b11[keep,]

dds.b11 <- estimateSizeFactors(dds.b11)
dds.b11 <- DESeq(dds.b11)
res.b11 <- results(dds.b11, contrast=c("condition","intracellular","control"), alpha=0.05)
summary(res.b11) # 21% up, 20% down
#res.b11 <- res.b11[order(res.b11$padj),]


mcols(res.b11)$description
lfc.b11 <- lfcShrink(dds.b11, coef="condition_intracellular_vs_control")
plotMA(lfc.b11)

lfc.b70 <- lfcShrink(dds.b70, coef="condition_intracellular_vs_control")
plotMA(lfc.b70)


# write to flat tables
write.table(res.b159, file="burk159.DEresults.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(res.b70, file="burk70.DEresults.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(res.b859, file="burk859.DEresults.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(res.b11, file="burk11.DEresults.txt", quote=F, sep="\t", row.names=T, col.names=T)

res.b159 <- read.table(file="burk159.DEresults.txt", sep="\t", h=T)
res.b70 <- read.table(file="burk70.DEresults.txt", sep="\t", h=T)
res.b859 <- read.table(file="burk859.DEresults.txt", sep="\t", h=T)
res.b11 <- read.table(file="burk11.DEresults.txt", sep="\t", h=T)


# run roary to make it simpler to compare across. 
# also used roary to port img go term annotations (from interproscan) to prokka annotations


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
elim.stat <- new("elimCount", testStatistic=GOFisherTest, name="Fisher (elim)")


#### AGRICOLARIS BA159 #### 
# read in go mappings
raw.map <- read.table("ba159.gene_association.img_topgo.txt", h=F)
colnames(raw.map) <- c("gene_id","go_id")

go_map.b159 <- by(raw.map$go_id, raw.map$gene_id, function(x) as.character(x))
str(head(go_map.b159))

# define gene lists
up <- gsub("_gene","",row.names(subset(res.b159, padj<=0.05 & log2FoldChange>0)))
down <- gsub("_gene","",row.names(subset(res.b159, padj<=0.05 & log2FoldChange<0)))
universe <- gsub("_gene","",row.names(subset(res.b159, is.na(padj)==FALSE)))

# define base data structure
geneList <- subset(res.b159, is.na(padj)==FALSE)[,6]
names(geneList) <- universe
str(geneList)

# set up test data structure
BP.b159data <- new("topGOdata", description="GO BP analysis of  All Genes", ontology="BP", allGenes=geneList , geneSel=diffGene, annot = annFUN.gene2GO, gene2GO=go_map.b159, nodeSize=10)

BP.b159up <- new("topGOdata", description="GO BP analysis of  Up Genes", ontology="BP", allGenes=geneList , geneSel=upGene, annot = annFUN.gene2GO, gene2GO=go_map.b159, nodeSize=10)

BP.b159down <- new("topGOdata", description="GO BP analysis of  Down Genes", ontology="BP", allGenes=geneList , geneSel=downGene, annot = annFUN.gene2GO, gene2GO=go_map.b159, nodeSize=10)

# apply test
all.resultElim <- getSigGroups(BP.b159data, elim.stat)
all.resultElim # 25 @ 0.05 sig terms, 5 @ 0.01
GenTable(BP.b159data, elim=all.resultElim, topNodes=40)

up.resultElim <- getSigGroups(BP.b159up, elim.stat)
up.resultElim # 5 sig terms @ 0.05, 4 @ 0.01
GenTable(BP.b159up, elim=up.resultElim, topNodes=30)

down.resultElim <- getSigGroups(BP.b159down, elim.stat)
down.resultElim # 25 sig terms @ 0.05, 12 @ 0.01
GenTable(BP.b159down, elim=down.resultElim, topNodes=50)


all.b159Res <- GenTable(BP.b159data, elim=all.resultElim, orderBy="elim", ranksOf="elim", topNodes=25, numChar=99)
up.b159Res <- GenTable(BP.b159up, elim=up.resultElim, orderBy="elim", ranksOf="elim", topNodes=5, numChar=99)
down.b159Res <- GenTable(BP.b159down, elim=down.resultElim, orderBy="elim", ranksOf="elim", topNodes=25, numChar=99)

# save results files for go-figure
write.table(all.b159Res, file="burk159.allreg.bt2.topGO.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(up.b159Res, file="burk159.upreg.bt2.topGO.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(down.b159Res, file="burk159.downreg.bt2.topGO.txt", quote=F, sep="\t", row.names=T, col.names=T)



#### AGRICOLARIS BA70 #### 
# read in go mappings
raw.map <- read.table("ba70.gene_association.img_topgo.txt", h=F)
colnames(raw.map) <- c("gene_id","go_id")

go_map.b70 <- by(raw.map$go_id, raw.map$gene_id, function(x) as.character(x))
str(head(go_map.b70))

# define gene lists
up <- gsub("_gene","",row.names(subset(res.b70, padj<=0.05 & log2FoldChange>0)))
down <- gsub("_gene","",row.names(subset(res.b70, padj<=0.05 & log2FoldChange<0)))
universe <- gsub("_gene","",row.names(subset(res.b70, is.na(padj)==FALSE)))

# define base data structure
geneList <- subset(res.b70, is.na(padj)==FALSE)[,6]
names(geneList) <- universe
str(geneList)

# set up test data structure
BP.b70data <- new("topGOdata", description="GO BP analysis of  All Genes", ontology="BP", allGenes=geneList , geneSel=diffGene, annot = annFUN.gene2GO, gene2GO=go_map.b70, nodeSize=10)

BP.b70up <- new("topGOdata", description="GO BP analysis of  All Genes", ontology="BP", allGenes=geneList , geneSel=upGene, annot = annFUN.gene2GO, gene2GO=go_map.b70, nodeSize=10)

BP.b70down <- new("topGOdata", description="GO BP analysis of  All Genes", ontology="BP", allGenes=geneList , geneSel=downGene, annot = annFUN.gene2GO, gene2GO=go_map.b70, nodeSize=10)

# apply test
all.resultElim <- getSigGroups(BP.b70data, elim.stat)
all.resultElim # 23 @ 0.05 sig terms, 10 @ 0.01
GenTable(BP.b70data, elim=all.resultElim, topNodes=50)

up.resultElim <- getSigGroups(BP.b70up, elim.stat)
up.resultElim # 5 sig terms @ 0.05, 4 @ 0.01
GenTable(BP.b70up, elim=up.resultElim, topNodes=30)

down.resultElim <- getSigGroups(BP.b70down, elim.stat)
down.resultElim # 39 sig terms @ 0.05, 15 @ 0.01
GenTable(BP.b70down, elim=down.resultElim, topNodes=50)

# get and save results
all.b70Res <- GenTable(BP.b70data, elim=all.resultElim, orderBy="elim", ranksOf="elim", topNodes=23, numChar=99)
up.b70Res <- GenTable(BP.b70up, elim=up.resultElim, orderBy="elim", ranksOf="elim", topNodes=5, numChar=99)
down.b70Res <- GenTable(BP.b70down, elim=down.resultElim, orderBy="elim", ranksOf="elim", topNodes=39, numChar=99)

write.table(all.b70Res, file="burk70.allreg.bt2.topGO.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(up.b70Res, file="burk70.upreg.bt2.topGO.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(down.b70Res, file="burk70.downreg.bt2.topGO.txt", quote=F, sep="\t", row.names=T, col.names=T)


#### BONNIEA BB859 #### 
# read in go mappings
raw.map <- read.table("bb859.gene_association.img_topgo.txt", h=F)
colnames(raw.map) <- c("gene_id","go_id")

go_map.b859 <- by(raw.map$go_id, raw.map$gene_id, function(x) as.character(x))
str(head(go_map.b859))

# define gene lists
up <- gsub("_gene","",row.names(subset(res.b859, padj<=0.05 & log2FoldChange>0)))
down <- gsub("_gene","",row.names(subset(res.b859, padj<=0.05 & log2FoldChange<0)))
universe <- gsub("_gene","",row.names(subset(res.b859, is.na(padj)==FALSE)))

# define base data structure
geneList <- subset(res.b859, is.na(padj)==FALSE)[,6]
names(geneList) <- universe
str(geneList)

# set up test data structure
BP.b859data <- new("topGOdata", description="GO BP analysis of  All Genes", ontology="BP", allGenes=geneList , geneSel=diffGene, annot = annFUN.gene2GO, gene2GO=go_map.b859, nodeSize=1)

BP.b859up <- new("topGOdata", description="GO BP analysis of  All Genes", ontology="BP", allGenes=geneList , geneSel=upGene, annot = annFUN.gene2GO, gene2GO=go_map.b859, nodeSize=1)

BP.b859down <- new("topGOdata", description="GO BP analysis of  All Genes", ontology="BP", allGenes=geneList , geneSel=downGene, annot = annFUN.gene2GO, gene2GO=go_map.b859, nodeSize=1)

# apply test
all.resultElim <- getSigGroups(BP.b859data, elim.stat)
all.resultElim # 0 @ 0.05 sig terms, 0 @ 0.01
GenTable(BP.b859data, elim=all.resultElim, topNodes=30)

up.resultElim <- getSigGroups(BP.b859up, elim.stat)
up.resultElim # 8 sig terms @ 0.05, 3 @ 0.01
GenTable(BP.b859up, elim=up.resultElim, topNodes=30)

down.resultElim <- getSigGroups(BP.b859down, elim.stat)
down.resultElim # 8 sig terms @ 0.05, 4 @ 0.01
GenTable(BP.b859down, elim=down.resultElim, topNodes=50)

# get and save results
#all.b859Res <- GenTable(BP.b859data, elim=all.resultElim, orderBy="elim", ranksOf="elim", topNodes=0, numChar=99)
up.b859Res <- GenTable(BP.b859up, elim=up.resultElim, orderBy="elim", ranksOf="elim", topNodes=8, numChar=99)
down.b859Res <- GenTable(BP.b859down, elim=down.resultElim, orderBy="elim", ranksOf="elim", topNodes=8, numChar=99)

#write.table(all.b859Res, file="burk859.allreg.bt2.topGO.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(up.b859Res, file="burk859.upreg.bt2.topGO.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(down.b859Res, file="burk859.downreg.bt2.topGO.txt", quote=F, sep="\t", row.names=T, col.names=T)


#### HAYLEYELLA BH11 #### 
# read in go mappings
raw.map <- read.table("bh11.gene_association.img_topgo.txt", h=F)
colnames(raw.map) <- c("gene_id","go_id")

go_map.b11 <- by(raw.map$go_id, raw.map$gene_id, function(x) as.character(x))
str(head(go_map.b11))

# define gene lists
up <- gsub("_gene","",row.names(subset(res.b11, padj<=0.05 & log2FoldChange>0)))
down <- gsub("_gene","",row.names(subset(res.b11, padj<=0.05 & log2FoldChange<0)))
universe <- gsub("_gene","",row.names(subset(res.b11, is.na(padj)==FALSE)))

# define base data structure
geneList <- subset(res.b11, is.na(padj)==FALSE)[,6]
names(geneList) <- universe
str(geneList)

# set up test data structure
BP.b11data <- new("topGOdata", description="GO BP analysis of  All Genes", ontology="BP", allGenes=geneList , geneSel=diffGene, annot = annFUN.gene2GO, gene2GO=go_map.b11, nodeSize=1)

BP.b11up <- new("topGOdata", description="GO BP analysis of  All Genes", ontology="BP", allGenes=geneList , geneSel=upGene, annot = annFUN.gene2GO, gene2GO=go_map.b11, nodeSize=1)

BP.b11down <- new("topGOdata", description="GO BP analysis of  All Genes", ontology="BP", allGenes=geneList , geneSel=downGene, annot = annFUN.gene2GO, gene2GO=go_map.b11, nodeSize=1)

# apply test
all.resultElim <- getSigGroups(BP.b11data, elim.stat)
all.resultElim # 7 @ 0.05 sig terms, 0 @ 0.01
GenTable(BP.b11data, elim=all.resultElim, topNodes=30)

up.resultElim <- getSigGroups(BP.b11up, elim.stat)
up.resultElim # 15 sig terms @ 0.05, 2 @ 0.01
GenTable(BP.b11up, elim=up.resultElim, topNodes=30)

down.resultElim <- getSigGroups(BP.b11down, elim.stat)
down.resultElim # 12 sig terms @ 0.05, 0 @ 0.01
GenTable(BP.b11down, elim=down.resultElim, topNodes=30)

# get and save results
all.b11Res <- GenTable(BP.b11data, elim=all.resultElim, orderBy="elim", ranksOf="elim", topNodes=7, numChar=99)
up.b11Res <- GenTable(BP.b11up, elim=up.resultElim, orderBy="elim", ranksOf="elim", topNodes=15, numChar=99)
down.b11Res <- GenTable(BP.b11down, elim=down.resultElim, orderBy="elim", ranksOf="elim", topNodes=12, numChar=99)

write.table(all.b11Res, file="burk11.allreg.bt2.topGO.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(up.b11Res, file="burk11.upreg.bt2.topGO.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(down.b11Res, file="burk11.downreg.bt2.topGO.txt", quote=F, sep="\t", row.names=T, col.names=T)

rm(all.b11Res, up.b11Res, down.b11Res, all.b159Res, up.b159Res, down.b159Res, all.b70Res, up.b70Res, down.b70Res, up.b859Res, down.b859Res)


#### COMPARISON ACROSS STRAINS ####
go.b159u <- read.table("b159.upreg.bt2.topGO.txt", h=T, sep="\t")
go.b70u <- read.table("b70.upreg.bt2.topGO.txt", h=T, sep="\t")
go.b859u <- read.table("b859.upreg.bt2.topGO.txt", h=T, sep="\t")
go.b11u <- read.table("b11.upreg.bt2.topGO.txt", h=T, sep="\t")

go.b159d <- read.table("b159.downreg.bt2.topGO.txt", h=T, sep="\t")
go.b70d <- read.table("b70.downreg.bt2.topGO.txt", h=T, sep="\t")
go.b859d <- read.table("b859.downreg.bt2.topGO.txt", h=T, sep="\t")
go.b11d <- read.table("b11.downreg.bt2.topGO.txt", h=T, sep="\t")


# combine GO data (enrichment true/false)
go.burk <- data.frame(matrix(NA, ncol=1, nrow=81))
row.names(go.burk) <- unique(sort(c(go.b159u$GO.ID, go.b70u$GO.ID, go.b859u$GO.ID, go.b11u$GO.ID, go.b159d$GO.ID, go.b70d$GO.ID, go.b859d$GO.ID, go.b11d$GO.ID)))
go.burk$up.b159 <- row.names(go.burk) %in% go.b159u$GO.ID
go.burk$up.b70 <- row.names(go.burk) %in% go.b70u$GO.ID
go.burk$up.b859 <- row.names(go.burk) %in% go.b859u$GO.ID
go.burk$up.b11 <- row.names(go.burk) %in% go.b11u$GO.ID
go.burk$down.b159 <- row.names(go.burk) %in% go.b159d$GO.ID
go.burk$down.b70 <- row.names(go.burk) %in% go.b70d$GO.ID
go.burk$down.b859 <- row.names(go.burk) %in% go.b859d$GO.ID
go.burk$down.b11 <- row.names(go.burk) %in% go.b11d$GO.ID
temp <- go.burk[,2:9]
go.burk <- temp

limma::vennDiagram(go.burk[,1:4]) 
limma::vennDiagram(go.burk[,1:2]) 
limma::vennDiagram(go.burk[,3:4]) 

limma::vennDiagram(go.burk[,5:8]) 
limma::vennDiagram(go.burk[,5:6]) 

# nonreduced only shared GO terms
write.table(subset(go.b159u, GO.ID %in% row.names(subset(go.burk, up.b159==1 & up.b70==1)))[,1:2], file="burkGO.nonreduced.up.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(subset(go.b159d, GO.ID %in% row.names(subset(go.burk, down.b159==1 & down.b70==1)))[,1:2], file="burkGO.nonreduced.down.txt", quote=F, sep="\t", row.names=T, col.names=T)

write.table(subset(go.b859u, GO.ID %in% row.names(subset(go.burk, up.b859==1 & up.b11==1)))[,1:2], file="burkGO.reduced.up.txt", quote=F, sep="\t", row.names=T, col.names=T)



############ ############# SECRETION SYSTEM AND EUKARYOTIC DOMAIN GENES ANALYSIS ############# #############
## manual hypergeometric test
# order: (gene list in target list, universe genes in target list, universe genes not in target list, gene list size) 
hg.test <- function(a,b,c,d) { min(1-cumsum(dhyper(0:(a-1),b,c,d))) }

genes.b159up <- gsub("_gene","",row.names(subset(res.b159, padj<=0.05 & log2FoldChange>0)))
genes.b70up <- gsub("_gene","",row.names(subset(res.b70, padj<=0.05 & log2FoldChange>0)))
genes.b11up <- gsub("_gene","",row.names(subset(res.b11, padj<=0.05 & log2FoldChange>0)))
genes.b859up <- gsub("_gene","",row.names(subset(res.b859, padj<=0.05 & log2FoldChange>0)))

genes.b159down <- gsub("_gene","",row.names(subset(res.b159, padj<=0.05 & log2FoldChange<0)))
genes.b70down <- gsub("_gene","",row.names(subset(res.b70, padj<=0.05 & log2FoldChange<0)))
genes.b11down <- gsub("_gene","",row.names(subset(res.b11, padj<=0.05 & log2FoldChange<0)))
genes.b859down <- gsub("_gene","",row.names(subset(res.b859, padj<=0.05 & log2FoldChange<0)))

universe.b159 <- gsub("_gene","",row.names(subset(res.b159, is.na(padj)==FALSE)))
universe.b70 <- gsub("_gene","",row.names(subset(res.b70, is.na(padj)==FALSE)))
universe.b11 <- gsub("_gene","",row.names(subset(res.b11, is.na(padj)==FALSE)))
universe.b859 <- gsub("_gene","",row.names(subset(res.b859, is.na(padj)==FALSE)))


# secretion system genes
b159.sec <- read.table("ba159.secretion.system.prokka.txt", h=F, sep="\t")
names(b159.sec) <- c("id","product")

b70.sec <- read.table("ba70.secretion.system.prokka.txt", h=F, sep="\t")
names(b70.sec) <- c("id","product")

b859.sec <- read.table("bb859.secretion.system.prokka.txt", h=F, sep="\t")
names(b859.sec) <- c("id","product")

b11.sec <- read.table("bh11.secretion.system.prokka.txt", h=F, sep="\t")
names(b11.sec) <- c("id","product")

# are secretion system genes found in DE genes? most upregulated, usually not downregulated
#b159.sec[b159.sec$gene %in% genes.b159down,2]
#b70.sec[b70.sec$gene %in% genes.b70down,2]
#b859.sec[b859.sec$gene %in% genes.b859down,2]
#b11.sec[b11.sec$gene %in% genes.b11down,2]


b159.sec[b159.sec$id %in% genes.b159up,] 
b70.sec[b70.sec$id %in% genes.b70up,]
b859.sec[b859.sec$id %in% genes.b859up,]
b11.sec[b11.sec$id %in% genes.b11up,]


# genes with eukaryotic protein domains
b159.euk <- read.table("ba159.eukaryotic.pfam.prokka.txt", h=F, sep="\t")
names(b159.euk) <- c("pfam","domain","id","product")

b70.euk <- read.table("ba70.eukaryotic.pfam.prokka.txt", h=F, sep="\t")
names(b70.euk) <- c("pfam","domain","id","product")

b859.euk <- read.table("bb859.eukaryotic.pfam.prokka.txt", h=F, sep="\t")
names(b859.euk) <- c("pfam","domain","id","product")

b11.euk <- read.table("bh11.eukaryotic.pfam.prokka.txt", h=F, sep="\t")
names(b11.euk) <- c("pfam","domain","id","product")

#b159.euk[b159.euk$gene %in% genes.b159down,2]
#b70.euk[b70.euk$gene %in% genes.b70down,2]
#b859.euk[b859.euk$gene %in% genes.b859down,2]
#b11.euk[b11.euk$gene %in% genes.b11down,2]

b159.euk[b159.euk$id %in% genes.b159up,] 
b70.euk[b70.euk$id %in% genes.b70up,]
b859.euk[b859.euk$id %in% genes.b859up,]
b11.euk[b11.euk$id %in% genes.b11up,]


## set up hyperG test order: (gene list in target list, universe genes in target list, universe genes not in target list, gene list size)
#b159
mylist.1 <- genes.b159up[genes.b159up %in% b159.sec$id]
mylist.2 <- genes.b159down[genes.b159down %in% b159.sec$id]
checklist.1 <- universe.b159[universe.b159 %in% b159.sec$id]

mylist.3 <- genes.b159up[genes.b159up %in% b159.euk$id]
mylist.4 <- genes.b159down[genes.b159down %in% b159.euk$id]
checklist.2 <- universe.b159[universe.b159 %in% b159.euk$id]

b159.hg.table <- data.frame(Term=c("Secretion system genes (up)","Secretion system genes (down)","Eukaryotic Pfam genes (up)","Eukaryotic Pfam genes (down)"),
                            ExpCount=c(length(genes.b159up)*length(checklist.1)/length(universe.b159), length(genes.b159down)*length(checklist.1)/length(universe.b159), length(genes.b159up)*length(checklist.2)/length(universe.b159),length(genes.b159down)*length(checklist.2)/length(universe.b159)),
                           ObsCount=c(length(mylist.1), length(mylist.2), length(mylist.3), length(mylist.4)),
                           SizeTerm=c(nrow(b159.sec), nrow(b159.sec), nrow(b159.euk), nrow(b159.euk)))

b159.hg.table$raw.p <- c(hg.test(length(mylist.1), length(checklist.1), length(universe.b159)-length(checklist.1), length(genes.b159up)),
                        hg.test(length(mylist.2), length(checklist.1), length(universe.b159)-length(checklist.1), length(genes.b159down)),
                        hg.test(length(mylist.3), length(checklist.2), length(universe.b159)-length(checklist.2), length(genes.b159up)),
                        hg.test(length(mylist.4), length(checklist.2), length(universe.b159)-length(checklist.2), length(genes.b159down)))
b159.hg.table$adj.p <- p.adjust(b159.hg.table$raw.p, method="BH", n=length(b159.hg.table$raw.p))


#b70
mylist.1 <- genes.b70up[genes.b70up %in% b70.sec$id]
mylist.2 <- genes.b70down[genes.b70down %in% b70.sec$id]
checklist.1 <- universe.b70[universe.b70 %in% b70.sec$id]

mylist.3 <- genes.b70up[genes.b70up %in% b70.euk$id]
mylist.4 <- genes.b70down[genes.b70down %in% b70.euk$id]
checklist.2 <- universe.b70[universe.b70 %in% b70.euk$id]

b70.hg.table <- data.frame(Term=c("Secretion system genes (up)","Secretion system genes (down)","Eukaryotic Pfam genes (up)","Eukaryotic Pfam genes (down)"),
                           ExpCount=c(length(genes.b70up)*length(checklist.1)/length(universe.b70), length(genes.b70down)*length(checklist.1)/length(universe.b70), length(genes.b70up)*length(checklist.2)/length(universe.b70),length(genes.b70down)*length(checklist.2)/length(universe.b70)),
                           ObsCount=c(length(mylist.1), length(mylist.2), length(mylist.3), length(mylist.4)),
                           SizeTerm=c(nrow(b70.sec), nrow(b70.sec), nrow(b70.euk), nrow(b70.euk)))

b70.hg.table$raw.p <- c(hg.test(length(mylist.1), length(checklist.1), length(universe.b70)-length(checklist.1), length(genes.b70up)),
                         hg.test(length(mylist.2), length(checklist.1), length(universe.b70)-length(checklist.1), length(genes.b70down)),
                         hg.test(length(mylist.3), length(checklist.2), length(universe.b70)-length(checklist.2), length(genes.b70up)),
                         hg.test(length(mylist.4), length(checklist.2), length(universe.b70)-length(checklist.2), length(genes.b70down)))
b70.hg.table$adj.p <- p.adjust(b70.hg.table$raw.p, method="BH", n=length(b70.hg.table$raw.p))


#b859
mylist.1 <- genes.b859up[genes.b859up %in% b859.sec$id]
mylist.2 <- genes.b859down[genes.b859down %in% b859.sec$id]
checklist.1 <- universe.b859[universe.b859 %in% b859.sec$id]

mylist.3 <- genes.b859up[genes.b859up %in% b859.euk$id]
mylist.4 <- genes.b859down[genes.b859down %in% b859.euk$id]
checklist.2 <- universe.b859[universe.b859 %in% b859.euk$id]

b859.hg.table <- data.frame(Term=c("Secretion system genes (up)","Secretion system genes (down)","Eukaryotic Pfam genes (up)","Eukaryotic Pfam genes (down)"),
                            ExpCount=c(length(genes.b859up)*length(checklist.1)/length(universe.b859), length(genes.b859down)*length(checklist.1)/length(universe.b859), length(genes.b859up)*length(checklist.2)/length(universe.b859),length(genes.b859down)*length(checklist.2)/length(universe.b859)),
                           ObsCount=c(length(mylist.1), length(mylist.2), length(mylist.3), length(mylist.4)),
                           SizeTerm=c(nrow(b859.sec), nrow(b859.sec), nrow(b859.euk), nrow(b859.euk)))

b859.hg.table$raw.p <- c(hg.test(length(mylist.1), length(checklist.1), length(universe.b859)-length(checklist.1), length(genes.b859up)),
                         hg.test(length(mylist.2), length(checklist.1), length(universe.b859)-length(checklist.1), length(genes.b859down)),
                         hg.test(length(mylist.3), length(checklist.2), length(universe.b859)-length(checklist.2), length(genes.b859up)),
                         hg.test(length(mylist.4), length(checklist.2), length(universe.b859)-length(checklist.2), length(genes.b859down)))
b859.hg.table$adj.p <- p.adjust(b859.hg.table$raw.p, method="BH", n=length(b859.hg.table$raw.p))


# b11 
mylist.1 <- genes.b11up[genes.b11up %in% b11.sec$id]
mylist.2 <- genes.b11down[genes.b11down %in% b11.sec$id]
checklist.1 <- universe.b11[universe.b11 %in% b11.sec$id]

mylist.3 <- genes.b11up[genes.b11up %in% b11.euk$id]
mylist.4 <- genes.b11down[genes.b11down %in% b11.euk$id]
checklist.2 <- universe.b11[universe.b11 %in% b11.euk$id]

b11.hg.table <- data.frame(Term=c("Secretion system genes (up)","Secretion system genes (down)","Eukaryotic Pfam genes (up)","Eukaryotic Pfam genes (down)"),
                           ExpCount=c(length(genes.b11up)*length(checklist.1)/length(universe.b11), length(genes.b11down)*length(checklist.1)/length(universe.b11), length(genes.b11up)*length(checklist.2)/length(universe.b11),length(genes.b11down)*length(checklist.2)/length(universe.b11)),
                           ObsCount=c(length(mylist.1), length(mylist.2), length(mylist.3), length(mylist.4)),
                           SizeTerm=c(nrow(b11.sec), nrow(b11.sec), nrow(b11.euk), nrow(b11.euk)))

b11.hg.table$raw.p <- c(hg.test(length(mylist.1), length(checklist.1), length(universe.b11)-length(checklist.1), length(genes.b11up)),
                        hg.test(length(mylist.2), length(checklist.1), length(universe.b11)-length(checklist.1), length(genes.b11down)),
                        hg.test(length(mylist.3), length(checklist.2), length(universe.b11)-length(checklist.2), length(genes.b11up)),
                        hg.test(length(mylist.4), length(checklist.2), length(universe.b11)-length(checklist.2), length(genes.b11down)))
b11.hg.table$adj.p <- p.adjust(b11.hg.table$raw.p, method="BH", n=length(b11.hg.table$raw.p))


b159.hg.table
b70.hg.table
b859.hg.table
b11.hg.table


rm(mylist.1, mylist.2, mylist.3, mylist.4, checklist.1, checklist.2)


############ ############# NETWORK ANALYSIS ############# #############
library(igraph)

# these look fundamentally different from host networks, with clear central nodes and other nodes fanning out
b159_for_igraph <- varianceStabilizingTransformation(dds.b159, blind=F)
b70_for_igraph <- varianceStabilizingTransformation(dds.b70, blind=F)
b11_for_igraph <- varianceStabilizingTransformation(dds.b11, blind=F)
b859_for_igraph <- varianceStabilizingTransformation(dds.b859, blind=F)
head(assay(b159_for_igraph),1)

xsig.b159 <- row.names(subset(res.b159, padj<=0.05 & abs(log2FoldChange)>=4))
xsig.b70 <- row.names(subset(res.b70, padj<=0.05 & abs(log2FoldChange)>=4))
xsig.b11 <- row.names(subset(res.b11, padj<=0.05 & abs(log2FoldChange)>=4))
xsig.b859 <- row.names(subset(res.b859, padj<=0.05 & abs(log2FoldChange)>=4))

b159_igraph <- assay(b159_for_igraph)[xsig.b159,]
b70_igraph <- assay(b70_for_igraph)[xsig.b70,]
b11_igraph <- assay(b11_for_igraph)[xsig.b11,]
b859_igraph <- assay(b859_for_igraph)[xsig.b859,]

row.names(b159_igraph) <- gsub("_gene","",row.names(b159_igraph))
row.names(b70_igraph) <- gsub("_gene","",row.names(b70_igraph))
row.names(b11_igraph) <- gsub("_gene","",row.names(b11_igraph))
row.names(b859_igraph) <- gsub("_gene","",row.names(b859_igraph))

rm(b159_for_igraph, b70_for_igraph, b11_for_igraph, b859_for_igraph)

#### AGRICOLARIS BA159 ####
g.b159 <- graph_from_adjacency_matrix(as.matrix(as.dist(cor(t(b159_igraph), method="pearson"))), mode="undirected", weighted=TRUE, diag=FALSE)
g.b159 <- simplify(g.b159, remove.multiple=T, remove.loops=T)

E(g.b159)[which(E(g.b159)$weight<0)]$color <- "slategrey"
E(g.b159)[which(E(g.b159)$weight>0)]$color <- "red"

E(g.b159)$weight <- abs(E(g.b159)$weight)
g.b159 <- delete_edges(g.b159, E(g.b159)[which(E(g.b159)$weight<0.8)])
g.b159 <- delete_vertices(g.b159, degree(g.b159)==0)

V(g.b159)$name <- V(g.b159)$name
V(g.b159)$vertex.frame.color <- "white"

# Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst.b159 <- mst(g.b159, algorithm="prim")

# clustering
mst.b159.communities <- cluster_edge_betweenness(mst.b159, weights=NULL)
modularity(mst.b159.communities)
mst.b159.clustering <- make_clusters(mst.b159, membership=mst.b159.communities$membership)

b159.commSummary <- data.frame(mst.b159.communities$names, mst.b159.communities$membership, mst.b159.communities$modularity)
colnames(b159.commSummary) <- c("Gene", "Community", "Modularity")
b159.commSummary[order(b159.commSummary$Community),]
degree(mst.b159)[degree(mst.b159)>4]

b159.list <- ifelse(degree(mst.b159)>4, V(mst.b159)$name, NA)



#### AGRICOLARIS BA70 ####
g.b70 <- graph_from_adjacency_matrix(as.matrix(as.dist(cor(t(b70_igraph), method="pearson"))), mode="undirected", weighted=TRUE, diag=FALSE)
g.b70 <- simplify(g.b70, remove.multiple=T, remove.loops=T)

E(g.b70)[which(E(g.b70)$weight<0)]$color <- "slategrey"
E(g.b70)[which(E(g.b70)$weight>0)]$color <- "red"

E(g.b70)$weight <- abs(E(g.b70)$weight)
g.b70 <- delete_edges(g.b70, E(g.b70)[which(E(g.b70)$weight<0.8)])
g.b70 <- delete_vertices(g.b70, degree(g.b70)==0)

V(g.b70)$name <- V(g.b70)$name
V(g.b70)$vertex.frame.color <- "white"

# Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst.b70 <- mst(g.b70, algorithm="prim")

# clustering
mst.b70.communities <- cluster_edge_betweenness(mst.b70, weights=NULL)
modularity(mst.b70.communities)
mst.b70.clustering <- make_clusters(mst.b70, membership=mst.b70.communities$membership)

b70.commSummary <- data.frame(mst.b70.communities$names, mst.b70.communities$membership, mst.b70.communities$modularity)
colnames(b70.commSummary) <- c("Gene", "Community", "Modularity")
b70.commSummary[order(b70.commSummary$Community),]
degree(mst.b70)[degree(mst.b70)>6]

b70.list <- ifelse(degree(mst.b70)>5, V(mst.b70)$name, NA)



#### BONNIEA BB859 ####
g.b859 <- graph_from_adjacency_matrix(as.matrix(as.dist(cor(t(b859_igraph), method="pearson"))), mode="undirected", weighted=TRUE, diag=FALSE)
g.b859 <- simplify(g.b859, remove.multiple=T, remove.loops=T)

E(g.b859)[which(E(g.b859)$weight<0)]$color <- "slategrey"
E(g.b859)[which(E(g.b859)$weight>0)]$color <- "red"

E(g.b859)$weight <- abs(E(g.b859)$weight)
g.b859 <- delete_edges(g.b859, E(g.b859)[which(E(g.b859)$weight<0.8)])
g.b859 <- delete_vertices(g.b859, degree(g.b859)==0)

V(g.b859)$name <- V(g.b859)$name
V(g.b859)$vertex.frame.color <- "white"

# Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst.b859 <- mst(g.b859, algorithm="prim")

# clustering
mst.b859.communities <- cluster_edge_betweenness(mst.b859, weights=NULL)
modularity(mst.b859.communities)
mst.b859.clustering <- make_clusters(mst.b859, membership=mst.b859.communities$membership)

b859.commSummary <- data.frame(mst.b859.communities$names, mst.b859.communities$membership, mst.b859.communities$modularity)
colnames(b859.commSummary) <- c("Gene", "Community", "Modularity")
b859.commSummary[order(b859.commSummary$Community),]
degree(mst.b859)[degree(mst.b859)>1]

b859.list <- ifelse(degree(mst.b859)>1, V(mst.b859)$name, NA)



#### HAYLEYELLA BH11 ####
g.b11 <- graph_from_adjacency_matrix(as.matrix(as.dist(cor(t(b11_igraph), method="pearson"))), mode="undirected", weighted=TRUE, diag=FALSE)
g.b11 <- simplify(g.b11, remove.multiple=T, remove.loops=T)

E(g.b11)[which(E(g.b11)$weight<0)]$color <- "slategrey"
E(g.b11)[which(E(g.b11)$weight>0)]$color <- "red"

E(g.b11)$weight <- abs(E(g.b11)$weight)
g.b11 <- delete_edges(g.b11, E(g.b11)[which(E(g.b11)$weight<0.8)])
g.b11 <- delete_vertices(g.b11, degree(g.b11)==0)

V(g.b11)$name <- V(g.b11)$name
V(g.b11)$vertex.frame.color <- "white"

# Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst.b11 <- mst(g.b11, algorithm="prim")

# clustering
mst.b11.communities <- cluster_edge_betweenness(mst.b11, weights=NULL)
modularity(mst.b11.communities)
mst.b11.clustering <- make_clusters(mst.b11, membership=mst.b11.communities$membership)

b11.commSummary <- data.frame(mst.b11.communities$names, mst.b11.communities$membership, mst.b11.communities$modularity)
colnames(b11.commSummary) <- c("Gene", "Community", "Modularity")
b11.commSummary[order(b11.commSummary$Community),]
degree(mst.b11)[degree(mst.b11)>1]

b11.list <- ifelse(degree(mst.b11)>1, V(mst.b11)$name, NA)


#### NETWORK GRAPHS ####
# Plot the tree object using Fruchterman-Reingold force-directed algorithm
plot(mst.b159, layout=layout_with_fr, vertex.label=b159.list,
     vertex.size=4, vertex.label.dist=-1, vertex.label.color="slategrey", 
     asp=FALSE, vertex.label.cex=0.75, edge.arrow.mode=0, main="b159")

# internal vertices become default color despite previous coding
graph.b159 <- plot(mst.b159.clustering, mst.b159, layout=layout_with_fr, vertex.label=b159.list,
                   vertex.size=4, vertex.label.dist=-1, vertex.label.color="red", 
                   asp=FALSE, vertex.label.cex=0.75, edge.arrow.mode=0, main="b159")

# Plot the tree object using Fruchterman-Reingold force-directed algorithm
plot(mst.b70, layout=layout_with_fr, vertex.label=b70.list,
     vertex.size=4, vertex.label.dist=-1, vertex.label.color="slategrey", 
     asp=FALSE, vertex.label.cex=0.75, edge.arrow.mode=0, main="b70")

# internal vertices become default color despite previous coding
graph.b70 <- plot(mst.b70.clustering, mst.b70, layout=layout_with_fr, vertex.label=b70.list,
                  vertex.size=4, vertex.label.dist=-1, vertex.label.color="red", 
                  asp=FALSE, vertex.label.cex=0.75, edge.arrow.mode=0, main="b70")

# Plot the tree object using Fruchterman-Reingold force-directed algorithm
plot(mst.b859, layout=layout_with_fr, vertex.label=b859.list,
     vertex.size=4, vertex.label.dist=-1, vertex.label.color="slategrey", 
     asp=FALSE, vertex.label.cex=0.75, edge.arrow.mode=0, main="b859")

# internal vertices become default color despite previous coding
graph.b859 <- plot(mst.b859.clustering, mst.b859, layout=layout_with_fr, vertex.label=b859.list,
                   vertex.size=4, vertex.label.dist=-1, vertex.label.color="red", 
                   asp=FALSE, vertex.label.cex=0.75, edge.arrow.mode=0, main="b859")

# Plot the tree object using Fruchterman-Reingold force-directed algorithm
plot(mst.b11, layout=layout_with_fr, vertex.label=b11.list,
     vertex.size=4, vertex.label.dist=-1, vertex.label.color="slategrey", 
     asp=FALSE, vertex.label.cex=0.75, edge.arrow.mode=0, main="b11")

# internal vertices become default color despite previous coding
graph.b11 <- plot(mst.b11.clustering, mst.b11, layout=layout_with_fr, vertex.label=b11.list,
                  vertex.size=4, vertex.label.dist=-1, vertex.label.color="red", 
                  asp=FALSE, vertex.label.cex=0.75, edge.arrow.mode=0, main="b11")



#### COMPARE GENES WITH HIGH VERTEX DEGREES ####
# are these descriptive or is there something to evaluate? 
# Hub score (higher score = more influence), only works for directed graphs
# Vertex degree (higher score = more influence)

# burk
sort(degree(mst.b159)[degree(mst.b159)>1], decreasing=T)
sort(degree(mst.b70)[degree(mst.b70)>1], decreasing=T)
sort(degree(mst.b859)[degree(mst.b859)>1], decreasing=T)
sort(degree(mst.b11)[degree(mst.b11)>1], decreasing=T)

names(degree(mst.b159)[degree(mst.b159)>1])
names(degree(mst.b70)[degree(mst.b70)>1])
names(degree(mst.b859)[degree(mst.b859)>1])
names(degree(mst.b11)[degree(mst.b11)>1])

# after much interpro, uniprot blast and alignment,
# only good one is probably PBONN_01851, ortholog in hayl is also sig DE

# check if interesting genes are actually DE
# yes for all these in agri
"PAGRI_01162_gene" %in% xsig.b159 #row.names(subset(res.b159, padj<=0.05))
"PA070_02950_gene" %in% xsig.b70 #row.names(subset(res.b70, padj<=0.05)) # degree > 1
"PAGRI_01176_gene" %in% xsig.b159 #row.names(subset(res.b159, padj<=0.05))
"PA070_02936_gene" %in% xsig.b70 #row.names(subset(res.b70, padj<=0.05)) # degree > 1

degree(mst.b159) # all low degree 1
degree(mst.b70) # all low degrees 1-2

# no for hayl
"PBONN_01851_gene" %in% xsig.b859 #row.names(subset(res.b859, padj<=0.05)) # degree > 1
"PHAYL_01420_gene" %in% xsig.b11 #row.names(subset(res.b11, padj<=0.05))

# no for bonn
"PBONN_01921_gene" %in% xsig.b859 #row.names(subset(res.b859, padj<=0.05))
"PHAYL_01354_gene" %in% xsig.b11 #row.names(subset(res.b11, padj<=0.05)) # degree > 1

# additional ones, but all low degree
"PAGRI_01202_gene" %in% xsig.b159 #row.names(subset(res.b159, padj<=0.05))
"PA070_02910_gene" %in% xsig.b70 #row.names(subset(res.b70, padj<=0.05))


############ ############# GENE NEIGHBORHOOD ANALYSIS ############# #############
# load libraries
library(geneviewer)
library(webshot2)

# read in gene annotation file
gv <- read.table("ss_struct_eff_geneviewer.txt", h=F) # order manually set to T3SS and effectors, T6SS and effectors
gene_cluster <- data.frame(gv)[,c(5,2,3,1)]
colnames(gene_cluster) <- c("name","start","end","class")
str(gene_cluster)
gene_cluster$cluster <- as.factor(c(rep(1,27), rep(2,27), rep(3,27), rep(4,27)))

# generate plot
cluster_plot <- GC_chart(gene_cluster, group = "class", cluster = "cluster", width = "800px", height="900px") %>%
  GC_clusterFooter(
    title = c("<i>P. agricolaris Ba159</i>", "<i>P. agricolaris Ba70</i>", "<i>P. bonniea Bb859</i>", "<i>P. hayleyella Bh11</i>"), 
    subtitle = c("1: 1,236,929 - 1,299,223", "1: 3,352,704 - 3,412,354", "1: 2,085,666 - 2,153,038", "1: 1,465,207 - 1,528,262"),
    align = "left",
    x = 50) %>%
  GC_labels("name", fontSize="8px", fontStyle="regular") %>% 
  GC_genes(marker="boxarrow", marker_size="small") %>% GC_scaleBar(y = 20) %>%
  GC_scale(scale_breaks = TRUE, scale_break_threshold=17, hidden = TRUE)


# save plot to temp.html file for download
htmlwidgets::saveWidget(cluster_plot, "temp.html", selfcontained = TRUE)

