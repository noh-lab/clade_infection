#### Visualize GO terms ####
#### new version, streamlined from previous method

# dicty geneid2go
dicty.map <- read.table("geneid2go.dictyBase.20230102.filter.topgo.txt", h=F)
colnames(dicty.map) <- c("gene_id","go_id")

goterm_dict <- read.table("gobp_names.2025.txt", h=F, sep="\t") # from geneontology.org
colnames(goterm_dict) <- c("term", "name")
goterm_dict$term <- gsub(pattern = "GO:", replacement = "GO_", goterm_dict$term)
str(goterm_dict)

# read in results tables from topGO
dicty.b159u <- read.table("dicty159.upreg.topGO.txt", h=T, sep="\t")
dicty.b70u <- read.table("dicty70.upreg.topGO.txt", h=T, sep="\t")
dicty.b859u <- read.table("dicty859.upreg.topGO.txt", h=T, sep="\t")
dicty.b11u <- read.table("dicty11.upreg.topGO.txt", h=T, sep="\t")

dicty.b159d <- read.table("dicty159.downreg.topGO.txt", h=T, sep="\t")
dicty.b70d <- read.table("dicty70.downreg.topGO.txt", h=T, sep="\t")
dicty.b859d <- read.table("dicty859.downreg.topGO.txt", h=T, sep="\t")
dicty.b11d <- read.table("dicty11.downreg.topGO.txt", h=T, sep="\t")

# only consider sig enriched GO terms with 0.01 threshold (table contains 0.05 threshold)
dicty.b159u.list <- subset(dicty.b159u, elim<=0.01)[,1]
dicty.b70u.list <- subset(dicty.b70u, elim<=0.01)[,1]
dicty.b11u.list <- subset(dicty.b11u, elim<=0.01)[,1]
dicty.b859u.list <- subset(dicty.b859u, elim<=0.01)[,1]

dicty.b159d.list <- subset(dicty.b159d, elim<=0.01)[,1]
dicty.b70d.list <- subset(dicty.b70d, elim<=0.01)[,1]
dicty.b11d.list <- subset(dicty.b11d, elim<=0.01)[,1]
dicty.b859d.list <- subset(dicty.b859d, elim<=0.01)[,1]

# make a list of GO ids, and make an appropriate version for filenames
ids <- dicty.b159u.list
ids_clean <- gsub(":", "_", ids)

# for each list, find gene ids that belong to each go term
for(i in 1:length(ids)) {
  # Match IDs against dicty.map
  matched <- subset(dicty.map, go_id==ids[i])[,1]
  filename <- paste0(ids_clean[i], ".txt")
  if (is.na(matched[i]) || length(matched[i]) == 0) {
    next   # skip this iteration
  }
  # Save gene ids belonging to each go id as individual text files
  write.table(matched, file=filename, sep="\t", row.names=F, col.names=F, quote=F)
}

ids <- dicty.b70u.list
ids_clean <- gsub(":", "_", ids)

for(i in 1:length(ids)) {
  matched <- subset(dicty.map, go_id==ids[i])[,1]
  filename <- paste0(ids_clean[i], ".txt")
  if (is.na(matched[i]) || length(matched[i]) == 0) {
    next   # skip this iteration
  }
  write.table(matched, file=filename, sep="\t", row.names=F, col.names=F, quote=F)
}

ids <- dicty.b859u.list
ids_clean <- gsub(":", "_", ids)

for(i in 1:length(ids)) {
  matched <- subset(dicty.map, go_id==ids[i])[,1]
  filename <- paste0(ids_clean[i], ".txt")
  if (is.na(matched[i]) || length(matched[i]) == 0) {
    next   # skip this iteration
  }
  write.table(matched, file=filename, sep="\t", row.names=F, col.names=F, quote=F)
}

ids <- dicty.b11u.list
ids_clean <- gsub(":", "_", ids)

for(i in 1:length(ids)) {
  matched <- subset(dicty.map, go_id==ids[i])[,1]
  filename <- paste0(ids_clean[i], ".txt")
  if (is.na(matched[i]) || length(matched[i]) == 0) {
    next   # skip this iteration
  }
  write.table(matched, file=filename, sep="\t", row.names=F, col.names=F, quote=F)
}

ids <- dicty.b159d.list
ids_clean <- gsub(":", "_", ids)

# for each list, find gene ids that belong to each go term
for(i in 1:length(ids)) {
  matched <- subset(dicty.map, go_id==ids[i])[,1]
  filename <- paste0(ids_clean[i], ".txt")
  if (is.na(matched[i]) || length(matched[i]) == 0) {
    next   # skip this iteration
  }
  write.table(matched, file=filename, sep="\t", row.names=F, col.names=F, quote=F)
}

ids <- dicty.b70d.list
ids_clean <- gsub(":", "_", ids)

for(i in 1:length(ids)) {
  matched <- subset(dicty.map, go_id==ids[i])[,1]
  filename <- paste0(ids_clean[i], ".txt")
  if (is.na(matched[i]) || length(matched[i]) == 0) {
    next   # skip this iteration
  }
  write.table(matched, file=filename, sep="\t", row.names=F, col.names=F, quote=F)
}

ids <- dicty.b859d.list
ids_clean <- gsub(":", "_", ids)

for(i in 1:length(ids)) {
  matched <- subset(dicty.map, go_id==ids[i])[,1]
  filename <- paste0(ids_clean[i], ".txt")
  if (is.na(matched[i]) || length(matched[i]) == 0) {
    next   # skip this iteration
  }
  write.table(matched, file=filename, sep="\t", row.names=F, col.names=F, quote=F)
}

ids <- dicty.b11d.list
ids_clean <- gsub(":", "_", ids)

for(i in 1:length(ids)) {
  matched <- subset(dicty.map, go_id==ids[i])[,1]
  filename <- paste0(ids_clean[i], ".txt")
  if (is.na(matched[i]) || length(matched[i]) == 0) {
    next   # skip this iteration
  }
  write.table(matched, file=filename, sep="\t", row.names=F, col.names=F, quote=F)
}

rm(list = ls(pattern = "^GO_"))

# scan in each text file, e.g. GO_0000028 <- scan("GO_0000028.txt", what="character")
source("terms_and_members_for_visualization.R")

# put together the lists for each sig enriched go id, and the gene ids that belong to it
go.b159u.list <- list(GO_0000028 = GO_0000028,
                GO_0002181 = GO_0002181,
                GO_0006096 = GO_0006096,
                GO_0006099 = GO_0006099,
                GO_0006412 = GO_0006412,
                GO_0006869 = GO_0006869,
                GO_0006909 = GO_0006909,
                GO_0006955 = GO_0006955,
                GO_0009253 = GO_0009253,
                GO_0015031 = GO_0015031,
                GO_0015986 = GO_0015986,
                GO_0019835 = GO_0019835,
                GO_0042742 = GO_0042742,
                GO_0055085 = GO_0055085,
                GO_1902600 = GO_1902600,
                GO_2001235 = GO_2001235)

go.b70u.list <- list(GO_0000028 = GO_0000028,
                      GO_0000462 = GO_0000462,
                      GO_0006096 = GO_0006096,
                      GO_0006412 = GO_0006412,
                      GO_0006414 = GO_0006414,
                      GO_0006417 = GO_0006417,
                      GO_0006606 = GO_0006606,
                      GO_0006730 = GO_0006730,
                      GO_0006955 = GO_0006955,
                      GO_0009253 = GO_0009253,
                     GO_0030490 = GO_0030490,
                     GO_0042254 = GO_0042254,
                      GO_0042273 = GO_0042273,
                      GO_0042742 = GO_0042742,
                      GO_2001235 = GO_2001235)

go.b859u.list <- list(GO_0000122 = GO_0000122,
                      GO_0006096 = GO_0006096,
                      GO_0006412 = GO_0006412,
                      GO_0006413 = GO_0006413,
                      GO_0006730 = GO_0006730,
                      GO_0006886 = GO_0006886,
                      GO_0009253 = GO_0009253,
                      GO_0016126 = GO_0016126,
                      GO_0019835 = GO_0019835,
                      GO_0042742 = GO_0042742,
                      GO_1902600 = GO_1902600)

go.b11u.list <- list(GO_0000079 = GO_0000079,
                      GO_0007010 = GO_0007010,
                      GO_0007155 = GO_0007155,
                      GO_0009253 = GO_0009253,
                      GO_0019835 = GO_0019835,
                      GO_0042742 = GO_0042742)

go.b159d.list <- list(GO_0000070 = GO_0000070,
                      GO_0000723 = GO_0000723,
                      GO_0006270 = GO_0006270,
                      GO_0006302 = GO_0006302,
                      GO_0006633 = GO_0006633,
                      GO_0006974 = GO_0006974,
                      GO_0007062 = GO_0007062,
                      GO_0007064 = GO_0007064,
                      GO_0007131 = GO_0007131,
                      GO_0007165 = GO_0007165,
                      GO_0008360 = GO_0008360,
                      GO_0019953 = GO_0019953,
                      GO_0030261 = GO_0030261,
                      GO_0043157 = GO_0043157,
                      GO_0043547 = GO_0043547,
                      GO_0051301 = GO_0051301)
                   
go.b70d.list <- list(GO_0000079 = GO_0000079,
                     GO_0000723 = GO_0000723,
                     GO_0007131 = GO_0007131,
                     GO_0007165 = GO_0007165,
                     GO_0007166 = GO_0007166,
                     GO_0007265 = GO_0007265,
                     GO_0008360 = GO_0008360,
                     GO_0016567 = GO_0016567,
                     GO_0019953 = GO_0019953,
                     GO_0043157 = GO_0043157,
                     GO_0048856 = GO_0048856)

go.b859d.list <- list(GO_0000281 = GO_0000281,
                      GO_0007163 = GO_0007163,
                      GO_0007265 = GO_0007265,
                      GO_0008360 = GO_0008360,
                      GO_0019953 = GO_0019953,
                      GO_0030865 = GO_0030865,
                      GO_0034599 = GO_0034599,
                      GO_0043652 = GO_0043652,
                      GO_0045454 = GO_0045454,
                      GO_0098869 = GO_0098869)

go.b11d.list <- list(GO_0006221 = GO_0006221,
                     GO_0007264 = GO_0007264,
                     GO_0008360 = GO_0008360,
                     GO_0030150 = GO_0030150,
                     GO_0030490 = GO_0030490,
                     GO_0031152 = GO_0031152,
                     GO_0031167 = GO_0031167,
                     GO_0031288 = GO_0031288,
                     GO_0042274 = GO_0042274,
                     GO_0043652 = GO_0043652)


# get fold change data for each go term from results
# this step makes it so that only sig DE genes are included in the mean foldchange
sig.strain159u <- subset(res.strain159, padj<=0.05 & log2FoldChange>0)
sig.strain70u <- subset(res.strain70, padj<=0.05 & log2FoldChange>0)
sig.strain859u <- subset(res.strain859, padj<=0.05 & log2FoldChange>0)
sig.strain11u <- subset(res.strain11, padj<=0.05 & log2FoldChange>0)

sig.strain159d <- subset(res.strain159, padj<=0.05 & log2FoldChange<0)
sig.strain70d <- subset(res.strain70, padj<=0.05 & log2FoldChange<0)
sig.strain859d <- subset(res.strain859, padj<=0.05 & log2FoldChange<0)
sig.strain11d <- subset(res.strain11, padj<=0.05 & log2FoldChange<0)

# get the expression level of each gene id for each go term
dicty.b159u.exp = list() 
for (name in names(go.b159u.list)) {
  ids <- go.b159u.list[[name]]
  # Skip if fewer than 3 elements
  if (length(ids) < 3) next
  subset_dt <- sig.strain159u[rownames(sig.strain159u) %in% ids, 2]
  dicty.b159u.exp[[name]] <- subset_dt
}  

dicty.b70u.exp = list() 
for (name in names(go.b70u.list)) {
  ids <- go.b70u.list[[name]]
  # Skip if fewer than 3 elements
  if (length(ids) < 3) next
  subset_dt <- sig.strain70u[rownames(sig.strain70u) %in% ids, 2]
  dicty.b70u.exp[[name]] <- subset_dt
}  

dicty.b859u.exp = list() 
for (name in names(go.b859u.list)) {
  ids <- go.b859u.list[[name]]
  # Skip if fewer than 3 elements
  if (length(ids) < 3) next
  subset_dt <- sig.strain859u[rownames(sig.strain859u) %in% ids, 2]
  dicty.b859u.exp[[name]] <- subset_dt
}  

dicty.b11u.exp = list() 
for (name in names(go.b11u.list)) {
  ids <- go.b11u.list[[name]]
  # Skip if fewer than 3 elements
  if (length(ids) < 3) next
  subset_dt <- sig.strain11u[rownames(sig.strain11u) %in% ids, 2]
  dicty.b11u.exp[[name]] <- subset_dt
}  

dicty.b159d.exp = list() 
for (name in names(go.b159d.list)) {
  ids <- go.b159d.list[[name]]
  # Skip if fewer than 3 elements
  if (length(ids) < 3) next
  subset_dt <- sig.strain159d[rownames(sig.strain159d) %in% ids, 2]
  dicty.b159d.exp[[name]] <- subset_dt
}  

dicty.b70d.exp = list() 
for (name in names(go.b70d.list)) {
  ids <- go.b70d.list[[name]]
  # Skip if fewer than 3 elements
  if (length(ids) < 3) next
  subset_dt <- sig.strain70d[rownames(sig.strain70d) %in% ids, 2]
  dicty.b70d.exp[[name]] <- subset_dt
}  

dicty.b859d.exp = list() 
for (name in names(go.b859d.list)) {
  ids <- go.b859d.list[[name]]
  # Skip if fewer than 3 elements
  if (length(ids) < 3) next
  subset_dt <- sig.strain859d[rownames(sig.strain859d) %in% ids, 2]
  dicty.b859d.exp[[name]] <- subset_dt
}  

dicty.b11d.exp = list() 
for (name in names(go.b11d.list)) {
  ids <- go.b11d.list[[name]]
  # Skip if fewer than 3 elements
  if (length(ids) < 3) next
  subset_dt <- sig.strain11d[rownames(sig.strain11d) %in% ids, 2]
  dicty.b11d.exp[[name]] <- subset_dt
}  

# save to file to add GO term descriptions
dicty.b159u.fig <- ldply(dicty.b159u.exp, data.frame)
names(dicty.b159u.fig) <- c("term","log2FC") 
write.table(dicty.b159u.fig, "dicty.b159u.fig.txt", quote=F, sep="\t", row.names = F)

dicty.b70u.fig <- ldply(dicty.b70u.exp, data.frame)
names(dicty.b70u.fig) <- c("term","log2FC") 
write.table(dicty.b70u.fig, "dicty.b70u.fig.txt", quote=F, sep="\t", row.names = F)

dicty.b859u.fig <- ldply(dicty.b859u.exp, data.frame)
names(dicty.b859u.fig) <- c("term","log2FC") 
write.table(dicty.b859u.fig, "dicty.b859u.fig.txt", quote=F, sep="\t", row.names = F)

dicty.b11u.fig <- ldply(dicty.b11u.exp, data.frame)
names(dicty.b11u.fig) <- c("term","log2FC") 
write.table(dicty.b11u.fig, "dicty.b11u.fig.txt", quote=F, sep="\t", row.names = F)

dicty.b159d.fig <- ldply(dicty.b159d.exp, data.frame)
names(dicty.b159d.fig) <- c("term","log2FC") 
write.table(dicty.b159d.fig, "dicty.b159d.fig.txt", quote=F, sep="\t", row.names = F)

dicty.b70d.fig <- ldply(dicty.b70d.exp, data.frame)
names(dicty.b70d.fig) <- c("term","log2FC") 
write.table(dicty.b70d.fig, "dicty.b70d.fig.txt", quote=F, sep="\t", row.names = F)

dicty.b859d.fig <- ldply(dicty.b859d.exp, data.frame)
names(dicty.b859d.fig) <- c("term","log2FC") 
write.table(dicty.b859d.fig, "dicty.b859d.fig.txt", quote=F, sep="\t", row.names = F)

dicty.b11d.fig <- ldply(dicty.b11d.exp, data.frame)
names(dicty.b11d.fig) <- c("term","log2FC") 
write.table(dicty.b11d.fig, "dicty.b11d.fig.txt", quote=F, sep="\t", row.names = F)


# add go term descriptions in unix 
# sed 's/GO_/GO:/g' dicty.b159u.fig.txt > dicty.b159u.fig.colon.txt
# awk -F "\t" 'NR==FNR {a[$1]=$2; next} ($1 in a){print $1,$2,a[$1]}' OFS="\t" gobp_names.2025.txt dicty.b159u.fig.colon.txt | sed 's/GO:/GO_/g' > dicty.b159u.fig.name.txt
# 
# sed 's/GO_/GO:/g' dicty.b70u.fig.txt > dicty.b70u.fig.colon.txt
# awk -F "\t" 'NR==FNR {a[$1]=$2; next} ($1 in a){print $1,$2,a[$1]}' OFS="\t" gobp_names.2025.txt dicty.b70u.fig.colon.txt | sed 's/GO:/GO_/g' > dicty.b70u.fig.name.txt
# 
# sed 's/GO_/GO:/g' dicty.b859u.fig.txt > dicty.b859u.fig.colon.txt
# awk -F "\t" 'NR==FNR {a[$1]=$2; next} ($1 in a){print $1,$2,a[$1]}' OFS="\t" gobp_names.2025.txt dicty.b859u.fig.colon.txt | sed 's/GO:/GO_/g' > dicty.b859u.fig.name.txt
# 
# sed 's/GO_/GO:/g' dicty.b11u.fig.txt > dicty.b11u.fig.colon.txt
# awk -F "\t" 'NR==FNR {a[$1]=$2; next} ($1 in a){print $1,$2,a[$1]}' OFS="\t" gobp_names.2025.txt dicty.b11u.fig.colon.txt | sed 's/GO:/GO_/g' > dicty.b11u.fig.name.txt
# 
# sed 's/GO_/GO:/g' dicty.b159d.fig.txt > dicty.b159d.fig.colon.txt
# awk -F "\t" 'NR==FNR {a[$1]=$2; next} ($1 in a){print $1,$2,a[$1]}' OFS="\t" gobp_names.2025.txt dicty.b159d.fig.colon.txt | sed 's/GO:/GO_/g' > dicty.b159d.fig.name.txt
# 
# sed 's/GO_/GO:/g' dicty.b70d.fig.txt > dicty.b70d.fig.colon.txt
# awk -F "\t" 'NR==FNR {a[$1]=$2; next} ($1 in a){print $1,$2,a[$1]}' OFS="\t" gobp_names.2025.txt dicty.b70d.fig.colon.txt | sed 's/GO:/GO_/g' > dicty.b70d.fig.name.txt
# 
# sed 's/GO_/GO:/g' dicty.b859d.fig.txt > dicty.b859d.fig.colon.txt
# awk -F "\t" 'NR==FNR {a[$1]=$2; next} ($1 in a){print $1,$2,a[$1]}' OFS="\t" gobp_names.2025.txt dicty.b859d.fig.colon.txt | sed 's/GO:/GO_/g' > dicty.b859d.fig.name.txt
# 
# sed 's/GO_/GO:/g' dicty.b11d.fig.txt > dicty.b11d.fig.colon.txt
# awk -F "\t" 'NR==FNR {a[$1]=$2; next} ($1 in a){print $1,$2,a[$1]}' OFS="\t" gobp_names.2025.txt dicty.b11d.fig.colon.txt | sed 's/GO:/GO_/g' > dicty.b11d.fig.name.txt


# visualize, but only terms with mean logFC >1 or <-1

dicty.b159u.fig <- read.table("dicty.b159u.fig.name.txt", h=F, sep="\t")
names(dicty.b159u.fig) <- c("term","log2FC","description") 
dicty.b159u.fig$term <- factor(dicty.b159u.fig$term)
temp <- aggregate(dicty.b159u.fig$log2FC, list(dicty.b159u.fig$term), mean)
temp[order(temp$x),] # order by increasing mean logFC
dicty.b159u.fig.levels <- temp[order(temp$x),1]
highest <- droplevels(factor(tail(dicty.b159u.fig.levels, 6), levels=dicty.b159u.fig.levels)) 

dicty.b159d.fig <- read.table("dicty.b159d.fig.name.txt", h=F, sep="\t")
names(dicty.b159d.fig) <- c("term","log2FC","description") 
dicty.b159d.fig$term <- factor(dicty.b159d.fig$term)
temp <- aggregate(dicty.b159d.fig$log2FC, list(dicty.b159d.fig$term), mean)
temp[order(temp$x),] # order by increasing mean logFC
dicty.b159d.fig.levels <- temp[order(temp$x),1]
lowest <- droplevels(factor(head(dicty.b159d.fig.levels, 1), levels=dicty.b159d.fig.levels))

dicty.b159.fig <- droplevels(rbind(subset(dicty.b159d.fig, term %in% lowest), subset(dicty.b159u.fig, term %in% highest)))
dicty.b159.fig.levels <- c(lowest, highest) # check that order is consistent
dicty.b159.fig$term <- factor(dicty.b159.fig$term, levels=dicty.b159.fig.levels)
dicty.b159.fig <- dicty.b159.fig[order(dicty.b159.fig$term),]
temp <- unique(dicty.b159.fig[,c(1,3)])
dicty.b159.label <- paste(gsub(pattern="_", replacement=":", temp$term), temp$description, sep="\n")

ggplot(dicty.b159.fig, aes(x=term, y=log2FC)) + geom_boxplot(fill="grey") + xlab("GO Biological Process") + manuscript_theme + coord_flip() + geom_point(shape=21, fill="white") + scale_x_discrete(labels=dicty.b159.label)
rm(temp, highest, lowest, dicty.b159u.fig, dicty.b159u.fig.levels, dicty.b159d.fig, dicty.b159d.fig.levels)


dicty.b70u.fig <- read.table("dicty.b70u.fig.name.txt", h=F, sep="\t")
names(dicty.b70u.fig) <- c("term","log2FC","description") 
dicty.b70u.fig$term <- factor(dicty.b70u.fig$term)
temp <- aggregate(dicty.b70u.fig$log2FC, list(dicty.b70u.fig$term), mean)
temp[order(temp$x),] # order by increasing mean logFC
dicty.b70u.fig.levels <- temp[order(temp$x),1]
highest <- droplevels(factor(tail(dicty.b70u.fig.levels, 4), levels=dicty.b70u.fig.levels))

dicty.b70d.fig <- read.table("dicty.b70d.fig.name.txt", h=F, sep="\t")
names(dicty.b70d.fig) <- c("term","log2FC","description") 
dicty.b70d.fig$term <- factor(dicty.b70d.fig$term)
temp <- aggregate(dicty.b70d.fig$log2FC, list(dicty.b70d.fig$term), mean)
temp[order(temp$x),] # order by increasing mean logFC
dicty.b70d.fig.levels <- temp[order(temp$x),1]
lowest <- droplevels(factor(head(dicty.b70d.fig.levels, 6), levels=dicty.b70d.fig.levels))

dicty.b70.fig <- droplevels(rbind(subset(dicty.b70d.fig, term %in% lowest), subset(dicty.b70u.fig, term %in% highest)))
dicty.b70.fig.levels <- c(lowest, highest) # check that order is consistent
dicty.b70.fig$term <- factor(dicty.b70.fig$term, levels=dicty.b70.fig.levels)
dicty.b70.fig <- dicty.b70.fig[order(dicty.b70.fig$term),]
temp <- unique(dicty.b70.fig[,c(1,3)])
dicty.b70.label <- paste(gsub(pattern="_", replacement=":", temp$term), temp$description, sep="\n")

ggplot(dicty.b70.fig, aes(x=term, y=log2FC)) + geom_boxplot(fill="grey") + xlab("GO Biological Process") + manuscript_theme + coord_flip() + geom_point(shape=21, fill="white") + scale_x_discrete(labels=dicty.b70.label) 
rm(temp, highest, lowest, dicty.b70u.fig, dicty.b70u.fig.levels, dicty.b70d.fig, dicty.b70d.fig.levels)


dicty.b859u.fig <- read.table("dicty.b859u.fig.name.txt", h=F, sep="\t")
names(dicty.b859u.fig) <- c("term","log2FC","description") 
dicty.b859u.fig$term <- factor(dicty.b859u.fig$term)
temp <- aggregate(dicty.b859u.fig$log2FC, list(dicty.b859u.fig$term), mean)
temp[order(temp$x),] # order by increasing mean logFC
dicty.b859u.fig.levels <- temp[order(temp$x),1]
highest <- droplevels(factor(tail(dicty.b859u.fig.levels, 3), levels=dicty.b859u.fig.levels))

dicty.b859d.fig <- read.table("dicty.b859d.fig.name.txt", h=F, sep="\t")
names(dicty.b859d.fig) <- c("term","log2FC","description") 
dicty.b859d.fig$term <- factor(dicty.b859d.fig$term)
temp <- aggregate(dicty.b859d.fig$log2FC, list(dicty.b859d.fig$term), mean)
temp[order(temp$x),] # order by increasing mean logFC
dicty.b859d.fig.levels <- temp[order(temp$x),1]
lowest <- droplevels(factor(head(dicty.b859d.fig.levels, 10), levels=dicty.b859d.fig.levels)) #[c(-2,-3,-4)]

dicty.b859.fig <- droplevels(rbind(subset(dicty.b859d.fig, term %in% lowest), subset(dicty.b859u.fig, term %in% highest)))
dicty.b859.fig.levels <- c(lowest, highest) # check that order is consistent
dicty.b859.fig$term <- factor(dicty.b859.fig$term, levels=dicty.b859.fig.levels)
dicty.b859.fig <- dicty.b859.fig[order(dicty.b859.fig$term),]
temp <- unique(dicty.b859.fig[,c(1,3)])
dicty.b859.label <- paste(gsub(pattern="_", replacement=":", temp$term), temp$description, sep="\n")

ggplot(dicty.b859.fig, aes(x=term, y=log2FC)) + geom_boxplot(fill="grey") + xlab("GO Biological Process") + manuscript_theme + coord_flip() + geom_point(shape=21, fill="white") + scale_x_discrete(labels=dicty.b859.label) 
rm(temp, highest, lowest, dicty.b859u.fig, dicty.b859u.fig.levels, dicty.b859d.fig, dicty.b859d.fig.levels)


dicty.b11u.fig <- read.table("dicty.b11u.fig.name.txt", h=F, sep="\t")
names(dicty.b11u.fig) <- c("term","log2FC","description") 
dicty.b11u.fig$term <- factor(dicty.b11u.fig$term)
temp <- aggregate(dicty.b11u.fig$log2FC, list(dicty.b11u.fig$term), mean)
temp[order(temp$x),] # order by increasing mean logFC
dicty.b11u.fig.levels <- temp[order(temp$x),1]
highest <- droplevels(factor(tail(dicty.b11u.fig.levels, 3), levels=dicty.b11u.fig.levels)) #[c(-1,-2)]

dicty.b11d.fig <- read.table("dicty.b11d.fig.name.txt", h=F, sep="\t")
names(dicty.b11d.fig) <- c("term","log2FC","description") 
dicty.b11d.fig$term <- factor(dicty.b11d.fig$term)
temp <- aggregate(dicty.b11d.fig$log2FC, list(dicty.b11d.fig$term), mean)
temp[order(temp$x),] # order by increasing mean logFC
dicty.b11d.fig.levels <- temp[order(temp$x),1]
lowest <- droplevels(factor(head(dicty.b11d.fig.levels, 5), levels=dicty.b11d.fig.levels)) #[c(-1,-4)]

dicty.b11.fig <- droplevels(rbind(subset(dicty.b11d.fig, term %in% lowest), subset(dicty.b11u.fig, term %in% highest)))
dicty.b11.fig.levels <- c(lowest, highest) # check that order is consistent
dicty.b11.fig$term <- factor(dicty.b11.fig$term, levels=dicty.b11.fig.levels)
dicty.b11.fig <- dicty.b11.fig[order(dicty.b11.fig$term),]
temp <- unique(dicty.b11.fig[,c(1,3)])
dicty.b11.label <- paste(gsub(pattern="_", replacement=":", temp$term), temp$description, sep="\n")

ggplot(dicty.b11.fig, aes(x=term, y=log2FC)) + geom_boxplot(fill="grey") + xlab("GO Biological Process") + manuscript_theme + coord_flip() + geom_point(shape=21, fill="white") + scale_x_discrete(labels=dicty.b11.label)
rm(temp, highest, lowest, dicty.b11u.fig, dicty.b11u.fig.levels, dicty.b11d.fig, dicty.b11d.fig.levels)



manuscript_theme = ggplot2::theme_bw() + ggplot2::theme(axis.line.x = element_line(color="black", linewidth = .7),axis.line.y = element_line(color="black", size = .7),text=element_text(size=16))

# 7
setEPS(width=10, height=8) # default is 7x7
postscript(file=paste("dicty.b159.GO",format(Sys.time(),"%Y%m%d"),"eps",sep="."))
ggplot(dicty.b159.fig, aes(x=term, y=log2FC)) + geom_boxplot(fill="#009999") + xlab("GO Biological Process") + manuscript_theme + coord_flip() + geom_point(shape=21, fill="white", position=position_jitter(0.1)) + scale_x_discrete(labels=dicty.b159.label) 
dev.off()

# 9
setEPS(width=10, height=10) # default is 7x7
postscript(file=paste("dicty.b70.GO",format(Sys.time(),"%Y%m%d"),"eps",sep="."))
ggplot(dicty.b70.fig, aes(x=term, y=log2FC)) + geom_boxplot(fill="#00CC00") + xlab("GO Biological Process") + manuscript_theme + coord_flip() + geom_point(shape=21, fill="white", position=position_jitter(0.1)) + scale_x_discrete(labels=dicty.b70.label) 
dev.off()

# 13
setEPS(width=10, height=14) # default is 7x7
postscript(file=paste("dicty.b859.GO",format(Sys.time(),"%Y%m%d"),"eps",sep="."))
ggplot(dicty.b859.fig, aes(x=term, y=log2FC)) + geom_boxplot(fill="#FFa200") + xlab("GO Biological Process") + manuscript_theme + coord_flip() + geom_point(shape=21, fill="white", position=position_jitter(0.1)) + scale_x_discrete(labels=dicty.b859.label) 
dev.off()

# 8
setEPS(width=10, height=9) # default is 7x7
postscript(file=paste("dicty.b11.GO",format(Sys.time(),"%Y%m%d"),"eps",sep="."))
ggplot(dicty.b11.fig, aes(x=term, y=log2FC)) + geom_boxplot(fill="#ff0000") + xlab("GO Biological Process") + manuscript_theme + coord_flip() + geom_point(shape=21, fill="white", position=position_jitter(0.1)) + scale_x_discrete(labels=dicty.b11.label)
dev.off()


