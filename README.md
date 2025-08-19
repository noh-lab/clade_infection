# Clade_infection
This repo contains the following files:

* `host_de_analysis_plus.*.R` - `R` code used to run statistical analyses using count files and generate figures
* `symbiont_de_analysis_plus.R` - `R` code used to run statistical analyses using count files and generate figures

In `genome_files/`:
* `*.fsa` - draft genome assembly for each _Paraburkholderia_ strain
* `*.fixed.gff` - gene annotation files for each _Paraburkholderia_ strain

In `supporting_files/`:
* `b*.gene_association.img_topgo.txt` - Processed GO annotation files, one for each _Paraburkholderia_ symbiont genome
* `b*.secretion.system.prokka.txt` - Table of symbiont genes with predicted function in bacterial secretion systems
* `b*.eukaryotic.pfam.prokka.txt` - Table of symbiont genes possessing eukaryotic protein domains
* `gene_association.dictybase.filter.gostats.txt` - Processed GO annotation file for _D. discoideum_ genome
* `phagocytosis_genes.tsv` - Table of _D. discoideum_ genes with known roles in phagocytosis, with columns in the following order: 
   * gene symbol
   * dictyBase gene id
   * decription of protein coded by gene
* `ss_struct_eff_geneviewer.txt` - Table of symbiont bacterial secretion system structural components and predicted effectors in order of genome coordinate

In `DE_results/`:
* `burk*.DEresults.txt` - Four `DESeq2` summary tables for each symbiont species or strain associated with hosts
* `dicty*.DEresults.txt` - Four `DESeq2` summary tables for hosts associated with each symbiont species or strain

In `GO_results/`:
* `burk*.*reg.bt2.topGO.txt` - Eight `topGO` summary tables for each symbiont species or strain associated with hosts; tables for up- and down-regulated groups of genes are separated
* `dicty*.*reg.topGO.txt` - Eight `topGO` summary tables for hosts associated with each symbiont species or strain; tables for up- and down-regulated groups of genes are separated
    
# Description of data processing and analysis
## I. Alignment to counts
### a. QC of raw sequencing reads
Files clipped `fastq_quality_filter`
```
for FILE in *R1*.fastq.gz; do
	INFILE="${FILE%.gz}"
	KEEP="${FILE%_S*}"
	echo "$FILE","$INFILE","$KEEP";
	#echo "Processing $FILE"
	#gunzip "$FILE"
	#cat "$FILE" | fastx_clipper -Q 33 -v -l 20 | fastq_quality_filter -Q 33 -v -q 30 -p 33 -o clipped/"$KEEP"_R1.clip.fastq
	#gzip "$FILE"
	#gzip clipped/"$KEEP"_R1.clip.fastq
done;
```
### b. Host
#### 1. Index reference `STAR`
```
STAR --runThreadN 12 --runMode genomeGenerate --genomeDir ./STAR_index --genomeFastaFiles dicty_chromosomal_num_nodup.fa --sjdbGTFfile dicty_chromosomal_num_fixed.gtf
```
#### 2. Align 
```
for FILE in /export/groups/snoh/shared/clade_infection/initial_data/qs*.fastq; do
	INFILE="${FILE/*data\/}"
	KEEP="${INFILE%.clip.fastq}"
	STRAIN="${KEEP%_*}"
	echo "$INFILE","$KEEP", "$STRAIN";
	echo "Processing $INFILE"
	STAR --runThreadN 12 --genomeDir /export/groups/snoh/shared/dicty_reference_20230116/STAR_index --readFilesIn "$FILE" --outFileNamePrefix dicty_1_star_aln/"$KEEP". --outSAMtype BAM Unsorted --outSAMattributes All
done;
```
#### 3. Sort and fix read groups with `samtools` and `picard` (repeat for each list)
```
list_1952_7="qs9_b11Aligned.out.bam qs9_b159Aligned.out.bam qs9_a1Aligned.out.bam"
list_1952_8="qs9_b859Aligned.out.bam qs9_b70Aligned.out.bam qs9_a2Aligned.out.bam"
list_1966_1="qs18_b11Aligned.out.bam qs18_b159Aligned.out.bam qs18_b1Aligned.out.bam"
list_1966_2="qs18_b859Aligned.out.bam qs18_b70Aligned.out.bam qs18_b2Aligned.out.bam"
list_1966_3="qs864_b11Aligned.out.bam qs864_b159Aligned.out.bam qs864_c1Aligned.out.bam"
list_1966_4="qs864_b859Aligned.out.bam qs864_b70Aligned.out.bam qs864_c2Aligned.out.bam"

for FILE in dicty_1_star_aln/${list_1966_4}; do
	INFILE="${FILE/*dicty_1_star_aln\/}"
	KEEP="${INFILE%Aligned.out.bam}"
	echo "$INFILE","$KEEP";
	echo "Processing $INFILE"
	samtools view -u -F 0x4 dicty_1_star_aln/"$INFILE" | samtools sort -O bam -T "$KEEP" - > dicty_1_star_aln/"$KEEP".sorted.bam
	picard AddOrReplaceReadGroups I=dicty_1_star_aln/"$KEEP".sorted.bam O=dicty_2_fix_rg/"$KEEP".star.rg.bam RGPU=1966_4 RGSM="$KEEP" RGLB="$KEEP" RGPL=illumina
done;
```

#### 4. Take primary alignment 
```
for FILE in dicty_2_fix_rg/q*.rg.bam; do
	INFILE="${FILE/*dicty_2_fix_rg\/}"
	KEEP="${INFILE%.star.rg.bam}"
	echo "$INFILE","$KEEP";
	echo "Processing $INFILE"
	samtools view -F 0x100 -b dicty_2_fix_rg/"$INFILE" > dicty_3_primary/"$KEEP".star.primary.bam
done;
```

#### 5. Convert primary alignment back to fastq (repeat for each list)
```
list_1952_7="qs9_b11.star.primary.bam qs9_b159.star.primary.bam qs9_a1.star.primary.bam"
list_1952_8="qs9_b859.star.primary.bam qs9_b70.star.primary.bam qs9_a2.star.primary.bam"
list_1966_1="qs18_b11.star.primary.bam qs18_b159.star.primary.bam qs18_b1.star.primary.bam"
list_1966_2="qs18_b859.star.primary.bam qs18_b70.star.primary.bam qs18_b2.star.primary.bam"
list_1966_3="qs864_b11.star.primary.bam qs864_b159.star.primary.bam qs864_c1.star.primary.bam"
list_1966_4="qs864_b859.star.primary.bam qs864_b70.star.primary.bam qs864_c2.star.primary.bam"

for FILE in dicty_3_primary/${list_1966_4}; do
	INFILE="${FILE/*dicty_3_primary\/}"
	KEEP="${INFILE%.star.primary.bam}"
	echo "$INFILE","$KEEP";
	echo "Processing $INFILE"
	bamToFastq -i dicty_3_primary/"$INFILE" -fq dicty_4_fq/"$KEEP".aln.fastq
done;
```

#### 6. Align to symbiont genome (repeat for each symbiont genome)
```
for FILE in dicty_4_fq/*b159*.fastq; do
	INFILE="${FILE/*dicty_4_fq\/}"
	KEEP="${INFILE%.aln.fastq}"
	echo "$INFILE","$KEEP";
	echo "Processing $INFILE"
	bowtie2 --no-unal -p 12 --sensitive -x /export/groups/snoh/shared/burk_reference_20230612/bt2_ba159 -U /export/groups/snoh/snoh/clade_infection/dicty_4_fq/"$INFILE" -S dicty_4_fq/"$KEEP".bt2.bam
done;
```

#### 7. Filter out symbiont-aligned reads and count 
(there are very few, usually less than 1000 reads)
```
for FILE in dicty_4_fq/*b859*.bt2.bam; do
	INFILE="${FILE/*dicty_4_fq\/}"
	KEEP="${INFILE%.bt2.bam}"
	echo "$INFILE","$KEEP";
	echo "Processing $INFILE"
	samtools view -F 0x100 -bq 20 "$FILE" | cut -f 1 |  sort | uniq > dicty_4_fq/"$KEEP".exclude
	picard FilterSamReads I=dicty_3_primary/"$KEEP".star.primary.bam O=dicty_5_star_count/"$KEEP".noburk.bam READ_LIST_FILE=dicty_4_fq/"$KEEP".exclude FILTER=excludeReadList
	samtools view dicty_5_star_count/"$KEEP".noburk.bam | htseq-count -s reverse -a 20 - /export/groups/snoh/shared/dicty_reference_20230116/dicty_chromosomal_num_fixed.gtf > dicty_5_star_count/"$KEEP".count
done;
```

### c. Symbiont
#### 1. Index reference `bowtie2` 
(repeat throughout for each genome)
```
bowtie2-build ../../shared/burk_reference_20230612/ba159.fsa ../../shared/burk_reference_20230612/bt2_ba159
```

#### 2. Align
```
for FILE in /export/groups/snoh/shared/clade_infection/initial_data/b859*.fastq; do
	INFILE="${FILE/*data\/}"
	KEEP="${INFILE%.clip.fastq}"
	#STRAIN="${KEEP%_*}"
	#echo "$INFILE","$KEEP"
	echo "Processing $INFILE"
	bowtie2 --no-unal -p 12 --sensitive -x /export/groups/snoh/shared/burk_reference_20230612/bt2_bb859 -U /export/groups/snoh/shared/clade_infection/initial_data/"$INFILE" -S burk_1_bt2_aln/"$KEEP".bt2.bam
done;
```
#### 3. Take primary alignment, sort and fix read groups
```
list_1952_7="b11_a.bt2.bam b11_qs9.bt2.bam b159_a.bt2.bam b159_qs9.bt2.bam"
list_1952_8="b70_a.bt2.bam b70_qs9.bt2.bam b859_a.bt2.bam b859_qs9.bt2.bam"
list_1966_1="b11_b.bt2.bam b11_qs18.bt2.bam b159_b.bt2.bam b159_qs18.bt2.bam"
list_1966_2="b70_b.bt2.bam b70_qs18.bt2.bam b859_b.bt2.bam b859_qs18.bt2.bam"
list_1966_3="b11_c.bt2.bam b11_qs864.bt2.bam b159_c.bt2.bam b159_qs864.bt2.bam"
list_1966_4="b70_c.bt2.bam b70_qs864.bt2.bam b859_c.bt2.bam b859_qs864.bt2.bam"

for FILE in burk_1_bt2_aln/${list_1966_4}; do
	INFILE="${FILE/*burk_1_bt2_aln\/}"
	KEEP="${INFILE%.bt2.bam}"
	#echo "$INFILE","$KEEP";
	echo "Processing $INFILE"
	samtools view -u -F 0x100 burk_1_bt2_aln/"$INFILE" | samtools sort -O bam -T "$KEEP" - > burk_2_bt2_rg/"$KEEP".sorted.bam
	picard AddOrReplaceReadGroups I=burk_2_bt2_rg/"$KEEP".sorted.bam O=burk_2_bt2_rg/"$KEEP".rg.bam RGPU=1966_4 RGSM="$KEEP" RGLB="$KEEP" RGPL=illumina
done;
```

#### 4. Convert to fastq and align to host genome
```
list_1952_7="b11_a.rg.bam b11_qs9.rg.bam b159_a.rg.bam b159_qs9.rg.bam"
list_1952_8="b70_a.rg.bam b70_qs9.rg.bam b859_a.rg.bam b859_qs9.rg.bam"
list_1966_1="b11_b.rg.bam b11_qs18.rg.bam b159_b.rg.bam b159_qs18.rg.bam"
list_1966_2="b70_b.rg.bam b70_qs18.rg.bam b859_b.rg.bam b859_qs18.rg.bam"
list_1966_3="b11_c.rg.bam b11_qs864.rg.bam b159_c.rg.bam b159_qs864.rg.bam"
list_1966_4="b70_c.rg.bam b70_qs864.rg.bam b859_c.rg.bam b859_qs864.rg.bam"

for FILE in burk_2_bt2_rg/${list_1966_4}; do
	INFILE="${FILE/*burk_2_bt2_rg\/}"
	KEEP="${INFILE%.rg.bam}"
	#echo "$INFILE","$KEEP";
	echo "Processing $INFILE"
	bamToFastq -i burk_2_bt2_rg/"$KEEP".rg.bam -fq burk_3_bt2_fq/"$KEEP".aln.fastq
done;

for FILE in burk_3_bt2_fq/b8*.fastq; do
	INFILE="${FILE/*burk_3_bt2_fq\/}"
	KEEP="${INFILE%.aln.fastq}"
	echo "$INFILE","$KEEP";
	echo "Processing $INFILE"
	STAR --runThreadN 12 --genomeDir /export/groups/snoh/shared/dicty_reference_20230116/STAR_index --readFilesIn "$FILE" --outFileNamePrefix burk_3_bt2_fq/"$KEEP". --outSAMtype BAM Unsorted
done;
```

#### 5. Filter out host-aligned reads
```
for FILE in burk_3_bt2_fq/b70*.Aligned.out.bam; do
	INFILE="${FILE/*burk_3_bt2_fq\/}"
	KEEP="${INFILE%.Aligned.out.bam}"
	echo "$INFILE","$KEEP";
	echo "Processing $INFILE"
	samtools view -F 0x100 -bq 20 "$FILE" | cut -f 1 |  sort | uniq > burk_3_bt2_fq/"$KEEP".exclude
	picard FilterSamReads I=burk_2_bt2_rg/"$KEEP".rg.bam O=burk_4_bt2_count/"$KEEP".nohost.bam READ_LIST_FILE=burk_3_bt2_fq/"$KEEP".exclude FILTER=excludeReadList
done;
```

#### 6. Count
```
for FILE in burk_4_bt2_count/b70*nohost.bam; do
	INFILE="${FILE/*burk_4_bt2_count\/}"
	KEEP="${INFILE%.nohost.bam}"
	echo "$INFILE","$KEEP";
	echo "Processing $INFILE"
	htseq-count -f bam -s reverse -a 20 -t gene -i ID burk_4_bt2_count/"$INFILE" /export/groups/snoh/shared/burk_reference_20230612/ba70.fixed.gff > burk_4_bt2_count/"$KEEP".bt2.count
done
```

## II. Statistical analysis
### a. Input files to `R`: 

#### 1. count files
```
# host
/export/groups/snoh/snoh/clade_infection/dicty_5_star_count/

# symbiont
/export/groups/snoh/snoh/clade_infection/burk_4_bt2_count/
```
#### 2. GO annotation files
```
# host
geneid2go.dictyBase.20230102.filter.topgo.txt

# symbiont
ba159.gene_association.img_topgo.txt
ba70.gene_association.img_topgo.txt
bb859.gene_association.img_topgo.txt
bh11.gene_association.img_topgo.txt
```

#### 3. Specific genes of interest
```
# host
phagocytosis_genes.tsv

# symbiont
b*secretion.system.prokka.txt
b*eukaryotic.pfam.prokka.txt
```

#### 4. Geneviewer file
```
ss_struct_eff_geneviewer.txt
```


### b. Analysis in `R`
#### 1. Differential expression analysis `DESeq2`
```
# host
host_de_analysis_plus.*.R

# symbiont
symbiont_de_analysis_plus.*.R
```
#### 2. Enrichment analysis `topGO` and Network analysis `igraph` 
Included in the R script above

#### 3. Gene neighborhood visualization `geneviewer`
Included in the R script above


### c. Results files
#### 1. Differential expression analysis
```
# host
dicty*.DEresults.txt

# symbiont
burk*.DEresults.txt
```

#### 2. Enrichment analysis
```
# host
dicty*.*reg.bt2.topGO.txt

# symbiont
burk*.*reg.bt2.topGO.txt
```
