# clade_infection
This repo contains

    alignment_to_count_summary.txt - Summary of pipeline to generate count tables from raw illumina reads
    dicty_de_analysis.R - R code used to run statistical analyses using count files and generate figures
    burk_de_analysis.R - R code used to run statistical analyses using count files and generate figures
    gene_association.dictybase.filter.gostats.txt - Processed GO annotation file for D. discoideum genome
    endo_phago_path.tsv - Table of D. discoideum genes with known roles in endocytosis and phagocytosis, with columns in the following order
        gene symbol
        dictyBase gene id
        decription of protein coded by gene
    clade_infection.de_*.txt - Three DESeq2 summary tables for the time series contrasts performed in the study
    clade_infection.GO_*.txt - Six GOstats summary tables for the time series contrasts performed in the study; tables for up- and down-regulated groups of genes are separated
