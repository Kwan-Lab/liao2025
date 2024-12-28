# Psilocybin Cortex snRNAseq

This repository contains analysis scripts for the manuscript "**Single-nucleus transcriptomics reveals time-dependent and cell-type-specific effects of psilocybin on gene expression**" by Liao et al, 2025. 

# Usage

## Reads (fastq) processing

Raw fastq files can be downloaded from the SRA Accession. 
Per cellranger v7.1.0 docs, each sample should be in its own directory entitled with the sample name (ie. `4797-CL-1/`)
Files split into different sequencing lanes should be concatenated.
(Note that the 4 samples requiring concatenation were all excluded from the analysis due to oversequencing).
Processing should follow the **Methods** section of the manuscript.
The provided `Snakefile` can be used as a template but was built specifically for the authors' file system.

## Cell-gene matrix processing

Loading the matrices, quality control and preprocessing, transformations, UMAP embeddings, and clustering can be found in the `notebooks` folder. This also includes differential expression analysis using the memento package.

## Figure generation

Code to generate figures appearing in the paper and supplement can be found in the `paper_figures` folder.

# Correspondance

Please refer to manuscript for correspondance details