1. Input files: fastq/ bam fragments
2. Bulk & single-cell
3. add peak calling (and choice)
 - whether provided with replicates, and how to take into account.
4. add whether need some intermediate files & downstream
5. investigate ENCODE standard



conda create -n bulk_atac_e2e python=3.9    
conda install -c bioconda -c conda-forge aioeasywebdav eido snakemake fastqc cutadapt trim-galore
mamba install -c bioconda -c conda-forge aioeasywebdav eido snakemake fastqc cutadapt trim-galore


1. fastqc
2. adapter trimming:
    -a. not know adapter seq: trim galore!
    -b. know adapter seq: cutadapt [40], AdapterRemoval v2 [41], Skewer [42], and trimmomatic [43] all requiring input of known adapter sequences
3. alignment: bwamem2 OR bowtie2
4. mark duplicates
5. peak calling

<!-- 
scRNA-seq
GWAS
pheWAS
variant call
eQTL
TCGA
1000 Genome project
GTEx -->