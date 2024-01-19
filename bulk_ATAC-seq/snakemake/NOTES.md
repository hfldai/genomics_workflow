1. Input files: fastq/ bam fragments
- name conventions:
    - FASTQ include replicates: xxx_RepX_1.fastq.gz, xxx_RepX_2.fastq.gz
    - FASTQ not include replicates: xxx_1.fastq.gz, xxx_2.fastq.gz
    - bam include replicates: xxx_RepX.bam ...
    - bam not include replicates: xxx.bam
2. Bulk & single-cell
3. peak calling: 
 - whether provided with replicates, and how to take into account.
4. add whether need some intermediate files & downstream
5. investigate ENCODE standard


2. fastqc
3. adapter trimming:
    -a. not know adapter seq: trim galore!
    -b. know adapter seq: cutadapt [40], AdapterRemoval v2 [41], Skewer [42], and trimmomatic [43] all requiring input of known adapter sequences
-[x] 4. alignment: bwamem2 OR bowtie2
-[x] 5. mark duplicates
5. peak calling
-[x] macs2
-[x] lanceotron
-[ ] HMMRATAC

<!-- 
scRNA-seq
GWAS
pheWAS
variant call
eQTL
TCGA
1000 Genome project
NCBI
Cosmic
ClinVar
UK Biobank
TCGA
dbSNP
OMIM
GTEx -->