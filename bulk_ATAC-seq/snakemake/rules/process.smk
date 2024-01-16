import os

### TODO: get PREFIX automatically
PREFIX = "ENCLB674HYQ_Rep1"
# "Nextera": CTGTCTCTTATACACATCT+AGATGTGTATAAGAGACAG
bwamem2idx_ext = [".fa.gz.0123", ".fa.gz.amb", ".fa.gz.ann", ".fa.gz.bwt.2bit.64", ".fa.gz.pac"]    
bowtie2idx_ext = [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]

reads_all, = glob_wildcards(os.path.join(config["general"]["input"], "{reads}.fastq.gz"))
if not reads_all:
    reads_all, = glob_wildcards(os.path.join(config["general"]["input"], "{reads}.fq.gz"))
reads_all = sorted(reads_all) # make sure _1 lists before _2
reads_all_after_trimming = [f"{x}_val_{i+1}" for i, x in enumerate(reads_all)] # TODO more automatic
postfix_fastqc = ["html", "zip"]

# module load FASTQC
# rule _01_fastqc_before_trimming_all:
#     input: expand(os.path.join(projdir, "01_fastqc/{reads}_fastqc.{postfix}"), reads=reads_all, postfix=postfix_fastqc)
#     resources: tmpdir=tmpdir

rule _01_fastqc_before_trimming_single:
    '''Run fastqc before adaptor trimming.'''
    input: os.path.join(config["general"]["input"], "{reads}.fastq.gz")
    output: 
        fastqc_out1 = os.path.join(projdir, "01_fastqc_before_adaptor_trimming/{reads}_fastqc.zip"),
        fastqc_out2 = os.path.join(projdir, "01_fastqc_before_adaptor_trimming/{reads}_fastqc.html")
    threads: 1
    resources: tmpdir=tmpdir
    params:
        outdir_fastqc = f"{projdir}/01_fastqc_before_adaptor_trimming",
        log = f"{logdir}/{proj}.log"
    shell: 
        """
            fastqc -t {threads} {input} -o {params.outdir_fastqc}
            echo "- rule _01_fastqc_before_trimming_single {input}" >> {params.log}
        """

rule _02_adaptor_trimming:
    '''Run trim-galore to remove adaptor
    TODO: add other options
    '''
    input:
        before_trimming_fastqc_all = expand(os.path.join(projdir, "01_fastqc_before_adaptor_trimming/{reads}_fastqc.{postfix}"), reads=reads_all, postfix=postfix_fastqc)
    output:
        trim_reads1 = os.path.join(projdir, f"02_adaptor_trimming/{PREFIX}_1_val_1.fq.gz"),
        trim_reads2 = os.path.join(projdir, f"02_adaptor_trimming/{PREFIX}_2_val_2.fq.gz")
    threads: workflow.cores
    resources: tmpdir=tmpdir
    params:
        outdir_adaptor_trimming = projdir + "/02_adaptor_trimming",
        read_prefix = PREFIX,
        log = f"{logdir}/{proj}.log",
        indir = config["general"]["input"]
    shell:
        """
        mkdir -p {params.outdir_adaptor_trimming}
        trim_galore --paired -a \" TAAGGCGA -a GTAAGGAG -a CTGTCTCTTATACACATCT\" -j {threads} {params.indir}/{params.read_prefix}*.gz -o {params.outdir_adaptor_trimming}
        echo "- rule _02_adaptor_trimming {input}" >> {params.log}
        """

rule _03_fastqc_after_trimming_all:
    input: expand(os.path.join(projdir, "03_fastqc_after_adaptor_trimming/{reads}_fastqc.{postfix}"), reads=reads_all_after_trimming, postfix=postfix_fastqc)
    resources: tmpdir=tmpdir

rule _03_fastqc_after_trimming_single:
    '''Run fastqc after adaptor trimming.'''
    input: os.path.join(projdir, "02_adaptor_trimming/{reads}.fq.gz")
    output: 
        fastqc_out1 = os.path.join(projdir, "03_fastqc_after_adaptor_trimming/{reads}_fastqc.zip"),
        fastqc_out2 = os.path.join(projdir, "03_fastqc_after_adaptor_trimming/{reads}_fastqc.html")
    threads: 1
    resources: tmpdir=tmpdir
    params:
        outdir_fastqc = projdir + "/03_fastqc_after_adaptor_trimming",
        log = f"{logdir}/{proj}.log"
    shell: 
        """
            fastqc -t {threads} {input} -o {params.outdir_fastqc}
            echo "- rule _03_fastqc_after_trimming_single {input}" >> {params.log}
        """

rule _04_alignment:
    '''Align reads to reference genome, using bowtie2 or bwamem2'''
    input:
        trim_reads1 = rules._02_adaptor_trimming.output.trim_reads1,
        trim_reads2 = rules._02_adaptor_trimming.output.trim_reads2,
        bowtie2idx = expand(config["general"]["data_download_to"] + f"/{build}/bowtie2_idx/{build}" + "{postfix}",  postfix=bowtie2idx_ext),
        bwamem2idx = expand(config["general"]["data_download_to"] + f"/{build}/bwamem2_idx/{build}" + "{postfix}",  postfix=bwamem2idx_ext),
        bwamem2_ref = config["general"]["data_download_to"] + f"/{build}/bwamem2_idx/{build}.fa.gz"
    output: f"{projdir}/04_alignment/{PREFIX}.bam"
    threads: workflow.cores
    resources: tmpdir=tmpdir
    params:
        build = build,
        outdir = f"{projdir}/04_alignment/",
        prefix = PREFIX,
        minq = 30,
        bam_tmp_bwamem2 = f"{tmpdir}/tmp.{PREFIX}.bwamem2.bam",
        bowtie2_indeces = config["general"]["data_download_to"] + f"/{build}/bowtie2_idx/{build}",
        bowtie2_xsize = config["alignment"]["bowtie2_xsize"],
        log = f"{logdir}/{proj}.log"
    shell:
        """
        if [ {aligner} == "bowtie2"  ];
        then
            echo "Align paired-end reads using bowtie2, and filter for properly mapped, aligned and high quality (MAPQ>={params.minq})..."
            bowtie2 --very-sensitive -X {params.bowtie2_xsize} -p {threads} -x {params.bowtie2_indeces} -1 {input.trim_reads1} -2 {input.trim_reads2} | samtools view -@ 8 -bh -f 2 -F 1804 -q {params.minq} -O BAM -o {output}
        elif [ {aligner} == "bwamem2" ];
        then
            echo "Align paired-end reads using bwa-mem2..."
            bwa-mem2 mem -t {threads} {input.bwamem2_ref} {input.trim_reads1} {input.trim_reads2} | samtools view -@ {threads} -bh > {params.bam_tmp_bwamem2}
            echo "Filter properly mapped, aligned and high quality (MAPQ>=20)..."
            samtools view -b -f 2 -F 1804 -q 20 -o {output} {params.bam_tmp_bwamem2}
            rm {params.bam_tmp_bwamem2}
        fi
        """
# samtools view -1 - 

rule all_alignment:
    input: 
        bam_sorted = f"{projdir}/04_alignment/{PREFIX}.sorted.bam",
        bam_sorted_index = f"{projdir}/04_alignment/{PREFIX}.sorted.bam.bai"
    resources: tmpdir=tmpdir

rule final_bam_coverage:
    input: "{}/05_bamCoverage/{}.sorted.bam.{}.bs{}.bw".format(projdir, PREFIX, config["bamCoverage"]["normalize"], config["bamCoverage"]["binSize"])
    resources: tmpdir=tmpdir

# ### Marke and Remove duplicates ###
# module load picard samtools

# INBAM="$1"

# OUTDIR=$(dirname $INBAM)"/1_duplicate_removal/"
# OUTBAM=$OUTDIR/$(basename "$INBAM" | sed 's/.bam$//g').dedup.bam
# METRICS=$OUTDIR/$(basename "$INBAM" | sed 's/.bam$//g').MarkDuplicates.metrics.txt

# mkdir -p $OUTDIR

# java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
# INPUT=$INBAM \
# OUTPUT=$OUTBAM \
# ASSUME_SORTED=true \
# REMOVE_DUPLICATES=true \
# METRICS_FILE=$METRICS \
# VALIDATION_STRINGENCY=LENIENT \
# TMP_DIR=/well/gerton/liangtid/tmp

# samtools index $OUTBAM
