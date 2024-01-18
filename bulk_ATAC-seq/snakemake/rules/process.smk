import os

input_type = config["general"]["input_type"].lower()
assert input_type in ["fastq", "bam"]
has_replicates = bool(config["general"]["has_replicates"])
# "Nextera": CTGTCTCTTATACACATCT+AGATGTGTATAAGAGACAG
bwamem2idx_ext = [".fa.gz.0123", ".fa.gz.amb", ".fa.gz.ann", ".fa.gz.bwt.2bit.64", ".fa.gz.pac"]    
bowtie2idx_ext = [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]

if input_type.lower() == "fastq":
    reads_all, = glob_wildcards(os.path.join(config["general"]["input"], "{reads}.fastq.gz"))
    if not reads_all:
        reads_all, = glob_wildcards(os.path.join(config["general"]["input"], "{reads}.fq.gz")) 
    reads_all = sorted(reads_all) # make sure _1 lists before _2
    reads_all_after_trimming = [f"{x}_val_{i}" for i, x in zip([1,2] * int(len(reads_all)/2), reads_all)] # TODO more automatic
    postfix_fastqc = ["html", "zip"]

postfix = f".{input_type}.gz" if input_type == "fastq" else f".{input_type}"
if has_replicates:
    prefix_all, reps_all, ids = glob_wildcards(os.path.join(config["general"]["input"], "{prefix}_Rep{rep}_{id}" + postfix))
    reps_proj = [f"_Rep{rep}" for rep in set(reps_all)]
else:
    prefix_all, = glob_wildcards(os.path.join(config["general"]["input"], "{prefix}" + postfix))
    reps_proj = [""]

assert len(set(prefix_all)) == 1 # make sure the same prefix for all reps
prefix_proj = prefix_all[0]


rule _01_fastqc_before_trimming_all:
    input: expand(os.path.join(projdir, "01_fastqc/{reads}_fastqc.{postfix}"), reads=reads_all, postfix=postfix_fastqc)
    resources: tmpdir=tmpdir

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
        trim_reads1 = os.path.join(projdir, "02_adaptor_trimming/{reads}_1_val_1.fq.gz"),
        trim_reads2 = os.path.join(projdir, "02_adaptor_trimming/{reads}_2_val_2.fq.gz")
    threads: workflow.cores
    resources: tmpdir=tmpdir
    params:
        outdir_adaptor_trimming = projdir + "/02_adaptor_trimming",
        read_prefix = "{reads}",
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
        fastqc_after_trimming_out1 = rules._03_fastqc_after_trimming_single.output.fastqc_out1,
        fastqc_after_trimming_out2 = rules._03_fastqc_after_trimming_single.output.fastqc_out2,
        bowtie2idx = expand(config["general"]["data_download_to"] + f"/{build}/bowtie2_idx/{build}" + "{postfix}",  postfix=bowtie2idx_ext),
        bwamem2idx = expand(config["general"]["data_download_to"] + f"/{build}/bwamem2_idx/{build}" + "{postfix}",  postfix=bwamem2idx_ext),
        bwamem2_ref = config["general"]["data_download_to"] + f"/{build}/bwamem2_idx/{build}.fa.gz"
    output: f"{projdir}/04_alignment/" + "{reads}.bam",
    threads: workflow.cores
    resources: tmpdir=tmpdir
    params:
        build = build,
        outdir = f"{projdir}/04_alignment/",
        minq = 30,
        bam_tmp_bwamem2 = f"{tmpdir}/tmp." + "{reads}.bwamem2.bam",
        bowtie2_indeces = config["general"]["data_download_to"] + f"/{build}/bowtie2_idx/{build}",
        bowtie2_xsize = config["alignment"]["bowtie2_xsize"],
        log = f"{logdir}/{proj}.log"
    wildcard_constraints:
        reads="^(sorted)"
    shell:
        """
        if [ {aligner} == "bowtie2" ];
        then
            echo "Align paired-end reads using bowtie2, and filter for properly mapped, aligned and high quality (MAPQ>={params.minq})..."
            bowtie2 --very-sensitive -X {params.bowtie2_xsize} -p {threads} -x {params.bowtie2_indeces} -1 {input.trim_reads1} -2 {input.trim_reads2} | samtools view -@ 8 -bh -f 2 -F 1804 -q {params.minq} -O BAM -o {output}
            echo "- rule _04_alignment ({aligner}) " >> {params.log}
        elif [ {aligner} == "bwamem2" ];
        then
            echo "Align paired-end reads using bwa-mem2..."
            bwa-mem2 mem -t {threads} {input.bwamem2_ref} {input.trim_reads1} {input.trim_reads2} | samtools view -@ {threads} -bh > {params.bam_tmp_bwamem2}
            echo "Filter properly mapped, aligned and high quality (MAPQ>=20)..."
            samtools view -b -f 2 -F 1804 -q 20 -o {output} {params.bam_tmp_bwamem2}
            rm {params.bam_tmp_bwamem2}
            echo "- rule _04_alignment ({aligner})" >> {params.log}
        else
            echo "Aligner not supported! Should be either bowtie2 or bwamem2"
            echo "- [FAILED] rule _04_alignment ({aligner}): aligner not supported." >> {params.log}
        fi
        """
# samtools view -1 - 

# merge and remove dup
rule _05_mark_remove_duplicates:
    input:
        bam = f"{projdir}/05_bam_merge/{prefix_proj}.sorted.bam"
    output:
        bam_markdup = f"{projdir}/05_bam_merge/{prefix_proj}.sorted.dedup.bam",
        markdup_metrics = f"{projdir}/05_bam_merge/{prefix_proj}.sorted.dedup.markduplicates.metrics.txt"
    threads: 1
    resources: tmpdir=tmpdir
    params:
        outdir = f"{projdir}/05_bam_merge/",
        log = f"{logdir}/{proj}.log"
    shell: 
        """
        mkdir -p {params.outdir}
        picard MarkDuplicates \
            INPUT={input.bam} \
            OUTPUT={output.bam_markdup} \
            ASSUME_SORTED=true \
            REMOVE_DUPLICATES=true \
            METRICS_FILE={output.markdup_metrics} \
            VALIDATION_STRINGENCY=LENIENT \
            TMP_DIR={tmpdir}
        """

rule _07_peak_calling_macs2:
    input:
        bam = f"{projdir}/05_bam_merge/{prefix_proj}.sorted.dedup.bam",
        bam_index = f"{projdir}/05_bam_merge/{prefix_proj}.sorted.dedup.bam.bai"
    output: f"{projdir}/07_peak_calling/.{prefix_proj}.macs.done"
    threads: workflow.cores
    resources: tmpdir=tmpdir
    params:
        outdir = f"{projdir}/07_peak_calling/",
        peak_caller = config["peak_calling"]["peak_caller"],
        q = config["peak_calling"]["macs_q_cutoff"],
        modes = "-nolambda --keep-dup all --call-summits", # -g hs
        blacklist = config["general"]["data_download_to"] + f"/{build}/{build}.blacklist.bed",
        log = f"{logdir}/{proj}.log"
    shell: 
        """
        out_raw={params.outdir}/{prefix_proj}.{params.peak_caller}.q{params.q}"_peaks.narrowPeak"
        out={params.outdir}/{prefix_proj}.{params.peak_caller}.q{params.q}".uniq.noblack.merge.bed"
        mkdir -p {params.outdir}
        if [ ! -f $out_raw ];
        then
            macs2 callpeak --format BAMPE {params.modes} -q {params.q} -t {input.bam} --outdir {params.outdir} --name {prefix_proj}.{params.peak_caller}.q{params.q}
        fi
        grep -E chr"[0-9]{{1,2}}\s|[XY]\s" $out_raw | cut -f 1,2,3 | uniq | bedtools intersect -a - -b {params.blacklist} -v | bedtools merge -i - > $out
        touch {output}
        """


rule _07_peak_calling_lanceotron:
    input: "{}/06_bamCoverage/{}.sorted.dedup.bam.{}.bs{}.bw".format(projdir, prefix_proj, "RPKM", config["bamCoverage"]["binSize"])
    output: f"{projdir}/07_peak_calling/.{prefix_proj}.lanceotron.done"
    threads: workflow.cores
    resources: tmpdir=tmpdir
    params:
        outdir = f"{projdir}/07_peak_calling",
        peak_caller = config["peak_calling"]["peak_caller"],
        lot_cutoff = config["peak_calling"]["lanceotron_callPeaks_cutoff"],
        lot_format = config["peak_calling"]["lanceotron_callPeaks_format"],
        lot_thr = config["peak_calling"]["lanceotron_callPeaks_threshold"],
        lot_window = config["peak_calling"]["lanceotron_callPeaks_window"],
        lot_skipheader = bool(config["peak_calling"]["lanceotron_callPeaks_skipheader"]),
        blacklist = config["general"]["data_download_to"] + f"/{build}/{build}.blacklist.bed",
        log = f"{logdir}/{proj}.log"
    shell: 
        """
        out_raw={params.outdir}/{prefix_proj}"_L-tron.bed"
        tmp_pos={params.outdir}/{prefix_proj}"_L-tron.tmp.pos"
        tmp_score={params.outdir}/{prefix_proj}"_L-tron.tmp.score"
        out={params.outdir}/{prefix_proj}"_L-tron.overall_peak_score.bed"
        out_postproc={params.outdir}/{prefix_proj}"_L-tron.overall_peak_score.noblack.bed"
        mkdir -p {params.outdir}
        mkdir -p {params.outdir}/lanceotron_processed

        # 1. peak calling by lanceotron. default cutoff 0.5
        echo "peak calling by lanceotron, cut off at {params.lot_cutoff}..."
        lot_command="lanceotron callPeaks {input} -c {params.lot_cutoff} --format {params.lot_format} -f {params.outdir} -t {params.lot_thr} -w {params.lot_window}"
        if [ ! -f $out_raw ];
        then
            if [ {params.lot_skipheader} == True ];
            then
                $lot_command --skipheader
            else
                $lot_command
            fi
        fi
        tail -n +2 $out_raw > "$out_raw".tmp
        mv "$out_raw".tmp $out_raw

        # 2. score each peak
        echo "score lanceotron peaks..."
        lanceotron scoreBed -b $out_raw -f {params.outdir} {input}

        # 3. bed -> bed6 (for plotting colors accord. peak score)
        cat $out_raw | grep "chr" | tail -n +2 | cut -f 1,2,3 | sed 's/$/\t./g' > $tmp_pos
        cat $out_raw | grep "chr" | tail -n +2 | cut -f 4 | sed 's/$/\t./g' > $tmp_score
        paste $tmp_pos $tmp_score > $out
        rm $tmp_pos $tmp_score

        # 4. postprocessing: rm blacklists & split by score
        echo "postprocess lanceotron peak calls..."
        # 4.1 select chromosomes & rm blacklists
        cat $out | grep "chr" | bedtools intersect -a - -b {params.blacklist} -v > $out_postproc
        # 4.2 split by overall_peak_score
        for score in 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 0.99 0.999 0.9999 0.99999 0.999999 0.9999999;
        do
            awk -v SCORE=$score '($5>=SCORE) {{print}}' $out_postproc > {params.outdir}/lanceotron_processed/"{prefix_proj}"_score_"$score".bed
        done
        touch {output}
        """

rule all_alignment:
    input: 
        bam = expand(f"{projdir}/04_alignment/{prefix_proj}" + "{rep}.sorted.bam", rep = reps_proj)
    resources: tmpdir=tmpdir

rule final_bam_coverage:
    input: "{}/06_bamCoverage/{}.sorted.dedup.bam.{}.bs{}.bw".format(projdir, prefix_proj, config["bamCoverage"]["normalize"], config["bamCoverage"]["binSize"])
    resources: tmpdir=tmpdir

rule final_peak_calling:
    input: f"{projdir}/07_peak_calling/.{prefix_proj}." + config["peak_calling"]["peak_caller"] + ".done"
