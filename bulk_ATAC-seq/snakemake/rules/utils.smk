import os

input_type = config["general"]["input_type"].lower()
assert input_type in ["fastq", "bam"]
has_replicates = bool(config["general"]["has_replicates"])
# "Nextera": CTGTCTCTTATACACATCT+AGATGTGTATAAGAGACAG
bwamem2idx_ext = [".fa.gz.0123", ".fa.gz.amb", ".fa.gz.ann", ".fa.gz.bwt.2bit.64", ".fa.gz.pac"]    
bowtie2idx_ext = [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]

rule _sort_bam:
    """
    Sorts a bam file by genome position coordinates.
    Order: chrom --> start --> end
    """
    input: 
        bam = "{some}.bam"
    output: 
        bam = "{some}.sorted.bam"
    threads: workflow.cores * 0.2
    resources: tmpdir=tmpdir
    params:
        some = "{some}",
        log = f"{logdir}/{proj}.log"
    shell: 
        """
            samtools sort -@ {threads} -l 1 -o {output.bam} {input.bam}
            echo "- rule _sort_bam {params.some}" >> {params.log}
        """

        
rule _index_bam: 
    """
    Indexes a bam file.
    Note: only a coordinate-sorted bam can be indexed
    """
    input: 
        bam = "{someindex}.bam"
    output:
        bam = "{someindex}.bam.bai"
    threads: workflow.cores * 0.2
    resources: tmpdir=tmpdir
    params:
        someindex = "{someindex}",
        log = f"{logdir}/{proj}.log"
    shell: 
        """
            samtools index -@ {threads} {input.bam}
            echo "- rule _index_bam {params.someindex}" >> {params.log}
        """

rule _bam_merge:
    input: expand(f"{projdir}/04_alignment/{prefix_proj}" + "{rep}.sorted.bam", rep = reps_proj)
    output: bam = f"{projdir}/05_bam_merge/{prefix_proj}.bam"
    threads: workflow.cores
    resources: tmpdir=tmpdir
    params: 
        outdir = f"{projdir}/05_bam_merge",
        has_replicates = bool(config["general"]["has_replicates"]),
        log = f"{logdir}/{proj}.log"
    shell:
        """
        inbams=$(ls -d {projdir}/04_alignment/* | grep "bam$")
        mkdir -p {params.outdir}
        if [ {params.has_replicates} == True ];
        then
            samtools merge -o {output.bam} {input} -@ {threads}
        else
            ln -sr {output.bam} $inbams
        fi
        """

rule _bamCoverage: 
    """
    Calculate the read coverage (.bw) for bam.
    """
    input: 
        bam = f"{projdir}/05_bam_merge/{prefix_proj}.sorted.dedup.bam", # TODO replace 04_alignment by something more general
        bai = f"{projdir}/05_bam_merge/{prefix_proj}.sorted.dedup.bam.bai"
    output:
        bw = f"{projdir}/06_bamCoverage/{prefix_proj}.sorted.dedup.bam" + ".{}.bs{}.bw".format(config["bamCoverage"]["normalize"], config["bamCoverage"]["binSize"])
    threads: workflow.cores
    resources: tmpdir=tmpdir
    params:
        normalize = config["bamCoverage"]["normalize"],
        binSize = config["bamCoverage"]["binSize"],
        outdir = f"{projdir}/06_bamCoverage/",
        log = f"{logdir}/{proj}.log"
    shell: 
        """
            mkdir -p {params.outdir}
            bamCoverage --bam {input.bam} -o {output.bw} -bs {params.binSize} --extendReads -p max --normalizeUsing {params.normalize} -v
            echo "- rule _bamCoverage for {input.bam} (-bs {params.binSize}, --normalizeUsing {params.normalize})" >> {params.log}
        """










