import os

rule _sort_bam:
    """
    Sorts a bam file by genome position coordinates.
    Order: chrom --> start --> end
    """
    input: 
        bam = "{somesort}.bam"
    output: 
        bam = "{somesort}.sorted.bam"
    threads: workflow.cores
    resources: tmpdir=tmpdir
    params:
        somesort = "{somesort}",
        log = f"{logdir}/{proj}.log"
    shell: 
        """
            samtools sort -@ {threads} -l 1 -o {output.bam} {input.bam}
            echo "- rule _sort_bam {params.somesort}" >> {params.log}
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
    threads: workflow.cores
    resources: tmpdir=tmpdir
    params:
        someindex = "{someindex}",
        log = f"{logdir}/{proj}.log"
    shell: 
        """
            samtools index -@ {threads} {input.bam}
            echo "- rule _index_bam {params.someindex}" >> {params.log}
        """

rule _bamCoverage: 
    """
    Calculate the read coverage (.bw) for bam.
    """
    input: 
        bam = "{some}.bam"
        bai = "{some}.bam.bai"
    output:
        bw = "{}/05_bamCoverage/".format(projdir) + os.path.basename("{some}.bam") + ".{}.bs{}.bw".format(config["bamCoverage"]["normalize"], config["bamCoverage"]["binSize"])
    threads: workflow.cores
    resources: tmpdir=tmpdir
    params:
        normalize = config["bamCoverage"]["normalize"],
        binSize = config["bamCoverage"]["binSize"],
        outdir = f"{projdir}/05_bamCoverage/",
        log = f"{logdir}/{proj}.log"
    shell: 
        """
            mkdir -p {params.outdir}
            bamCoverage --bam {input.bam} -o {output.bw} -bs {params.binSize} --extendReads -p max --normalizeUsing {params.normalize} -v
            echo "- rule _bamCoverage for {input.bam} (-bs {params.binSize}, --normalizeUsing {params.normalize})" >> {params.log}
        """


# ### Merge bam replicates ###
# BAMDIR="$1"
# COMMON="$2" # common pattern to for bams
# OUTBAM="$3"

# # RANDPREF=$(date +'%Y%m%dT%T' | sed 's/://g')
# DIR_OUTBAM=$(dirname $OUTBAM)
# mkdir -p $DIR_OUTBAM

# BAMS=$(ls -d $BAMDIR/*  | grep ".bam$" | grep $COMMON )
# samtools merge $OUTBAM $BAMS -@ 16
# samtools index $OUTBAM






