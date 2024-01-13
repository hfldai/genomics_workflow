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

