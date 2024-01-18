### referred to avocato ###

# if aligner == "bwa":
#     ext =  [".fa.gz.amb", ".fa.gz.ann", ".fa.gz.bwt", ".fa.gz.pac", ".fa.gz.sa"]
# else:
bwamem2idx_ext = [".fa.gz.0123", ".fa.gz.amb", ".fa.gz.ann", ".fa.gz.bwt.2bit.64", ".fa.gz.pac"]    
bowtie2idx_ext = [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]

"""
Downloads the reference genome of the corresponding build."""
rule _download_reference:
    output: 
        fagz = config["general"]["data_download_to"] + f"/{build}/{build}.fa.gz",
        bwamem2_ls_fagz = config["general"]["data_download_to"] + f"/{build}/bwamem2_idx/{build}.fa.gz",
        blacklist = config["general"]["data_download_to"] + f"/{build}/{build}.blacklist.bed"
    resources: tmpdir = tmpdir
    params:
        path_download = os.path.join(config["general"]["data_download_to"], build),
        path_download_bwamem2 = config["general"]["data_download_to"] + f"/{build}/bwamem2_idx",
        tmp_blacklist = config["general"]["data_download_to"] + f"/{build}/{build}.blacklist.bed.gz",
        log = f"{logdir}/{proj}.log"
    shell:
        """
            build_lower=$(echo {build} | tr '[:upper:]' '[:lower:]')
            if [ $build_lower == "hg19"  ] || [ $build_lower == "grch37"  ];
            then
                download_pref="hg19"
                url_blacklist_proj="ENCFF001TDO"
            elif [ $build_lower == "hg38"  ] || [ $build_lower == "grch38"  ];
            then
                download_pref="hg38"
                url_blacklist_proj="ENCFF356LFX"
            fi
            mkdir -p {params.path_download}
            mkdir -p {params.path_download_bwamem2}
            if [ ! -f {output.fagz} ];
            then
                wget https://hgdownload.cse.ucsc.edu/goldenpath/"$download_pref"/bigZips/"$download_pref".fa.gz -O {output.fagz}
            fi
            ln -sr {output.fagz} {output.bwamem2_ls_fagz}
            wget -O {params.tmp_blacklist} https://www.encodeproject.org/files/"$url_blacklist_proj"/@@download/"$url_blacklist_proj".bed.gz
            zcat {params.tmp_blacklist} > {output.blacklist}
            rm {params.tmp_blacklist}
            echo "- rule _download_reference {output}" >> {params.log}
        """
# wget {public_path}/"$download_pref"/"$download_pref".fa -O {params.path_to_download}/"$download_pref".fa

# """
# Downloads the ref gene annotations for 
# a given genome build. Used for the MotifSearch.
# Currently (as is true for the whole pipeline)
# only support hg19 or hg38."""
# rule _download_genes:
#     output: 
#         genes = f"runs/common/ref_gene_{build}.txt"
#     params: 
#         genes_tmp = f"runs/common/ref_gene_{build}.txt.gz",
#         log = "runs/"+config["general"]["analysis_name"]+"/log/avocato.log"
#     resources:
#         #tmpdir = config["general"]["tmpdir"]
#         tmpdir = "runs/"+config["general"]["analysis_name"]+"/tmp"
#     run:
#         os.system(f"wget {public_path}/{build}/ref_gene_{build}.txt -O {output.genes}")
#         file = open(params.log, "a")
#         file.write("- rule _download_genes\n")
#         file.close

# """
# Downloads the chrom sizes for a given genome build.
# Used in generating coverage track from a fragments file.
# """
# rule _download_chrom_sizes:
#     output: 
#         chrom_sizes = f"runs/common/{build}.chrom.sizes"
#     resources:
#         #tmpdir = config["general"]["tmpdir"]
#         tmpdir = "runs/"+config["general"]["analysis_name"]+"/tmp"
#     params:
#         log = "runs/"+config["general"]["analysis_name"]+"/log/avocato.log"
#     run:
#         os.system(f"wget {public_path}/{build}/{build}.chrom.sizes -O {output.chrom_sizes}")
#         file = open(params.log, "a")
#         file.write("- rule _download_chrom_sizes\n")
#         file.close

rule _download_bwamem2idx_all:
    input:
        expand(config["general"]["data_download_to"] + f"/{build}/bwamem2_idx/{build}" + "{postfix}",  postfix=bwamem2idx_ext)
    resources: tmpdir=tmpdir
    params: log = f"{logdir}/{proj}.log"

rule _download_bwamem2idx_single:
    output: config["general"]["data_download_to"] + f"/{build}/bwamem2_idx/{build}" + "{postfix}"
    threads: 1
    resources: tmpdir=tmpdir
    params:
        pf = "{postfix}",
        path_download = config["general"]["data_download_to"] + f"/{build}/bwamem2_idx",
        log = f"{logdir}/{proj}.log"
    shell:
        """
            build_lower=$(echo {build} | tr '[:upper:]' '[:lower:]')
            if [ $build_lower == "hg19"  ] || [ $build_lower == "grch37"  ];
            then
                download_pref="hg19"
            elif [ $build_lower == "hg38"  ] || [ $build_lower == "grch38"  ];
            then
                download_pref="hg38"
            fi
            mkdir -p {params.path_download}
            if [ ! -f {output} ];
            then
                wget {public_path}/"$download_pref"/bwamem2/"$download_pref"{params.pf} -O {output}
            fi
            echo "- rule _download_index_bwamem2_single {output}" >> {params.log}
        """
        

rule _download_bowtie2idx_all:
    input:
        expand(config["general"]["data_download_to"] + f"/{build}/bowtie2_idx/{build}" + "{postfix}",  postfix=bowtie2idx_ext)
    resources: tmpdir=tmpdir
    params: log = f"{logdir}/{proj}.log"


rule _download_bowtie2idx_single:
    output: config["general"]["data_download_to"] + f"/{build}/bowtie2_idx/{build}" + "{postfix}"
    threads: 1
    resources: tmpdir=tmpdir
    params:
        pf = "{postfix}",
        path_download = config["general"]["data_download_to"] + f"/{build}/bowtie2_idx",
        log = f"{logdir}/{proj}.log"
    shell:
        """
            build_lower=$(echo {build} | tr '[:upper:]' '[:lower:]')
            if [ $build_lower == "hg19"  ] || [ $build_lower == "grch37"  ];
            then
                download_pref="GRCh37"
            elif [ $build_lower == "hg38"  ] || [ $build_lower == "grch38"  ];
            then
                download_pref="grch38_1kgmaj"
            fi
            mkdir -p {params.path_download}
            if [ ! -f {output} ];
            then
                wget https://genome-idx.s3.amazonaws.com/bt/"$download_pref"{params.pf} -O {output}
            fi
            echo "- rule _download_index_bowtie_single {output}" >> {params.log}
        """

# https://genome-idx.s3.amazonaws.com/bt/GRCh37.rev.2.bt2