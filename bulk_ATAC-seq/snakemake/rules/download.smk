### referred to avocato ###

if aligner == "bwamem2":
    ext =  [".fa.gz.amb", ".fa.gz.ann", ".fa.gz.bwt", ".fa.gz.pac", ".fa.gz.sa"]
else:
    ext =  [".fa.gz.0123", ".fa.gz.amb", ".fa.gz.ann", ".fa.gz.bwt.2bit.64", ".fa.gz.pac"]    

"""
Downloads the reference genome of the corresponding build."""
# rule _download_reference:
#     output: 
#         fas = multiext(f"runs/common/{build}", ".fa", ".fa.gz")
#     resources: tmpdir = tmpdir
#     params: log = f"{logdir}/{proj}.log"
#     run:
#         os.system(f"wget {public_path}/{build}/{build}.fa -O runs/common/{build}.fa")
#         os.system(f"wget {public_path}/{build}/{build}.fa -O runs/common/{build}.fa.gz")
#         file = open(params.log, "a")
#         file.write("- rule _download_reference\n")
#         file.close

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
        
# rule _download_index_bwa:
#     output:
#         indeces = multiext(f"runs/common/{build}", ext[0], ext[1], ext[2], ext[3], ext[4])
#     resources:
#         #tmpdir = config["general"]["tmpdir"]
#         tmpdir = "runs/"+config["general"]["analysis_name"]+"/tmp"
#     params:
#         log = "runs/"+config["general"]["analysis_name"]+"/log/avocato.log"
#     run: 
#         os.system(f"wget {public_path}/{build}/{check_aligner}/{build}%s -O runs/common/{build}%s"%(ext[0], ext[0]))
#         os.system(f"wget {public_path}/{build}/{check_aligner}/{build}%s -O runs/common/{build}%s"%(ext[1], ext[1]))
#         os.system(f"wget {public_path}/{build}/{check_aligner}/{build}%s -O runs/common/{build}%s"%(ext[2], ext[2]))
#         os.system(f"wget {public_path}/{build}/{check_aligner}/{build}%s -O runs/common/{build}%s"%(ext[3], ext[3]))
#         os.system(f"wget {public_path}/{build}/{check_aligner}/{build}%s -O runs/common/{build}%s"%(ext[4], ext[4]))
#         file = open(params.log, "a")
#         file.write("- rule _download_index_bwa\n")
#         file.close
        

rule _download_bowtie2idx_all:
    input:
        expand(config["general"]["data_download_to"] + f"/{build}/bowtie2idx/{build}" + "{postfix}",  postfix=[".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"])
    resources: tmpdir=tmpdir
    params: log = f"{logdir}/{proj}.log"


rule _download_bowtie2idx_single:
    output: config["general"]["data_download_to"] + f"/{build}/bowtie2idx/{build}" + "{postfix}"
    threads: 1
    resources: tmpdir=tmpdir
    params:
        pf = "{postfix}",
        path_download = config["general"]["data_download_to"] + f"/{build}/bowtie2idx",
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
            wget https://genome-idx.s3.amazonaws.com/bt/"$download_pref"{params.pf} -O {output}
            echo "- rule _download_index_bowtie_single {output}" >> {params.log}
        """

# https://genome-idx.s3.amazonaws.com/bt/GRCh37.rev.2.bt2