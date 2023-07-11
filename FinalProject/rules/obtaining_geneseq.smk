
## extracting gene list

rule create_gene_list:
    input:
        "results/featurecounts/allsamples.featureCounts"
    output:
        "results/gene_sequences/gene_list.tsv"
    # get rid of the header, saving gene name, "chr/contig", gene start and gene end
    shell:
        "tail -n +3 {input} | cut -f 1-4 -d $'\t'   > {output}"

# sort each bam file according to the genomic position of the reads
rule sort_bam:
    input:
        "results/hisat2/mapped/{sample}.bam"
    output:
        "results/hisat2/mapped/{sample}_sorted.bam"
    conda:
        "../envs/env.yaml"
    shell:
        "samtools sort {input} > {output}"

## adds index to the sorted bam files (positions of the reads in the Bam file)
rule index_bam:
    input:
         "results/hisat2/mapped/{sample}_sorted.bam"
    output:
         "results/hisat2/mapped/{sample}_sorted.bam.bai"
    conda:
        "../envs/env.yaml"
    shell:
        "samtools index {input}"

# Obtaining gene sequences 

rule samtools_faidx:
    input:
        #"{sample}.fa",
        genome=config["ref"] if os.path.exists(config["ref"]) else "results/assembly/pilon/"+ wgs_name + ".fasta"
    output:
        #"{sample}.fa.fai",
        ref_path + ".fai" if os.path.exists(config["ref"]) else "results/assembly/pilon/"+ wgs_name + ".fasta.fai" #TODO check proper naming
    log:
        "logs/samtools/reference_faidx.log",
    params:
        extra="",  # optional params string
    wrapper:
        "v2.2.0/bio/samtools/faidx"

# Obtain one region to run consensus rule
#regions = pd.read_csv("results/gene_sequences/gene_list.tsv", sep='\t')
#regions["Gene"] = regions.iloc[:,1].astype(str)+":"+regions.iloc[:,2].astype(str)+"-"+regions.iloc[:,3].astype(str)
#gene_regions = list(regions.iloc[:,1].astype(str)+":"+regions.iloc[:,2].astype(str)+"-"+regions.iloc[:,3].astype(str))

r#ule samtools_consensus:
 #   input:
 #       "results/hisat2/mapped/{sample}_sorted.bam",
 #       region = "{gene_region}"
 #   output:
 #      "results/consensus/{sample}/{sample}_{gene_region}.fasta",
 #   conda:
 #       "../envs/env.yaml"
 #   log:
 #       "logs/consensus/{sample}.log",
 #   params:
 #       #region="k141_147_pilon:286-486",
 #   threads: 8
 #   shell:
 #       "samtools consensus -r {input.region} -f fasta {input[0]} -d 10 --show-del yes -a --show-ins yes"
