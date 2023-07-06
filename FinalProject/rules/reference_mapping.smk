
## building STAR index

rule build_index:
    input:
        ref=config["ref"],
        gtf=config["ref_anno"]
    output:
        directory("reference/genome_idx"),
        "reference/genome_idx/Genome",
        "reference/genome_idx/genomeParameters.txt",
        "reference/genome_idx/sjdbInfo.txt",
        "reference/genome_idx/sjdbList.fromGTF.out.tab",

    threads:
        config["software"]["star"]["indx_threads"]
    conda:
        "../envs/env.yaml"
    log:
        "logs/star/build_idx.log"
    shell:
        "STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output[0]} --genomeFastaFiles {input.ref} --sjdbGTFfile {input.gtf} &> {log}"


def get_fq1_for_star(wildcards):
	if config["use_trimmed"]:
		return "results/trimmed/"+ wildcards.sample + "_1.fastq"
	else:
		return samples.loc[wildcards.sample][0]

def get_fq2_for_star(wildcards):
	if config["use_trimmed"]:
		return  "results/trimmed/"+ wildcards.sample + "_2.fastq"
	else:
		return samples.loc[wildcards.sample][1]


rule star_pe_multi:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fq1=get_fq1_for_star,
        # paired end reads needs to be ordered so each item in the two lists match
        fq2=get_fq2_for_star,  #optional
        # path to STAR reference genome index
        idx="reference/genome_idx",
    output:
        # see STAR manual for additional output files
        aln="results/star/pe/{sample}/pe_aligned.sam",
        log="logs/pe/{sample}/Log.out",
        sj="results/star/pe/{sample}/SJ.out.tab",
    log:
        "logs/star/pe/{sample}.log",
    params:
        # optional parameters
        extra=config["software"]["star"]["align_extra"],

    threads: config["software"]["star"]["align_threads"]
    wrapper:
        "v2.0.0/bio/star/align"

rule feature_counts:
    input:
        sam=expand("results/star/pe/{sample}/pe_aligned.sam" , sample=list(samples.index)), # list of sam or bam files
        annotation=config["ref_anno"],
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    output:
        multiext("results/featurecounts/count",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts")
    threads:
        1
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="-O --fracOverlap 0.2 -p"
    log:
        "logs/featurecounts/featurecounts.log"
    wrapper:
        "0.72.0/bio/subread/featurecounts"