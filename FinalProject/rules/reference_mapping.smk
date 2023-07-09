import os.path

############################################ QC & Preprocessing ################################################################
def get_rawread (wildcards):
	if wildcards.index == "1":
		return( samples.loc[wildcards.sample][0])
	else:
		return (samples.loc[wildcards.sample][1])

rule fastqc_raw:
	input:
		get_rawread
	output:
		html="results/qc/RNA_seq/raw/{sample}_{index}_fastqc.html",
		zip="results/qc/RNA_seq/raw/{sample}_{index}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename

	log:
		"logs/fastqc/{sample}_{index}.log"

	wrapper:
		"v1.31.1/bio/fastqc"

rule cutadapt:
	input:
		[get_fastq_pair1, get_fastq_pair2]
	output:
		fastq1="results/trimmed/RNA_seq/{sample}_1.fastq",
		fastq2="results/trimmed/RNA_seq/{sample}_2.fastq",
		qc="results/qc/RNA_seq/trimmed/{sample}.qc.txt"

	params:
		adapters= config["software"]["cutadapt"]["adapters"],
		extra= config["software"]["cutadapt"]["extra"]
	log:
		"logs/cutadapt/{sample}.log"
	conda:
		"../envs/env.yaml"
	threads:
		config["software"]["cutadapt"]["threads"]
	wrapper:
		"v1.31.1/bio/cutadapt/pe"

rule fastqc_trimmed:
	input:
		"results/trimmed/RNA_seq/{sample}_{index}.fastq"
	output:
		html="results/qc/RNA_seq/trimmed/{sample}_{index}_fastqc.html",
		zip="results/qc/RNA_seq/trimmed/{sample}_{index}.fastqc.zip"  # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename

	log:
		"logs/fastqc_trimmed/{sample}_{index}.log"
	conda:
		"../envs/env.yaml"
	wrapper:
		"v1.31.1/bio/fastqc"

################################################### Reference Mapping & Gene Counting ###################################################################

## building STAR index
# for de nov assembly transfer gff to gtf
rule convert_to_gtf:
    input:
        "results/assembly/annotation/" + wgs_name + ".gff"
    output:
        "results/assembly/annotation/" + wgs_name + ".gtf"
    log:
        "logs/gtf_convertion/"+ wgs_name + ".log"

    conda:
        "../envs/env.yaml"
    shell:
        "gffread {input} -T -o {output} &>{log}"


rule build_index:
    input:
        ref=config["ref"] if os.path.exists(config["ref"]) else "results/assembly/pilon/"+ wgs_name + ".fasta",
        gtf=config["ref_anno"] if os.path.exists(config["ref_anno"]) else "results/assembly/annotation/" + wgs_name + ".gtf"
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
		return "results/trimmed/RNA_seq/"+ wildcards.sample + "_1.fastq"
	else:
		return samples.loc[wildcards.sample][0]

def get_fq2_for_star(wildcards):
	if config["use_trimmed"]:
		return  "results/trimmed/RNA_seq/"+ wildcards.sample + "_2.fastq"
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
        annotation=config["ref_anno"] if os.path.exists(config["ref_anno"]) else "results/assembly/annotation/" + wgs_name + ".gtf"
        # optional input
        # chr_names="",           # implicitly sets the -A flag
        # fasta="genome.fasta"      # implicitly sets the -G flag
    output:
        multiext("results/featurecounts/allsamples",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts")
    threads:
        config["software"]["feature_counts"]["threads"]
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="-O --fracOverlap 0.2 -p"
    log:
        "logs/featurecounts/featurecounts.log"
    wrapper:
        "0.72.0/bio/subread/featurecounts"


################### Mapping QC & MultiQC ###########################################

rule RNAseq_samtobam:
	input:
		"results/star/pe/{sample}/pe_aligned.sam"
	output:
		"results/star/pe/{sample}/pe_aligned.bam"
	log:
		"logs/samtools/samtobam/{sample}.log"
	conda:
		"../envs/env.yaml"
	threads:
		config["software"]["samtools"]["threads"]
	shell:
		"samtools view -@ {threads} -Sb {input} > {output} 2> {log}"

rule qualimap:
    input:
        # BAM aligned, splicing-aware, to reference genome
        bam="results/star/pe/{sample}/pe_aligned.bam",
        # GTF containing transcript, gene, and exon data
        gtf= config["ref_anno"] if os.path.exists(config["ref_anno"]) else "results/assembly/annotation/" + wgs_name + ".gtf"
    output:
        directory("results/qc/RNA_seq/mapping/{sample}") 
    log:
        "logs/qualimap/{sample}.log"
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    wrapper:
        "v2.2.0/bio/qualimap/rnaseq"


rule RNA_Seq_multiqc:
    input:
        [expand("results/qc/RNA_seq/raw/{sample}_{index}_fastqc.html", index=numbers, sample = list(samples.index)),
         expand("results/qc/RNA_seq/trimmed/{sample}_{index}_fastqc.html",sample=list(samples.index),index= numbers),
         expand("results/qc/RNA_seq/mapping/{sample}", sample=list(samples.index))] if config["use_trimmed"] else
        [expand("results/qc/RNA_seq/raw/{sample}_{index}_fastqc.html", index=numbers, sample = list(samples.index)), expand("results/qc/RNA_seq/mapping/{sample}", sample=list(samples.index))]
    output:
        "results/qc/RNA_seq/multiqc_report.html"
    log:
        "logs/multiqc/RNA_seq_multiqc.log"
    conda:
        "../envs/env.yaml"
    shell:
        "multiqc results/qc/RNA_seq  -o results/qc/RNA_seq &> {log}"