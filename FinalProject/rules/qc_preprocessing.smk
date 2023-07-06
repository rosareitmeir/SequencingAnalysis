
def get_rawread (wildcards):
	if wildcards.index == "1":
		return( samples.loc[wildcards.sample][0])
	else:
		return (samples.loc[wildcards.sample][1])

rule fastqc_raw:
	input:
		get_rawread
	output:
		html="results/qc/raw/{sample}_{index}_fastqc.html",
		zip="results/qc/raw/{sample}_{index}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename

	log:
		"logs/fastqc/{sample}_{index}.log"

	wrapper:
		"v1.31.1/bio/fastqc"

rule cutadapt:
	input:
		[get_fastq_pair1, get_fastq_pair2]
	output:
		fastq1="results/trimmed/{sample}_1.fastq",
		fastq2="results/trimmed/{sample}_2.fastq",
		qc="results/qc/trimmed/{sample}.qc.txt"

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
		"results/trimmed/{sample}_{index}.fastq"
	output:
		html="results/qc/trimmed/{sample}_{index}_fastqc.html",
		zip="results/qc/trimmed/{sample}_{index}.fastqc.zip"  # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename

	log:
		"logs/fastqc_trimmed/{sample}_{index}.log"
	conda:
		"../envs/env.yaml"
	wrapper:
		"v1.31.1/bio/fastqc"


rule run_multiqc:
	input:
		[expand("results/qc/raw/{sample}_{index}_fastqc.html", index=numbers, sample = list(samples.index)),
		expand("results/qc/trimmed/{sample}_{index}_fastqc.html",sample=list(samples.index),index= numbers)] if config["use_trimmed"] else
		expand("results/qc/raw/{sample}_{index}_fastqc.html", index=numbers, sample = list(samples.index))
	output:
		"results/qc/multiqc_report.html"
	log:
		"logs/multiqc/fasta_multiqc.log"
	conda:
		"../envs/env.yaml"
	shell:
		"multiqc results/qc  -o results/qc &> {log}"


