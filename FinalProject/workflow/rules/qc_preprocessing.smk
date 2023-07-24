
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



rule run_multiqc:
	input:
		[expand("results/qc/RNA_seq/raw/{sample}_{index}_fastqc.html", index=numbers, sample = list(samples.index)),
		expand("results/qc/RNA_seq/trimmed/{sample}_{index}_fastqc.html",sample=list(samples.index),index= numbers)] if config["use_trimmed"] else
		expand("results/qc/RNA_seq/raw/{sample}_{index}_fastqc.html", index=numbers, sample = list(samples.index))
	output:
		"results/qc/RNA_seq/multiqc_report.html"
	log:
		"logs/multiqc/RNA_seq_multiqc.log"
	conda:
		"../envs/env.yaml"
	shell:
		"multiqc results/qc/RNA_seq  -o results/qc/RNA_seq &> {log}"


