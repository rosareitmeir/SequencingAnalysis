
## Kraken2
rule kraken_classification:
	input:
		fq1= wgs_path + "_1.fastq.gz",
		fq2= wgs_path + "_2.fastq.gz"
	output:
		"results/kraken2/"+ wgs_name + "_output.txt",
		"results/kraken2/" + wgs_name + "_report.txt"
	params:
		threads= config["software"]["kraken2"]["threads"],
		db=config["software"]["kraken2"]["database"]
	log:
		"logs/kraken2/" + wgs_name + "_classification.log"
	conda:
		"../envs/env.yaml"
	shell:
		"kraken2 --use-names  --db {params.db}  --report {output[1]} --output {output[0]} --threads {params.threads}   --paired {input.fq1} {input.fq2}  &> {log}"

rule screen:
	input:
		"results/kraken2/"+ wgs_name + "_output.txt"
	output:
		"results/qc/kraken2/multiqc_report.html"
	log:
		"logs/multiqc/kraken2.log"
	conda:
		"../envs/env.yaml"
	shell:
		"multiqc results/kraken2  -o results/qc/kraken2 &> {log}"

# Decontamination with bowtie2 + samtools
# Using flag -f 12 on samtools to ensure all target reads, in this case read unmapped mate unmapped, with respect to the contamination reference, are kept.

rule cont_bowtie_index:
	input:
		reference=config["contaminations"],
	output:
		expand(cont_path + "/" + cont_name + ".{index}.bt2",  index=range(1,5)),
		expand(cont_path + "/" + cont_name + ".rev.{index}.bt2", index=range(1,3))
	log:
		"logs/bowtie2/cont_build_index/contaminations.log"
	conda:
		"../envs/env.yaml"
	params:
		out_dir=cont_path + "/" + cont_name
	shell:
		"bowtie2-build -f {input.reference} {params.out_dir} &> {log}"

rule cont_mapping:
	input:
		["results/trimmed/wgs/"+ wgs_name + "_1.fastq", 
		"results/trimmed/wgs/"+ wgs_name + "_2.fastq"] if config["WGS_trimmed"] else
		[wgs_path + "_1.fastq.gz",
		wgs_path + "_2.fastq.gz"],
		expand(cont_path + "/" + cont_name + ".{index}.bt2",  index=range(1,5)),
		expand(cont_path + "/" + cont_name + ".rev.{index}.bt2", index=range(1,3))
	output:
		"results/decontamination/mapping/sam/" + wgs_name + ".sam"
	log:
		"logs/bowtie2/cont_mapping/" + wgs_name + ".log"
	conda:
		"../envs/env.yaml"
	threads:
		config["software"]["bowtie2"]["threads"]
	params:
		ref_idx=cont_path + "/" + cont_name
	shell:
		"bowtie2 -p {threads} -x {params.ref_idx} -1 {input[0]} -2 {input[1]} -S {output} &> {log}"

# Extract target reads -> unmapped to contaminations, keep reads where none of the two reads in a pair map

rule cont_read_filter:
	input:
		"results/decontamination/mapping/sam/" + wgs_name + ".sam"
	output:
		"results/decontamination/mapping/sam/" + wgs_name + "_decon.sam"
	log:
		"logs/samtools/cont_filter/" + wgs_name + "_decon.log"
	conda:
		"../envs/env.yaml"
	threads:
		config["software"]["samtools"]["threads"]
	shell:
		"samtools view -f 12 {input} > {output} 2> {log}"

rule decon_sam_to_bam:
	input:
		"results/decontamination/mapping/sam/" + wgs_name + "_decon.sam"
	output:
		"results/decontamination/mapping/bam/" + wgs_name + "_decon.bam"
	log:
		"logs/samtools/decon_samtobam/" + wgs_name + "_decon.log"
	conda:
		"../envs/env.yaml"
	threads:
		config["software"]["samtools"]["threads"]
	shell:
		"samtools view -@ {threads} -Sb {input} > {output} 2> {log}"

rule decon_sort_bam:
	input:
		"results/decontamination/mapping/bam/" + wgs_name + "_decon.bam"
	output:
		"results/decontamination/mapping/bam/" + wgs_name + "_decon_sorted.bam"
	log:
		"logs/samtools/decon_sort/" + wgs_name + "_decon_sort.log"
	conda:
		"../envs/env.yaml"
	threads:
		config["software"]["samtools"]["threads"]
	shell:
		"samtools sort -@ {threads} {input} -o {output} &> {log}"

rule decon_bam_to_fastq:
	input:
		"results/decontamination/mapping/bam/" + wgs_name + "_decon_sorted.bam"
	output:
		"results/decontamination/" + wgs_name + "_1.fastq.gz",
		"results/decontamination/" + wgs_name + "_2.fastq.gz"
	log:
		"logs/samtools/sam_to_fastq/" + wgs_name + "_bamtofq.log"
	conda:
		"../envs/env.yaml"
	threads:
		config["software"]["samtools"]["threads"]
	shell:
		"samtools fastq -@ {threads} {input} -1 {output[0]} -2 {output[1]} -0 /dev/null -s /dev/null -n &> {log}"

# Quality control after decontamination

rule kraken_decon_classification:
	input:
		decon1="results/decontamination/" + wgs_name + "_1.fastq.gz",
		decon2="results/decontamination/" + wgs_name + "_2.fastq.gz"
	output:
		"results/kraken2/decon/" + wgs_name + "_decon_output.txt",
		"results/kraken2/decon/" + wgs_name + "_decon_report.txt"
	params:
		threads= config["software"]["kraken2"]["threads"],
		db=config["software"]["kraken2"]["database"]
	log:
		"logs/kraken2/decon/" + wgs_name + "_classification.log"
	conda:
		"../envs/env.yaml"
	shell:
		"kraken2 --use-names  --db {params.db}  --report {output[1]} --output {output[0]} --threads {params.threads}   --paired {input.decon1} {input.decon2}  &> {log}"


rule screen_after_decon:
	input:
		"results/kraken2/decon/" + wgs_name + "_decon_report.txt"
	output:
		"results/qc/kraken2/decon/multiqc_report.html"
	log:
		"logs/multiqc/multiqc_decon.log"
	conda:
		"../envs/env.yaml"
	shell:
		"multiqc results/kraken2  -o results/qc/kraken2/decon &> {log}"
