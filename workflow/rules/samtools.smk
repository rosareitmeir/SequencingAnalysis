rule samtobam:
	input:
		"../results/sam/{sample}.sam"
	output:
		"../results/bam/{sample}.bam"
	conda:
		"../envs/env.yaml"
	threads:
		config["software"]["samtools"]["threads"]
	shell:
		"samtools view -@ {threads} -Sb {input} > {output}"

rule sort:
	input:
		"../results/bam/{sample}.bam"
	output:
		"../results/bam_sorted/{sample}_sorted.bam"
	conda:
		"../envs/env.yaml"
	threads:
		config["software"]["samtools"]["threads"]
	shell:
		"samtools sort -@ {threads} {input} -o {output}"

rule index:
	input:
		"../results/bam_sorted/{sample}_sorted.bam"
	output:
		"../results/bam_sorted/{sample}_sorted.bam.bai"
	conda:
		"../envs/env.yaml"
	threads:
		config["software"]["samtools"]["threads"]
	shell:
		"samtools index -@ {threads} -b {input}"

rule stats:
	input:
		"../results/bam_sorted/{sample}_sorted.bam",
		"../results/bam_sorted/{sample}_sorted.bam.bai"
	output:
		"../results/stats/{sample}.stats"
	conda:
		"../envs/env.yaml"
	shell:
		"samtools idxstats {input[0]} > {output}"