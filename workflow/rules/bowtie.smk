rule bowtie_index:
	input:
		reference=config["ref"],
	params:
		idx_prefix=config["idx_prefix"]
	output:
		expand("{idxprefix}.{index}.bt2", idxprefix=config["idx_prefix"], index=range(1,5)),
		expand("{idxprefix}.rev.{index}.bt2", idxprefix=config["idx_prefix"], index=range(1,3))
	conda:
		#path relative to rules folder
		"../envs/env.yaml"
	shell:
		"bowtie2-build -f {input} {params.idx_prefix}"

rule bowtie_map:
	input:	
		get_fastq_pair1,
		get_fastq_pair2,
		expand("{idxprefix}.{index}.bt2", idxprefix=config["idx_prefix"], index=range(1,5)),
		expand("{idxprefix}.rev.{index}.bt2", idxprefix=config["idx_prefix"], index=range(1,3))
	params:
		idx_prefix = expand("{idxprefix}", idxprefix=config["idx_prefix"]),
		bowtie2_opts = config["software"]["bowtie2"]["map_opts"]
	output:
		"../results/sam/{sample}.sam"
	conda:
		"../envs/env.yaml"
	threads: 
		config["software"]["bowtie2"]["threads"]
	shell:
		"bowtie2 -p {threads} {params.bowtie2_opts} -x {params.idx_prefix} -1 {input[0]} -2 {input[1]} -S {output}"