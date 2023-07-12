
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




## Decontamination with Bowtie
rule bowtie_index:
    input:
        reference=get_contamination,
    output:
        expand("contaminations/{{cont}}.{index}.bt2",  index=range(1,5)),
        expand("contaminations/{{cont}}.rev.{index}.bt2", index=range(1,3))

    log:
        "logs/bowtie2/build_index/{cont}.log"
    conda:
        "../envs/env.yaml"
    shell:
        "bowtie2-build -f {input.reference} contaminations/{wildcards.cont} &> {log}"

