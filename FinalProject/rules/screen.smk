
## Kraken2
rule kraken_classification:
    input:
        fq1= get_fastq_pair1,
        fq2= get_fastq_pair2
    output:
        "results/kraken2/{sample}_output.txt",
        "results/kraken2/{sample}_report.txt"

    params:
        threads= config["software"]["kraken2"]["threads"],
        db=config["software"]["kraken2"]["database"]

    log:
        "logs/kraken2/{sample}_classification.log"
    conda:
        "../envs/env.yaml"
    shell:
        "kraken2 --use-names  --db {params.db}  --report {output[1]} --output {output[0]} --threads {params.threads}   --paired {input.fq1} {input.fq2}  &> {log}"



rule screen:
    input:
        expand("results/kraken2/{sample}_report.txt", sample= list(samples.index))
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


def get_fqs_for_decontamination(wildcards):
    if config["use_trimmed_data"]:
            return ["results/trimmed/"+ wildcards.sample + "_1.fastq", "results/trimmed/"+ wildcards.sample + "_2.fastq"]
    else:
        return [ samples.loc[wildcards.sample][0], samples.loc[wildcards.sample][1]]


rule decontamination:
    input:
        fqs= get_fqs_for_decontamination,
        index=expand("contaminations/{cont}.{index}.bt2", cont=list(cont_table.index), index=range(1, 5)),
        revindex=expand("contaminations/{cont}.rev.{index}.bt2", cont= list(cont_table.index), index=range(1,3))

    output:
        decon1="results/decontamination/{sample}_1.fastq.gz",
        decon2="results/decontamination/{sample}_2.fastq.gz"
    params:
       bowtie_indeces = list(cont_table.index),
       threads=config["software"]["bowtie2"]["threads"]
    log:
        "logs/bowtie2/{sample}_decontamination.log"
    conda:
        "../envs/env.yaml"
    shell:
        """
        fq1="{input.fqs[0]}";
        fq2="{input.fqs[1]}";
        list="{params.bowtie_indeces}";
        for i in $list; do 
                temp="results/decontamination/{wildcards.sample}_${{i}}";
                bowtie2 -p {params.threads}  -x "contaminations/${{i}}" -1 $fq1 -2 $fq2 --un-conc-gz $temp &>> {log} ;
                fq1="${{temp}}.1";
                fq2="${{temp}}.2";
            done
        ## renaming finished decontaminated file from temp 
        mv $fq1 {output.decon1};
        mv $fq2 {output.decon2}; 
        ## remove temp files
        rm results/decontamination/{wildcards.sample}_*.1;
        rm results/decontamination/{wildcards.sample}_*.2;
         
                """

rule kraken_decon_classification:
    input:
        decon1="results/decontamination/{sample}_1.fastq.gz",
        decon2="results/decontamination/{sample}_2.fastq.gz"
    output:
        "results/kraken2/decon/{sample}_decon_output.txt",
        "results/kraken2/decon/{sample}_decon_report.txt"

    params:
        threads= config["software"]["kraken2"]["threads"],
        db=config["software"]["kraken2"]["database"]

    log:
        "logs/kraken2/decon/{sample}_classification.log"
    conda:
        "../envs/env.yaml"
    shell:
        "kraken2 --use-names  --db {params.db}  --report {output[1]} --output {output[0]} --threads {params.threads}   --paired {input.decon1} {input.decon2}  &> {log}"



rule screen_after_decon:
    input:
        expand("results/kraken2/decon/{sample}_decon_report.txt", sample= list(samples.index))
    output:
        "results/qc/kraken2/decon/multiqc_report.html"
    log:
        "logs/multiqc/multiqc_decon.log"
    conda:
        "../envs/env.yaml"
    shell:
        "multiqc results/kraken2  -o results/qc/kraken2/decon &> {log}"
