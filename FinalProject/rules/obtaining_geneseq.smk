
## extracting gene list

rule create_gene_list:
    input:
        "results/featurecounts/count.featureCounts"
    output:
        "results/gene_sequences/gene_list.tsv"
    # get rid of the header, saving gene name, "chr/contig", gene start and gene end
    shell:
        "tail -n +3 {input} | cut -f 1-4 -d $'\t'   > {output}"

## create bam from the STAR sam file
rule sam_to_bam:
    input:
        "results/star/pe/{sample}/pe_aligned.sam"

    output:
        "results/star/pe/{sample}/pe_aligned.bam"

    conda:
         "../envs/env.yaml"

    shell:
        "samtools view -b {input} > {output}"


# sort each bam file according to the genomic position of the reads
rule sort_bam:
    input:
        "results/star/pe/{sample}/pe_aligned.bam"
    output:
        "results/star/pe/{sample}/pe_aligned_sorted.bam"

    conda:
        "../envs/env.yaml"

    shell:
        "samtools sort {input} > {output}"

## adds index to the sorted bam files (positions of the reads in the Bam file)
rule index_bam:
    input:
         "results/star/pe/{sample}/pe_aligned_sorted.bam"
    output:
         "results/star/pe/{sample}/pe_aligned_sorted.bam.bai"
    conda:
        "../envs/env.yaml"
    shell:
        "samtools index {input}"


rule obtain_gene_sequences:
    input:
        bam= "results/star/pe/{sample}/pe_aligned_sorted.bam",
        bam_index= "results/star/pe/{sample}/pe_aligned_sorted.bam.bai",
        gene_list="results/gene_sequences/gene_list.tsv"

    output:
        out="results/gene_sequences/{sample}.tsv"

    params:
        min_coverage= config["software"]["OGS"]["min_coverage"]

    conda:
        "../envs/env.yaml"

    script:
        "../scripts/OGS.py"
