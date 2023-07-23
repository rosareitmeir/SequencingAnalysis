## extracting gene list
rule create_gene_bed:
    input:
        "results/featurecounts/allsamples.featureCounts"
    output:
        "results/consensus/gene_list.bed"
    # get rid of the header, saving gene name, "chr/contig", gene start and gene end
    log:
        "logs/gene_bed/create_gene_bed.log"
    params:
        expression= r"""{{print $2 "\t" $3 "\t" $4 "\t" $1}}"""
    shell:
        "tail -n +3 {input} | awk  '{params.expression}'  > {output}"


# Sort each bam file according to the genomic position of the reads

rule sort_bam:
    input:
        "results/hisat2/mapped/{sample}" + ".bam"
    output:
        "results/hisat2/mapped/{sample}"+ ".sorted.bam"
    conda:
        "../envs/env.yaml"
    log:
        "logs/samtools/{sample}_sorting.log"
    shell:
        "samtools sort {input} > {output} 2> {log}"


## adds index to the sorted bam files (positions of the reads in the Bam file)

rule index_bam:
    input:
         "results/hisat2/mapped/{sample}.sorted.bam"
    output:
         "results/hisat2/mapped/{sample}.sorted.bam.bai"
    conda:
        "../envs/env.yaml"
    log:
        "logs/samtools/{sample}_indexing.log"
    shell:
        "samtools index {input} 2>> {log}"


# Obtaining gene sequences 

rule samtools_faidx:
    input:
        #"{sample}.fa",
        genome=config["ref"] if os.path.exists(config["ref"]) else "results/assembly/pilon/"+ wgs_name + ".fasta"
    output:
        #"{sample}.fa.fai",
        ref_path + ".fai" if os.path.exists(config["ref"]) else "results/assembly/pilon/"+ wgs_name + ".fasta.fai" 
    log:
        "logs/samtools/reference_faidx.log",
    params:
        extra="",  # optional params string
    wrapper:
        "v2.2.0/bio/samtools/faidx"


rule samtools_whole_consensus:
    input:
        bam="results/hisat2/mapped/{sample}.sorted.bam",
        bam_idx="results/hisat2/mapped/{sample}.sorted.bam.bai",
        #gene_list="results/gene_sequences/gene_list.tsv"
    output:
        "results/consensus/{sample}/{sample}_WholeConsensus.fasta",
    conda:
        "../envs/env.yaml"
    log:
        "logs/obtain_gene_seqs/samtools_consensus/{sample}.log",
    params:
        #d=config["software"]["consensus"]["call_d"],
        d=config["software"]["samtools"]["consensus"]["min_depth"]
    shell:
        "samtools consensus  -f fasta {input.bam} -d {params.d} -o {output} --show-del yes -a --show-ins yes 2>> {log}"


rule extract_gene_sequences:
    input:
        gene_bed = "results/consensus/gene_list.bed",
        consensus= "results/consensus/{sample}/{sample}_WholeConsensus.fasta"
    output:
        "results/consensus/{sample}/{sample}_gene_seqs.fasta",
    conda:
        "../envs/env.yaml"
    log:
        "logs/obtain_gene_seqs/bedtools_extractgenes/{sample}.log"
    shell:
        "bedtools getfasta -name -fi {input.consensus} -bed {input.gene_bed} -fo {output} 2> {log}"


## filter out genes that do not have sufficient coverage for each sample

rule filter_out_genes:
    input:
        "results/consensus/{sample}/{sample}_gene_seqs.fasta"
    output:
        "results/consensus/{sample}/{sample}_minCov_genes.txt",
        "results/consensus/{sample}/{sample}_minCov_ConsensusSeqs.fasta"
    log:
        "logs/obtain_gene_seqs/filter_genes/{sample}.log"
    params:
        min_coverage = config["software"]["FilterGenes"]["min_coverage"]
    script:
        "../scripts/FilterGenes.py"

## get gene intersection of all samples with sufficient coverage

rule create_subset:
    input:
        expand("results/consensus/{sample}/{sample}_minCov_genes.txt", sample=list(samples.index))
    output:
        "results/consensus/Subset_Genes.txt"
    shell:
        """
        common_lines=$(cat {input[0]})
        for file in {input}; do 
            common_lines=$(comm -12 <(sort <(echo "$common_lines")) <(sort <(cat "$file")))
        done 
        echo "$common_lines" > {output}
        """

## preparation for MSA -> merge the sequences from the genes in intersection

rule create_merged_sequences:
    input:
        fastas=expand("results/consensus/{sample}/{sample}_minCov_ConsensusSeqs.fasta", sample=list(samples.index)),
        genes="results/consensus/Subset_Genes.txt"
    output:
        "results/consensus/MSAinput.fasta"
    script:
        "../scripts/MergeSequences.py"
