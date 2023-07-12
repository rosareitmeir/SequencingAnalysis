## extracting gene list
rule create_gene_list:
    input:
        "results/featurecounts/allsamples.featureCounts"
    output:
        "results/gene_sequences/gene_list.tsv"
    # get rid of the header, saving gene name, "chr/contig", gene start and gene end
    shell:
        "tail -n +3 {input} | cut -f 1-4 -d $'\t'   > {output}"

# sort each bam file according to the genomic position of the reads
rule sort_bam:
    input:
        "results/hisat2/mapped/{sample}.bam"
    output:
        "results/hisat2/mapped/{sample}_sorted.bam"
    conda:
        "../envs/env.yaml"
    shell:
        "samtools sort {input} > {output}"

## adds index to the sorted bam files (positions of the reads in the Bam file)
rule index_bam:
    input:
         "results/hisat2/mapped/{sample}_sorted.bam"
    output:
         "results/hisat2/mapped/{sample}_sorted.bam.bai"
    conda:
        "../envs/env.yaml"
    log:
        "logs/samtools/{sample}indexing.log"
    shell:
        "samtools index {input} 2>> {log}"

# Obtaining gene sequences 

rule samtools_faidx:
    input:
        #"{sample}.fa",
        genome=config["ref"] if os.path.exists(config["ref"]) else "results/assembly/pilon/"+ wgs_name + ".fasta"
    output:
        #"{sample}.fa.fai",
        ref_path + ".fai" if os.path.exists(config["ref"]) else "results/assembly/pilon/"+ wgs_name + ".fasta.fai" #TODO check proper naming
    log:
        "logs/samtools/reference_faidx.log",
    params:
        extra="",  # optional params string
    wrapper:
        "v2.2.0/bio/samtools/faidx"

# Obtain one region to run consensus rule
#regions = pd.read_csv("results/gene_sequences/gene_list.tsv", sep='\t')
#regions["Gene"] = regions.iloc[:,1].astype(str)+":"+regions.iloc[:,2].astype(str)+"-"+regions.iloc[:,3].astype(str)
#gene_regions = list(regions.iloc[:,1].astype(str)+":"+regions.iloc[:,2].astype(str)+"-"+regions.iloc[:,3].astype(str))

rule samtools_consensus:
    input:
        bam="results/hisat2/mapped/{sample}_sorted.bam",
        bam_idx="results/hisat2/mapped/{sample}_sorted.bam.bai",
        gene_list="results/gene_sequences/gene_list.tsv"
    output:
       "results/consensus/{sample}/{sample}_ConsensusSeqs.fasta",
    conda:
        "../envs/env.yaml"
    log:
        "logs/obtain_gene_seqs/samtools_consensus/{sample}.log",
    params:
       d=config["software"]["consensus"]["call_d"],
    shell:
        """
        while read a b c d; do 
            seq=$(samtools consensus -r "$b:$c-$d" -f fasta {input.bam} -d {params.d} --show-del yes -a --show-ins yes 2>> {log} | sed 1d ) ;
            echo -e ">$a $b $c $d\n$seq" >> {output} ;
            done < {input.gene_list}
            """
        ## a: gene name , b: chr/contig, c: gene start , d: gene end


## filter out genes that do not have sufficient coverage for each sample
rule filter_out_genes:
    input:
        "results/consensus/{sample}/{sample}_ConsensusSeqs.fasta"

    output:
        "results/consensus/{sample}/{sample}_minCov_genes.txt",
        "results/consensus/{sample}/{sample}_minCov_ConsensusSeqs.fasta"


    log:
        "logs/obtain_gene_seqs/filter_genes/{sample}.log"

    params:
        min_coverage = config["software"]["consensus"]["min_coverage"]

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





































