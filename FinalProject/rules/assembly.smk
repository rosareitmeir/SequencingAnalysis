
## TODO: put to the normalize rule
def get_files_for_assembly(wildcards):
    if not cont_table.empty :
        return ["results/decontamination/"+ wildcards.sample +"_1.fastq.gz","results/decontamination/"+ wildcards.sample +"_2.fastq.gz"]
    elif (config["use_trimmed_data"]):
        return ["results/trimmed/"+ wildcards.sample + "_1.fastq", "results/trimmed/"+ wildcards.sample + "_1.fastq"]
    else:
        return [ samples.loc[wildcards.sample][0], samples.loc[wildcards.sample][1]]

### TODO chaning the input to the normalized one
rule assembly:
    input:
        get_files_for_assembly
    params:
        output_dir = "results/assembly/megahit/{sample}",
        out_prefix = "{sample}"
    output:
        "results/assembly/megahit/{sample}/checkpoints.txt",
        "results/assembly/megahit/{sample}/done",
        directory("results/assembly/megahit/{sample}/intermediate_contigs"),
        "results/assembly/megahit/{sample}/options.json",
        "results/assembly/megahit/{sample}/{sample}.contigs.fa",
        "results/assembly/megahit/{sample}/{sample}.log"
    log: 
        "logs/megahit/{sample}.log"
    conda:
        "../envs/env.yaml"
    threads:
        config["software"]["megahit"]["threads"]
    shell:
        """
        rm -rf {params.output_dir}
        megahit -1 {input[0]} -2 {input[1]} -t {threads} -o {params.output_dir} --out-prefix {params.out_prefix} &> {log}
        """

## BLAST the assemblies
## TODO ANA : adding those rules
rule build_database:
    input:
        config["ref"]
    output:
        multiext("results/blastn/database/db", ".ndb", ".nhr", ".nin",".not", ".nsq", ".ntf",".nto" )
    log:
        "blast_build_db.log"
    conda:
        "../envs/env.yaml"

    shell:
        "makeblastdb -in {input} -dbtype nucl -out results/blastn/database/db &> {log} "

## TODO ANA is the input of the query correct ( output of the assembly)
rule blastn_assemblies:
    input:
        query = "results/assembly/megahit/{sample}/{sample}.contigs.fa",
        blastdb=multiext("results/blastn/database/db", ".ndb", ".nhr", ".nin",".not", ".nsq", ".ntf",".nto" )
    output:
        "results/blastn/{sample}.blast.txt"
    log:
        "logs/{sample}.blast.log"
    threads:
        config["software"]["blastn"]["threads"]
    params:
        # Usable options and specifiers for the different output formats are listed here:
        # https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/blast/blastn.html.
        format="6 qseqid sseqid evalue pident length mismatch gapopen qstart qend sstart send",
        extra=config["software"]["blastn"]["extra"]
    conda:
        "../envs/env.yaml"
    shell:
        """blastn -query {input.query} -db results/blastn/database/db -out {output} -outfmt '{params.format}' {params.extra} &> {log}
           echo "qseqid\tsseqid\tevalue\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend" | cat - {output}> results/blastn/temp_file && mv results/blastn/temp_file {output}
        """

## TODO ANA adding the rule and the two defintions
rule find_best_reference:
    input:
        expand("results/blastn/{sample}.blast.txt", sample=list(samples.index))
    output:
        "results/blastn/best_reference.tsv"
    params:
        criteria= config["software"]["blastn"]["criteria"],
        multi_fasta_file=config["ref"]
    conda:
        "../envs/env.yaml"
    script:
        "../scripts/choose_best_reference.py"


##
def get_best_reference(wildcards, path):
    sample=wildcards.sample
    ref_table = pd.read_csv(path, sep='\t', header=None, names=["sample", "ref"])
    refID= ref_table.loc[ref_table["sample"] == sample, "ref"].iloc[0]
    return("reference/"+ refID + ".fasta")

def get_refID(wildcards, path):
    sample = wildcards.sample
    ref_table = pd.read_csv(path,sep='\t',header=None,names=["sample", "ref"])
    refID = ref_table.loc[ref_table["sample"] == sample, "ref"].iloc[0]
    return(refID)


## TODO: also adding this one, deleting old version of the scaffolding and deleting all the merge rules
rule scaffolding:
    input:
        "results/assembly/megahit/{sample}/{sample}.contigs.fa",
        "results/blastn/best_reference.tsv",

    params:
        output_dir = "results/assembly/ragtag/{sample}",
        ref= lambda wildcards: get_best_reference(wildcards,"results/blastn/best_reference.tsv")
    output:
        "results/assembly/ragtag/{sample}/ragtag.scaffold.asm.paf",
        "results/assembly/ragtag/{sample}/ragtag.scaffold.asm.paf.log",
        "results/assembly/ragtag/{sample}/ragtag.scaffold.err",
        "results/assembly/ragtag/{sample}/ragtag.scaffold.fasta"

    log:
        "logs/ragtag/{sample}/{sample}.log"
    conda:
        "../envs/env.yaml"
    threads:
        config["software"]["ragtag"]["threads"]
    shell:
        "ragtag.py scaffold {params.ref} {input[0]} -o {params.output_dir} -u -t {threads} &> {log}"

## TODO ANA: replace old rule with this new one
## since the output of the scaffolding merge may contain multiple scaffolds, we reduce it to one scaffold -> the longest and mapped scaffold
rule extract_max_scaffold_fasta:
    input:
        "results/assembly/ragtag/{sample}/ragtag.scaffold.fasta",
        "results/blastn/best_reference.tsv",

    output:
        "results/assembly/ragtag/{sample}/max_scaffold.fasta",

    params:
        sample="{sample}",
        ref=lambda wildcards: get_refID(wildcards, "results/blastn/best_reference.tsv")
    conda:
        "../envs/env.yaml"
    script:
        "../scripts/helper_maxseq.py"



## TODO ANA nothing changed here

#Generate multifasta file for clustalo input
rule get_infile:
    input:
        expand("results/assembly/ragtag/{sample}/max_scaffold.fasta", sample=list(samples.index))
    output:
        "results/assembly/clustalo_infile.txt"
    shell:
        "cat {input} > results/assembly/clustalo_infile.txt"

rule msa_clustalo:
    input:
        "results/assembly/clustalo_infile.txt"
    output:
       "results/msa/max_scaffolds.msa.fa"
    params:
        extra=""
    log:
        "logs/clustalo/test/msa.log"
    threads:
        config["software"]["clustalo"]["threads"]
    wrapper:
        "v1.31.1/bio/clustalo"

