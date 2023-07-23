################### Multiple Sequence Alignment ###########################################

rule msa_clustalo:
	input:
		"results/consensus/MSAinput.fasta"
	output:
		"results/msa/MSA.msa.fa"
	params:
		extra=""
	log:
		"logs/clustalo/msa.log"
	threads:
		config["software"]["clustalo"]["threads"]
	wrapper:
		"v1.31.1/bio/clustalo"

################### Phylogenetic Tree Generation ###########################################

rule fasttree:
    input:
        alignment="results/msa/MSA.msa.fa" # Input alignment file
    output:
        tree="results/fasttree/MSA.nwk",  # Output tree file
    log:
        "logs/fasttree/fasttree.log",
    params:
        extra="",
    wrapper:
        "v1.31.1/bio/fasttree"

################### Phylogenetic Tree Visualization ###########################################

rule toytree:
    input:
        "results/fasttree/MSA.nwk"
    output:
        "results/toytree/tree_plot.pdf"
    log:
        "logs/toytree/alignment.log"
    conda:
        "../envs/env.yaml"
    script:
        "../scripts/TreeVisualization.py"
