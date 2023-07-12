### rule for MSA




rule fasttree:
    input:
        alignment="results/msa/max_scaffolds.msa.fa" # Input alignment file
    output:
        tree="results/fasttree/max_scaffolds.nwk",  # Output tree file
    log:
        "logs/fasttree/fasttree.log",
    conda:
        "../envs/env.yaml"
    params:
        extra="",
    wrapper:
        "v1.31.1/bio/fasttree"


rule toytree:
    input:
        "results/fasttree/max_scaffolds.nwk"
    output:
        "results/toytree/tree_plot.pdf"
    log:
        "logs/toytree/alignment.log"
    conda:
        "../envs/env.yaml"
    script:
        "../scripts/tree_view.py"
