
rule differential_gene_expression:
        input:
            counts= "results/featurecounts/allsamples.featureCounts",
            group_assignment= config["group_assignments"]
        output:
            table= "results/DGEAnalysis/DEG_results.tsv",
            volcano="results/DGEAnalysis/volcano_DESeq2.pdf",
            heatmap = "results/DGEAnalysis/heatmap_DESeq2.pdf",
            venn= "results/DGEAnalysis/venn.pdf"
        params:
            alpha = config["software"]["DEG"]["alpha"]
        log:
            "logs/DGE_Rscript.log"
        script:
            "../scripts/DGE_analysis.R"
