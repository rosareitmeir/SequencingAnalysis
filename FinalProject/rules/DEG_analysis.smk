
rule differential_gene_expression:
        input:
            counts= "example/allsamples.featureCounts" ,# results/featurecounts/allsamples.featureCounts",
            group_assignment= config["group_assignment"]
        output:
            table= "results/DGEAnalysis/DEG_results.tsv",
            volcano="results/DGEAnalysis/volcano.pdf",
            heatmap = "results/DGEAnalysis/heatmap.pdf",
            ma= "results/DGEAnalysis/MA.pdf",
            pca="results/DGEAnalysis/PCA.pdf"
        params:
            alpha = config["software"]["DEG"]["alpha"],
            topn = config["software"]["DEG"]["topn"]
        log:
            "logs/DGE_Rscript.log"
        script:
            "../scripts/DGE_analysis.R"
