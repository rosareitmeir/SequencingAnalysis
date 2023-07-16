### Differential Gene Expression Analysis
library(dplyr)
library(DESeq2)
library(ggplot2)

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Passing argumentd from command line 
# 1 feature counts input 
args = commandArgs(trailingOnly = TRUE)

featurecount_path = snakemake@input[['counts']]
group_assignment_path = snakemake@input[['group_assignment']]

output_table=snakemake@output[["table"]]
output_volcano=snakemake@output[["volcano"]]
output_heatmap=snakemake@output[["heatmap"]]
output_venn=snakemake@output[["venn"]]

                              
alpha= snakemake@params[["alpha"]]



#### Import of Data 
# just deleting the first columns of the feature counts matrix 
# "/home/rosa/SequencingAnalysis/FinalProject/results/featurecounts/allsamples.featureCounts"
count_data <- read.csv2(featurecount_path, 
                     header=T, sep="\t", blank.lines.skip = T, skip=1 )

row.names(count_data) <- count_data$Geneid
count_data <- count_data[, -c(1:6)]

colnames(count_data)<- gsub("results.hisat2.mapped.", "", colnames(count_data)) 
colnames(count_data)<- gsub(".bam", "", colnames(count_data)) 


###### to be deleted 
# "/home/rosa/SequencingAnalysis/FinalProject/example/count.featureCounts"

colnames(count_data)<- gsub("results.star.pe.", "", colnames(count_data)) 
colnames(count_data)<- gsub(".pe_aligned.sam", "", colnames(count_data)) 

cond <- data.frame(sample=colnames(count_data), group=c("n", "t", "n", "t","n", "t"))

condition <- cond

## importing the condition table 
# example 
condition <- read.csv2(group_assignment_path, sep="\t", row.names = NULL)
condition$group <- factor(condition$group)
#dataframe(sample=colnames(example)) %>% left_join(condition, by ="sample") %>% pull(sample)
#sort it in the order of the count matrix
condition <- condition[match(colnames(example), condition$sample),]

condition_vec <- condition$group



#### DESeq2
 dds <- DESeqDataSetFromMatrix(countData =count_data,
                                 colData= condition, ~ group) 
dds <- DESeq(dds)
  
res<- results(dds,alpha=alpha)

# for volcano plot and also nice dataframe for output 
deseq_df<- data.frame(log2FoldChange=res$log2FoldChange, padj= res$padj,pval=res$pvalue, genes=rownames(res))
deseq_df<- deseq_df %>%
  mutate(gene_type = case_when(log2FoldChange > 0 & padj < 0.05 ~ "upregulated",
                               log2FoldChange < 0 & padj < 0.05 ~ "downregulated",
                               TRUE ~ "not significant"))   

### edgeR

library(edgeR)

edge_list <- DGEList(counts=as.matrix(count_data), group=factor(group))

# filter out genes with low counts
keep <- filterByExpr(y=edge_list)
edge_list <- edge_list[keep, , keep.lib.sizes=FALSE]

# normalization / scaling 
edge_list <- calcNormFactors(object=edge_list)

# dispersion estimation
edge_list <- estimateDisp(y=edge_list )

# finally: DGE calling 

edge_res <- exactTest(object=edge_list)

# get adjusted p values 

edge_res <- data.frame(topTags(object=edge_res, n = "Inf"))

edge_df<- data_frame(log2FoldChange=edge_res$logFC, padj= edge_res$FDR,pval=edge_res$PValue, genes=rownames(edge_res))
edge_df <- edge_df %>%
  mutate(gene_type = case_when(log2FoldChange > 0 & padj < 0.05 ~ "upregulated",
                               log2FoldChange < 0 & padj < 0.05 ~ "downregulated",
                               TRUE ~ "not significant"))   


### limma
library(limma)

limma.d <-  DGEList(counts=as.matrix(count_data), group=factor(group))
limma.d <- calcNormFactors(limma.d)

# filter out genes with low counts 
drop <- which(apply(cpm(limma.d), 1, max) < 1)
limma.d <- limma.d[-drop,] 

mm <- model.matrix(~ group)
voom <- voom(limma.d, mm)

f1 = lmFit(voom, mm)
ef1 = eBayes(f1)

limma_res <- topTable(ef1, n=Inf)

limma_df<- data_frame(log2FoldChange=limma_res$logFC, padj= limma_res$adj.P.Val,pval=limma_res$P.Value, genes=rownames(limma_res))
limma_df <- limma_df %>%
  mutate(gene_type = case_when(log2FoldChange > 0 & padj < 0.05 ~ "upregulated",
                               log2FoldChange < 0 & padj < 0.05 ~ "downregulated",
                               TRUE ~ "not significant"))   

### Visualization

## Venn Diagram and Intersection of the three methods 
library(ggvenn)
limma_signif <- filter(limma_df, gene_type != "not significant") %>% pull(genes)
edge_signif <- filter(edge_df, gene_type != "not significant") %>% pull(genes)
deseq_signif <- filter(deseq_df, gene_type != "not significant") %>% pull(genes)


ggvenn(list(limma=limma_signif, edgeR=edge_signif, DESeq2=deseq_signif),
     
       fill_color  = c( "#E69F00", "#56B4E9", "#999999"))

intersection <- deseq_df %>% filter ( genes %in% intersect(intersect(limma_signif, edge_signif), deseq_signif) )

ggsave(output_venn, plot =  hm_plt, device = 'pdf')


## create volcano plot 


create_volcano <- function(df, path){
  
  return(ggplot(data=df, aes(x=log2FoldChange, y=-log10(padj) , col= gene_type  ) )+ 
           geom_point() + 
           theme(legend.title= element_blank()) +
           xlab("log2 Fold Change") + 
           ylab("-log10( adj. p-value)"))
  
}


volcano_plt <- create_volcano(deseq_df)
ggsave(output_volcano, plot =  volcano_plt, device = 'pdf')



## heatmap 

library("pheatmap")
top=50

select <-  deseq_df %>%filter (gene_type != "not significant") %>% top_n(top, abs(log2FoldChange)) %>% pull(genes)

ntd<-normTransform(dds)

df <- data.frame(group=colData(dds)[,c("group")])
rownames(df) <- colnames(dds)
hm_plt <- pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

ggsave(output_heatmap, plot =  hm_plt, device = 'pdf')


### TODO  table output 
write.csv2(deseq_df, file= output_table, quote= F, sep = "\t", col.names =T, row.names=F)




 