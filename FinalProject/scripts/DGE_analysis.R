library(dplyr)
library(DESeq2)
library(ggplot2)

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Passing argumentd from command line 
featurecount_path = snakemake@input[['counts']]
group_assignment_path = snakemake@input[['group_assignment']]

output_table=snakemake@output[["table"]]
output_volcano=snakemake@output[["volcano"]]
output_heatmap=snakemake@output[["heatmap"]]
#output_venn=snakemake@output[["venn"]]


alpha= snakemake@params[["alpha"]]



#### Import of Data 
# just deleting the first columns of the feature counts matrix 
#featurecount_path= "/home/rosa/SequencingAnalysis/FinalProject/example/count.featureCounts"
count_data <- read.csv2(featurecount_path, 
                        header=T, sep="\t", blank.lines.skip = T, skip=1 )

row.names(count_data) <- count_data$Geneid
count_data <- count_data[, -c(1:6)]

colnames(count_data)<- gsub("results.hisat2.mapped.", "", colnames(count_data)) 
colnames(count_data)<- gsub(".bam", "", colnames(count_data)) 


## importing the condition table

condition <- read.csv2(group_assignment_path, sep="\t", row.names = NULL)
condition$group <- factor(condition$group)
#dataframe(sample=colnames(example)) %>% left_join(condition, by ="sample") %>% pull(sample)
#sort it in the order of the count matrix
condition <- condition[match(colnames(count_data), condition$sample),]

condition_vec <- factor(condition$group)



#### DESeq2
dds <- DESeqDataSetFromMatrix(countData =count_data,
                              colData= condition, ~ group) 
dds <- DESeq(dds)

res<- results(dds,alpha=alpha)

# for volcano plot and also nice dataframe for output 
deseq_df<- data.frame(GeneID=rownames(res), log2FoldChange=res$log2FoldChange, padj= res$padj, pval=res$pvalue)
deseq_df<- deseq_df %>%
  mutate(gene_type = case_when(log2FoldChange > 0 & padj < 0.05 ~ "upregulated",
                               log2FoldChange < 0 & padj < 0.05 ~ "downregulated",
                               TRUE ~ "not significant"))   

create_volcano <- function(df, method, intersection, path){

  return(ggplot(data=df, aes(x=log2FoldChange, y=-log10(padj)  , col = gene_type) )+ 
           geom_point() + 
           theme(legend.title= element_blank()) +
           xlab("log2 Fold Change") + 
           ylab("-log10( adj. p-value)"))
  
  
}


volcano_plt <- create_volcano(deseq_df)
ggsave(output_volcano, plot =  volcano_plt, device = 'pdf')

## heatmap plot
library("pheatmap")

top= snakemake@params[["topn"]]

select <-  deseq_df %>%filter (gene_type != "not significant") %>% top_n(top, abs(log2FoldChange)) %>% pull(genes)

ntd<-normTransform(dds)

df <- data.frame(group=colData(dds)[,c("group")])
rownames(df) <- colnames(dds)
hm_plt <- pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

ggsave(output_heatmap, plot =  hm_plt, device = 'pdf')









## table output
write.csv2(deseq_df, file= output_table, quote= F, sep = "\t", col.names =T, row.names=F)

