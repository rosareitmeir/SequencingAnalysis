log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(dplyr)
library(DESeq2)
library(ggplot2)
library(pheatmap)


## Passing argumentd from command line 
featurecount_path = snakemake@input[['counts']]
group_assignment_path = snakemake@input[['group_assignment']]

output_table=snakemake@output[["table"]]
output_volcano=snakemake@output[["volcano"]]
output_heatmap=snakemake@output[["heatmap"]]
output_ma=snakemake@output[["ma"]]
output_pca=snakemake@output[["pca"]]

alpha= snakemake@params[["alpha"]]
top= snakemake@params[["topn"]]


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

## pre filtering 
total_count_threshold <- length(condition_vec)

# Calculate the total count per gene
total_counts <- rowSums(count_data)

# Filter out genes with low counts
filtered_counts <- count_data[total_counts >= total_count_threshold, ]


## running deseq2

dds <- DESeqDataSetFromMatrix(countData =filtered_counts,
                              colData= condition, ~ group) 
dds <- DESeq(dds)

res<- results(dds,alpha=alpha)

#  dataframe for table and plotting 
deseq_df<- data.frame(GeneID=rownames(res), log2FoldChange=res$log2FoldChange, padj= res$padj, pval=res$pvalue)
deseq_df<- deseq_df %>%
  mutate(gene_type = case_when(log2FoldChange > 0 & padj < 0.05 ~ "upregulated",
                               log2FoldChange < 0 & padj < 0.05 ~ "downregulated",
                               TRUE ~ "not significant"))   


## volcano plot 

volcano_plt <- ggplot(data=deseq_df, aes(x=log2FoldChange, y= -log10(padj)  , col = gene_type) )+ 
           geom_point() + 
           theme(legend.title= element_blank()) +
           xlab("log2 Fold Change") + 
           ylab("-log10( adj. p-value)") +
          xlim(-5,5)+
          ylim(0,20)


library(EnhancedVolcano)

volcano <- EnhancedVolcano(deseq_df,
                           lab = NA ,
                           x = "log2FoldChange",
                           y = 'padj',
                           pCutoff = 0.05,
                           FCcutoff = 2,
                          # xlim = xlim,
                          # ylim=ylim,
                           title="", subtitle="")


ggsave(output_volcano, plot =volcano, device = 'pdf')


## ma plot 
pdf(output_ma)  
plot_ma <- plotMA(res)
print(plot_ma)
dev.off()


## heatmap plot

select <-  deseq_df %>%filter (gene_type != "not significant") %>% top_n(top, abs(log2FoldChange)) %>% pull(GeneID)

ntd<-normTransform(dds)

df <- data.frame(group=colData(dds)[,c("group")])
rownames(df) <- colnames(dds)

pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, filename = output_heatmap)

## PCA 
vsd<- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c("group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca_plt <- ggplot(pcaData, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 


ggsave(output_pca, plot =pca_plt, device = 'pdf')



## table output
# quote= F, sep = "\t", col.names =T, row.names=F
write.csv2(deseq_df, file= output_table, row.names = F, sep="\t")

