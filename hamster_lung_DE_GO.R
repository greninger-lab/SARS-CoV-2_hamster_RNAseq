#normalize raw counts of lung, infected and uninfected with SARS-CoV-2
#4 replicates combined into single tube

setwd("~/Desktop")

library(tidyverse)
library(DESeq2)
library(pheatmap)
library(org.Mm.eg.db)
library(clusterProfiler)
library(viridis)


metadata_file <- file.choose()   
metadata <- read_csv(metadata_file)

lung_meta <- filter(metadata, tissue == "lung")
lung_meta <- column_to_rownames(lung_meta, "sample") 

counts_file <- file.choose() 
counts <- read_csv(counts_file)
counts <- column_to_rownames(counts, "X1") 

lung_counts <- counts[, rownames(lung_meta)]

#lung
#normalization
dds_lung <- DESeqDataSetFromMatrix(countData = lung_counts, colData = lung_meta, design= ~ Infection)
keep <- (rowSums(counts(dds_lung)) >= nrow(lung_meta)) #pre-filter to remove any genes without avg 1 count per sample
dds_lung <- dds_lung[keep,]
dds_lung <- estimateSizeFactors(dds_lung)
norm_counts_lung <- as.data.frame(counts(dds_lung, normalized=TRUE))
colnames(norm_counts_lung) <- c("infected_treated", "infected_untreated", "uninfected")
write.csv(norm_counts_lung, "~/Desktop/hamster_norm_counts_lung.csv")

norm_counts_lung_calcs <- norm_counts_lung + 1
norm_counts_lung_calcs$treated_vs_uninfected_ratio <- norm_counts_lung_calcs$infected_treated / norm_counts_lung_calcs$uninfected
norm_counts_lung_calcs$untreated_vs_uninfected_ratio <- norm_counts_lung_calcs$infected_untreated / norm_counts_lung_calcs$uninfected
norm_counts_lung_calcs$treated_vs_untreated_ratio <- norm_counts_lung_calcs$infected_treated / norm_counts_lung_calcs$infected_untreated

write.csv(norm_counts_lung_calcs, "~/Desktop/lung_counts_and_ratios.csv")

#heatmap 
norm_counts_corr_lung <- filter(norm_counts_lung, uninfected > 5) + 1
FC_lung <- norm_counts_corr_lung / norm_counts_corr_lung$uninfected
l2FC_lung <- log2(FC_lung)
l2FC_lung <- rownames_to_column(l2FC_lung)
most_var_lung <- arrange(l2FC_lung, desc(abs(infected_untreated)))
most_var_lung <- most_var_lung[1:50,]
rownames(most_var_lung) <- most_var_lung$rowname
most_var_lung <- most_var_lung[,-(c(1:2))] #removal of unrelated sample
colnames(most_var_lung) <- c("Infected", "Uninfected")
df_lung <- lung_meta[, c("tissue", "Infection")]  
df_lung <- df_lung[-1, -1, drop=FALSE] #removal of unrelated sample
rownames(df_lung) <- c("Infected", "Uninfected")

alt_lung_names <- read_csv("most_var_lung.csv")
rownames(most_var_lung) <- alt_lung_names$manual_gene_name

png(filename="~/Desktop/hamster_lung_heatmap.png", width=6, height=8, units="in", res=300)
print(pheatmap(most_var_lung, annotation_col=df_lung, show_colnames = F, cluster_cols=FALSE)) #creates dendrogram for samples and genes. if cluster_cols=FALSE it does not cluster
dev.off()



#Pathway Analysis
#restricting search space to genes with > 5 counts in uninfected (to match heatmap), > abs(2) l2fc (ie > 4 fold change)
#Create results by calculating l2fc of infected_untreated / uninfected

lung_res <- rownames_to_column(norm_counts_lung)
lung_res$l2fc <- log((lung_res$infected_untreated / lung_res$uninfected), base = 2)

sig_lung_genes <- filter(lung_res, abs(l2fc) > 2, infected_untreated > 1, uninfected > 5)

mm <- org.Mm.eg.db
entrez <- AnnotationDbi::select(mm, 
                                keys = lung_res$rowname,
                                columns = c("ENTREZID", "SYMBOL"),
                                keytype = "SYMBOL")
entrez <- distinct(entrez, SYMBOL, .keep_all = TRUE) #removes duplicate entries if there was multiple mapping
lung_res$entrez <- entrez$ENTREZID

write.csv(lung_res, "~/Desktop/lung_results.csv", row.names = FALSE)

lung_res <- filter(lung_res, entrez != "NA") #remove entries where there was not mapping to entrez id. this should mostly affect poorly characterized genes 
lung_res <- filter(lung_res, infected_untreated > 1, uninfected > 5)


geneList <- lung_res$l2fc     
names(geneList) = as.character(lung_res$entrez)   
gene <- dplyr::filter(lung_res, abs(l2fc) > 2)

write.csv(gene, "~/Desktop/lung_sig_genes.csv", row.names = FALSE)

gene <- gene$entrez

#GO enrichment
ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",        #BP, CC, or MF
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.1,
                readable      = TRUE)
go_df <- as.data.frame(ego@result)
#remove duplicate categories, keeping one with lowest pvalue
go_df <- distinct(go_df, geneID, .keep_all = TRUE)
go_df_20 <- go_df[1:20,]

des <- go_df_20$Description
des <- rev(des)
go_df_20$Description <- factor(go_df_20$Description, levels = go_df_20$Description[order(-(go_df_20$pvalue))])

#GO Figure Lung 
png(filename="~/Desktop/lung_GO_plot_l2fc2_BP.png", width=8, height=6, units="in", res=300)
ggplot(go_df_20, aes(Description, Count, fill=p.adjust)) + 
  geom_col() +
  coord_flip() +
  scale_x_discrete(labels = str_wrap(des, width = 60)) +
  scale_fill_viridis() +
  ylab("Number Enriched") +
  xlab("GO Term") +
  theme(
    legend.position = "right",
    text = element_text(size = 10),
    axis.title.x = element_text(size = 10), 
    axis.title.y = element_text(size = 10))
dev.off()
