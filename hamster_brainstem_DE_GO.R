#normalize raw counts of brain, infected and uninfected with SARS-CoV-2
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

brain_meta <- filter(metadata, tissue == "brain stem")
brain_meta <- column_to_rownames(brain_meta, "sample")

counts_file <- file.choose() 
counts <- read_csv(counts_file)
counts <- column_to_rownames(counts, "X1") 

brain_counts <- counts[, rownames(brain_meta)]

#brain
#normalization
dds_brain <- DESeqDataSetFromMatrix(countData = brain_counts, colData = brain_meta, design= ~ Infection) 
keep_brain <- (rowSums(counts(dds_brain)) >= nrow(brain_meta))  #pre-filter to remove any genes without avg 1 count per sample
dds_brain <- dds_brain[keep_brain,]
dds_brain <- estimateSizeFactors(dds_brain)
norm_counts_brain <- as.data.frame(counts(dds_brain, normalized=TRUE))
colnames(norm_counts_brain) <- c("infected_treated", "infected_untreated", "uninfected")
write.csv(norm_counts_brain, "~/Desktop/hamster_norm_counts_brain.csv")

norm_counts_brain_calcs <- norm_counts_brain + 1
norm_counts_brain_calcs$treated_vs_uninfected_ratio <- norm_counts_brain_calcs$infected_treated / norm_counts_brain_calcs$uninfected
norm_counts_brain_calcs$untreated_vs_uninfected_ratio <- norm_counts_brain_calcs$infected_untreated / norm_counts_brain_calcs$uninfected
norm_counts_brain_calcs$treated_vs_untreated_ratio <- norm_counts_brain_calcs$infected_treated / norm_counts_brain_calcs$infected_untreated

write.csv(norm_counts_brain_calcs, "~/Desktop/brain_counts_and_ratios.csv")

#heatmap 
norm_counts_corr_brain <- filter(norm_counts_brain, uninfected > 1) + 1
FC_brain <- norm_counts_corr_brain / norm_counts_corr_brain$uninfected
l2FC_brain <- log2(FC_brain)
l2FC_brain <- rownames_to_column(l2FC_brain)
most_var_brain <- arrange(l2FC_brain, desc(abs(infected_untreated)))
most_var_brain <- most_var_brain[1:50,]
rownames(most_var_brain) <- most_var_brain$rowname
most_var_brain <- most_var_brain[,-(c(1:2))] #removal of unrelated sample
colnames(most_var_brain) <- c("Infected", "Uninfected")
df_brain <- brain_meta[, c("tissue", "Infection")]  
df_brain <- df_brain[-1, -1, drop=FALSE] #removal of unrelated sample
rownames(df_brain) <- c("Infected", "Uninfected")

alt_brain_names <- read_csv("most_var_brain.csv")
rownames(most_var_brain) <- alt_brain_names$manual_gene_name

png(filename="~/Desktop/hamster_brain_heatmap.png", width=6, height=8, units="in", res=300)
print(pheatmap(most_var_brain, annotation_col=df_brain, show_colnames = F, cluster_cols=FALSE)) #creates dendrogram for samples and genes. if cluster_cols=FALSE it does not cluster
dev.off()

#Pathway Analysis
#restricting search space to genes with > 5 counts in uninfected (to match heatmap), > abs(2) l2fc (ie > 4 fold change)
#Create results by calculating l2fc of infected_untreated / uninfected

brain_res <- rownames_to_column(norm_counts_brain)
brain_res$l2fc <- log((brain_res$infected_untreated / brain_res$uninfected), base = 2)

sig_brain_genes <- filter(brain_res, abs(l2fc) > 2, infected_untreated > 1, uninfected > 1)

mm <- org.Mm.eg.db
entrez <- AnnotationDbi::select(mm, 
                                keys = brain_res$rowname,
                                columns = c("ENTREZID", "SYMBOL"),
                                keytype = "SYMBOL")
entrez <- distinct(entrez, SYMBOL, .keep_all = TRUE) #removes duplicate entries if there was multiple mapping
brain_res$entrez <- entrez$ENTREZID

write.csv(brain_res, "~/Desktop/brain_results.csv", row.names = FALSE)

brain_res <- filter(brain_res, entrez != "NA") #remove entries where there was not mapping to entrez id. this should mostly affect poorly characterized genes 
brain_res <- filter(brain_res, infected_untreated > 1, uninfected > 1)


geneList <- brain_res$l2fc     
names(geneList) = as.character(brain_res$entrez)   
gene <- dplyr::filter(brain_res, abs(l2fc) > 2)

write.csv(gene, "~/Desktop/brain_sig_genes.csv", row.names = FALSE)

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

#GO Figure brain 
png(filename="~/Desktop/brain_GO_plot_l2fc2_BP.png", width=8, height=6, units="in", res=300)
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
