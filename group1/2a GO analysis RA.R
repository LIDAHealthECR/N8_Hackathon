library(clusterProfiler)
library(org.Hs.eg.db)


# Download Data -----------------------------------------------------------
RA_DEGs <- read.csv("Data/RA_limma_sig_results.csv")
RA_all <- read.csv("Files/RA_batch_corrected.csv")


# Prep data -----------------------------------------------------------

up_RA <- RA_DEGs[RA_DEGs$logFC > 0, ]
down_RA <- RA_DEGs[RA_DEGs$logFC < 0, ]

# Convert to ENTREZID
# All DEGs
genes <- RA_DEGs$X

# Convert
gene_ids <- bitr(genes,
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)


# Identify background genes

all_genes <- RA_all$geneNames

bg_ids_RA <- bitr(all_genes,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)

# analyse RA data -----------------------------------------------------------


up_ids_RA <- bitr(up_RA$X, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

ego_up_RA <- enrichGO(
  gene = up_ids_RA$ENTREZID,
  universe = bg_ids_RA$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP"
)


down_ids_RA <- bitr(down_RA$X, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

ego_down_RA <- enrichGO(
  gene = down_ids_RA$ENTREZID,
  universe = bg_ids_RA$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP"
)



library(enrichplot)


# Create dotplots
p_up <- dotplot(ego_up_RA, showCategory = 15) +
  ggtitle("Upregulated GO terms") +
  theme(plot.title = element_text(hjust = 0.5))

p_down <- dotplot(ego_down_RA, showCategory = 15) +
  ggtitle("Downregulated GO terms") +
  theme(plot.title = element_text(hjust = 0.5))

# Combine
p_up + p_down



p_up_bar <- barplot(ego_up_RA, showCategory = 10) +
  ggtitle("Upregulated")

p_down_bar <- barplot(ego_down_RA, showCategory = 10) +
  ggtitle("Downregulated")

p_up_bar + p_down_bar


library(clusterProfiler)

# Prepare gene lists
gene_list <- list(
  Up   = up_ids$ENTREZID,
  Down = down_ids$ENTREZID
)

compare <- compareCluster(
  geneCluster = gene_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  ont = "BP"
)

dotplot(compare)


library(patchwork)

combined_dot <- p_up + p_down

ggsave(
  filename = "Figures/RA_GO_dotplot_up_vs_down.png",
  plot = combined_dot,
  width = 12,
  height = 6,
  dpi = 300
)


combined_bar <- p_up_bar + p_down_bar

ggsave(
  filename = "Figures/RA_GO_barplot_up_vs_down.png",
  plot = combined_bar,
  width = 12,
  height = 6,
  dpi = 300
)

compare_plot <- dotplot(compare)

ggsave(
  filename = "Figures/RA_GO_compareCluster.png",
  plot = compare_plot,
  width = 8,
  height = 6,
  dpi = 300
)
