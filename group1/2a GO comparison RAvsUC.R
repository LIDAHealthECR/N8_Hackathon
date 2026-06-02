gene_list <- list(
  UC_up   = up_ids_UC$ENTREZID,
  UC_down = down_ids_UC$ENTREZID,
  RA_up   = up_ids_RA$ENTREZID,
  RA_down = down_ids_RA$ENTREZID
)

up_ids_UC <- bitr(up_UC$X, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
down_ids_UC <- bitr(down_UC$X, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

up_ids_RA <- bitr(up_RA$X, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
down_ids_RA <- bitr(down_RA$X, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

compare_full <- compareCluster(
  geneCluster = gene_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  ont = "BP"
)


p_compare_full <- dotplot(compare_full, showCategory = 10) +
  ggtitle("GO comparison: RA vs UC (Top 10 per group)")



ggsave(
  "Figures/RA_UC_GO_full_comparison.png",
  p_compare_full,
  width = 10,
  height = 20,
  dpi = 300
)
