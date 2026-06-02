# Download Data -----------------------------------------------------------
Common_DEGs <- read.csv("Data/RA_UC_conserved_DEGs.csv")



common_up <- Common_DEGs[
  Common_DEGs$logFC_UC > 0 & Common_DEGs$logFC_RA > 0,
]

common_down <- Common_DEGs[
  Common_DEGs$logFC_UC < 0 & Common_DEGs$logFC_RA < 0,
]



common_up_ids <- bitr(common_up$X,
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)

ego_common_up <- enrichGO(
  gene      = common_up_ids$ENTREZID,
  universe  = bg_ids_RA$ENTREZID,   # optional but good
  OrgDb     = org.Hs.eg.db,
  ont       = "BP",
  readable  = TRUE
)


common_down_ids <- bitr(common_down$X,
                        fromType = "SYMBOL",
                        toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db)

ego_common_down <- enrichGO(
  gene      = common_down_ids$ENTREZID,
  universe  = bg_ids_RA$ENTREZID,
  OrgDb     = org.Hs.eg.db,
  ont       = "BP",
  readable  = TRUE
)



p_up <- dotplot(ego_common_up, showCategory = 15) +
  ggtitle("Shared UP (RA + UC)")

p_down <- dotplot(ego_common_down, showCategory = 10) +
  ggtitle("Shared DOWN (RA + UC)")

p_up + p_down


ggsave("Figures/shared_GO_up_down.png",
       p_up + p_down, width = 10, height = 8, dpi = 300)

