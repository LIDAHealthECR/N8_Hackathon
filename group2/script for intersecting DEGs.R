
library(tidyverse)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db) 
library(AnnotationDbi)

# differential gene expression -------------------------------------------------

RA_DEGs <- read.csv(paste0("C:/Users/hazel/OneDrive/Documents/hackathon/Differentially expressed genes for rheumatoid arthritis.csv"))
UC_DEGs <- read.csv(paste0("C:/Users/hazel/OneDrive/Documents/hackathon/Differentially expressed genes for ulcerative colitis.csv"))

colnames(UC_DEGs)
RA_DEGs_filtered <- RA_DEGs %>%
  filter(abs(logFC) > 1 & adj.P.Val < 0.05)

UC_DEGs_filtered <- UC_DEGs %>%
  filter(abs(logFC) > 1 & adj.P.Val < 0.05)

RA_DEGs_vector <- RA_DEGs_filtered$id
UC_DEGs_vector <- UC_DEGs_filtered$id

intersect_DEGs <- intersect(RA_DEGs_vector,UC_DEGs_vector)
length(intersect_DEGs)

intersect_DEGs_df <- left_join(RA_DEGs,UC_DEGs, by = "id") %>%
  filter(id %in% intersect_DEGs) 

intersect_DEGs_df_filtered <- intersect_DEGs_df %>%
  filter((logFC.x > 0 & logFC.y > 0) | (logFC.x < 0 & logFC.y < 0))
# What genes do we have that are not found in the paper?
# removed conflicting regulations, and got 88 intersecting genes

# Data formatting all sig DEGs for ORA 

intersect_DEGs_df_up <- intersect_DEGs_df_filtered %>%
  filter(logFC.y > 1)

intersect_DEGs_df_down <- intersect_DEGs_df_filtered %>%
  filter(logFC.y < 1)

length(intersect_DEGs_df_down$id)



      # ENTREZID annotation mapping
ids <- bitr(intersect_DEGs_df_up$id, fromType = "SYMBOL", toType = "ENTREZID", OrgDb= org.Hs.eg.db)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
df2 = intersect_DEGs_df_up[intersect_DEGs_df_up$id %in% dedup_ids$SYMBOL,]
df2$ENTREZID = dedup_ids$ENTREZID
ENTREZID_gene_list <- df2$logFC.y
names(ENTREZID_gene_list) <- df2$ENTREZID
ENTREZID_gene_list_sorted <- sort(ENTREZID_gene_list, decreasing = TRUE)

# ORA KEGG ~~~~~~~~~~~~~

url <- "https://rest.kegg.jp/link/hsa/pathway"
tryCatch(
  { readLines(url) },
  error = function(e) { print(e) }
)

# ORA KEGG
ora_kegg_down <- enrichKEGG(gene = names(ENTREZID_gene_list_sorted),
                       organism = "hsa",
                       minGSSize = 10,
                       maxGSSize = 500,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       keyType = "kegg")


# plot
dotplot(ora_kegg_down,  showCategory= 10)

# ORA GO ~~~~~~~~~~~

ora_go_up <- enrichGO(
  gene = names(ENTREZID_gene_list_sorted),  
  OrgDb = org.Hs.eg.db,  
  keyType = "ENTREZID", 
  ont = "BP",  
  minGSSize = 10,  
  maxGSSize = 500,  
  pvalueCutoff = 0.05,  
  pAdjustMethod = "BH"  
)

head(ora_go_up@result$geneID)

dotplot(ora_go_up, showCategory= 10,  x = "Count")

ora_go_up_filtered <- ora_go_up@result %>%
  filter(Description == "chemotaxis")%>%
  pull(geneID)

ora_go_up_geneids <- strsplit(ora_go_up_filtered, "/")
ora_go_up_immune <- c("79931", "4057",  "931" ,  "3500",  "2633" , "54900" ,"973",   "971",   "3537",
  "1380" ,"930",   "1236",  "3507",  "952",   "28639", "3495")

ora_go_up_chemotaxis <- c("6373" , "4283" , "6372"  ,"10563", "3627",  "6362" , "6363",  
                          "7852",  "7474" , "3689" ,
                           "1236" , "11151" ,"5880" , "7130",  "64218" ,"10507", "8600" , "7941" )

gene_symbols_immune <- mapIds(org.Hs.eg.db, 
                       keys = ora_go_up_immune, 
                       column = "SYMBOL", 
                       keytype = "ENTREZID", 
                       multiVals = "first")

gene_symbols_overlap <- intersect(gene_symbols_immune, gene_symbols_chemotaxis)
gene_symbols_overlap
writeClipboard(intersect_DEGs_df_up$id)


# downregulation of metabolism pathways 
# upregulation of immune response, cell movement 
