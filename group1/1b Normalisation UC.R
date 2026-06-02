library("tidyverse")
library("dplyr")
library("ggplot2")
library("ggpubr")
library("tidyr")
library("readr")
library("readxl")
library("limma")
library("sva")



# Download Data -----------------------------------------------------------
UC_raw <- read.csv("Data/Ulcerative Colitis dataset before normalization.normalization..csv")

# Batch correct -----------------------------------------------------------

# Set gene names as rownames
rownames(UC_raw) <- UC_raw$geneNames

# Remove geneNames column
edata <- UC_raw[, -1]

# Convert to matrix
edata <- as.matrix(edata)

sample_names <- colnames(edata)

# Condition
condition <- ifelse(grepl("Control", sample_names), "Control", "Treat")

# Batch (dataset)
batch <- sub("_(GSM.*)", "", sample_names)   # extracts GSE IDs

# Convert to factors
condition <- factor(condition)
batch <- factor(batch)

# Build matrices
mod  <- model.matrix(~ condition)
mod0 <- model.matrix(~ 1)


# ComBat correction
edata_combat <- ComBat(dat = edata,
                       batch = batch,
                       mod = mod)

# Then run SVA
n.sv <- num.sv(edata_combat, mod, method = "leek")


# PCA

get_pca_df <- function(mat, batch, condition) {
  pca <- prcomp(t(mat), scale. = TRUE)
  
  df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    batch = batch,
    condition = condition
  )
  
  return(df)
}

# Before correction
pca_before <- get_pca_df(edata, batch, condition)

# After correction
pca_after  <- get_pca_df(edata_combat, batch, condition)

library(ggplot2)
library(gridExtra)

# Colours for batches
batch_colors <- c(
  "GSE36807" = "#E41A1C",  # red
  "GSE87473" = "#377EB8",  # blue
  "GSE92415" = "#FF7F00"   # orange
)

# Shapes for condition
shape_values <- c(
  "Control" = 16,  # circle
  "Treat"   = 17   # triangle
)

p1 <- ggplot(pca_before, aes(PC1, PC2, color = batch, shape = condition)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = batch_colors) +
  scale_shape_manual(values = shape_values) +
  ggtitle("Before batch correction") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(colour = "black", fill = NA)
  )

p2 <- ggplot(pca_after, aes(PC1, PC2, color = batch, shape = condition)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = batch_colors) +
  scale_shape_manual(values = shape_values) +
  ggtitle("After batch correction") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(colour = "black", fill = NA)
  )

combined_plot <- gridExtra::grid.arrange(p1, p2, ncol = 2)


# Save last plot
ggsave("Figures/UC_PCA_batch_correction.png",
       plot = combined_plot,
       width = 15,
       height = 5,
       dpi = 300)


# Save batch corrected csv file

# Convert to data frame
edata_out <- as.data.frame(edata_combat)

# Add gene names as a column
edata_out$geneNames <- rownames(edata_out)

# Move geneNames to first column
edata_out <- edata_out[, c("geneNames", colnames(edata_combat))]

write.csv(edata_out,
          file = "Files/UC_batch_corrected.csv",
          row.names = FALSE)
