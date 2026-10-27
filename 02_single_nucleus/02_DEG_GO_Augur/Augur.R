library(Seurat)
library(Augur)
library(dplyr)
library(ggplot2)
library(tidyr)

# ==========================
# 1. Add cell-type annotations
# ==========================
data.all <- readRDS('snRNA.rds')
new.cluster.ids <- c(
  'Pa','VT','Pa','Pa','Ov','PC','Pa','Pa','PC','Pa',
  'Pa','Pa','Pa','Ep','Ep','PC','VT','RT','Ov','RT',
  'Pa','VT','Ov','Ov','VT','VT','PC','VT','VT','VT','VT'
)
names(new.cluster.ids) <- levels(data.all)

data.all <- RenameIdents(data.all, new.cluster.ids)
data.all$celltype <- Idents(data.all)

# ==========================
# 2. Create MO versus C+P labels
# ==========================
data.all$label_MO_vs_CP <- ifelse(data.all$orig.ident == "MO", "MO",
                                  ifelse(data.all$orig.ident %in% c("C","P"), "C_P", NA))

# Exclude U, which is labelled NA
data.sub <- subset(data.all, !is.na(label_MO_vs_CP))

table(data.sub$label_MO_vs_CP)
# Cell counts for the MO and C_P groups

# ==========================
# 3. Run Augur
# ==========================
augur <- calculate_auc(data.sub,
                       label_col = "label_MO_vs_CP",
                       cell_type_col = "celltype",
                       n_threads = 70)

head(augur$AUC)

# ==========================
# 4. Plot AUC values
# ==========================
plot_df <- augur$AUC %>%
  arrange(desc(auc)) %>%
  mutate(cell_type = factor(cell_type, levels = unique(cell_type)))

p <- ggplot(plot_df, aes(x = cell_type, y = auc, fill = cell_type)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = round(auc,2)), vjust = -0.5, size = 3) +
  scale_fill_brewer(palette = 'Set3') +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Cell Type") + ylab("AUC") +
  ggtitle("Augur: MO vs C+P")

print(p)

# Save results
ggsave("Augur_MO_vs_CP.pdf", p, width = 7, height = 5)
save.image(file="Augur_MO_vs_CP.Rdata")
