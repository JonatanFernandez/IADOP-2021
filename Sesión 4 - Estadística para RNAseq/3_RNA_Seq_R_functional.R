library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(ggplot2)
library(cowplot)

data(geneList)

head(geneList)

gene <- names(geneList)[abs(geneList) > 2]

GO_overrepresented <- enrichGO(gene          = gene,
              universe      = names(geneList),
              OrgDb         = org.Hs.eg.db,
              ont           = "MF",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.01,
              qvalueCutoff  = 0.05,
              readable      = TRUE)

GO_GSEA <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "MF",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

p1 <- dotplot(GO_overrepresented, showCategory=30) + ggtitle("Dotplot for GO enrichment")
p2 <- dotplot(GO_GSEA, showCategory=30) + ggtitle("Dotplot for GSEA")
plot_grid(p1, p2, ncol=2)

