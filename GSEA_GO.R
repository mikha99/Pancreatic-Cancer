organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
library(clusterProfiler)
#Load DESeq2 results file as input file
result <- read.table("/Users/animikha/Desktop/PINE2/DESeq2_TCGAGTEX_2.txt", sep="\t", header = TRUE, row.names=1)
#Check top 10 lines in results 
head(result)
# Convert gene IDs for gseGO function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(rownames(result), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb=organism)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
res_dedup = result[dedup_ids$ENSEMBL,]
# Create a new column in df2 with the corresponding ENTREZ IDs
res_dedup$SYMBOL = dedup_ids$SYMBOL
#Extract the fold changes from the results
foldchanges <- res_dedup$log2FoldChange

# Name fold change vector with ENTREZ ids
names(foldchanges) <- res_dedup$SYMBOL


# omit any NA values 
foldchanges<-na.omit(foldchanges)
# sort the list in decreasing order (required for clusterProfiler)
foldchanges = sort(foldchanges, decreasing = TRUE)

#check top 5 values
head(foldchanges)

#run GSEA for GO respository
gse <- gseGO(geneList=foldchanges, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")
# Save GSEA GO result
write.csv(gse,"TCGA_GTEX_GSEA.csv")
install.packages("repr")
library(repr)
# Change plot size to 8 inches wide
options(repr.plot.width=10)
require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
install.packages("ggnewscale")
# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(gse,font.size=4, categorySize="geneNum", foldChange=foldchanges,max.overlaps=50,node_label="category")
# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(gse,font.size=4, categorySize="geneNum", foldChange=foldchanges,max.overlaps=50,node_label="gene")
###upsetplot

install.packages("ggupset")
library(ggupset)
library(ggplot2)
enrichplot::upsetplot(gse, n = 10) + theme(axis.text.y = element_text(size = 9))
#ridgeplot
install.packages("ggridges")
library(ggridges)
ridgeplot(gse) + labs(x = "Enrichment distribution") + 
  theme(text = element_text(size = 8),
        axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 8))

## heatplot
heatplot(gse, showCategory = 10, foldChange=foldchanges) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=9)) + 
  theme(legend.position="bottom")
###
