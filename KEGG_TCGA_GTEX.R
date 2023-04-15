if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")

BiocManager::install("gage")
library(gage)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("gageData")
library(gageData)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pathview")
library(pathview)
#Load DESeq2 results file as input file
result <- read.table("/Users/animikha/Desktop/PINE2/DESeq2_TCGAGTEX_2.txt", sep="\t", header = TRUE, row.names=1)
#Check top 10 lines in results 
head(result)
#list of KEGG pathways genes
data(kegg.sets.hs)
#Index of numbers of pathway
data(sigmet.idx.hs)
#Extract cleaner gene set of important pathways
kegg.sets.hs <- kegg.sets.hs[sigmet.idx.hs]


BiocManager::install("clusterProfiler")
library(clusterProfiler)
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(rownames(result), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)

# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
res_dedup = result[dedup_ids$ENSEMBL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
res_dedup$ENTREZID = dedup_ids$ENTREZID

#Extract the fold changes from the results
foldchanges <- res_dedup$log2FoldChange

# Name fold change vector with ENTREZ ids
names(foldchanges) <- res_dedup$ENTREZID

# omit any NA values 
foldchanges<-na.omit(foldchanges)

# sort the list in decreasing order (required for clusterProfiler)
foldchanges = sort(foldchanges, decreasing = TRUE)

#check top 5 values
head(foldchanges)

#Check top 10 values
head(foldchanges)
#Get all KEGG pathways associated with genes
kegg_result <- gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
#Write results of all pathways in a file
write.table(kegg_result, file="all_pathways_TCGAGTEX.txt", sep="\t", row.names=FALSE,quote=FALSE)
#Check attributes
attributes(kegg_result)
#Check upregulated pathways
head(kegg_result$greater)
#check down regulated pathways
head(kegg_result$less)
#Extract all upregulated pathways data
Upreg_pathways <- data.frame(Pathway_id=rownames(kegg_result$greater), kegg_result$greater)
#write all upregulated pathways into a file
write.table(Upreg_pathways, file="upregulated_pathways.txt", sep="\t", row.names=F, quote=FALSE)
#Extract all downregulated pathways data
Downreg_pathways <- data.frame(Pathway_id=rownames(kegg_result$less), kegg_result$less)
#write all downregulated pathways into a file
write.table(Downreg_pathways, file="Downregulated_pathways.txt", sep="\t", row.names=F, quote=FALSE)
#Visualize one upregulated pathway
pathview(gene.data=foldchanges, pathway.id="hsa04514")
#Obtain a different PDF based output of the same data
pathview(gene.data=foldchanges, pathway.id="hsa04514", kegg.native=FALSE)
#Extract the top 5 upregulated pathways 
keggrespathways_up <- rownames(kegg_result$greater)[1:5]
#Extract the IDs part of each string
keggresids_up <- substr(keggrespathways_up, start=1, stop=8)
head(keggresids_up)
#write Ids in a file
write.table(keggresids_up, file="top5_upreg_pathways_ids.txt", sep="\t", row.names=F, quote=FALSE)
#Draw plots for top 5 upregulated pathways
pathview(gene.data=foldchanges, pathway.id=keggresids_up, species="hsa")
#Extract top 5 down regulated pathways 
keggrespathways_down <- rownames(kegg_result$less)[1:5]
#Extract the IDs 
keggresids_down = substr(keggrespathways_down, start=1, stop=8)
head(keggresids_down)
#write Ids in a file
write.table(keggresids_down, file="top5_downreg_pathways_ids.txt", sep="\t", row.names=F, quote=FALSE)
#Draw plots for top 5 downregulated pathways
pathview(gene.data=foldchanges, pathway.id=keggresids_down, species="hsa")
