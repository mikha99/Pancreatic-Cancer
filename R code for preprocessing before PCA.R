#### loading library

### Setting working directory
### write dor table is for saving in the workspace 

library(tidyverse)
library(dplyr)

setwd("~/Desktop/PINE2")

## loading tumor data

tumor_data<- read.table("TCGA-PAAD.htseq_counts modified del version.txt", header=TRUE , sep="\t")

## verity tumor data
tumor_data [1:10 , 1:10 ]
## loading normal data

normal_data <- read.table("gene_reads_pancreas_control_modified.txt" , header=TRUE , sep="\t")
### veryfying the normal data
normal_data[1:10,1:5]
### we need to merge the data
merged_data <-inner_join(tumor_data, normal_data, by=c("Ensembl_ID" = "Name"))
### veryfying merged data
merged_data[1:5,1:5]
### saving merged data in wordking directory
write.table(merged_data, file="tcgatumorgtexcontrolcount.txt", sep="\t", quote = FALSE,  row.names = FALSE, col.names = TRUE)


### count the tumor and normal samples
names(merged_data)
#### annotation for running PCA
### make annotation row for tumor and normal samples
Group<-c("Group",rep("Tumor",182),rep("Normal",328))
annotated_merged_data<- rbind(Group, merged_data)
#### save the merged data

write.table(annotated_merged_data,file="tcgatumorgtexcontrolcountannotated.txt",sep="\t", quote = FALSE,  row.names = FALSE, col.names = TRUE)




           