# PIPELINE FOR TRANSCRIPTOMIC DATASET ANALYSIS ISSUED FROM MULTISEP PRE-TREATMENT
# DESCRIPTION OF DATASET ORIGIN:
# A precursor dataset was obtained from: 
# The MMRF (Multiple Myeloma Research Foundation) 
# CoMMpass (Clinical Outcomes in Multiple Myeloma to Personal Assessment) data resource (https://themmrf.org/; Shaw, 2021; Settino, 2021), 
# IA21 data section, “Expression Estimates – Gene Based” rubric: ‘MMRF_CoMMpass_IA21_salmon_geneUnstranded_tpm.tsv’. 
# It was processed using the MultiSEp package in R. (Wappett et al., 2021; https://www.overton-lab.uk/app/synlegg/#tuid)
# Using the MultiSEp algorithm, RNA-seq data can be statistically processed, decomposed, and restructured to identify cell clusters formed by different gene expressions. 
# (MultiSEp Analysis was run in the Ian Overton Laboratory by PhD student Adeline McKie)
#===============================================================
# Upload Libraries.
#===============================================================
library(dplyr)
library(tidyverse)
library(tidyr)
library(stringr)

#================================================================
# Uploading Dataset.
#================================================================
# Initial dataset was saved in `pattern` variable, the cleaned original dataset was saved in `pattern_new` variable.
pattern <- read.table (file="../gdr/Combined_Binomial63million_padjusted_qvalLt0.05.txt", head = TRUE)

# Removing rows containing NA values.
pattern_new <- na.omit(pattern)

#=================================================================
# Preliminary Dataset Analysis.
#=================================================================
# Explore categories and total number of table cells (Cluster Combination) types in different types of tables (Total Cluster Combinations).
categories <- unique(pattern_new$Cluster_Combination)
numberOfCategories <- length(categories)

# Explore categories and total number of Total Cluster Combinations types.
categories_1 <- unique(pattern_new$TotalCluster_Combinations) #number of table types
numberOfCategories_1 <- length(categories_1)

# Explore counts of Cluster Combination types for each Total Cluster Combinations table type.
# Example of code for cluster combination count “by hand”.
pattern_2_2_4_1x2 <- pattern_new %>% filter(Gene1_ClusterCount == 2, Gene2_ClusterCount == 2, TotalCluster_Combinations == 4, Cluster_Combination == "1_x_2", New_q_value <= 0.05)

# Code for automatization of Cluster Combination type count for each Total Cluster Combinations table type.
cl_comb <- pattern_new %>% group_by(Gene1_ClusterCount, Gene2_ClusterCount, TotalCluster_Combinations, Cluster_Combination) %>% summarise(number = n())

#====================================================================================================================================
# Standartization of Total Cluster Combinations table types. 
#=====================================================================================================================================
# Transformation of the initial dataset into a new dataset where all Total Cluster Combination table types (n x m) have n ≤ m.
# Example: Table types with form 3x2, 4x2, 4x3, 5x2, 5x4 transformed accordingly into 2x3, 2x4, 3x4, 2x5, 4x5. 
# Result: The standardization step reduce the number of table types and make analysis more efficient.
# Method: Because data was swapped in all columns reflecting gene pair co-expression levels (Gene1, Gene2, Gene1_ClusterCount, Gene2_ClusterCount, Cluster_Combination), new dataset preserve this information correctly.

# STEP 1. **Setting the Condition**:
# Create in the dataset an additional column with TRUE/FALSE conditions which then indicate where data from Gene1/Gene2, Gene1_ClusterCount/Gene2_ClusterCount, first_cluster_id/second_cluster_id columns should be swapped.

pattern_new$condition <- pattern_new$Gene1_ClusterCount > pattern_new$Gene2_ClusterCount

# STEP 2. **Splitting the Cluster_Combination Coordinates**: create two separate indices so that they could be treated individually.
# Separate Cluster_Combination column in two columns, `first_cluster_id` and `second_cluster_id`, to swap the cluster ids in cluster combinations.

pattern_new <- pattern_new %>% separate_wider_delim(Cluster_Combination, delim = "_x_", names = c("first_cluster_id", "second_cluster_id"), cols_remove = FALSE)

# STEP 3. **Swap the Variables**: Data for Genes, Gene_ClusterCount and Cluster_Combination were swapped between the columns according the estimated condition.
# Data Integrity Protection: 
# To bypass the simultaneous variable overwriting problem inherent to sequential column transposition, 
# the pipeline temporarily routes data through explicit staging variables (`gene1_new`, `gene1_cc_new`, `ClusterComb_new`). 
# This layout prevents memory collisions and structural data loss during the matrix transposition loop.
# Generation of New Columns:
# For storage of information from the `Gene2` column, new column `gene1_new` is created.
# For storage of information from the `Gene2_ClusterCount` column, new column `gene1_cc_new column` is created. 
# For storage of indices from `Cluster_Combination` column, new column `ClusterComb_new` is created to switch the data between `first_cluster_id` and `second_cluster_id` columns.
# Swap the Data Based on Condition:
# If the condition is TRUE, the value of `Gene2` column is copied into a new temporary column `gene1_new`. 
# If the condition is TRUE, the value of `Gene2_ClusterCount` column is copied into a new temporary column `gene1_cc_new`. 
# For the `Cluster_Combination` column indices, if the condition is TRUE, `second_cluster_id` is copied to the `ClusterComb_new` column. 
# Create a new variable for the transformed initial dataset with new columns.

sw_pattern <- pattern_new %>% mutate(gene1_new = if_else(condition == TRUE, Gene2, Gene1), 
                                      gene1_cc_new = if_else(condition == TRUE, Gene2_ClusterCount, Gene1_ClusterCount), 
                                      Gene2 = if_else(condition == TRUE, Gene1, Gene2), 
                                      Gene2_ClusterCount = if_else(condition == TRUE, Gene1_ClusterCount, Gene2_ClusterCount), 
                                      ClusterComb_new = if_else(condition == TRUE, second_cluster_id, first_cluster_id), 
                                      second_cluster_id = if_else(condition == TRUE, first_cluster_id, second_cluster_id))

# STEP 4. **Cleaning and Renaming Columns:** Assign names to newly created columns containing swapped information from `Gene2`, `Gene2_ClusterCount` and `second_cluster_id` columns and delete old columns.

sw_pattern <- sw_pattern %>% select(-Gene1, Gene1 = gene1_new, -Gene1_ClusterCount, Gene1_ClusterCount = gene1_cc_new, -Cluster_Combination, -condition, -first_cluster_id, first_cluster_id = ClusterComb_new)

# STEP 5. Unite two columns with Cluster Combination indices in one column Cluster_Combination as it was in the initial dataset. 

sw_pattern <- sw_pattern %>% unite("Cluster_Combination", c(first_cluster_id, second_cluster_id), sep = "_x_")

# STEP 6. **Reorganization of Column Positions:** Reorder column positions according to the initial dataset order.

sw_pattern <- sw_pattern %>% relocate(Gene1)
sw_pattern <- sw_pattern %>% relocate(Gene1_ClusterCount, .after = Gene2)
sw_pattern <- sw_pattern %>% relocate(Cluster_Combination, .after = TotalCluster_Combinations)

# Save new dataset in a file. 
write.csv(sw_pattern, row.names=TRUE, file="/home/svetlana/switch/Switched_Pattern.txt")

# Explore counts of Cluster Combination types for each Total Cluster Combinations table type in new dataset with standardized Total Cluster Combinations table types.
# Automatized version. 
scf_count <- sw_pattern %>% group_by(Gene1_ClusterCount, Gene2_ClusterCount, TotalCluster_Combinations) %>% summarise(count=n()) %>% arrange(TotalCluster_Combinations)

#==================================================================================================================================================================================
# Generating Combination Dataset (Table with Unique Gene Pairs' Co-Expression Table Characteristics, or Co-Expression Patterns)
#==================================================================================================================================================================================
# Generation of the Pattern table with combinations of Cluster Combinations corresponding to each gene pair.
# Each gene pair can have only one type of combination of Cluster Combinations in only one table type. 
# Combination of Cluster Combinations can contain single or multiple Cluster Combinations (e.g. single: 1_x_1, multiple: 1_x_1.2_x_1.2_x_2).
# These patterns are the unified transcriptomic signature inventory individual for a certain set of gene pairs.

# Creating combinations table
Comb <- sw_pattern %>% select(Gene1, Gene2, Gene1_ClusterCount, Gene2_ClusterCount, TotalCluster_Combinations, Cluster_Combination) %>% distinct() 
                %>% group_by(Gene1, Gene2, Gene1_ClusterCount, Gene2_ClusterCount, TotalCluster_Combinations) 
                %>% arrange(Cluster_Combination) 
                %>% summarise(Result_Comb = str_c(Cluster_Combination, collapse = ".")) 
                %>% arrange(nchar(Result_Comb), Result_Comb) 
                %>% ungroup()

# Creating combination tables for each table type by filtering Combination Dataset by TotalCluster_Combination type (example for 4-cell table).
comb4 <- Comb %>% filter(TotalCluster_Combinations == 4) %>% arrange(nchar(Result_Comb))

# Combination count for each table type (example for 4-cell table). 
# Shows both existing combinations and their number.
comb4_n <- comb4 %>% group_by(Result_Comb) %>% summarise(count=n()) %>% arrange(desc(count))
length(comb4_n)

#============================================================================================================================================================
# Creating gene pairs and single gene list
#============================================================================================================================================================
# Based on the combination count, the most numerous combinations for each table type were chosen to create gene lists for subsequent enrichment analysis.
genes_pairs_4 <- comb4 %>% filter(Result_Comb == "2_x_1")
length(genes_pairs_4$Result_Comb)
unique_gp_4 <- genes_pairs_4 %>% distinct(Gene1, Gene2)
length(unique_gp_4$Gene1)
length(unique_gp_4$Gene2)
gene_list_4 <- unique_gp_4 %>% pivot_longer(c(Gene1, Gene2), values_to = "Genes") %>% select(-name) %>% distinct(Genes) %>% arrange(Genes)








