# Preliminary Dataset Analysis.
> categories <- unique(pattern$Cluster_Combination) #number of table cell positions types in different types of tables
> numberOfCategories <- length(categories)
> categories
> numberOfCategories
> categories_1 <- unique(pattern$TotalCluster_Combinations) #number of table types
> numberOfCategories_1 <- length(categories_1)
> categories_1
> numberOfCategories_1
> print(sort(categories_1))
> print(sort(categories)) #more simple visualisation of the cell table type list

> scf_count <- sw_pattern %>% group_by(Gene1_ClusterCount, Gene2_ClusterCount, TotalCluster_Combinations) %>% summarise(count=n()) %>% arrange(TotalCluster_Combinations)
> scf_count %>% print(width=Inf)

# Cluster Combination count.
# Example of code for cluster combination count “by hand”.
> pattern_2_2_4_1x2 <- pattern %>% filter(Gene1_ClusterCount == 2, Gene2_ClusterCount == 2, TotalCluster_Combinations == 4, Cluster_Combination == "1_x_2", New_q_value <= 0.05)> print(nrow(pattern_2_2_4_1x2))

# Code for automatization.
> cl_comb <- pattern_new %>% group_by(Gene1_ClusterCount, Gene2_ClusterCount, TotalCluster_Combinations, Cluster_Combination) %>% summarise(number = n())

# Table type standardization
> pattern <- read.table (file="../gdr/Combined_Binomial63million_padjusted_qvalLt0.05.txt", head = TRUE)
> pattern_new <- na.omit(pattern)
> library(dplyr)
> library(tidyr)
> library(stringr)
> pattern_new <- pattern_new %>% separate_wider_delim(Cluster_Combination, delim = "_x_", names = c("first_cluster_id", "second_cluster_id"), cols_remove = FALSE)
> pattern_new <- pattern_new %>% mutate(gene1_new = if_else(condition == TRUE, Gene2, Gene1), gene1_cc_new = if_else(condition == TRUE, Gene2_ClusterCount, Gene1_ClusterCount), Gene2 = if_else(condition == TRUE, Gene1, Gene2), Gene2_ClusterCount = if_else(condition == TRUE, Gene1_ClusterCount, Gene2_ClusterCount), ClusterComb_new = if_else(condition == TRUE, second_cluster_id, first_cluster_id), second_cluster_id = if_else(condition == TRUE, first_cluster_id, second_cluster_id))
> pattern_new <- pattern_new %>% select(-Gene1, Gene1 = gene1_new, -Gene1_ClusterCount, Gene1_ClusterCount = gene1_cc_new, -Cluster_Combination, -condition, -first_cluster_id, first_cluster_id = ClusterComb_new)
> pattern_new <- pattern_new %>% unite("Cluster_Combination", c(first_cluster_id, second_cluster_id), sep = "_x_")
> pattern_new <- pattern_new %>% relocate(Gene1)
> pattern_new <- pattern_new %>% relocate(Gene1_ClusterCount, .after = Gene2)
> pattern_new <- pattern_new %>% relocate(Cluster_Combination, .after = TotalCluster_Combinations)
> Switched_Pattern <- write.csv(pattern_new, row.names=TRUE, file="/home/svetlana/switch/Switched_Pattern.txt")
> pattern_new %>% print(width = Inf)

# Creating combinations table
> Comb <- sw_pat %>% select(Gene1, Gene2, Gene1_ClusterCount, Gene2_ClusterCount, TotalCluster_Combinations, Cluster_Combination) %>% distinct() %>% group_by(Gene1, Gene2, Gene1_ClusterCount, Gene2_ClusterCount, TotalCluster_Combinations) %>% arrange(Cluster_Combination) %>% summarise(Result_Comb = str_c(Cluster_Combination, collapse = ".")) %>% arrange(nchar(Result_Comb), Result_Comb) %>% ungroup()
> Comb %>% print(width=Inf, 1:5,)

# Creating combination tables for each table type (example for 4-cell table).
> comb4 <- Comb %>% filter(TotalCluster_Combinations == 4) %>% arrange(nchar(Result_Comb))

# Combination count for each table type (example for 4-cell table).
> comb4_n <- comb4 %>% group_by(Result_Comb) %>% summarise(count=n()) %>% arrange(desc(count))
> length(comb4_n)
> comb4_n

# Creating gene pairs and single gene list
> genes_pairs_4 <- comb4 %>% filter(Result_Comb == "2_x_1")
> genes_pairs_4[1:3,]
> length(genes_pairs_4$Result_Comb)
> unique_gp_4 <- genes_pairs_4 %>% distinct(Gene1, Gene2)
> unique_gp_4[1:3,]
> length(unique_gp_4$Gene1)
> length(unique_gp_4$Gene2)
> gene_list_4 <- unique_gp_4 %>% pivot_longer(c(Gene1, Gene2), values_to = "Genes") %>% select(-name) %>% distinct(Genes) %>% arrange(Genes)
> head(gene_list_4)
> gene_list_4 %>% print(width=Inf, 1:3)






