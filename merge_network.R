df_m=read.table(file="p_net-edges_4448min.txt")
dim(df_m)
df_m=df_m[,-(4:6)]
colnames(df_m)=df_m[1,]
df_m=df_m[-1,]

#df1=read.csv(file = "p_network_4448.csv")
df2=read.csv(file="STRING_network_ID4_4448.csv")


df2_filtered <- df2[df2$fromNode %in% gene_names & df2$toNode %in% gene_names, ]

df2_dup <- df2[!(df2$fromNode %in% gene_names & df2$toNode %in% gene_names), ]
#df2_dup1 <- df2[!(df2$fromNode %in% gene_names), ]
#df2_dup2 <- df2[!(df2$toNode %in% gene_names), ]

length(unique(df2_filtered$fromNode))
length(unique(df2_filtered$toNode))


library(dplyr)
###########The two data frames df1 and df2 are merged according to their common "fromNode" and "toNode" columns, and the weights are calculated and processed after merging.
merged_df <- merge(df_m, df2_filtered, by.x = c("fromNode", "toNode"), by.y = c("fromNode", "toNode"), all.x = TRUE) %>%
  mutate(weight.x = as.numeric(weight.x),
         weight.y = as.numeric(weight.y),
         weight = ifelse(!is.na(weight.y), weight.x + weight.y, weight.x))#df1,df2 on the contrary 


########
merged_1=merged_df[!is.na(merged_df$weight.y),]###on the contrary

##############
merged_df2 <- merge(df_m, df2_filtered, by.x = c("fromNode", "toNode"), by.y = c("toNode", "fromNode"), all.x = TRUE) %>%
  mutate(weight.x = as.numeric(weight.x),
         weight.y = as.numeric(weight.y),
         weight = ifelse(!is.na(weight.y), weight.x + weight.y, weight.x))

merged_df3 <- merge(df2_filtered, df_m, by.x = c("fromNode", "toNode"), by.y = c("toNode", "fromNode"), all.x = TRUE) %>%
  mutate(weight.x = as.numeric(weight.x),
         weight.y = as.numeric(weight.y),
         weight = ifelse(!is.na(weight.y), weight.x + weight.y, weight.x)) 


merged_1_2=merged_df2[!is.na(merged_df2$weight.y),]##df1，negative is the same
merged_2=merged_df3[!is.na(merged_df3$weight.y),]##df2, negative is the same
#######
merged_df_rows<- anti_join(merged_df, merged_1_2, by = c("fromNode", "toNode"))#In df1, positive is equally correlated, minus negative is equally correlated
########
merged_all <- bind_rows(merged_1, merged_2)#df2, same crossing
remaining_rows_df2 <- anti_join(df2_filtered, merged_all, by = c("fromNode", "toNode"))#df2 minus cross identical
########
merged_df_rows <- bind_rows(merged_df_rows, remaining_rows_df2)#

#write.csv(merged_df_rows,file="p_merged_df_rows.csv")
#write.csv(merged_1_2,file="p_merged_df_add.csv")
#################
rm(ls=merged)

merged_df_rows=merged_df_rows[,-(3:4)]
merged_1_2=merged_1_2[,-(3:4)]
all_net=rbind(merged_df_rows,merged_1_2)

all_net$weight=all_net$weight/2
write.csv(all_net,file="p_net_all_4448min.csv",row.names = F)
