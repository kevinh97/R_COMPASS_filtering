library("tidyverse") 
library("mergeutils")

combined_tbl <- read.csv(file = '/Users/chanfreaulab/Documents/combined_tbl.csv')

rep_names <- levels(combined_tbl$sample_name)

for(i in rep_names) {
print(i)
assign(i, subset(x = combined_tbl, combined_tbl$sample_name == i))
}

#Trying out a new way to merge rows
#kh_merge_test <- combined_tbl %>% 
 # group_by(intron_ID) %>% 
#  arrange_all(.by_group = TRUE)
#  summarise_all(funs(select(-intron_ID, everything())))

RRP6_1_YPD_truncated <- RRP6_1_YPD_[,c(-3:-7)]

kh_merge_test <- Reduce(function(x,y) merge(x,y,by="intron_ID",all=TRUE) ,list(RRP6_1_YPD_, RRP6_2YPD_, RRP6_3YPD_))
kh_merge_test1 <- Reduce(function(x,y) merge(x,y,by="intron_ID",all=TRUE) ,list(UPF1_1_YPD_, UPF1_2YPD_, UPF1_3YPD_))
kh_merge_test2 <- Reduce(function(x,y) merge(x,y,by="intron_ID",all=TRUE) ,list(U_R_1YPD_, U_R_2YPD_, U_R_3YPD_))
kh_merge_test_final <- Reduce(function(x,y) merge(x,y,by="intron_ID",all=TRUE) ,list(kh_merge_test, kh_merge_test1))
kh_merge_test_final_trim <- kh_merge_test_final %>% 
  group_by(intron_ID) %>% 
  summarise_all(funs(trimws(paste(., collapse = ''))))
kh_merge_test_final <- Reduce(function(x,y) merge(x,y,by="intron_ID",all=TRUE) ,list(kh_merge_test_final_trim, kh_merge_test2))

df_list <- list(RRP6_1_YPD_, RRP6_2YPD_, RRP6_3YPD_, U_R_1YPD_, U_R_2YPD_, U_R_3YPD_)

#Remove unnecessary columns

drop_test <- c("X.x.x","seqnames.x.x","start.x.x", "end.x.x", "width.x.x", "strand.x.x", 
               "annotated_junction.x.x", "chrom.x.x", "RNA_strand.x.x", "ann_5SS.x.x", "ann_3SS.x.x", "canonical_5SS.x.x",
               "intron_size.x.x", "poly_U_count.x.x", "poly_Y_count.x.x", "X.x.y", "seqnames.x.y", "start.x.y", "end.x.y",
               "width.x.y", "strand.x.y", "annotated_junction.x.y", "chrom.x.y", "RNA_strand.x.y", "ann_5SS.x.y", "ann_3SS.x.y",
               "canonical_5SS.x.y", "intron_size.x.y", "poly_U_count.x.y", "poly_Y_count.x.y", "US_5SS_1.x.x", "US_5SS_2.x.x", "US_5SS_3.x.x",
               "DS_3SS_1.x.x", "DS_3SS_2.x.x", "DS_3SS_3.x.x", "US_5SS_1.x.y", "US_5SS_2.x.y", "US_5SS_3.x.y", "DS_3SS_1.x.y",
               "DS_3SS_2.x.y", "DS_3SS_3.x.y", "US_5SS_10nt.x.x", remove_column$)
remove_column <- kh_merge_test_final[,!(names(kh_merge_test_final) %in% drop_test)]
remove_column <- remove_column[, c(-3:-29)]
#, U_R_1YPD_, U_R_2YPD_, U_R_3YPD_, UPF1_1_YPD_, UPF1_2YPD_, UPF1_3YPD_

test_df <- multimerge(df_list, by = "intron_ID", all = TRUE ,suffixes = rep_names)

#THIS IS OK BUT NEEDS TO BE FIXED IT ONLY DOES 3 DATAFRAMES SO FAR
merged <- merge(RRP6_1_YPD_, RRP6_2YPD_, by = "intron_ID", suffixes = c("RRP6_1_YPD_", "RRP6_2YPD_"), all.x = TRUE, all.y= TRUE)
merged <- merge(merged, U_R_1YPD_,by = "intron_ID", all.x = TRUE, all.y= TRUE)

#THESE 2 LINES MAY NOT BE NEEDED 
merged$FAnS <- as.numeric(merged$FAnS)
merged$FAnSRRP6_1_YPD_ <- as.numeric(merged$FAnSRRP6_1_YPD_)

merged[is.na(merged)] <- 0 #THIS IS BAD DO NOT DO ONLY DO THIS FOR COLUMNS WHERE NA ACTUALLY MEANS 0

upf1rrp6_rrp6_plot <- ggplot(data = merged)+
  geom_point(aes(x = FAnSRRP6_1_YPD_, y = FAnS, color = COMPASS_counts))+
  ggtitle(label = "rrp6∆ vs rrp6∆upf1∆")+
  geom_abline(slope=1, intercept=0) + 
  xlab("RRP6∆ FAnS") +
  ylab("UPF1∆RRP6∆ FAnS")+
  scale_x_log10()+
  scale_y_log10()+
  scale_color_gradient(low = "yellow3", high = "green4", trans = "log2")

print(upf1rrp6_rrp6_plot)

