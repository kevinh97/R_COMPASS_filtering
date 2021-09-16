junction_df_summarised$seqnames <- paste("chr", junction_df_summarised$seqnames)
junction_df_summarised$seqnames <- str_replace_all(junction_df_summarised$seqnames, fixed(" "), "")
junction_df_summarised$intron_ID <- gsub(" ", "_", junction_df_summarised$intron_ID)
junction_df_summarised$intron_ID <- paste("chr", junction_df_summarised$intron_ID)
junction_df_summarised$intron_ID <- str_replace_all(junction_df_summarised$intron_ID, fixed(" "), "")
previous_SJ_specific_study$intron_ID <- gsub(" ", "_", previous_SJ_specific_study$intron_ID)
previous_SJ$intron_ID <- gsub(" ", "_", previous_SJ$intron_ID)
junction_df_filtered$intron_ID <- gsub(" ", "_", junction_df_filtered$intron_ID)
junction_df_filtered$intron_ID <- paste("chr", junction_df_filtered$intron_ID)
junction_df_filtered$intron_ID <- str_replace_all(junction_df_filtered$intron_ID, fixed(" "), "")
ann_introns_overlap_df$intron_ID <- gsub(" ", "_", ann_introns_overlap_df$intron_ID)
ann_introns_overlap_df$intron_ID <- paste("chr", ann_introns_overlap_df$intron_ID)
ann_introns_overlap_df$intron_ID <- str_replace_all(ann_introns_overlap_df$intron_ID, fixed(" "), "")
comparison <- compare(ann_introns_overlap_df$intron_ID,previous_SJ_specific_study$intron_ID,allowAll=TRUE)
as.data.frame(comparison)
df <- data.frame(matrix(unlist(comparison), byrow=TRUE),stringsAsFactors=FALSE)
new_sj <- subset(ann_introns_overlap_df$intron_ID, ann_introns_overlap_df$intron_ID == previous_SJ$intron_ID)
new_sj <- as.data.frame(new_sj)
previous_SJ$same <- ifelse(previous_SJ$intron_ID %in% ann_introns_overlap_df$intron_ID,"true","false")
ann_introns_overlap_df$blah <- ifelse(ann_introns_overlap_df$intron_ID %in% previous_SJ$intron_ID,"true","false")
z <- subset(ann_introns_overlap_df, ann_introns_overlap_df$annotated_junction == "FALSE")
new_sj <- subset(z, z$same == "false")
previous_SJ_with_COMPASS$intron_ID <- paste(previous_SJ_with_COMPASS$intron_ID, previous_SJ_with_COMPASS$strand)
previous_SJ_with_COMPASS$intron_ID <- gsub(" ", "_", previous_SJ_with_COMPASS$intron_ID)
ann_introns_overlap_df$compare_to_compass <- ifelse(ann_introns_overlap_df$intron_ID %in% previous_SJ_with_COMPASS$intron_ID,"true","false")
new_sj_compass_count_filtered <- subset(new_sj, new_sj$COMPASS_counts > 20)

combined_tbl$intron_ID <- gsub(" ", "_", combined_tbl$intron_ID)
combined_tbl$intron_ID <- paste("chr", combined_tbl$intron_ID)
combined_tbl$intron_ID <- str_replace_all(combined_tbl$intron_ID, fixed(" "), "")
combined_tbl$previously_reported <- ifelse(combined_tbl$intron_ID %in% previous_SJ$intron_ID,"true","false")
compass_new_introns <- subset(combined_tbl, combined_tbl$previously_reported == 'false')
compass_nonoverlapping_introns <- subset(compass_new_introns, compass_new_introns$annotated_junction == "FALSE")

combined_tbl %>% group_by(sample_name) %>% tally()
write.csv(combined_tbl, file = "/Users/chanfreaulab/Documents/combined_tbl.csv")
combined_tbl <- read.csv(file = '/Users/chanfreaulab/Documents/combined_tbl.csv')

#created mean FAnS by grouping the replicates based on sample_prefix (rrp6, upf1, or u_r)
combined_tbl$FAnS_mean <- combined_tbl %>%
  group_by(sample_prefix) %>%
  summarise(across(FAnS, mean, na.rm = TRUE))

test_df <- combined_tbl %>%
  group_by(intron_ID)

combine_replicates <- combined_tbl %>% 
  group_by(intron_ID, gene, sample_prefix) %>% 
  summarize(Ave = mean(FAnS), SE = sd(FAnS)/sqrt(n())) %>% 
  ungroup()

ggplot(data = combine_replicates,aes(x = intron_ID, y = Ave, fill = sample_prefix)) + geom_point(position = position_dodge(width = 1), stat = "identity")

ggplot(data = combine_replicates,aes(x = intron_ID, y = Ave)) + geom_point(aes(colour = sample_prefix))

#subset the data into separate frames based on which knockout they are
rrp6 <- subset(combine_replicates, combine_replicates$sample_prefix == 'RRP6')
upf1 <- subset(combine_replicates, combine_replicates$sample_prefix == 'UPF1')
u_r <- subset(combine_replicates, combine_replicates$sample_prefix == 'U_R')

#removed spaces from intron_ID
rrp6$intron_ID <- gsub(" ", "_", rrp6$intron_ID)
u_r$intron_ID <- gsub(" ", "_", u_r$intron_ID)
upf1$intron_ID <- gsub(" ", "_", upf1$intron_ID)

#are the different introns found in the same dataframes?
u_r$found_in_rrp6 <- ifelse(u_r$intron_ID %in% rrp6$intron_ID,"true","false")
u_r$found_in_upf1 <- ifelse(u_r$intron_ID %in% upf1$intron_ID,"true","false")
rrp6$found_in_rrp6 <- ifelse(rrp6$intron_ID %in% u_r$intron_ID,"true","false")
upf1$found_in_upf1 <- ifelse(upf1$intron_ID %in% u_r$intron_ID,"true","false")
subset_upf1 <- subset(upf1, upf1$found_in_upf1=='true')
subset_rrp6 <- subset(rrp6, rrp6$found_in_rrp6=='true')
u_r_rrp6_common <- subset(u_r, u_r$found_in_rrp6=='true')

u_r_upf1_common <- subset(u_r, u_r$found_in_upf1=='true')

#Plots
ggplot() +
  geom_point(data=u_r_upf1_common, aes(x=subset_upf1$Ave, y=Ave)) + 
  geom_abline(slope=1, intercept=0) + 
  scale_x_log10() + 
  scale_y_log10() + 
  labs(title = "FAnS comparision of UPF1∆RRP6∆ vs UPF1∆") +
  xlab("UPF1∆ FAnS") +
  ylab("UPF1∆RRP6∆ FAnS")

ggplot() 
  + geom_point(data=u_r_rrp6_common, aes(x=subset_rrp6$Ave, y=Ave)) + ylim(0, 1) + xlim(0, 1) + geom_abline(slope=1, intercept=0)

ggplot() + geom_point(data=u_r_upf1_common, aes(x=subset_upf1$Ave, y=Ave)) + ylim(0, 0.15) + xlim(0, 0.15) + geom_abline(slope=1, intercept=0)
ggplot() + geom_point(data=u_r_rrp6_common, aes(x=subset_rrp6$Ave, y=Ave)) + ylim(0, 0.15) + xlim(0, 0.15) + geom_abline(slope=1, intercept=0)

#list.data <- list(data.frame(x=rrp6$intron_ID, y=rrp6$Ave),data.frame(x=u_r$intron_ID, y=u_r$Ave))
#ggplot(data = do.call(rbind, list.data), ) + 
#  geom_point(aes(x, y))

#combine_replicates <- aggregate(combined_tbl$FAnS,by=list(name=combined_tbl$sample_prefix,sample_name=combined_tbl$sample_name,intron_ID=combined_tbl$intron_ID,FAnS=combined_tbl$FAnS,name=combined_tbl$Name,gene=combined_tbl$gene),data=combined_tbl,FUN=mean)
