library(tidyverse)
library(vwr)
library(plyranges)
library(Biostrings)
library(readxl)

#bind_rows is not optimal since it can append unnecessary rows if not run from source

combined_junctions = tibble()
DIR <- '/Users/chanfreaulab/Documents/'

# last_name_first_author <- 'Xu'
# SRR_PREFIX <- 'SRR72'

last_name_first_author <- "He" 
SRR_PREFIX <- "x_" 
files <- Sys.glob(paths = paste0(DIR, SRR_PREFIX, "*COMPASS_splice_junctions_with_sequence_info.tsv")) # list.files(path = DIR, pattern = "SRR72*_COMPASS_splice_junctions_with_sequence_info.tsv" )
files

for (file in files) { print(file) 
  junction_df <- read_tsv(file)
  combined_junctions <- bind_rows(combined_junctions, junction_df)
}

#last_name_first_author <- 'Talkish'
#SRR_PREFIX <- 'SRR50'
# sample_key <- read_excel('/Volumes/SPxDrive/COMPASS/splicing_datasets/Talkish_et_al_hungry_spliceosome_Ares/Talkish_sample_key.xlsx')
#files <- Sys.glob(paths = paste0(DIR, SRR_PREFIX, "*_COMPASS_splice_junctions_with_sequence_info.tsv")) # list.files(path = DIR, pattern = "SRR72*_COMPASS_splice_junctions_with_sequence_info.tsv" )
#files

#for (file in files) { print(file) 
#  junction_df <- read_tsv(file)
#  combined_junctions <- bind_rows(combined_junctions, junction_df)
#}

#last_name_first_author <- 'Aslanzadeh'
#SRR_PREFIX <- 'SRR55'
# sample_key <- read_excel('/Volumes/SPxDrive/COMPASS/splicing_datasets/Aslanzadeh_et_al/Aslanzadeh_sample_key.xlsx')
#files <- Sys.glob(paths = paste0(DIR, SRR_PREFIX, "*_COMPASS_splice_junctions_with_sequence_info.tsv")) # list.files(path = DIR, pattern = "SRR72*_COMPASS_splice_junctions_with_sequence_info.tsv" )
#files

#for (file in files) { print(file) 
#  junction_df <- read_tsv(file)
#  combined_junctions <- bind_rows(combined_junctions, junction_df)
#}

last_name_first_author <- 'He_Kevin'
#sample_key <- read_excel('/Volumes/GoogleDrive/My Drive/Scripts/COMPASS/combined_sample_key.xlsx') # '/Volumes/SPxDrive/COMPASS/splicing_datasets/combined_sample_key.xlsx')

annotated_introns_filename = '/Users/chanfreaulab/Documents/saccharomyces_cerevisiae_R64-2-1_20150113_introns_no_chr.tsv' # '/Volumes/SPxDrive/COMPASS/S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113_introns.tsv'
#annotated_introns <- read_tsv(annotated_introns_filename) %>% filter(Name != 'YBR090C_intron') ## this intron is duplicated, as it is also a 5'UTR intron for NH6PB
annotated_introns <- read_tsv(annotated_introns_filename)

junction_df_test <- combined_junctions %>% mutate(start = amb_start + 1, end = amb_stop + 1, strand = RNA_strand) %>% 
  left_join(annotated_introns) %>% 
  mutate(seqnames = chrom,
         intron_size = amb_stop - amb_start + 1,
         intron_ID = paste(chrom, start, end, strand),
         ann_5SS = if_else(RNA_strand == '+',
                           start %in% annotated_introns$start,
                           end %in% annotated_introns$end), 
         ann_3SS = if_else(RNA_strand == '+',
                           end %in% annotated_introns$end,
                           start %in% annotated_introns$start) ,
         annotated_junction = ann_5SS & ann_3SS & !is.na(intron_type),
         amb_junc_US = (substr(US_5SS_10nt, 10, 10) == substr(US_3SS_10nt, 10, 10) ),
         amb_junc_DS = (substr(DS_5SS_10nt, 1, 1) == substr(DS_3SS_10nt, 1, 1) ) ,
         edit_dist_US_junction = levenshtein.distance(US_5SS_10nt, US_3SS_10nt),
         edit_dist_DS_junction = levenshtein.distance(DS_5SS_10nt, DS_3SS_10nt), 
         hamming_dist_US_junction = hamming.distance(US_5SS_10nt, US_3SS_10nt),
         hamming_dist_DS_junction = hamming.distance(DS_5SS_10nt, DS_3SS_10nt), 
         fraction_reads_with_junction_proximal_mismatch = mismatch_near_junction_counts / COMPASS_counts,
         five_SS_2nt = substr(fiveSS_seq, 1, 2), 
         three_SS_2nt = substr(threeSS_seq, 2, 3),
         consensus_branchpoint = 'TACTAAC',
         branchpoint_consensus_hamming_dist = hamming.distance(best_BP_sequence, consensus_branchpoint),
         averaged_fiveSS_threeSS_unspliced_reads = (fiveSS_unspliced_reads + threeSS_unspliced_reads) / 2,
         spliced_to_unspliced_ratio = COMPASS_counts / averaged_fiveSS_threeSS_unspliced_reads ,
  ) %>% 
  select(-index, -amb_start, -amb_stop, -intron_length, -five_SS, -three_SS) %>%
  relocate(intron_ID, annotated_junction, chrom, seqnames, start, end, RNA_strand:intron_size, poly_U_count:DS_3SS_10nt, fiveSS_seq, threeSS_seq, fiveSS_type:hamming_dist_DS_junction, five_SS_2nt: branchpoint_consensus_hamming_dist,
           sample_name:threeSS_unspliced_reads, COMPASS_counts:median_alignment_score, upstream_5SS, downstream_5SS, upstream_3SS, downstream_3SS, ) #%>%
  #left_join(sample_key)

junction_df_test %>% write_tsv( paste0(DIR, last_name_first_author, "_COMPASS_splice_junctions_junction_df.tsv")  )

#last_name_first_author <- 'He_Kevin'
#junction_df <- read_tsv( paste0(DIR, last_name_first_author, "_COMPASS_splice_junctions_junction_df.tsv") )

# this summarize step groups together all of the samples, and summarizes the evidence for each junction
junction_df_summarised <- junction_df_test %>%
  relocate(intron_ID:branchpoint_consensus_hamming_dist, unique_US_perfect_matches, unique_DS_perfect_matches, 
           max_US_perfect_matches, max_DS_perfect_matches, best_alignment_score, fraction_reads_with_junction_proximal_mismatch, mean_alignment_score) %>%
  group_by_at(.vars = vars(intron_ID:branchpoint_consensus_hamming_dist) ) %>% 
  summarise(across(mean_alignment_score, mean), 
            across(best_alignment_score, min),
            across(unique_US_perfect_matches:max_DS_perfect_matches, max),
            across(total_reads:perfect_gapped_alignment_counts, sum),
            across(mean_five_SS_Q_score:mean_three_SS_Q_score, mean),
  )  %>% 
  mutate(fraction_reads_with_junction_proximal_mismatch = mismatch_near_junction_counts / COMPASS_counts,
         averaged_fiveSS_threeSS_unspliced_reads = (fiveSS_unspliced_reads + threeSS_unspliced_reads) / 2,
         spliced_to_unspliced_ratio = COMPASS_counts / averaged_fiveSS_threeSS_unspliced_reads,
  ) 

junction_df_summarised %>% write_tsv( paste0(DIR, last_name_first_author, "_COMPASS_splice_junctions_junction_df_summarised.tsv")  )

junction_df_summarised <- read_tsv( paste0(DIR, last_name_first_author, "_COMPASS_splice_junctions_junction_df_summarised.tsv") )

filtered_junction_df_summarised <- junction_df_summarised %>% filter(chrom != 'chrmt', intron_size > 40, 
                                                                     mean_alignment_score < 2, best_alignment_score == 0,
                                                                     unique_US_perfect_matches > 2, unique_DS_perfect_matches > 2,
                                                                    # branchpoint_consensus_hamming_dist < 3, YTRAY_branchpoint,
                                                                     max_DS_perfect_matches > 10, max_US_perfect_matches > 8, 
                                                                     fraction_reads_with_junction_proximal_mismatch < 0.2, 
                                                                     edit_dist_US_junction > 3, edit_dist_DS_junction > 3,
                                                                     num_amb_junctions < 10, # this is the length of ambiguous bp around the 5'SS and 3'SS
                                                                     !(chrom == 'chrXII' & start > 450000 & start < 470000), ## rRNA
                                                                     !(chrom == 'chrII' & start > 680500 & start < 682000), ## LSR1  ## potentially follow up on these, perhaps some are real?
                                                                     !(chrom == 'chrXV' & start > 721800 & start < 722800), ## HIS3
                                                                    !(chrom == 'chrIV' & start > 461700 & start < 462400), ## TRP1
                                                                    !(chrom == 'chrV' & start > 115500 & start < 117200), ## URA3
                                                                    !(chrom == 'chrIII' & start > 91200 & start < 92600), ## LEU2
                                                                    !(chrom == 'chrII' & start > 276000 & start < 276300), ## deletion in GAL1,10 in fast/slow Pol II mutants?
                                                                    # chrII 276125 276179
                                                                    # CT      AC  junctions appear to be artifacts on strand switching of abundant GT-AG junctions?
                                                                    ! (adj_5SS == 'CT' & adj_3SS == 'AC' & ann_5SS & ann_3SS),
                                                                     # homopolymer filter
                                                                     !str_detect(US_5SS_10nt, strrep('A',8) ),
                                                                     !str_detect(US_5SS_10nt, strrep('T',8) ),
                                                                     !str_detect(US_5SS_10nt, strrep('G',8) ),
                                                                     !str_detect(US_5SS_10nt, strrep('C',8) ),
)

junction_df_filtered <- junction_df_test %>% filter(intron_ID %in% filtered_junction_df_summarised$intron_ID)

# The vast majority of junctions are removed by the combined action of the filter steps
junction_df_summarised # 214,417
filtered_junction_df_summarised # 3,538 

# What is the impact of each filter?
junction_df_summarised %>% filter(chrom == 'chrmt') # 1,941
junction_df_summarised %>% filter(intron_size <= 40) # 43,450 
junction_df_summarised %>% filter(mean_alignment_score >= 2) # 142,556
junction_df_summarised %>% filter(best_alignment_score > 0) # 159,227
junction_df_summarised %>% filter(unique_US_perfect_matches <= 2) # 198,058
junction_df_summarised %>% filter(unique_DS_perfect_matches <= 2) # 187,249
junction_df_summarised %>% filter(max_DS_perfect_matches <= 10) # 91,855
junction_df_summarised %>% filter(max_US_perfect_matches <= 8) # 122,306 for 10, 116,930 for 8
junction_df_summarised %>% filter(fraction_reads_with_junction_proximal_mismatch >= 0.2) # 141,554
junction_df_summarised %>% filter(edit_dist_US_junction <= 3) # 23,379
junction_df_summarised %>% filter(edit_dist_DS_junction <= 3) # 26,187
junction_df_summarised %>% filter(num_amb_junctions >= 10) # 7,224
junction_df_summarised %>% filter(chrom == 'chrXII' & start > 450000 & start < 470000) ## rRNA # 3,144   
junction_df_summarised %>% filter(chrom == 'chrII' & start > 680500 & start < 682000) ## LSR1 # 1,543
junction_df_summarised %>% filter(chrom == 'chrXV' & start > 721800 & start < 722800) ## HIS3 # 24
junction_df_summarised %>% filter(chrom == 'chrIV' & start > 461700 & start < 462400) ## TRP1 # 15
junction_df_summarised %>% filter(chrom == 'chrV' & start > 115500 & start < 117200) ## URA3 # 45
junction_df_summarised %>% filter(chrom == 'chrIII' & start > 91200 & start < 92600) ## LEU2 # 24
junction_df_summarised %>% filter(chrom == 'chrII' & start > 276000 & start < 276300) ## deletion in GAL1,10 in fast/slow Pol II mutants? # 4
junction_df_summarised %>% filter(adj_5SS == 'CT' & adj_3SS == 'AC' & ann_5SS & ann_3SS) ## CT-AC  junctions opposite strand of annotated introns # 24
junction_df_summarised %>% filter(str_detect(US_5SS_10nt, strrep('A',8) )) # 1,374
junction_df_summarised %>% filter(str_detect(US_5SS_10nt, strrep('T',8) )) # 4,721 
junction_df_summarised %>% filter(str_detect(US_5SS_10nt, strrep('G',8) )) # 1
junction_df_summarised %>% filter(str_detect(US_5SS_10nt, strrep('C',8) )) # 0

# check why annotated junctions do not agree with annotated splice sites
filtered_junction_df_summarised %>% filter(ann_5SS, ann_3SS, !annotated_junction) %>% select(adj_5SS, adj_3SS, intron_ID)

# how many unique junctions passing the global filters are found in each sample by replicate
junction_df_filtered %>% group_by(sample_info, replicate) %>% tally()

unique(junction_df_filtered$sample_info)
# need to add zeros for junctions with zero counts in specific samples

## OVERLAP WITH ANNOTATED INTRONS ##
################################################################################################################################
library(plyranges)

annotated_introns %>% group_by(intron_type) %>% tally()

annotated_introns_with_junction_counts <- annotated_introns %>% left_join( junction_df_filtered %>% ungroup() %>%
                                                                             select(chrom, start, end, strand, BP_3SS_dist, sample_name, COMPASS_counts) ) %>%
                                                                             dplyr::rename(annotated_COMPASS_counts = COMPASS_counts, 
                                                                                           annotated_start = start,
                                                                                           annotated_end = end,
                                                                                           annotated_BP_3SS_dist = BP_3SS_dist) 

ann_introns <- annotated_introns %>% 
  dplyr::rename(seqnames = chrom) %>% 
  mutate(annotated_start = start, annotated_end = end) %>% as_granges()

junction_df_filtered_granges <- junction_df_filtered %>% 
  select(-type, -Name, -intron_type, -gene) %>%
  mutate(seqnames = chrom) %>% as_granges()

ann_introns_overlap_df <- join_overlap_left(junction_df_filtered_granges, ann_introns) %>% as_tibble()
combined_tbl <- ann_introns_overlap_df %>% left_join(annotated_introns_with_junction_counts)

combined_tbl$intron_ID

previous_SJ <-readxl::read_excel('~/Downloads/previous_studies_SJ.xlsx') %>% 
  #  '/Volumes/SPxDrive/COMPASS/splicing_datasets/previous_studies_SJ.xlsx') %>% 
  mutate(start = if_else(source == "Talkish et al.", start + 1, start)) %>% 
  mutate(start = if_else(source == "Gould et al.", start + 2, start)) %>%
  mutate(intron_ID = paste(seqnames, start, end, strand))
last_name_first_author = "Talkish"
previous_SJ_specific_study <- previous_SJ %>% filter(source == paste(last_name_first_author, 'et al.') )

previous_SJ$intron_ID

combined_tbl <- combined_tbl %>% mutate(FAnS = COMPASS_counts / annotated_COMPASS_counts, 
                                        spliced_to_unspliced_ratio_each_sample = COMPASS_counts / (fiveSS_unspliced_reads + threeSS_unspliced_reads)/2,
                                        overlapping_RPG_annotated_intron = str_detect(gene, 'RPL') | str_detect(gene, 'RPS'),
                                        overlapping_annotated_intron = !is.na(Name),
                                        in_previous_datasets = intron_ID %in% previous_SJ$intron_ID,
                                        in_previous_SJ_specific_study = intron_ID %in% previous_SJ_specific_study$intron_ID,
                                        ) %>%
  replace_na(list(overlapping_RPG_annotated_intron = 'not_overlapping_intron')) 

combined_tbl %>% group_by(overlapping_RPG_annotated_intron) %>% tally()

junction_df_filtered_summarized_by_replicates <- combined_tbl  %>% relocate(sample_info, gene, annotated_start, annotated_end, annotated_BP_3SS_dist, in_previous_datasets, in_previous_SJ_specific_study, overlapping_RPG_annotated_intron, overlapping_annotated_intron) %>%
  group_by_at(.vars = vars(sample_info:branchpoint_consensus_hamming_dist) )  %>%
  summarise(across(spliced_to_unspliced_ratio, list(mean = mean, min = min, max = max)),
            across(FAnS, list(mean = mean,  min = min, max = max) ),
  ) 

junction_df_filtered_summarized_by_replicates %>% group_by(overlapping_RPG_annotated_intron) %>% tally()

junction_df_filtered_summarized_by_replicates %>% write_tsv( paste0(DIR, last_name_first_author, "_COMPASS_splice_junctions_junction_df_filtered_summarised_by_replicates.tsv")  )
junction_df_filtered_summarized_by_replicates <- read_tsv( paste0(DIR, last_name_first_author, "_COMPASS_splice_junctions_junction_df_filtered_summarised_by_replicates.tsv") )


junction_df_filtered_summarized_by_replicates_wider <- junction_df_filtered_summarized_by_replicates %>%
  pivot_wider(names_from = sample_info, values_from = c(spliced_to_unspliced_ratio_mean, 
                                                        spliced_to_unspliced_ratio_min, 
                                                        spliced_to_unspliced_ratio_max,
                                                        FAnS_mean,
                                                        FAnS_min,
                                                        FAnS_max
                                                        ), values_fill = 0) 

# check how many filtered junctions overlap known introns
junction_df_filtered_summarized_by_replicates_wider %>% group_by(overlapping_RPG_annotated_intron) %>% tally()

# check how many annotated vs unannotated junctions are reported in a specific previous study 
junction_df_filtered_summarized_by_replicates_wider %>% group_by(annotated_junction, in_previous_SJ_specific_study, ) %>% tally()

# or in any of the previous splicing datasets
junction_df_filtered_summarized_by_replicates_wider %>% group_by(annotated_junction, in_previous_datasets ) %>% tally()
## TODO check in previous datasets

# previous_SJ
previous_SJ %>% group_by(source) %>% tally()

# spot checks
################################################################################################################################
combined_tbl$in_previous_SJ_specific_study

combined_tbl %>% filter(adj_5SS != 'GT', FAnS > 0.1 | is.na(FAnS), !in_previous_datasets, min_score_unannotated_intron_junctions_disagree/COMPASS_counts < .1) %>% 
  arrange(-spliced_to_unspliced_ratio_each_sample) %>% # potential_5SS, potential_3SS, 
  select(intron_ID, adj_5SS, adj_3SS, bbmap_counts, HISAT2_noncanonical_counts, HISAT2_default_counts, FAnS,
         spliced_to_unspliced_ratio_each_sample, sample_info ) %>%
  print(n = 100)


combined_tbl  %>% filter(!annotated_junction, !in_previous_datasets) %>% 
  arrange(-FAnS) %>% 
  select(sample_info, gene, seqnames, start, end, adj_5SS, adj_3SS, FAnS, COMPASS_counts, bbmap_counts, STAR_default_counts, HISAT2_default_counts) %>% print(n = 50)

combined_tbl  %>% filter(gene == 'APE2', adj_3SS != 'AG') %>% 
  arrange(-FAnS) %>% 
  select(sample_info, gene, seqnames, start, end, adj_5SS, adj_3SS, FAnS, COMPASS_counts, bbmap_counts, STAR_default_counts, HISAT2_default_counts) %>% print(n = 50)

## Aslanzadeh analysis
################################################################################################################################
COUNTS_FILTER <- 1

options(scipen=10000)

junction_df_filtered_summarized_by_replicates_wider %>% filter(!annotated_junction) %>%
  group_by(overlapping_annotated_intron, overlapping_RPG_annotated_intron, in_previous_SJ_specific_study) %>%
  tally()

junction_df_filtered_summarized_by_replicates_wider <- junction_df_filtered_summarized_by_replicates_wider %>% mutate(
  annotation_for_plot = if_else(in_previous_SJ_specific_study, 'in Aslanzadeh et al.', overlapping_RPG_annotated_intron),
  alt_ann_3SS_dist = if_else(strand == '+', end - annotated_end, annotated_start - start),
  FAnS_fast_vs_normal = FAnS_mean_fast / FAnS_mean_normal,
  FAnS_slow_vs_normal = FAnS_mean_slow / FAnS_mean_normal,
)

junction_df_filtered_summarized_by_replicates_wider %>% 
  filter(!annotated_junction, 
         (spliced_to_unspliced_ratio_min_fast > 0) | (spliced_to_unspliced_ratio_min_normal > 0)) %>%
  group_by(annotation_for_plot) %>% tally()

junction_df_filtered_summarized_by_replicates_wider$overlapping_RPG_annotated_intron

ggplot(junction_df_filtered_summarized_by_replicates_wider %>% 
         filter(!annotated_junction, (spliced_to_unspliced_ratio_min_slow > 0) | (spliced_to_unspliced_ratio_min_normal > 0),
                              adj_3SS %in% c( 'TG', 'GG', 'AT', 'AC', 'AG', 'AA', 'CG'), #,  'AA' ),
) , 
aes(x = spliced_to_unspliced_ratio_mean_normal, y = spliced_to_unspliced_ratio_mean_slow, 
    color = overlapping_RPG_annotated_intron) ) + 
  geom_errorbar(aes(ymin = spliced_to_unspliced_ratio_min_slow,
                    ymax = spliced_to_unspliced_ratio_max_slow), size = 0.5) +
  geom_errorbarh(aes(xmin = spliced_to_unspliced_ratio_min_normal,
                     xmax = spliced_to_unspliced_ratio_max_normal), size = 0.5) +
  geom_jitter(alpha = .5) + scale_x_log10( ) + scale_y_log10( ) + geom_abline(slope = 1) + 
  facet_wrap( ~ adj_3SS, nrow = 1) + theme_Publication()
ggsave(paste0(DIR, 'Aslanzadeh_slow_vs_normal_3SS_by_overlap_annotated.svg'), width = 15, height = 4)


ggplot(junction_df_filtered_summarized_by_replicates_wider %>% 
         filter(!annotated_junction, (spliced_to_unspliced_ratio_min_fast > 0) | (spliced_to_unspliced_ratio_min_normal > 0),
                adj_3SS %in% c( 'TG', 'GG', 'AT', 'AC', 'AG', 'AA', 'CG'), #,  'AA' ),
         ) , 
       aes(x = spliced_to_unspliced_ratio_mean_normal, y = spliced_to_unspliced_ratio_mean_fast, 
           color = overlapping_RPG_annotated_intron) ) + 
  geom_errorbar(aes(ymin = spliced_to_unspliced_ratio_min_fast,
                    ymax = spliced_to_unspliced_ratio_max_fast), size = 0.5) +
  geom_errorbarh(aes(xmin = spliced_to_unspliced_ratio_min_normal,
                     xmax = spliced_to_unspliced_ratio_max_normal), size = 0.5) +
  geom_jitter(alpha = .5) + scale_x_log10( ) + scale_y_log10( ) + geom_abline(slope = 1) + 
  facet_wrap( ~ adj_3SS, nrow = 1) + theme_Publication()

ggsave(paste0(DIR, 'Aslanzadeh_fast_vs_normal_3SS_by_overlap_annotated.svg'), width = 15, height = 4)

ggplot(junction_df_filtered_summarized_by_replicates_wider %>% 
         filter(!annotated_junction, (spliced_to_unspliced_ratio_min_fast > 0) | (spliced_to_unspliced_ratio_min_slow > 0),
                adj_3SS %in% c( 'TG', 'GG', 'AT', 'AC', 'AG', 'AA', 'CG'), #,  'AA' ),
         ) , 
       aes(x = spliced_to_unspliced_ratio_mean_slow, y = spliced_to_unspliced_ratio_mean_fast, 
           color = overlapping_RPG_annotated_intron) ) + 
  geom_errorbar(aes(ymin = spliced_to_unspliced_ratio_min_fast,
                    ymax = spliced_to_unspliced_ratio_max_fast), size = 0.5) +
  geom_errorbarh(aes(xmin = spliced_to_unspliced_ratio_min_slow,
                     xmax = spliced_to_unspliced_ratio_max_slow), size = 0.5) +
  geom_jitter(alpha = .5) + scale_x_log10( ) + scale_y_log10( ) + geom_abline(slope = 1) + 
  facet_wrap( ~ adj_3SS, nrow = 1) + theme_Publication()
ggsave(paste0(DIR, 'Aslanzadeh_fast_vs_slow_3SS_by_overlap_annotated.svg'), width = 15, height = 4)

ggplot(junction_df_filtered_summarized_by_replicates_wider %>% 
         filter(!annotated_junction, (spliced_to_unspliced_ratio_min_upf1 > 0) | (spliced_to_unspliced_ratio_min_upf1_prp18 > 0),
                adj_3SS %in% c( 'TG', 'GG', 'AT', 'AC', 'AG', 'AA', 'CG'), #,  'AA' ),
         ) , 
       aes(x = spliced_to_unspliced_ratio_mean_upf1, y = spliced_to_unspliced_ratio_mean_upf1_prp18, 
           color = overlapping_RPG_annotated_intron) ) + 
  geom_errorbar(aes(ymin = spliced_to_unspliced_ratio_min_upf1_prp18,
                    ymax = spliced_to_unspliced_ratio_max_upf1_prp18), size = 0.5) +
  geom_errorbarh(aes(xmin = spliced_to_unspliced_ratio_min_upf1,
                     xmax = spliced_to_unspliced_ratio_max_upf1), size = 0.5) +
  geom_jitter(alpha = .5) + scale_x_log10( ) + scale_y_log10( ) + geom_abline(slope = 1) + 
  facet_wrap( ~ adj_3SS, nrow = 1) + theme_Publication()
ggsave(paste0(DIR, 'Aslanzadeh_upf1_vs_upf1_prp18_3SS_by_overlap_annotated.svg'), width = 15, height = 4)

ggplot(junction_df_filtered_summarized_by_replicates_wider %>% 
         filter(!annotated_junction, (FAnS_min_upf1 > 0) | (FAnS_min_upf1_prp18 > 0),
                adj_3SS %in% c( 'TG', 'GG', 'AT', 'AC', 'AG', 'AA', 'CG'), #,  'AA' ),
         ) , 
       aes(x = FAnS_mean_upf1, y = FAnS_mean_upf1_prp18, 
           color = overlapping_RPG_annotated_intron) ) + 
  geom_errorbar(aes(ymin = FAnS_min_upf1_prp18,
                    ymax = FAnS_max_upf1_prp18), size = 0.5) +
  geom_errorbarh(aes(xmin = FAnS_min_upf1,
                     xmax = FAnS_max_upf1), size = 0.5) +
  geom_jitter(alpha = .5) + scale_x_log10( ) + scale_y_log10( ) + geom_abline(slope = 1) + 
  facet_wrap( ~ adj_3SS, nrow = 1) + theme_Publication()
ggsave(paste0(DIR, 'Aslanzadeh_upf1_vs_upf1_prp18_FAnS_3SS_by_overlap_annotated.svg'), width = 15, height = 4)

## Talkish analysis
################################################################################################################################
junction_df_filtered_summarized_by_replicates_wider %>% 
  filter(spliced_to_unspliced_ratio_min_rapamycin_60_min > .01,
         spliced_to_unspliced_ratio_max_time_zero < .01,
         adj_3SS == 'AC',
  )  %>% ungroup() %>%
  select(gene, seqnames, start, end, annotated_junction, ann_5SS, ann_3SS, fiveSS_seq, adj_3SS, spliced_to_unspliced_ratio_min_rapamycin_60_min, spliced_to_unspliced_ratio_max_time_zero, FAnS_mean_time_zero, FAnS_mean_rapamycin_60_min) %>% print(n = 50)

junction_df_filtered_summarized_by_replicates_wider$in_previous_SJ_specific_study
junction_df_filtered_summarized_by_replicates_wider$in_previous_SJ_specific_study

junction_df_filtered_summarized_by_replicates_wider <- junction_df_filtered_summarized_by_replicates_wider %>% 
  mutate(
  annotation_for_plot = if_else(in_previous_SJ_specific_study, 'in Talkish et al.', overlapping_RPG_annotated_intron),
  alt_ann_3SS_dist = if_else(strand == '+', end - annotated_end, annotated_start - start),
  FAnS_rapamycin_vs_time_zero = FAnS_mean_rapamycin_60_min / FAnS_mean_time_zero,
  spliced_to_unspliced_rapamycin_vs_time_zero = spliced_to_unspliced_ratio_mean_rapamycin_60_min / spliced_to_unspliced_ratio_mean_time_zero,
)

junction_df_filtered_summarized_by_replicates_wider %>% group_by(annotation_for_plot) %>% tally()

junction_df_filtered_summarized_by_replicates_wider %>% group_by(overlapping_RPG_annotated_intron) %>% tally()

junction_df_filtered_summarized_by_replicates_wider$alt_ann_3SS_dist

annotated_spliced_to_unspliced_ratio <- junction_df_filtered_summarized_by_replicates_wider %>% ungroup() %>%
  filter(annotated_junction) %>% 
  select(gene, annotated_start, annotated_end, spliced_to_unspliced_rapamycin_vs_time_zero, FAnS_rapamycin_vs_time_zero) %>% 
  dplyr::rename(annotated_spliced_to_unspliced_rapamycin_vs_time_zero = spliced_to_unspliced_rapamycin_vs_time_zero,
                annotated_FAnS_rapamycin_vs_time_zero = FAnS_rapamycin_vs_time_zero)

junction_df_filtered_summarized_by_replicates_wider <- junction_df_filtered_summarized_by_replicates_wider %>% left_join(annotated_spliced_to_unspliced_ratio)

junction_df_filtered_summarized_by_replicates_wider <- junction_df_filtered_summarized_by_replicates_wider %>% 
  mutate(spliced_to_unspliced_in_rapamycin_vs_annotated = spliced_to_unspliced_rapamycin_vs_time_zero / annotated_spliced_to_unspliced_rapamycin_vs_time_zero )

junction_df_filtered_summarized_by_replicates_wider %>% 
  filter(!annotated_junction, 
         (spliced_to_unspliced_ratio_min_rapamycin_60_min > 0) | (spliced_to_unspliced_ratio_max_time_zero > 0)) %>%
  group_by(annotation_for_plot) %>% tally()

junction_df_filtered_summarized_by_replicates_wider %>% ungroup() %>% filter(annotated_junction, threeSS_type == 'AC') %>%
  select(gene,  strand, RNA_strand, threeSS_type, fiveSS_seq, threeSS_seq) %>% unique()

junction_df_filtered_summarized_by_replicates_wider$overlapping_RPG_annotated_intron = factor(junction_df_filtered_summarized_by_replicates_wider$overlapping_RPG_annotated_intron, levels = c('TRUE', 'FALSE', 'not_overlapping_intron') )

junction_df_filtered_summarized_by_replicates_wider %>% group_by(overlapping_RPG_annotated_intron) %>% tally()

junction_df_filtered_summarized_by_replicates_wider$annotated_spliced_to_unspliced_rapamycin_vs_time_zero
junction_df_filtered_summarized_by_replicates_wider$annotated_junction

# ggplot(junction_df_filtered_summarized_by_replicates_wider %>%
#          filter(overlapping_annotated_intron,
#                 (spliced_to_unspliced_ratio_min_rapamycin_60_min > 0) | (spliced_to_unspliced_ratio_min_time_zero > 0),
#                 # adj_3SS %in% c( 'TG', 'GG', 'AT', 'AC', 'AG', 'AA', 'CG'), #,  'AA' ),
#          ) ,
#        aes(x = annotated_spliced_to_unspliced_rapamycin_vs_time_zero, y = spliced_to_unspliced_rapamycin_vs_time_zero,
#            color = overlapping_RPG_annotated_intron) ) +
#   # geom_errorbar(aes(ymin = FAnS_min_rapamycin_60_min,
#   #                   ymax = FAnS_max_rapamycin_60_min), size = 0.5) +
#   # geom_errorbarh(aes(xmin = FAnS_min_time_zero,
#   #                    xmax = FAnS_max_time_zero), size = 0.5) +
#   geom_point(alpha = .2, position = position_jitterdodge()) +
#   scale_y_log10( ) + scale_x_log10() +
#   theme_Publication() + facet_wrap(~ annotated_junction)
# ggsave(paste0(DIR, 'Talkish_FAnS_3SS_boxplots_overlap_annotated.svg'), width = 12, height = 4)

junction_df_filtered_summarized_by_replicates_wider %>% 
  filter(overlapping_annotated_intron, 
         (spliced_to_unspliced_ratio_min_rapamycin_60_min > 0) | (spliced_to_unspliced_ratio_min_time_zero > 0) )

## RPG vs non RPG on alt vs ann splicing FAnS
ggplot(junction_df_filtered_summarized_by_replicates_wider %>%
         filter(overlapping_annotated_intron, (spliced_to_unspliced_ratio_min_rapamycin_60_min > 0) | (spliced_to_unspliced_ratio_min_time_zero > 0),
                # adj_3SS %in% c( 'TG', 'GG', 'AT', 'AC', 'AG', 'AA', 'CG'), #,  'AA' ),
         ) ,
       aes(x = threeSS_type, y = FAnS_rapamycin_vs_time_zero,
           color = overlapping_RPG_annotated_intron) ) +
  # geom_errorbar(aes(ymin = FAnS_min_rapamycin_60_min,
  #                   ymax = FAnS_max_rapamycin_60_min), size = 0.5) +
  # geom_errorbarh(aes(xmin = FAnS_min_time_zero,
  #                    xmax = FAnS_max_time_zero), size = 0.5) +
  geom_boxplot() +
  geom_jitter(alpha = .2, position = position_jitterdodge()) +
  scale_y_log10( ) +
 theme_Publication() + facet_wrap(~ annotated_junction)
# ggsave(paste0(DIR, 'Talkish_FAnS_3SS_boxplots_overlap_annotated.svg'), width = 12, height = 4)

ggplot(junction_df_filtered_summarized_by_replicates_wider %>% 
         filter(!annotated_junction, (FAnS_min_time_zero > 0) | (FAnS_min_rapamycin_60_min > 0),
                adj_3SS %in% c( 'TG', 'GG', 'AT', 'AC', 'AG', 'AA', 'CG'), #,  'AA' ),
         ) , 
       aes(x = FAnS_mean_time_zero, y = FAnS_mean_rapamycin_60_min, 
           color = overlapping_RPG_annotated_intron) ) + 
  # geom_errorbar(aes(ymin = spliced_to_unspliced_ratio_min_rapamycin_60_min,
  #                   ymax = spliced_to_unspliced_ratio_max_rapamycin_60_min), size = 0.5) +
  # geom_errorbarh(aes(xmin = spliced_to_unspliced_ratio_min_time_zero,
  #                    xmax = spliced_to_unspliced_ratio_max_time_zero), size = 0.5) +
  geom_jitter(alpha = .5) + scale_x_log10( ) + scale_y_log10( ) + geom_abline(slope = 1) + 
  facet_grid(ann_3SS ~  fiveSS_type) + theme_Publication()
# ggsave(paste0(DIR, 'Talkish_spliced_to_unspliced_3SS_by_overlap_annotated_only.svg'), width = 9, height = 7)

ggplot(junction_df_filtered_summarized_by_replicates_wider %>% 
       filter(!annotated_junction, (spliced_to_unspliced_ratio_min_rapamycin_60_min > 0) | (spliced_to_unspliced_ratio_min_time_zero > 0),
                # adj_3SS %in% c( 'TG', 'GG', 'AT', 'AC', 'AG', 'AA', 'CG'), #,  'AA' ),
         ) , 
       aes(x = spliced_to_unspliced_ratio_mean_time_zero, y = spliced_to_unspliced_ratio_mean_rapamycin_60_min, 
           color = overlapping_RPG_annotated_intron) ) + 
  # geom_errorbar(aes(ymin = spliced_to_unspliced_ratio_min_rapamycin_60_min,
  #                   ymax = spliced_to_unspliced_ratio_max_rapamycin_60_min), size = 0.5) +
  # geom_errorbarh(aes(xmin = spliced_to_unspliced_ratio_min_time_zero,
  #                    xmax = spliced_to_unspliced_ratio_max_time_zero), size = 0.5) +
  geom_jitter(alpha = .5) + scale_x_log10( ) + scale_y_log10( ) + geom_abline(slope = 1) + 
  facet_grid(ann_5SS ~ threeSS_type) + theme_Publication()
# ggsave(paste0(DIR, 'Talkish_spliced_to_unspliced_3SS_by_overlap_annotated.svg'), width = 15, height = 4)


junction_df_filtered_summarized_by_replicates_wider$alt_ann_3SS_dist

junction_df_filtered_summarized_by_replicates_wider %>% ungroup() %>%
  filter(!annotated_junction, # abs(alt_ann_3SS_dist) > 10,
         (spliced_to_unspliced_ratio_min_rapamycin_60_min > 0) | (spliced_to_unspliced_ratio_min_time_zero > 0) ) %>%
  arrange(-FAnS_mean_rapamycin_60_min) %>%
  select(gene, strand, chrom, start, end, fiveSS_type, threeSS_type, annotated_start, annotated_end, alt_ann_3SS_dist,
         fiveSS_seq, threeSS_seq, spliced_to_unspliced_ratio_mean_rapamycin_60_min,
         spliced_to_unspliced_ratio_mean_time_zero,
         spliced_to_unspliced_rapamycin_vs_time_zero,
         FAnS_mean_rapamycin_60_min,
         FAnS_mean_time_zero,
         FAnS_rapamycin_vs_time_zero,
         annotated_spliced_to_unspliced_rapamycin_vs_time_zero,
         spliced_to_unspliced_in_rapamycin_vs_annotated,
         in_previous_datasets, in_previous_SJ_specific_study, overlapping_RPG_annotated_intron,
         overlapping_annotated_intron) %>%
  mutate_if(is.numeric, round, digits = 3) %>%
  write_tsv(paste0(DIR, 'Talkish_alt_introns_in_rapamycin.tsv') )


ggsave(paste0(DIR, 'Talkish_spliced_over_unspliced_3SS_boxplots_overlap_annotated.svg'), width = 12, height = 4)

## RPG vs non RPG on alt vs ann splicing FAnS
ggplot(junction_df_filtered_summarized_by_replicates_wider %>% 
         filter(! annotated_junction, (spliced_to_unspliced_ratio_min_rapamycin_60_min > 0) | (spliced_to_unspliced_ratio_min_time_zero > 0),
                adj_3SS %in% c( 'TG', 'GG', 'AT', 'AC', 'AG', 'AA', 'CG'), #,  'AA' ),
         ) , 
       aes(x = spliced_to_unspliced_ratio_mean_time_zero, y = spliced_to_unspliced_ratio_mean_rapamycin_60_min, 
           color = annotation_for_plot) ) + 
  # geom_errorbar(aes(ymin = spliced_to_unspliced_ratio_min_rapamycin_60_min,
  #                   ymax = spliced_to_unspliced_ratio_max_rapamycin_60_min), size = 0.5) +
  # geom_errorbarh(aes(xmin = spliced_to_unspliced_ratio_min_time_zero,
  #                    xmax = spliced_to_unspliced_ratio_max_time_zero), size = 0.5) +
  geom_jitter(alpha = .5) + scale_x_log10( ) + scale_y_log10( ) + geom_abline(slope = 1) + 
  facet_grid( ~ adj_3SS) # + theme_Publication()
ggsave(paste0(DIR, 'Talkish_spliced_to_unspliced_3SS_by_overlap_annotated.svg'), width = 15, height = 4)




ggsave(paste0(DIR, 'Aslanzadeh_3SS_by_overlap_annotated.svg'), width = 15, height = 4)
################################################################################################################################


ggplot(junction_df_filtered_summarized_by_replicates_wider %>% 
         filter(!annotated_junction, !ann_3SS, (spliced_to_unspliced_ratio_min_fast > 0) | (spliced_to_unspliced_ratio_min_normal > 0),
              #  adj_3SS %in% c( 'TG', 'GG', 'AT', 'AC', 'AG', 'AA', 'CG'), #,  'AA' ),
         ) , 
       aes(x = alt_ann_3SS_dist, y = FAnS_slow_vs_normal, 
           color = annotation_for_plot) ) + 
  # geom_errorbar(aes(ymin = spliced_to_unspliced_ratio_min_fast,
  #                   ymax = spliced_to_unspliced_ratio_max_fast), size = 0.5) + 
  # geom_errorbarh(aes(xmin = spliced_to_unspliced_ratio_min_normal,
  #                    xmax = spliced_to_unspliced_ratio_max_normal), size = 0.5) +
  geom_jitter(alpha = .5) + scale_y_log10( ) + xlim(-30, 30) +
  theme_Publication() #

splice_site_tally <- sample_counts_junction_df %>% group_by(adj_5SS, adj_3SS, sample_info) %>% tally(wt = COMPASS_counts)


ggplot(splice_site_tally %>% filter(adj_3SS != 'AG'),  
       aes(x = adj_3SS, y = n, fill = sample_info)) +
  geom_col(position = position_dodge()) + theme_Publication() +
  scale_fill_brewer(type = 'qual', drop=FALSE)

ggplot(splice_site_tally %>%
         filter(adj_5SS != 'GT', # !amb_junc_US, !amb_junc_DS,
                # !str_detect(aligner, 'default') 
         ),  
       aes(x = adj_5SS, y = n, fill = sample_info)) +
  geom_col(position = position_dodge()) + theme_Publication() +
  scale_fill_brewer(type = 'qual', drop=FALSE)

## CHECK VARIOUS FILTERS ##
################################################################################################################################

junction_df_summarised %>% filter(!annotated_junction, spliced_to_unspliced_ratio > 100)

junction_df_summarised %>% filter(annotated_junction, fraction_reads_with_junction_proximal_mismatch > 0.5) %>%
  ungroup() %>% 
  select(chrom, start, end, Name, intron_type, COMPASS_counts, fiveSS_unspliced_reads, threeSS_unspliced_reads ) %>%
  print(n = 100)
# select(chrom, start, end) %>% write_tsv( paste0(DIR, 'ann_junc_with_mismatches.txt'))
filtered_junction_df_summarised %>% filter(annotated_junction | (!annotated_junction & spliced_to_unspliced_ratio < 10 ) )

ggplot(filtered_junction_df_summarised, aes(x = fraction_reads_with_junction_proximal_mismatch)) + geom_histogram() + facet_wrap(~ annotated_junction)

filtered_junction_df_summarised %>% filter(!amb_junc_US, !amb_junc_DS) %>% 
  group_by(adj_5SS) %>% tally(sort = TRUE) %>% print(n = 50)

filtered_junction_df_summarised %>% ungroup() %>%
  filter(start == 548632) %>% select(chrom, start, end, strand, intron_size, edit_dist_US_junction, edit_dist_DS_junction, contains('10nt'))

filtered_junction_df_summarised %>% ungroup() %>% group_by(five_SS_2nt, amb_junc_US | amb_junc_DS) %>% tally() %>% print(n =100)

filtered_junction_df_summarised %>% ungroup() %>% filter(five_SS_2nt == 'CA', !amb_junc_US, !amb_junc_DS) %>% 
  arrange(desc(COMPASS_counts)) %>% select(five_SS_2nt, three_SS_2nt, chrom, start, end,
                                           intron_size, perfect_gapped_alignment, COMPASS_counts, bbmap_counts) %>%
  print(n = 70)

filtered_junction_df_summarised %>% ungroup() %>% filter(five_SS_2nt == 'AT', !amb_junc_US, !amb_junc_DS,
                                                         
) %>% 
  arrange(desc(COMPASS_counts)) %>% select(five_SS_2nt, three_SS_2nt, chrom, start, end, strand, intron_size,fraction_reads_with_junction_proximal_mismatch,  COMPASS_counts) %>%
  print(n = 70)

## plots
library(cowplot)
library(ggExtra)
library("RColorBrewer")
library(RColorBrewer)
myColors <- brewer.pal(6,"Set1")
aligner_lst <- factor(c("COMPASS", "bbmap",  "HISAT2 default", "HISAT2 noncanonical", "STAR default", "STAR noncanonical"))
#Create a custom color scale
names(myColors) <- levels(aligner_lst)
colScale <- scale_colour_manual(name = aligner_lst, values = myColors)

filtered_junction_df_summarised %>% filter(bbmap_counts / COMPASS_counts < .1, 
                                           min_score_unannotated_intron_junctions_disagree / COMPASS_counts < 0.01, 
                                           max_DS_perfect_matches > max_US_perfect_matches) %>% ungroup() %>%
  arrange(-COMPASS_counts) %>%
  select(intron_ID, adj_5SS, adj_3SS, ann_5SS, ann_3SS, COMPASS_counts, bbmap_counts, HISAT2_noncanonical_counts) %>%
  print(n = 100)

filtered_junction_df_summarised_longer <- filtered_junction_df_summarised %>% 
  pivot_longer(cols = COMPASS_counts:HISAT2_noncanonical_counts, values_to = "total_counts", names_to = "aligner") %>% 
  mutate(aligner = str_remove(aligner, "_counts")) %>% 
  mutate(aligner = str_replace(aligner, "_", " ")) 

unique(junction_df$sample_info)

junction_df_longer <- junction_df %>% 
  pivot_longer(cols = COMPASS_counts:HISAT2_noncanonical_counts, 
               values_to = "total_counts", names_to = "aligner" ) %>% 
  mutate(aligner = str_remove(aligner, "_counts")) %>% 
  mutate(aligner = str_replace(aligner, "_", " ")) 

combined_tbl %>% filter(gene == 'BIG1')

UBC12 <- combined_tbl  %>% filter(gene.y == 'UBC12') %>% 
  arrange(end) %>% 
  select(gene.y, seqnames, start, end, adj_5SS, threeSS_seq, adj_3SS, COMPASS_counts, bbmap_counts, 
         STAR_default_counts, STAR_noncanonical_counts, HISAT2_default_counts, HISAT2_noncanonical_counts)

UBC12 <- junction_df_longer %>% filter(amb_start == 744153, amb_stop >= 744228, amb_stop <= 744289)
UBC12$aligner <- factor(UBC12$aligner, levels = factor(aligner_lst) )

combined_tbl  %>% filter(gene == 'BIG1') %>% 
  arrange(-FAnS) %>% 
  select(sample_name, gene, FAnS, seqnames, start, end, adj_5SS, threeSS_seq, adj_3SS, COMPASS_counts, bbmap_counts, 
         STAR_default_counts, STAR_noncanonical_counts, HISAT2_default_counts, HISAT2_noncanonical_counts)

filtered_junction_df_summarised_longer_by_gene <- filtered_junction_df_summarised_longer %>% 
  left_join(combined_tbl %>% select(intron_ID, gene) %>% dplyr::rename(overlapping_intron = gene) %>% unique())

unique(junction_df_filtered_summarized_by_replicates$sample_info)
junction_df_filtered_summarized_by_replicates$sample_info = factor(junction_df_filtered_summarized_by_replicates$sample_info, levels = c("time_zero", "rapamycin_60_min"))


filtered_junction_df_summarised_longer_by_gene %>% filter(overlapping_intron == gene_to_check) %>% write_tsv(paste0(DIR, gene_to_check, '.tsv'))
filtered_junction_df_summarised_longer_by_gene$aligner <- factor(filtered_junction_df_summarised_longer_by_gene$aligner, levels = factor(aligner_lst) )

gene_chrom_df <- junction_df_filtered_summarized_by_replicates %>% ungroup() %>% select(gene, chrom, strand) %>% unique()
chrom_lst <- gene_chrom_df$chrom
names(chrom_lst) <- gene_chrom_df$gene
library(ggrepel)
strand_lst <- gene_chrom_df$strand
names(strand_lst) <- gene_chrom_df$gene


# MOB2, YOS1, SEC14,  ,  IWR1, RPS26A YBL092W_five_prime_UTR_intron MAF1, YGL189C_five_prime_UTR_intron
# YGL189C_five_prime_UTR_intron
# positive strand: SPT14, UBC12, RPS10A, RPL37A
# negative strand: YCL002C, REC107, YSC84, PHO85, RAD14, BIG1, TUB3, YGL189C_five_prime_UTR_intron
gene_to_check <- 'UBC12'
strand_lst[[gene_to_check]]

## for positive strand
## which aligners get which junctions

ggplot(filtered_junction_df_summarised_longer_by_gene %>% 
         filter(overlapping_intron == gene_to_check), 
       aes(x = interaction(start, threeSS_seq), y = total_counts, fill = aligner)) + geom_col(position = 'dodge') + # geom_col(position = position_dodge(preserve = "single")) +
  colScale + scale_y_log10() +
  theme_Publication() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_fill_brewer(type = 'qual', drop=FALSE) + facet_wrap(~  interaction('5′SS', chrom_lst[[gene_to_check]], start, fiveSS_seq, sep = ' ')  ) +
  ggtitle(paste(gene_to_check, 'alternative splicing'))
ggsave( paste0( DIR, gene_to_check, '_', last_name_first_author, '_optimal_alignments_by_aligner.svg'), width = 9, height = 7.5)

## spliced_to_unspliced_ratio_mean
p1 <- ggplot(junction_df_filtered_summarized_by_replicates %>% filter(gene == gene_to_check), 
       aes(x = end, y = spliced_to_unspliced_ratio_mean, fill = sample_info)) + geom_col(position = position_dodge(preserve = "single")) +
  colScale + scale_y_log10() +
  geom_errorbar(aes(ymin = spliced_to_unspliced_ratio_min,
                    ymax = spliced_to_unspliced_ratio_max), position='dodge', size = 0.5) +
  theme_Publication() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  geom_text_repel(min.segment.length = 10, data = junction_df_filtered_summarized_by_replicates %>% filter(sample_info == 'rapamycin_60_min', gene == gene_to_check),
            aes(x = end, y = spliced_to_unspliced_ratio_max, label = threeSS_seq), vjust = -1.3) +
  scale_fill_brewer(type = 'qual', drop=FALSE)  + facet_wrap(~ interaction('5′SS', chrom_lst[[gene_to_check]], start, fiveSS_seq, sep = ' ') ) +
  labs(x = paste('3′SS coordinate(', chrom_lst[[gene_to_check]], ')'), 
       title = parse(text = paste0("italic('", gene_to_check, "')", "~alternative~splicing~efficiency\n(", last_name_first_author, "~et~al.)")))

## FAnS_mean
p2 <- ggplot(junction_df_filtered_summarized_by_replicates %>% filter(gene == gene_to_check), 
       aes(x = end, y = FAnS_mean, fill = sample_info)) + geom_col(position = position_dodge(preserve = "single")) +
  colScale + scale_y_log10() +
  geom_errorbar(aes(ymin = FAnS_min,
                    ymax = FAnS_max), position='dodge', size = 0.5) +
  theme_Publication() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  geom_text(data = junction_df_filtered_summarized_by_replicates %>% filter(sample_info == 'time_zero', gene == gene_to_check),
                                                                        aes(x = end, y = FAnS_max, label = threeSS_seq), vjust = -1) +
  scale_fill_brewer(type = 'qual', drop=FALSE)  + facet_wrap(~ interaction('5′SS', chrom_lst[[gene_to_check]], start, fiveSS_seq, sep = ' ') ) +
  labs(x = paste('3′SS coordinate(', chrom_lst[[gene_to_check]], ')'), 
       title = parse(text = paste0("italic('", gene_to_check, "')", "~alternative~splicing~FAnS\n(", last_name_first_author, "~et~al.)")))

plot_grid(p1, p2, labels = "auto")
ggsave( paste0( DIR, gene_to_check, '_', last_name_first_author, '_FAnS_by_alt_SS.pdf'), width = 15, height = 6)

## for negative strand
## which aligners get which junctions
ggplot(filtered_junction_df_summarised_longer_by_gene %>% 
         filter(overlapping_intron == gene_to_check), 
       aes(x = interaction(end, threeSS_seq), y = total_counts, fill = aligner)) + geom_col(position = 'dodge') + # geom_col(position = position_dodge(preserve = "single")) +
  colScale + scale_y_log10() +
  theme_Publication() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_fill_brewer(type = 'qual', drop=FALSE) + facet_wrap(~ interaction(start, fiveSS_seq) ) +
  ggtitle(paste(gene_to_check, 'alternative splicing'))
ggsave( paste0( DIR, gene_to_check, '_', last_name_first_author, '_optimal_alignments_by_aligner.svg'), width = 9, height = 7.5)

## spliced_to_unspliced_ratio_mean
p1 <- ggplot(junction_df_filtered_summarized_by_replicates %>% filter(gene == gene_to_check), 
       aes(x = start, y = spliced_to_unspliced_ratio_mean, fill = sample_info)) + geom_col(position = position_dodge(preserve = "single")) +
  colScale + 
  geom_errorbar(aes(ymin = spliced_to_unspliced_ratio_min,
                    ymax = spliced_to_unspliced_ratio_max), position='dodge', size = 0.5) +
  theme_Publication() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  geom_text_repel(min.segment.length = 10, data = junction_df_filtered_summarized_by_replicates %>% filter(sample_info == 'rapamycin_60_min', gene == gene_to_check),
                  aes(x = start, y = spliced_to_unspliced_ratio_max, label = threeSS_seq), vjust = -1.3) +
  scale_fill_brewer(type = 'qual', drop=FALSE)  + facet_wrap(~ interaction('5′SS', chrom_lst[[gene_to_check]], end, fiveSS_seq, sep = ' ') ) +
  labs(x = paste('3′SS coordinate(', chrom_lst[[gene_to_check]], ')'), 
       title = parse(text = paste0("italic('", gene_to_check, "')", "~alternative~splicing~efficiency\n(", last_name_first_author, "~et~al.)")))

## FAnS_mean
p2 <- ggplot(junction_df_filtered_summarized_by_replicates %>% filter(gene == gene_to_check), 
       aes(x = start, y = FAnS_mean, fill = sample_info)) + geom_col(position = position_dodge(preserve = "single")) +
  colScale + 
  geom_errorbar(aes(ymin = FAnS_min,
                    ymax = FAnS_max), position='dodge', size = 0.5) +
  theme_Publication() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  geom_text(data = junction_df_filtered_summarized_by_replicates %>% filter(sample_info == 'time_zero', gene == gene_to_check),
            aes(x = start, y = FAnS_max, label = threeSS_seq), vjust = -1) +
  scale_fill_brewer(type = 'qual', drop=FALSE)  + facet_wrap(~ interaction('5′SS', chrom_lst[[gene_to_check]], end, fiveSS_seq, sep = ' ') ) +
  labs(x = paste('3′SS coordinate(', chrom_lst[[gene_to_check]], ')'), 
       title = parse(text = paste0("italic('", gene_to_check, "')", "~alternative~splicing~FAnS\n(", last_name_first_author, "~et~al.)")))
dev.off()
plot_grid(p1, p2, labels = "auto")
ggsave( paste0( DIR, gene_to_check, '_', last_name_first_author, '_FAnS_by_alt_SS.pdf'), width = 9, height = 6)

##  minus
## if splice sites are far apart
## spliced_to_unspliced_ratio_mean
p1 <- ggplot(junction_df_filtered_summarized_by_replicates %>% filter(gene == gene_to_check), 
             aes(x = reorder(factor(start), start), y = spliced_to_unspliced_ratio_mean, fill = sample_info)) + geom_col(position = position_dodge(preserve = "single")) +
  colScale + 
  geom_errorbar(aes(ymin = spliced_to_unspliced_ratio_min,
                    ymax = spliced_to_unspliced_ratio_max), position='dodge', size = 0.5) +
  theme_Publication() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  geom_text_repel(min.segment.length = 10, data = junction_df_filtered_summarized_by_replicates %>% filter(sample_info == 'rapamycin_60_min', gene == gene_to_check),
                  aes(x = reorder(factor(start), start), y = spliced_to_unspliced_ratio_max, label = threeSS_seq), vjust = -1.3) +
  scale_fill_brewer(type = 'qual', drop=FALSE)  + facet_wrap(~ interaction('5′SS', chrom_lst[[gene_to_check]], end, fiveSS_seq, sep = ' ') ) +
  labs(x = paste('3′SS coordinate(', chrom_lst[[gene_to_check]], ')'), 
       title = parse(text = paste0("italic('", gene_to_check, "')", "~alternative~splicing~efficiency\n(", last_name_first_author, "~et~al.)")))

## FAnS_mean
p2 <- ggplot(junction_df_filtered_summarized_by_replicates %>% filter(gene == gene_to_check), 
             aes(x = reorder(factor(start), start), y = FAnS_mean, fill = sample_info)) + geom_col(position = position_dodge(preserve = "single")) +
  colScale + 
  geom_errorbar(aes(ymin = FAnS_min,
                    ymax = FAnS_max), position='dodge', size = 0.5) +
  theme_Publication() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  geom_text(data = junction_df_filtered_summarized_by_replicates %>% filter(sample_info == 'time_zero', gene == gene_to_check),
            aes(x = reorder(factor(start), start), y = FAnS_max, label = threeSS_seq), vjust = -1) +
  scale_fill_brewer(type = 'qual', drop=FALSE)  + facet_wrap(~ interaction('5′SS', chrom_lst[[gene_to_check]], end, fiveSS_seq, sep = ' ') ) +
  labs(x = paste('3′SS coordinate(', chrom_lst[[gene_to_check]], ')'), 
       title = parse(text = paste0("italic('", gene_to_check, "')", "~alternative~splicing~FAnS\n(", last_name_first_author, "~et~al.)")))
p2
plot_grid(p1, p2, labels = "auto")
ggsave( paste0( DIR, gene_to_check, '_', last_name_first_author, '_FAnS_by_alt_SS.pdf'), width = 18, height = 6)


### plus
## spliced_to_unspliced_ratio_mean
p1 <- ggplot(junction_df_filtered_summarized_by_replicates %>% filter(gene == gene_to_check), 
             aes(x = reorder(factor(end), end), y = spliced_to_unspliced_ratio_mean, fill = sample_info)) + geom_col(position = position_dodge(preserve = "single")) +
  colScale + 
  geom_errorbar(aes(ymin = spliced_to_unspliced_ratio_min,
                    ymax = spliced_to_unspliced_ratio_max), position='dodge', size = 0.5) +
  theme_Publication() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  geom_text_repel(min.segment.length = 10, data = junction_df_filtered_summarized_by_replicates %>% filter(sample_info == 'rapamycin_60_min', gene == gene_to_check),
                  aes(x = reorder(factor(end), end), y = spliced_to_unspliced_ratio_max, label = threeSS_seq), vjust = -1.3) +
  scale_fill_brewer(type = 'qual', drop=FALSE)  + facet_wrap(~ interaction('5′SS', chrom_lst[[gene_to_check]], start, fiveSS_seq, sep = ' ') ) +
  labs(x = paste('3′SS coordinate(', chrom_lst[[gene_to_check]], ')'), 
       title = parse(text = paste0("italic('", gene_to_check, "')", "~alternative~splicing~efficiency\n(", last_name_first_author, "~et~al.)")))

## FAnS_mean
p2 <- ggplot(junction_df_filtered_summarized_by_replicates %>% filter(gene == gene_to_check), 
             aes(x = reorder(factor(end), end), y = FAnS_mean, fill = sample_info)) + geom_col(position = position_dodge(preserve = "single")) +
  colScale + 
  geom_errorbar(aes(ymin = FAnS_min,
                    ymax = FAnS_max), position='dodge', size = 0.5) +
  theme_Publication() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  geom_text(data = junction_df_filtered_summarized_by_replicates %>% filter(sample_info == 'time_zero', gene == gene_to_check),
            aes(x = reorder(factor(end), end), y = FAnS_max, label = threeSS_seq), vjust = -1) +
  scale_fill_brewer(type = 'qual', drop=FALSE)  + facet_wrap(~ interaction('5′SS', chrom_lst[[gene_to_check]], start, fiveSS_seq, sep = ' ') ) +
  labs(x = paste('3′SS coordinate(', chrom_lst[[gene_to_check]], ')'), 
       title = parse(text = paste0("italic('", gene_to_check, "')", "~alternative~splicing~FAnS\n(", last_name_first_author, "~et~al.)")))

plot_grid(p1, p2, labels = "auto")
ggsave( paste0( DIR, gene_to_check, '_', last_name_first_author, '_FAnS_by_alt_SS.pdf'), width = 9, height = 6)


counts_tally <- junction_df_longer %>% group_by(aligner, annotated_junction) %>% tally(wt = total_counts)
counts_tally$aligner <- factor(counts_tally$aligner, levels = factor(aligner_lst) )

counts_tally <- filtered_junction_df_summarised_longer %>% group_by(aligner, annotated_junction) %>% tally(wt = total_counts)
counts_tally$aligner <- factor(counts_tally$aligner, levels = factor(aligner_lst) )

ggplot(counts_tally,
       aes(x = aligner, y = n, fill = aligner)) + geom_col(position = 'dodge') +
  colScale +
  theme_Publication() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  facet_wrap(~ factor(annotated_junction, levels = c(TRUE, FALSE)), scales = "free_y", nrow = 1) +
  scale_fill_brewer(type = 'qual', drop=FALSE) 

ggsave( paste0( DIR, 'optimal_alignments_by_aligner.pdf'), width = 6, height = 6)


filtered_junction_df_summarised_longer$aligner <- factor(filtered_junction_df_summarised_longer$aligner ,
                                                         levels = factor( aligner_lst ) )

## 5SS sites
p1 <- ggplot(filtered_junction_df_summarised_longer %>%
               filter(total_counts != 0, adj_5SS == 'GT', !annotated_junction, # !amb_junc_US, !amb_junc_DS,
                      # !str_detect(aligner, 'default') 
               ),  
             aes(x = forcats::fct_infreq(adj_5SS), fill = aligner)) +
  geom_bar(position = position_dodge(preserve = "single")) + theme_Publication() +
  scale_fill_brewer(type = 'qual', drop=FALSE)
#p1
# ggsave( paste0( DIR, 'nonGT_5SS_optimal_alignments_by_aligner.pdf'), width = 8, height = 6) 

p2 <- ggplot(filtered_junction_df_summarised_longer %>%
               filter(total_counts != 0, adj_5SS != 'GT',  # !amb_junc_US, !amb_junc_DS,
                      #  !str_detect(aligner, 'default') 
               ),  
             aes(x = forcats::fct_infreq(adj_5SS), fill = aligner)) +
  geom_bar(position = position_dodge(preserve = "single")) + theme_Publication() + 
  scale_fill_brewer(type = 'qual', drop=FALSE)

# ggsave( paste0( DIR, 'GT_5SS_optimal_alignments_by_aligner.pdf'), width = 2, height = 6)

## 3SS sites
p3 <- ggplot(filtered_junction_df_summarised_longer %>%
               filter(total_counts != 0, adj_3SS == 'AG', !annotated_junction, # !amb_junc_US, !amb_junc_DS,
                      #  !str_detect(aligner, 'default') 
               ),  
             aes(x = forcats::fct_infreq(adj_3SS), fill = aligner)) +
  geom_bar(position = position_dodge(preserve = "single")) + theme_Publication()  +
  scale_fill_brewer(type = 'qual', drop=FALSE)
#p3
# ggsave( paste0( DIR, 'nonAG_3SS_optimal_alignments_by_aligner.pdf'), width = 8, height = 6)

p4 <- ggplot(filtered_junction_df_summarised_longer %>%
               filter(total_counts != 0, adj_3SS != 'AG', # !amb_junc_US, !amb_junc_DS,
                      # !str_detect(aligner, 'default') 
               ),  
             aes(x = forcats::fct_infreq(adj_3SS), fill = aligner)) +
  geom_bar(position = position_dodge(preserve = "single")) + theme_Publication() +
  scale_fill_brewer(type = 'qual', drop=FALSE)

# ggsave( paste0( DIR, 'AG_3SS_optimal_alignments_by_aligner.pdf'), width = 2, height = 6)

plot_grid(p1, p2, p3, p4, labels = 'auto', rel_widths = c(1, 4, 1, 4))
ggsave( paste0( DIR, last_name_first_author, '_optimal_alignments_by_aligner_by_splice_site.pdf'), width = 8, height = 8)
