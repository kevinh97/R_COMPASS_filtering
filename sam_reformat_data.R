combined_tbl <- read.csv("./combined_tbl.csv")

subset_combined_tbl <- function(combined_compass_tbl) {

  samples_vector <<- levels(combined_compass_tbl$sample_name)
  
  master_intron_ID_vector <<- levels(combined_compass_tbl$intron_ID)
  
  master_df <<- data.frame(row.names = master_intron_ID_vector)
  
  list_tbl <<- vector(mode = "list", length = length(samples_vector))
  
  #Subset each sample to an individual df and combinded those dfs into a list
    for (i in 1:length(samples_vector)) {
      subset_tbl_name <- paste(samples_vector[i], "tbl", sep = "")
      print(subset_tbl_name)
      
      subset_tbl <- subset(combined_compass_tbl, combined_compass_tbl$sample_name == samples_vector[i])
      
      list_tbl[[i]] <<- subset_tbl
    }
    for (i in 1:nrow(master_df)) {
      for (j in 1:length(list_tbl)){
        matched_row <- list_tbl[[j]][list_tbl[[j]]["intron_ID"] == rownames(master_df)[i]]
        if (length(matched_row) == 103) {
          #print(paste(rownames(master_df)[i], "in df", j))
          master_df[i,1:51] <<- matched_row[1:51]
          # uniques columns [52:101]
          break
        } else if (length(matched_row) == 0) {
          #print(paste(rownames(master_df)[i], "not in df", j))
        } else {
          # print(
          #   paste(
          #   "intron_ID shows up", 
          #   (length(matched_row)/103), 
          #   "times in df",
          #   j,
          #   "SKIPPING")
          #       )
        }
      }
    }
  }

subset_combined_tbl(combined_tbl)

#rownames(master_df)[17] is in list_tbl[[1]]["intron_ID"]
#list_tbl[[1]][list_tbl[[1]]["intron_ID"] == rownames(master_df)[17]]