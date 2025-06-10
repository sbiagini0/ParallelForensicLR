############################################################
# Function: Computes family-wise likelihood ratios (LRs) for a large set of
#            Persons of Interest (POIs) against multiple reference pedigrees.
#            First filters POIs by an exclusions list (with `computeExclusions()`), 
#            then processes them in parallel and in manageable chunks to reduce memory overhead.  
#            Each LR is returned in the form  "<LR_total> (<nMarkers>) / <LR_mismatch> (<nMarkers_mismatch>)".
# pm: A named list of POI pedigrees.
# am: A named list of reference family pedigrees.
# exclusions_list: Named list (by FAMID) whose elements are vectors of POI
#                  IDs allowed to be compared with each family.
# ncores: Integer; number of CPU cores to use for the parallel backend
#         (default: total cores ??? 1).
# chunk_size: Integer; number of POIs to process per parallel batch
#             (default: 500).  Adjust to balance speed vs. memory.
# Returns: A data frame where each row corresponds to a POI and contains:
#          . POI     - POI identifier  
#          . Sex     - Sex inferred via `getSexID()`  
#          . <FAMID> - One column per family, holding the formatted LR result
#                      string or NA if the POI was excluded for that family.
#          The columns "POI" and "Sex" are always first, followed by the
#          families in the order of `am`.  A timing message (minutes) is
#          printed to the console on completion.
ParallelForensicLR <- function(pm, am, exclusions_list, ncores = parallel::detectCores() - 1, chunk_size = 500) {
  if (length(pm) == 0 || length(am) == 0) {
    stop("Must supply at least one POI (pm) and one pedigree (am).")
  }
  
  start_time <- Sys.time()
  
  poi_valid <- intersect(names(pm), unique(unlist(exclusions_list)))
  pm_reduced <- pm[poi_valid]
  all_ids <- names(pm_reduced)
  chunks <- split(all_ids, ceiling(seq_along(all_ids) / chunk_size))
  
  df_output_all <- list()
  
  for (chunk in chunks) {
    pm_chunk <- pm_reduced[chunk]
    
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    clusterExport(cl, varlist = c("calculate_LR", "am", "exclusions_list"), envir = environment())
    clusterEvalQ(cl, {
      library(pedtools)
      library(forrel)
    })
    
    results <- foreach(poi_id = names(pm_chunk), .combine  = rbind, .inorder  = TRUE) %dopar% {
      poi <- pm_chunk[[poi_id]]
      
      df_temp <- setNames(rep(NA_character_, length(am)), names(am))
      for (j in seq_along(am)) {
        fam_id <- names(am)[j]
        ref <- am[[j]]
        
        if (poi_id %in% exclusions_list[[fam_id]]) {
          
          res <- calculate_LR(ref, poi)
          
          df_temp[fam_id] <- paste0(
            res[[1]], " (", res[[2]], ") / ",
            ifelse(is.na(res[[3]]), "NA", res[[3]]),
            " (", ifelse(is.na(res[[4]]), "NA", res[[4]]), ")"
          )
        }
      }
      rm(poi, res)
      gc()
      
      df_temp <- as.data.frame(t(df_temp), stringsAsFactors = FALSE)
      
      rownames(df_temp) <- poi_id
      df_temp
    }
    
    stopCluster(cl)
    
    df_chunk <- as.data.frame(results)
    df_chunk$POI <- rownames(df_chunk)
    df_chunk$Sex <- getSexID(pm_chunk)$Sex
    df_output_all[[length(df_output_all) + 1]] <- df_chunk
  }
  
  df_output <- do.call(rbind, df_output_all)
  df_output <- df_output[, c("POI", "Sex", setdiff(names(df_output), c("POI", "Sex")))]
  
  elapsed_min <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  message(sprintf("Elapsed time: %.3f minutes", elapsed_min))
  
  return(df_output)
}
