#' Run parallel LR calculation for a set of POIs (pm) against pedigrees (am)
#'
#' @param pm A list of POIs (pedigree objects) named by ID
#' @param am A list of family pedigrees named by ID
#' @param ncores Number of parallel workers (default = detectCores() - 1)
#' @return A data.frame of LR values (rows = pm, cols = am)
#' @export
ParallelForensicLR <- function(pm, am, ncores = parallel::detectCores() - 1) {
  if (length(pm) == 0 || length(am) == 0) {
    stop("Must supply at least one POI (pm) and one pedigree (am).")
  }
  
  # start timer
  start_time <- Sys.time()
  
  # set up cluster
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  clusterExport(cl, varlist = c("calculate_LR", "combine_dataframes"))
  parallel::clusterEvalQ(cl, {
    library(pedtools)
    library(forrel)
  })
  
  # perform LR calculations in parallel
  results <- foreach(poi = pm, .combine = rbind) %dopar% {
    # one‐row data.frame, columns = names(am)
    df_temp <- setNames(
      as.data.frame(matrix(nrow = 1, ncol = length(am)), stringsAsFactors = FALSE),
      names(am)
    )
    for (j in seq_along(am)) {
      df_temp[1, j] <- calculate_LR(am[[j]], poi)
    }
    df_temp
  }
  
  # tear down
  stopCluster(cl)
  
  # assemble output
  df_output <- dplyr::bind_rows(results)
  rownames(df_output) <- names(pm)

  output = combine_dataframes(pm, df_output)
  
  # report timing (minutes only)
  elapsed_min <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  message(sprintf("Elapsed time: %.3f minutes", elapsed_min))
  
  return(output)
}
