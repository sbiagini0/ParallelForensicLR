############################################################
# File: R/helpers.R
############################################################

# Function: Automatically installs (if missing) and loads the R packages specified in the input list `pkgs`.
# pkgs: A character vector of package names to be checked, installed if necessary, and loaded.
# Returns: Invisibly returns the result of loading the packages; primarily used for its side effect of ensuring
#          required packages are available in the current R session.
load_packages <- function(pkgs) {
  installed_packages <- pkgs %in% rownames(installed.packages())
  if (any(!installed_packages)) {
    install.packages(pkgs[!installed_packages])
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

############################################################
# Function to read a Familias file and extract necessary information
# file: Path to the .fam file containing family pedigrees and profiles
# Returns: A list with processed family pedigrees and profiles for analysis
read_famfile <- function(file) {
  tryCatch({
    pedFamilias::readFam(file, useDVI = T, simplify1 = T, deduplicate = T, prefixAdded = "EXTRAA", verbose = F)
  }, error = function(e) {
    message("Error reading Familias file: ", e$message)
    NULL
  })
}

############################################################
# Function to obtain the sex (M/F) of each individual in the pedigree
# pm: List of pedigrees with individual data
# Returns: A dataframe with individual IDs and their respective sex (M or F)
getSexID <- function(pm) {
  getSx <- function(s) ifelse(s == 1, "M", ifelse(s == 2, "F", NA))
  map_dfr(pm, ~ data.frame(ID = .x$ID, Sex = getSx(getSex(.x, ids = .x$ID, named = F))))
}

############################################################
# Function to check if a pedigree is connected or contains disconnected parts
# mpi: List of pedigrees (input)
# Returns: Logical TRUE if pedigree is connected, FALSE if there are disconnected parts
connectedPed <- function(mpi) {
  #Pedigree not connected
  result <- c()
  for (i in 2:length(mpi)){
    pedigree <- is.pedList(mpi[[i]])
    if(pedigree == F){
      res = paste0(names(mpi[i]), " - ", pedigree)
      result <- c(result, res)
    }
    return(result)
  }
}

############################################################
# Function to format pedigree data in mpi by removing unwanted list, reordering pedigrees, and sorting genotypes
# mpi: List containing pedigree data, including family pedigrees and unidentified persons
# Returns: A formatted mpi list with standardized pedigree and genotype ordering
pedFormat <- function(mpi) {
  #Deleting "_comp2"
  mpi[2:length(mpi)] <- map(mpi[2:length(mpi)], ~.x[['_comp1']])
  
  #Ordering
  for (i in 2:length(mpi)) {
    mpi[[i]] = reorderPed(mpi[[i]])
    mpi[[i]] = sortGenotypes(mpi[[i]])
  }
  
  for (j in seq_along(mpi$`Unidentified persons`)) {
    mpi$`Unidentified persons`[[j]] = sortGenotypes(mpi$`Unidentified persons`[[j]])
  }
  
  return(mpi)
}

############################################################
# Function to assign unique family IDs in mpi for both families and unidentified persons
# mpi: List containing pedigree data with families and unidentified persons
# Returns: A mpi list with updated family and individual identifiers
pedNames <- function(mpi) {
  #AM
  for (i in 2:length(mpi)){
    famid(mpi[[i]]) <- names(mpi[i])
  }
  
  #PM
  for (i in names(mpi[1]$`Unidentified persons`)) {
    famid(mpi$`Unidentified persons`[[i]]) = mpi$`Unidentified persons`[[i]]$ID
    new_name <- mpi$`Unidentified persons`[[i]]$ID
    names(mpi$`Unidentified persons`)[names(mpi$`Unidentified persons`) == i] <- new_name
  }
  return(mpi)
}

############################################################
# Function to perform a Mendelian consistency check on pedigrees, removing inconsistencies
# mpi: List of pedigrees for Mendelian check
# Returns: A mpi list with pedigrees checked and inconsistencies removed
mendelian <- function(mpi){
  for (i in 2:length(mpi)) {
    print(names(mpi)[i])
    mpi[[i]] <- mendelianCheck(mpi[[i]], remove = T, verbose = T)
  }
  return(mpi)
}

############################################################
# Function to remove the mutation model from pedigrees and unidentified persons by setting mutation rate to 0
# mpi: List containing pedigrees and unidentified persons
# Returns: A mpi list with mutation models removed (rate set to 0)
removeMutModel <- function(mpi){
  for (i in 2:length(mpi)){
    mpi[[i]] = setMutmod(mpi[[i]], model = "equal", rate = 0, update = T)
  }
  
  for (i in 1:length(mpi$'Unidentified persons')){
    mpi[[1]][[i]] = setMutmod(mpi[[1]][[i]], model = "equal", rate = 0, update = T)
  }
  return(mpi)
}

############################################################
# Function to set a mutation model for pedigrees and unidentified persons with customizable parameters
# mpi: List containing pedigrees and unidentified persons
# model: Mutation model to use ("stepwise" or "equal")
# rate: Mutation rate, as a list for "stepwise" (e.g., list(female = 0.01, male = 0.02)) or a single value for "equal"
# range: Range parameter for the "stepwise" mutation model (ignored if model is "equal")
# rate2: Secondary rate parameter for the "stepwise" mutation model (ignored if model is "equal")
# Returns: A mpi list with mutation models set according to specified parameters
setMutModel <- function(mpi, model = "stepwise", rate = list(female = 0.01, male = 0.02), range = 0.1, rate2 = 1e-6) {
  for (i in 2:length(mpi)) {
    if (model == "stepwise") {
      mpi[[i]] = setMutmod(mpi[[i]], model = model, rate = rate, range = range, rate2 = rate2, update = T)
    } else if (model == "equal") {
      mpi[[i]] = setMutmod(mpi[[i]], model = model, rate = rate, update = T)
    }
  }
  
  for (i in 1:length(mpi$'Unidentified persons')) {
    if (model == "stepwise") {
      mpi$`Unidentified persons`[[i]] = setMutmod(mpi$`Unidentified persons`[[i]], model = model, rate = rate, range = range, rate2 = rate2, update = T)
    } else if (model == "equal") {
      mpi$`Unidentified persons`[[i]] = setMutmod(mpi$`Unidentified persons`[[i]], model = model, rate = rate, update = T)
    }
  }
  
  return(mpi)
}

############################################################
# Function: Computes an exclusion matrix between a list of Persons of Interest (POIs) and a list of reference families.
#           The function runs in parallel to improve performance and determines which POIs have an acceptable number
#           of inconsistencies (??? maxIncomp) with each reference family.
# pm: A named list of POI pedigrees.
# am: A list of reference family pedigrees.
# maxIncomp: Maximum number of allowed marker inconsistencies for a POI to be considered compatible with a family (default: 3).
# ncores: Number of CPU cores to use for parallel processing (default: all available cores minus one).
# Returns: A named list where each element corresponds to a family (by FAMID) and contains the vector of POI names that
#          are not excluded from that family based on the defined inconsistency threshold.
computeExclusions <- function(pm, am, maxIncomp = 3, ncores = detectCores() - 1) {
  # Initialize parallel backend
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  clusterEvalQ(cl, {
    library(pedtools)
    library(forrel)
    library(dvir)
  })
  
  # Compute exclusion matrix in parallel
  results <- foreach(poi = pm, .combine = rbind, .inorder  = TRUE) %dopar% {
    df_temp <- data.frame(matrix(nrow = 1, ncol = 0))
    colnames(df_temp) <- character(0)
    
    for (j in seq_along(am)) {
      ref <- am[[j]]
      
      intersec_markers <- intersect(pedtools::name(ref), pedtools::name(poi))
      ref_sel <- pedtools::selectMarkers(ref, markers = intersec_markers)
      poi_sel <- pedtools::selectMarkers(poi, markers = intersec_markers)
      
      dviData <- dvir::dviData(poi_sel, ref_sel, missing = "Missing person", generatePairings = TRUE)
      exclusion_matrix <- dvir::exclusionMatrix(dviData, pairings = NULL, ignoreSex = TRUE)
      colnames(exclusion_matrix) <- ref$FAMID
      
      df_temp <- cbind(df_temp, exclusion_matrix)
    }
    
    df_temp
    
  }
  
  stopCluster(cl)
  
  # Post-process exclusion matrix
  exclusions <- bind_rows(results)
  exclusions$POI <- names(pm)
  exclusions <- exclusions[, c("POI", setdiff(names(exclusions), "POI"))]
  
  fams <- sapply(am, `[[`, "FAMID")
  exclusions_list <- setNames(
    lapply(fams, function(fam_id) {
      vals <- exclusions[[fam_id]]
      sel  <- !is.na(vals) & vals >= 0 & vals <= maxIncomp
      exclusions$POI[sel]
    }),
    fams
  )
  
  return(exclusions_list)
}

############################################################
# Function: Computes the kinship Likelihood Ratio (LR) between a reference
#           pedigree (with a placeholder "Missing person") and a Person of
#           Interest (POI) profile.  It reports the overall LR, the number of
#           informative markers, and-after removing any Mendelian-incompatible
#           markers-the "mismatch" LR and its marker count.
# ref: A pedigree object for the reference family that includes the
#      placeholder individual "Missing person".
# poi: A pedigree object representing the Person of Interest.
# Returns: A list of four elements  
#          1. LR_total               - Total LR across all shared markers  
#          2. nMarkers               - Count of informative overlapping markers  
#          3. LR_mismatch            - Total LR after excluding incompatible
#                                      markers (NA if none were incompatible)  
#          4. nMarkers_mismatch      - Count of informative markers after
#                                      exclusion (NA if none were incompatible)
calculate_LR <- function(ref, poi) {
  intersec_markers <- intersect(pedtools::name(ref), pedtools::name(poi))
  ref <- pedtools::selectMarkers(ref, markers = intersec_markers)
  poi <- pedtools::selectMarkers(poi, markers = intersec_markers)
  
  lr = forrel::missingPersonLR(
    reference = ref,
    poi = poi,
    missing = "Missing person",
    markers = intersec_markers,
    verbose = F
  )
  
  lr_per_marker_filtered <- lr[["LRperMarker"]][intersec_markers]
  
  overlapping_markers <- sum(abs(lr_per_marker_filtered - 1) > 1e-6)
  
  lr_total <- lr[["LRtotal"]][["H1:H2"]]
  lr_total <- signif(lr_total, 5)
  
  exclusions_markers <- forrel::findExclusions(x = ref, id = "Missing person", candidate = poi)
  
  if (length(exclusions_markers) > 0) {
    
    dif = setdiff(intersec_markers, exclusions_markers)
    
    lr_missmatch = forrel::missingPersonLR(
      reference = ref,
      poi = poi,
      missing = "Missing person",
      markers = dif,
      verbose = F
    )
    
    lr_per_marker_filtered <- lr_missmatch[["LRperMarker"]][dif]
    
    overlapping_markers_missmatch <- sum(abs(lr_per_marker_filtered - 1) > 1e-6)
    
    lr_missmatch <- lr_missmatch[["LRtotal"]][["H1:H2"]]
    
    lr_missmatch <- signif(lr_missmatch, 5)
    
  } else {
    lr_missmatch <- NA
    overlapping_markers_missmatch <- NA
  }
  
  return(list(lr_total, overlapping_markers, lr_missmatch, overlapping_markers_missmatch))
}

# Function to combine lr_aux with id_sex and rearrange columns
# pm: List of pedigrees with individual data for sex determination
# lr_aux: Dataframe containing LR values for each POI and family
# Returns: A reorganized dataframe with POI, Sex, and family columns
combine_dataframes <- function(pm, lr_aux) {
  id_sex <- getSexID(pm)
  lr_aux <- cbind(lr_aux, id_sex)
  lr_aux$POI <- rownames(lr_aux)
  rownames(lr_aux) <- NULL
  
  lr_aux <- lr_aux[, c("POI", "Sex", colnames(lr_aux)[1:(ncol(lr_aux) - 2)])]
  return(lr_aux)
}

############################################################
# Function to export the combined dataframe as an .xlsx file, sorted by maximum LR value
# df: The dataframe containing combined and cleaned data for export
# filename: Name base for the output .xlsx file
# Output: Saves an .xlsx file with sorted results for analysis
export_to_csv <- function(df, filename) {
  # Remove ID column if present
  if ("ID" %in% colnames(df)) {
    df <- df[, colnames(df) != "ID"]
  }
  
  # Function to extract the numeric part of an LR_total cell
  extract_lr_total <- function(x) {
    suppressWarnings(as.numeric(sub(" .*", "", x)))
  }
  
  # Select LR columns (from the third column onward)
  lr_values <- df[, 3:ncol(df), drop = FALSE]
  lr_numeric <- apply(lr_values, 2, extract_lr_total)
  
  # If there's only one LR column, ensure lr_numeric is a matrix
  if (is.null(dim(lr_numeric))) {
    lr_numeric <- matrix(lr_numeric, ncol = 1)
    colnames(lr_numeric) <- colnames(lr_values)
  }
  
  # Compute the maximum LR per row for sorting
  max_lr <- matrixStats::rowMaxs(lr_numeric, na.rm = TRUE)
  sorted_df <- df[order(max_lr, decreasing = TRUE), ]
  
  # Export to CSV
  name_csv <- paste0(filename, "_ParallelForensicLR.csv")
  write.table(sorted_df,
              file = name_csv,
              sep = ";",
              row.names = FALSE,
              col.names = TRUE,
              quote = FALSE)
  
  message("Export completed: ", name_csv)
}
