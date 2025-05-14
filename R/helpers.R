############################################################
# File: R/helpers.R
############################################################

load_packages <- function(pkgs) {
  installed_packages <- pkgs %in% rownames(installed.packages())
  if (any(!installed_packages)) {
    install.packages(pkgs[!installed_packages])
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

# Function to read a Familias file and extract necessary information
# file: Path to the .fam file containing family pedigrees and profiles
# Returns: A list with processed family pedigrees and profiles for analysis
read_famfile <- function(file) {
  tryCatch({
    pedFamilias::readFam(file, useDVI = T, simplify1 = T, deduplicate = T, prefixAdded = "EXTRA", verbose = F)
  }, error = function(e) {
    message("Error reading Familias file: ", e$message)
    NULL
  })
}

# Function to obtain the sex (M/F) of each individual in the pedigree
# pm: List of pedigrees with individual data
# Returns: A dataframe with individual IDs and their respective sex (M or F)
getSexID <- function(pm) {
  getSx <- function(s) ifelse(s == 1, "M", ifelse(s == 2, "F", NA))
  map_dfr(pm, ~ data.frame(ID = .x$ID, Sex = getSx(getSex(.x, ids = .x$ID, named = F))))
}

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

# Function to set a mutation model for pedigrees and unidentified persons with customizable parameters
# mpi: List containing pedigrees and unidentified persons
# model: Mutation model to use ("stepwise" or "equal")
# rate: Mutation rate, as a list for "stepwise" (e.g., list(female = 0.01, male = 0.02)) or a single value for "equal"
# range: Range parameter for the "stepwise" mutation model (ignored if model is "equal")
# rate2: Secondary rate parameter for the "stepwise" mutation model (ignored if model is "equal")
# Returns: A mpi list with mutation models set according to specified parameters
settingMutModel <- function(mpi, model = "stepwise", rate = list(female = 0.01, male = 0.02), range = 0.1, rate2 = 1e-6) {
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

# Function to calculate the Likelihood Ratio (LR) between a POI and a family pedigree
# poi: Profile of the Person of Interest (POI)
# pedigree: Family pedigree data for comparison
# Returns: Calculated LR value for the POI and pedigree
calculate_LR <- function(ref, poi) {
  intersec_markers <- intersect(pedtools::name(ref), pedtools::name(poi))
  ref <- pedtools::selectMarkers(ref, markers = intersec_markers)
  poi <- pedtools::selectMarkers(poi, markers = intersec_markers)
  
  forrel::missingPersonLR(
    reference = ref,
    poi = poi,
    missing = "Missing person",
    markers = intersec_markers
  )[["LRtotal"]]
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

# Function to export the combined dataframe as an .xlsx file, sorted by maximum LR value
# df: The dataframe containing combined and cleaned data for export
# filename: Name base for the output .xlsx file
# Output: Saves an .xlsx file with sorted results for analysis
export_to_xlsx <- function(df, filename) {
  if ("ID" %in% colnames(df)) {
    df <- df[, colnames(df) != "ID"]
  }
  
  df[, 3:ncol(df)] <- lapply(df[, 3:ncol(df)], function(x) as.numeric(as.character(x)))
  
  maxdf <- rowMaxs(as.matrix(df[, 3:ncol(df)]), na.rm = T)
  sorted_df <- df[order(maxdf, decreasing = T),]
  
  name_xlsx <- paste0(filename, "_ParallelForensicLR.xlsx")
  writexl::write_xlsx(sorted_df, name_xlsx)
  
  message("Export completed: ", name_xlsx)
}