# functions to get plate and data into data frame

# LIBRARIES

require(XLConnect)

# FUNCTIONS

# function to convert incucyte data from data frame to numeric matrix, with row and column names
df_to_t_matrix <- function(dat) {
  # transpose matrix
  colnames <- dat[,2]
  rownames <- colnames(dat)[3:ncol(dat)]
  df_t <- t(dat[,c(3:ncol(dat))])
  # make matrix numeric
  mat <- matrix(as.numeric(unlist(df_t)), nrow=nrow(df_t))
  colnames(mat) <- colnames
  rownames(mat) <- rownames
  return(mat)
}

# function to parse compound key data
# @param key_filename -- file name for compound kay
parse_compound_key <- function(key_filename) {
    key_wb <- loadWorkbook(key_filename$datapath)
    key_data_raw <- readWorksheet(key_wb, sheet = getSheets(key_wb))
    key_data <- lapply(key_data_raw, function(key_dat) { return(as.character(unlist(data.frame(t(key_dat))))) } )
    return(key_data)
}

# function to parse compound info data
# @param compound_info_filename -- file name for compound info
parse_compound_info <- function(compound_info_filename) {
  compound_wb <- loadWorkbook(compound_info_filename$datapath)
  input_compound_info <- readWorksheet(compound_wb, sheet = getSheets(compound_wb))
  colnames(input_compound_info)[which(colnames(input_compound_info) == "Target.class..11Mar15.")] <- "target"
  colnames(input_compound_info)[which(colnames(input_compound_info) == "Pathway")] <- "Pathway"
  colnames(input_compound_info)[which(colnames(input_compound_info) == "Product.Name")] <- "Compound"
  return(input_compound_info)
}

# function to parse incucyte data
# @param incucyte_data_filname -- file name for incucyte info
# @param key_data -- compound key 
parse_incucyte_data <- function(incucyte_data_filename, key_data) {
  incucyte_wb <- loadWorkbook(incucyte_data_filename$datapath)
  incucyte_data_raw <- readWorksheet(incucyte_wb, sheet = getSheets(incucyte_wb))
  last_empty_row <- which(incucyte_data_raw[[1]][,1] == "Date Time")
  incucyte_data <- lapply(incucyte_data_raw, function(dat) { return(dat[-(1:last_empty_row),]) })
  
  # plate numbers corresponding to each sheet of the incucyte data
  plate_numbers_in_incucyte <- as.numeric(gsub(".*(\\d+)", "\\1", getSheets(incucyte_wb)))
  
  # for each sheet of the incucyte data file, add correct key names as column names
  for (i in 1:length(incucyte_data)) { 
    colnames(incucyte_data[[i]])[3:length(colnames(incucyte_data[[i]]))] <- key_data[[plate_numbers_in_incucyte[i]]]  
  }
  
  # convert to transposed numeric matrix
  input_matrices <- lapply(incucyte_data, df_to_t_matrix)
  names(input_matrices) <- gsub("(.*)(\\d+)", "\\1:\\2", names(input_matrices))

  return(input_matrices)
}

# function to get plate positions
# @param incucyte_data_filname -- file name for incucyte info
get_plate_positions <- function(incucyte_data_filename) {
  incucyte_wb <- loadWorkbook(incucyte_data_filename$datapath)
  incucyte_data_raw <- readWorksheet(incucyte_wb, sheet = getSheets(incucyte_wb))
  last_empty_row <- which(incucyte_data_raw[[1]][,1] == "Date Time")
  plate_positions <- lapply(incucyte_data_raw, function(dat) { return(unlist(dat[last_empty_row,c(3:ncol(dat))])) })
  
  return(plate_positions)
}

