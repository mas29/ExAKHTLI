# Script to convert images for a selected compound from TIFF to JPEG

# Function to get the folder structure of the archive directory for year/month/days levels
# param archive_dir -- archive directory
get_experiment_days <- function(archive_dir) {
  # Experiment years/months
  expt_yrs_months <- list.files(path = archive_dir)
  expt_days <- vector()
  expt_hrs_mins <- list()

  # For each year/month, get the experiment days
  for (i in 1:length(expt_yrs_months)) {
    new_expt_days_list <- list.files(path = paste(archive_dir, "/", expt_yrs_months[i], sep=""))
    new_expt_days <- vector()
    for(j in 1:length(new_expt_days_list)) {
      new_expt_days <- c(new_expt_days, paste(expt_yrs_months[i], "/", new_expt_days_list[j], sep=""))
    }
    expt_days <- c(expt_days, new_expt_days)
  }

  return(expt_days)
}

# Function to get the folder structure of the archive directory for hours/mins level
# param archive_dir -- archive directory
get_experiment_hours_and_mins <- function(expt_days, archive_dir) {
  # For each year/month/day, get the hours/mins for the image data
  expt_hrs_mins <- list()
  for (j in 1:length(expt_days)) {
    expt_new_hrs_mins <- list.files(path = paste(archive_dir,expt_days[j],sep=""))
    expt_hrs_mins[[j]] <- expt_new_hrs_mins
  }
  expt_hrs_mins <- as.list(setNames(expt_hrs_mins, expt_days))
  return(expt_hrs_mins)
}

# Function to get plate numbers as they appear in the archive (e.g. folders 282, 283, 284, 285, 286  correspond to plates 1, 2, 3, 4, 5)
# param archive_dir -- archive directory
get_plate_nums_in_archive_dir <- function(expt_days, expt_hrs_mins, archive_dir) {
  plate_nums_etc <- list.files(path = paste(archive_dir,expt_days[1],"/",expt_hrs_mins[[1]][1],sep=""))
  plate_nums <- plate_nums_etc[grep("^[[:digit:]]*$", plate_nums_etc)]
  return(plate_nums)
}

# Function to get the images corresponding to the compound of interest
# param compound -- compound of interest
# param archive_dir -- archive directory
get_images <- function(compound, archive_dir, data_wide, image_types, time_elapsed) {

  # Get info for this compound
  expt_days <- get_experiment_days(archive_dir) # Archive structure
  expt_hrs_mins <- get_experiment_hours_and_mins(expt_days, archive_dir) # Archive structure
  plate_nums <- get_plate_nums_in_archive_dir(expt_days, expt_hrs_mins, archive_dir) # Get names of plates as they exist in the archive
  position <- data_wide[data_wide$Compound == compound,]$Position[1] # Position of compound in plate
  plate <- data_wide[data_wide$Compound == compound,]$Plate[1] # Plate of compound

  # Check if images for this compound has already been loaded
  if (!file.exists(paste("www/Plate",plate, "_Position", position, "_Image", image_types[1], "_t_", as.character(time_elapsed[1]), ".jpeg", sep=""))){

    # Set suffix for the different image types
    suffixes <- list()
    for (i in 1:length(image_types)) {
      new_suffix <- paste("-1-",image_types[i],".tif",sep="")
      suffixes[[i]] <- new_suffix
    }
    suffixes <- as.list(setNames(suffixes, image_types))

    # Remove all files in www directory of shiny "explore" app
    if(length(list.files("www/")) != 0) {
      files_to_delete <- paste("www/", list.files("www/"), sep="")
      file.remove(files_to_delete)
    }

    # Now convert and add the images to the www directory of shiny "explore" app
    count <- 1 # To keep track of which timepoint we're at
    for (i in 1:length(expt_days)) { # For each experiment day
      for (j in 1:length(expt_hrs_mins[[i]])) { # For each hour/minute of image capture in that experiment day
        plate_name_in_archive <- plate_nums[as.numeric(plate)]

        for (x in 1:length(image_types)) {
          image_dir <- paste(archive_dir,expt_days[i],"/",expt_hrs_mins[[i]][j],"/",plate_name_in_archive,"/",toupper(position),suffixes[[x]],sep="")
          img <- suppressWarnings(readTIFF(image_dir, native=TRUE))
          print((img*1.2)[1:3])
          print((img)[1:3])
          writeJPEG(img, target = paste("www/Plate",plate, "_Position", position, "_Image", image_types[x], "_t_", as.character(time_elapsed[count]), ".jpeg", sep=""), quality = 1)
        }

        count <- count + 1
      }
    }
  }
}

