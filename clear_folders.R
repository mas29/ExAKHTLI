#' Clear Folders
#'
#' Clears folders in media folder of Shiny app
#' @export
clear_folders <- function() {
  # Remove all files in www directory of shiny "explore" app
  if(length(list.files("shiny_scripts/explore/www/")) != 0) {
    files_to_delete <- paste("shiny_scripts/explore/www/", list.files("shiny_scripts/explore/www/"), sep="")
    file.remove(files_to_delete)
  }

  # Remove all files in Output directory
  if(length(list.files("Output/")) != 0) {
    files_to_delete <- paste("Output/", list.files("Output/"), sep="")
    file.remove(files_to_delete)
  }
}

