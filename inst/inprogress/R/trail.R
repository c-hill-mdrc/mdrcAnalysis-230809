trail <- function () {
  # The purpose of this function is to create a data frame with session information
  # This information can be joined to extracts for documentation purposes
  # Data include:
  #      - filepath (if script run in RStudio)
  #      - username 
  #      - datetime with timezone 
  #        (Note: timezone will be EST because that's where MDRC Servers Run)
  #      - R Version
  #      - packages and version numbers
  
  # Create named vector with each package and version
  pkg_vrs1 <- sapply(sessionInfo()[[7]], function(x) x$Version)
  
  # Flatten vector to single string with package names and versions
  pkg_vrs2 <- paste(names(pkg_vrs1), pkg_vrs1, sep = " ", collapse = ", ")
  
  # Get file path
  ## Only works in RStudio, otherwise returns blank
  t <- try(rstudioapi::getActiveDocumentContext()$path)
  
  # Create data frame of trail information
  df <- data.frame(filepath = t[[1]],
                   user     = Sys.info()[[7]],
                   datetime = format(Sys.time(), format = "%F %R %Z"),
                   RVersion = sessionInfo()[[1]]$version.string,
                   Packages = pkg_vrs2
                   )
  
  return(df)
}

# Test of data frame
# This can be assigned to object
# traildf <- trail()
# This can be joined to another data frame
# dplyr::full_join(iris, df, by = character())
