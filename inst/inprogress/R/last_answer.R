#' Identifies frequency each variable is the "last" value for a row
#' @param .df A tibble to find the rowwise last non-NA
#' @return Tibble displaying frequency table of each original variable comparing
#' amount of times a row has the variable as the last value
#' @export
last_distribution <- function(.df){
  # Ensuring that .df is a tibble with no row names
  df <- tibble::tibble(.df)
  # Looping through each row and identifying last column answered
  output <-   purrr::pmap_dfr(df, function(...){
    row_list <- list(...)
    # Removing missing values
    row_list_cleaned <- row_list %>%
      purrr::discard(~is.na(.))
    # returning last name remaining
    last_answer <- tail(names(row_list_cleaned), 1)
    return(tibble::tibble("Last Answer" = last_answer))
  }) %>%
    # Transforming into frequency table
    dplyr::count(`Last Answer`)

  return(output)
}
