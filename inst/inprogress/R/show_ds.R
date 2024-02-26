#' Data summary
#'
#' Generates data summary of input dataframe. Always prints and returns summary
#' as a tibble
#'
#' @param .df Input datframe to generate summary on
#' @return Dataframe displaying summary row for each column of original data.
#' Displays variable name, datatype, number of distinct values, number of
#' missing values, number of non-missing values, univariate statistics if
#' numeric, and a small sample of values
#' @export
show_ds <- function(.df) {
  get_ds <- function(vector, column_name) {
    if (is.numeric(vector)) {
      summary_list <- list(
        Variable = column_name,
        DataType = typeof(vector),
        DistinctValues = dplyr::n_distinct(vector, na.rm = F),
        NAValues = sum(is.na(vector), na.rm = T),
        ValidValues = sum(!is.na(vector), na.rm = T),
        Mean = round(mean(vector, na.rm = T), 3),
        Sum = round(sum(vector, na.rm = T), 3),
        Min = as.character(round(min(vector, na.rm = T), 3)),
        Max = as.character(round(max(vector, na.rm = T), 3)),
        Obs = length(vector),
        SampleValues = stringr::str_trunc(paste(vector, collapse = ","),
                                          20, "right"))
    } else if (is.factor(vector)) {
      summary_list <- list(
        Variable = column_name,
        DataType = "factor",
        DistinctValues = dplyr::n_distinct(vector, na.rm = F),
        NAValues = sum(is.na(vector), na.rm = T),
        ValidValues = sum(!is.na(vector), na.rm = T),
        Mean = NA_real_,
        Sum = NA_real_,
        Min = as.character(head(levels(vector), 1)),
        Max = as.character(tail(levels(vector), 1)),
        Obs = length(vector),
        SampleValues = stringr::str_trunc(paste(vector, collapse = ","),
                                          20, "right"))
      } else {
      summary_list <- list(
        Variable = column_name,
        DataType = if_else(typeof(vector) == "double",
                           class(vector)[1], typeof(vector)),
        DistinctValues = dplyr::n_distinct(vector, na.rm = F),
        NAValues = sum(is.na(vector), na.rm = T),
        ValidValues = sum(!is.na(vector), na.rm = T),
        Mean = NA_real_,
        Sum = NA_real_,
        Min = as.character(min(vector, na.rm = T)),
        Max = as.character(max(vector, na.rm = T)),
        Obs = length(vector),
        SampleValues = stringr::str_trunc(paste(vector, collapse = ","),
                                          20, "right"))
    }
    return(summary_list)
  }

  output_df <- purrr::imap_dfr(.df, ~get_ds(vector = .x, column_name = .y))
  print(output_df)
  invisible(output_df)
}
