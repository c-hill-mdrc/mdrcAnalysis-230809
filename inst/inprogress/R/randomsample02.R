#' Generate Random Sample
#'
#' Creates a random sample of given size or proportion
#'
#' @param .df Source dataset
#' @param .select_vars Optional character vector, Variables to include in sample
#'  dataframe Defaults to all
#' @param .grouping_vars Optional character vector,will take sample of
#'  specified size on a per group basis. Defaults to not gorup
#' @param .SampSize Numeric, size of sample. If between 0 and 1, will treat it
#' as a proprortion and the sample size will be the closest integer to that
#' proportion of the total number of rows in the dataframe
#' @param .Replacement Logical, Sample take with replacment?
#' @param .Flag Logical, Sample has a column added called FlagVar set equal to 1
#' @returns Sample dataframe
#' @export
randomsample02 <- function(.df,
                           .select_vars,
                           .grouping_vars = NULL,
                           .SampSize = 50,
                           .seed = NULL,
                           .Replacement = FALSE,
                           .Flag = FALSE) {
  # Setting seed
  set.seed(.seed)
  # If sample size is greater than 1, then leave as is
  # If sample size is less than 1, then calculate as rounded integer proportion
  # of number of rows in input dataframe
  .SampSize <- ifelse(.SampSize > 1,
                      .SampSize,
                      round((.SampSize) * nrow(.df))
  )

  .df <- .df %>%
    # Group by .grouping.vars if exists
    dplyr::group_by(dplyr::across( {{.grouping_vars}} )) %>%
    # Take sample of desired sample size, with replacement if explicitly set
    dplyr::slice_sample(n = .SampSize, replace = .Replacement) %>%
    dplyr::ungroup()
  # If .Flag is T, then add a Flag column, otherwise set column with name
  # FlagVar and set value equal 1
  if (.Flag) {
    .df <- .df %>%
      dplyr::mutate(FlagVar = 1)
  }

  # If .select.vars are not explicitly set, keep all columns
  # otherwise select only explicitly set variables
  if (!missing(.select_vars)) {
    .df <- .df %>%
      dplyr::select({{.select_vars}}, dplyr::any_of("FlagVar"))
  }
  return(.df)
}




