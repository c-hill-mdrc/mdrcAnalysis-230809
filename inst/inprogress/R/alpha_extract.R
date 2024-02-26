#' Cronbach Alphas Extract
#'
#' Generates Cronbach alphas. Prints reliabiltiy analysis and reliability
#' if an item is dropped into a single dataframe and excel output if file name
#' and path is outputted
#'
#' @param .df Required, dataframe
#' @param .items Optional, Items to create alphas on, defaults to all variables
#' @param .labels Optional, Character vector of labels matched with items. Must
#' be of equivalent length
#' @param .use Name of correlation method to use, defaults to pairwise
#' @param .keys Optional, Character vector of items to reverse key
#' @param .file_path Optional, Path to output excel file
#' @param .file_name Optional, Name of output excel file
#' @param ... Other arguments to pass to psych::alpha
#' @returns Tibble of total alphas in first row and drop alphas after
#' @export
alpha_extract <- function(.df,
                          .items,
                          .labels,
                          .use = c("pairwise",
                                   "everything",
                                   "all.obs",
                                   "complete.obs",
                                   "na.or.complete",
                                   "pairwise.complete.obs"),
                          .keys = NA,
                          .file_path,
                          .file_name,
                          ...) {
  # Explicitly calculate maximal number of categories to pass to psych::alpha
  .use <- match.arg(.use)
  # Selecting items if provided, otherwise assuming that user passed subsetxc
  if (!missing(.items)) {
    .df <- .df %>%
      dplyr::select(dplyr::all_of(.items))
  }

  if (!missing(.labels)) {
    names(.labels) <- names(.df)
  }
  maximal_ncat <- purrr::map_int(.df, ~dplyr::n_distinct(.)) %>%
    max()
  # Calculating alphas
  alphas <- .df %>%
    psych::alpha(keys = .keys, use = .use, max = maximal_ncat, ...)

  # Extracting from alphas list
  ## Total Alphas
  total_alphas <- alphas$total %>%
    tibble::as_tibble() %>%
    dplyr::select(raw_alpha, std_alpha = std.alpha, mean, sd) %>%
    dplyr::mutate(items = "Total")
  ## Drop Alphas
  drop_alphas <- alphas$alpha.drop %>%
    tibble::as_tibble(rownames = "items") %>%
    dplyr::mutate(items = names(.df)) %>%
    dplyr::select(items, raw_alpha, std_alpha = raw_alpha)
  ## Univariate statistics
  item_univariates <- alphas$item.stats %>%
    tibble::as_tibble(rownames = "items") %>%
    dplyr::select(items, mean, sd, n)
  drop_alphas <- drop_alphas %>%
    full_join(item_univariates, by = "items")

  # Combining output
  alpha_output <- total_alphas %>%
    dplyr::bind_rows(drop_alphas)
  # Attaching labels if possible
  if (!missing(.labels)) {
    alpha_output <- alpha_output %>%
      mutate(Label = c("", .labels), .after = 1)
  }
  if (!missing(.file_path) & !missing(.file_name)) {
    create_excel(.file_path, .file_name, alpha_output)
  }

  return(alpha_output)
}
