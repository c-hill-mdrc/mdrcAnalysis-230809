#' Searches for missing proxies across a dataset
#' @description Missing data is a frequent occurrence; however, missing data
#' frequently appears in forms more complex than just an NA. Often descriptions
#' of the cause of the missing data are present which are represented in the
#' data with a string or number. Once some of these proxies are identified,
#' they can be searched through across a dataset to see other potential missings
#' @param .data Dataframe
#' @param  .factordb Optional, factor dataframe listing all factors
#' @param .missing_proxy Vector of missing proxy values
#'
missing_check <- function(.data, .factordb = NULL, .missing_proxy = c(77, 88, 99)){

  check_full <- purrr::imap_dfr(.data, function(.column, .colname){
    # Has actual NA value
    missing_actual <- any(is.na(.column))

    missing_tibble <- dplyr::tibble(
      ColName = .colname,
      UniqueVal = dplyr::n_distinct(column, na.rm = T),
      Actual_NA = missing_actual
    )

    if(!is.null(.factordb)){
      if(.colname %in% .factordb$Question){
        # Take factordb that matches question and has values in data
        factordb_subset <- .factordb %>%
          dplyr::filter(Question == .colname) %>%
          dplyr::filter(NumVal %in% unique(.column))
        # Take factordb where numval matches proxy
        missing_proxied <- factordb_subset %>%
          dplyr::filter(NumVal %in% .missing_proxy) %>%
          dplyr::pull(Val)
          # Take factordb where string matches common missing words
        missing_factors <- .factordb %>%
          dplyr::filter(Question == .colname) %>%
          filter(NumVal %in% unique(.column)) %>%
          filter(stringr::str_detect(
            Val,
            # Looks for strings commonly used to prepresent missing data
            stringr::regex(
              "n\\/{0,1}a|missing|not answered|skipped|routing error",
              ignore_case = T)
          )) %>%
          dplyr::pull(Val)
        missing_tibble <- missing_tibble %>%
          dplyr::mutate(
            NA_proxy = paste(missing_proxied, collapse = ", "),
            NA_factors = paste(missing_factors, collapse = ", ")
          )
      } else{
        # Displaying all .missing_proxy that appear in actual data if not
        # retrievable from factors
        missing_proxied <- .missing_proxy[match(.column, .missing_proxy)] %>%
          unique() %>%
          sort()
        missing_tibble <- missing_tibble %>%
          dplyr::mutate(NA_proxy = paste(missing_proxied, collapse = ", "))
      }
    } else{
      # Displaying all .missing_proxy that appear in actual data if not
      # retrievable from factors
      missing_proxied <- .missing_proxy[match(.column, .missing_proxy)] %>%
        unique() %>%
        sort()
      missing_tibble <- missing_tibble %>%
        dplyr::mutate(NA_proxy = paste(missing_proxied, collapse = ", "))
    }
    return(missing_tibble)
  })

  # Creating filtered check where problem spotted
  names_check <- check_full %>%
    purrr::discard(. %in% c("ColName", "UniqueVal", "Actual_NA"))
  check_full <- check_full %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(names_check),
                                ~dplyr::if_else(. == "", NA_character_, .)))
  check_filtered <- check_full %>%
    dplyr::filter(dplyr::if_any(dplyr::all_of(names_check), ~!is.na(.)))
  print(check_filtered)
  invisible(list(
    "Missing Check" = check_full,
    "Missing Filtered" = check_filtered
  ))
}
