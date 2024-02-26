#' Transforms a categorical variable of n-levels into n-1 dummy variables
#' Notably distinguished from the functions available in the fastDummies package
#' as it is designed for use with a factor dataframe and label dataframe for
#' metadata storage (See metaSQL for details)
#' If these metadata structures are not present in the function call will mostly
#' replicate the behavior found in fastDummies
#' @import magrittr
#' @import rlang
#' @param df A tibble containing the dummy we are attempting to coerce OR a list of
#' the df, factordb, and labeldb. Used for easy piping
#' @param dummy_vector A character vector containing the name of the categorical variable
#' we are attempting to coerce
#' @param factordb Optional dataframe storing original numeric and verbose factors
#' @param labeldb Optional dataframe storing verbose labels
#' @param na_drop Will not create NA dummy variable option if TRUE
#' @param original_drop Will drop original categorical variable if TRUE
#' @return List of 3 dataframes. First element is the original dataframe but with
#' dummy variable columns created based on the parameters specified. Second element
#' is a new factor database based on the metadata provided, used to append to
#' pre-existing factors database. Third element is a new label database based on
#' the metadata provided, used to append to pre-existing factors database. New metadata
#' dataframe will be empty tibbles if no metadata was provided/could be produced. 4th
#' element is a vector of new variables
dummy_replace <- function(.df,
                          .dummy_vector,
                          .factordb = NULL,
                          .labeldb = NULL,
                          .na_drop = T,
                          .original_drop = T){
  # Storing original names
  original_colnames <- names(.df)
  # Creating empty vector to store named vector for new dummied variables
  new_dummies.vector <- character()

  # Creating empty metadata structures if none provided
  if(is.null(factordb)){
    factordb <- dplyr::tibble(
      "Question" = character(),
      "NumVal" = integer(),
      "CharVal" = character(),
      "Val" = character()
    )
  }
  if(is.null(labeldb)){
    labeldb <- dplyr::tibble(
      "variable_name" = character(),
      "variable_label" = character()
    ) }
  # Looping through all dummy_vars in dummy_vector
  # Adding dummy columns according to the specifications and adjusting df according
  # to specifications above
  # Also adjusting factor and label db
  for(dummy_var in dummy_vector){
    # Store the variable to coerce as dummy as a vector
    original.vector <- df[, dummy_var]
    # # Storing all unique values, excluding NA if na_drop = T
    # if(na_drop == T){
    #   original_unique.vector <- unique(original.vector)
    #   original_unique.vector <- original_unique.vector[!is.na(original_unique.vector)]
    # } else{ original_unique.vector <- unique(original_vector) }

    original_unique.vector <- sort(unique(original.vector), na.last = T)
    if(na_drop){
      original_unique.vector <- original_unique.vector[!is.na(original_unique.vector)]
    }
    # Creating dummy columns based on original vector by looping over all unique
    dummied.columns <- purrr::map_dfc(original_unique.vector, function(dummy_val){
      # Special case for creating "NA" dummy variable
      if(is.na(dummy_val)){
        dummied.vector <- dplyr::if_else(is.na(original.vector), 1L, 0L)
      } else{ dummied.vector <- dplyr::case_when(
        # If equal to current dummy value, return 1
        original.vector == dummy_val ~ 1L,
        # If not equal but implicitly not missing, return 0
        original.vector != dummy_val ~ 0L,
        # If missing and don't want missing column, return NA_integer_
        is.na(original.vector) & na_drop == T ~ NA_integer_,
        # If missing and want missing column, return 0
        is.na(original.vector) & na_drop == F ~ 0L
      )}

      # Return as named tibble column
      dummied.column <- dplyr::tibble(dummied.vector)
      names(dummied.column) <- paste0(dummy_var, "_", dummy_val)
      return(dummied.column)
    })

    # Attaching dummied columns to df
    dummied.tibble_unordered <- dplyr::bind_cols(df, dummied.columns)

    # Attach dummied columns to original in original column location
    dummy_var.index <- which(names(df) == dummy_var)
    # Collect names to the left of this index, inclusive to preserve order
    left.names <- names(df)[1:dummy_var.index]
    # Adding the dummied names
    new_left.names <- c(left.names, names(dummied.columns))

    # Reordering final.tibble_unordered based on names
    dummied.tibble <- dplyr::relocate(dummied.tibble_unordered,
                                    # Using !!! to unpack the vector to individual args
                                    !!!new_left.names)
    if(original_drop == T){
      dummied.tibble <- dummied.tibble[-dummy_var.index]
    }

    df <- dummied.tibble

    # Appending new dummied columns to new_dummies.vector
    append_dummies.vector <- names(dummied.columns)
    names(append_dummies.vector) <- rep_len(dummy_var, length(append_dummies.vector))
    new_dummies.vector <- c(new_dummies.vector, append_dummies.vector)
    # Creating new factors, treating all as simple yes no binary
    new_factordb <- purrr::map_dfr(
      names(dummied.columns), function(name){
        dplyr::tibble("Question" = name,
                      "NumVal" = 0:1,
                      "CharVal" = c("No", "Yes"),
                      "Val" = c("0: No", "1: Yes"))
      })

    # Creating new labels if possible
    if(dummy_var %in% labeldb$variable_name){
      original_label <- labeldb$variable_label[labeldb$variable_name == dummy_var]
      label_component <- rep_along(original_unique.vector, original_label)
    } else{
      label_component <- rep_along(original_unique.vector, dummy_var)
    }

    if(dummy_var %in% factordb$Question){
      factordb_filtered <- dplyr::filter(factordb, Question == dummy_var)
      factor_component <- purrr::map_chr(original_unique.vector, function(numval){
        dplyr::case_when(
          numval %in% factordb_filtered$NumVal ~
            paste0("_", factordb_filtered$Val[factordb_filtered$NumVal == numval]),
          is.na(numval) ~ "_NA: Missing",
          TRUE ~ paste0("_", numval)
        )[1]  } )
    } else{
      # Otherwise use only numeric value
      factor_component <- paste0("_", original_unique.vector)
    }
    new_labeldb <- dplyr::tibble(
      "variable_name" = names(dummied.columns),
      "variable_label" = paste0(label_component, factor_component) )
      labeldb <- dplyr::bind_rows(labeldb, new_labeldb)
      factordb <- dplyr::bind_rows(factordb, new_factordb)
  }

  # Return as a list
  # # Dataframe with new dummied variables
  # # Updated factor db
  # # Updated label db
  final.list <- list(
    "dummied_df" = df,
    "new_labeldb" = labeldb,
    "new_factordb" = factordb,
    "new_dummies" = new_dummies.vector
  )
  return(final.list)
}
