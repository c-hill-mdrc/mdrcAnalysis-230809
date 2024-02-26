#' Reading Qualtrics
#'
#' Not to be used for subsequent imports,
#' intended only for first import of data source
#'
#' Outputs a list of 3 elements
#' First element is label database
#' Second element is factor database
#' Third element is numeric data, with numbers changed if desired
#'
#' @import magrittr
#' @param .numeric_path Path to Numeric File
#' @param .character_path Path to Character File
#' @param .schema String for schema name
#' @param .table_name String for table name
#' @param .data_name String giving "dataset name" If not given, constructed from
#' .schema and .table_name
#' @param .row_skip Number of rows to skip in data before data begins
#' @param .factorvars Optional argument, vector of factor variable names
#'  If not included, factor variables taken from the variable names from
#' @param .numcoerce Optional argument, logical value indicating whether numeric
#'  should be coerced from numbers in the character value if present.
#' @return List containing metadata tibbles and coerced numeric data
#' @export
initial_read <- function(.numeric_path,
                         .character_path,
                         .schema = NA_character_,
                         .table_name = NA_character_,
                         .data_name = NULL,
                         .row_skip = 2,
                         .factorvars = NULL,
                         .numcoerce = T) {
  # Constructing dataname based on other arguments if not provided
  data_name <- dplyr::if_else(is.null(.data_name),
                              paste0("[", .schema, "].", "[", .table_name, "]"),
                              .data_name)
  # Creating factor database
  meta_list <- read_meta(.numeric_path,
                         .character_path,
                         .schema,
                         .table_name,
                         .data_name = data_name,
                         .row_skip, .factorvars,
                         .numcoerce)
  factor_db <- meta_list$FactorDatabase
  # Create labels and coerced numeric data from character data
  final_data <- read_data(.character_path,
                          .factordb = factor_db,
                          .row_skip,
                          .numcoerce)

  # Adding factordb to final exported list
  final_list <- c(meta_list, "NumericData" = list(final_data))
  return(final_list)
}

#' Qualtrics Data
#'
#' Extract character data  from Qualtrics "Character"\code{.csv} output
#' Character data is coerced to the numeric factor based on the factordb
#'
#' This function loads a Qualtrics formated character \code{.csv}
#' into a list of 2 tibbles.
#' The tibble in index 1 contains the numeric data as coded by Qualtrics.
#' It may seem silly to load the character data to create numeric data but some
#' numeric values are changed from their original in Qualtrics
#'
#' @import magrittr
#' @param .character_path Path to Character File
#' @param .data_name String giving "dataset name" If not given, constructed from
#' .schema_name and .table_name
#' @param .row_skip Number of rows to skip in data before data begins
#' @param .factorvars Optional argument, vector of factor variable names
#' If not included, factor variables taken from the variable names in
#' \code{factordb}
#' @param .numcoerce Optional argument indicating if common qualtrics string
#' manipulations should take place or if data should be stored as is.
#' Should match option creating metadata
#' @return A list containing character data coerced to factor equivalent
#' @export
read_data <- function(.character_path,
                      .factordb = NULL,
                      .row_skip = 2,
                      .numcoerce = T) {
  # Reading in character data, assuming either csv or xlsx file
  if(stringr::str_detect(.character_path, "\\.csv$")){
    character_data <- utils::read.csv(file = .character_path, check.names = F)
  } else{character_data <- openxlsx::read.xlsx(xlsxFile = .character_path,
                                               check.names = F)}
  Numeric <- tail(character_data, -(.row_skip))  %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character))

  # # Identical function exists in Create_Factors.R
  # Necessary for consistency if numcoerce is TRUE in QualtricsFactorsLib() call
  # used to generate .factordb
  CharReplace <- function(x){
    x <- stringr::str_replace_all(x, "\\{e://Field/", "\\{FIELD ") %>%
      stringr::str_replace_all( "\\$", "")
  }

  if (.numcoerce == T) {
    # Coerce all empty string entries to NA regardless of numcoerce,
    # Modify character string to match character string modifications in
    # QualtricsFactorsLib
    Numeric <- Numeric %>%
      mutate(across(where(is.character), ~na_if(.,"")),
             across(where(is.character), ~CharReplace(.x)))
  } else {
    Numeric <- Numeric %>%
      mutate(across(where(is.character), ~na_if(.,"")))
  }
  # Rewrite all character factors as their numeric version

  # Vector of factors
  if (!is.null(.factordb)) {
    factor_vars <- names(Numeric)[names(Numeric) %in% .factordb$Question]
    Numeric <- Numeric %>%
      dplyr::mutate(
        dplyr::across(
          all_of(factor_vars),
          ~as.character(
            factor(.,
                   levels = unique(.factordb$CharVal[.factordb$Question ==
                                                       dplyr::cur_column()]),
                   labels = unique(.factordb$NumVal[.factordb$Question ==
                                                      dplyr::cur_column()])
            ))))
  }
  return(Numeric)
}


#' Qualtircs Metadata
#'
#' Outputs a list of 2 dataframes
#' First dataframe matches factors in the character datafile with numbers in the
#' numeric file, with numeric coercion if on
#' Second dataframe matches variable names and variable labels
#' Requires both the character and numeric \code{.csv} Qualtrics output files.
#'
#' @import magrittr
#' @param .numeric_path Path to Numeric File
#' @param .character_path Path to Character File
#' @param .schema_name String for schema name
#' @param .table_name String for table name
#' @param .data_name String giving "dataset name" If not given, constructed from
#' .schema_name and .table_name
#' @param .row_skip Number of rows to skip in data before data begins
#' @param .factorvars Optional argument, vector of factor variable names
#' @param .numcoerce Optional argument, logical value indicating whether numeric
#'  should be coerced from numbers in the character value if present.
#' @return List of factor tibble and label tibble
#' @export
read_meta <- function(.numeric_path,
                      .character_path,
                      .schema_name = NA_character_,
                      .table_name = NA_character_,
                      .data_name = NULL,
                      .row_skip = 2,
                      .factorvars = NULL,
                      .numcoerce = T){
  # Reading in numeric data, assuming either csv or xlsx file
  if(stringr::str_detect(.numeric_path, "\\.csv$")){
    numeric_data <- utils::read.csv(file = .numeric_path, check.names = F)
  } else{numeric_data <- openxlsx::read.xlsx(xlsxFile = .numeric_path, check.names = F)}

  # Reading in character data, assuming either csv or xlsx file
  if(stringr::str_detect(.character_path, "\\.csv$")){
    character_data <- utils::read.csv(file = .character_path, check.names = F)
  } else{character_data <- openxlsx::read.xlsx(xlsxFile = .character_path, check.names = F)}

  data_name <- dplyr::if_else(is.null(.data_name),
                              paste0("[", .schema_name, "].[", .table_name, "]"),
                              .data_name
  )

  # Creating label database
  labeldb <- dplyr::tibble(
    variable_name = names(character_data),
    variable_label = purrr::map_chr(character_data, ~as.character(.[1])),
    Dataset = data_name,
    Schema_Name = .schema_name,
    Table_Name = .table_name
  ) %>%
    dplyr::arrange(as.numeric(stringr::str_extract(variable_name, "\\d+")))

  #Used in Character Format df creation
  # You can add any other consistent string manipulations you wish here
  # IMPORTANT: ANY CHANGES HERE MUST BE MADE TO THE CHARACTER FILE WHEN
  # RUNING read_numeric(), OTHERWISE, YOU WILL LOSE DATA
  CharReplace<- function(x){
    x <- stringr::str_replace_all(x, "\\{e://Field/", "\\{FIELD ") %>%
      stringr::str_replace_all( "\\$", "")
  }

  # Skipping lines on files according to specified
  numeric_data <- tail(numeric_data, -(.row_skip))
  character_data <- tail(character_data, -(.row_skip))

  # Only working with factorvars if defined
  if(!is.null(.factorvars)){
    numeric_data <- numeric_data %>%
      dplyr::select(dplyr::all_of(.factorvars))
    character_data <- character_data %>>%
      dplyr::select(dplyr::all_of(.factorvars))
  } else{
    # If not explicitly set, working with questions that start with "Q"
    numeric_data <- numeric_data %>%
      dplyr::select(dplyr::matches("^Q\\d"))
    character_data <- character_data %>%
      dplyr::select(dplyr::all_of(names(numeric_data)))
  }

  # Coerce all empty string entries to NA
  numeric_data <- dplyr::mutate(
    numeric_data,
    dplyr::across(where(is.character),
                  ~dplyr::if_else(. == "", NA_character_, .)))
  character_data <- dplyr::mutate(
    character_data,
    dplyr::across(where(is.character),
                  ~dplyr::if_else(. == "", NA_character_, .)))

  # Only working with columns where numeric_data is only a number
  numeric_data <- numeric_data %>%
    dplyr::select(where(~all(is.numeric(.)|stringr::str_detect(., "^\\d+\\.*\\d*$"),
                             na.rm = T
    )))
  character_data <- dplyr::select(
    character_data, dplyr::all_of(names(numeric_data))
  )
  # Creating list of lists
  # Each sublist is 2 tibbles, first is original if modified, empty tibble if not
  # Second is modified factors tibble
  factors.list <- list(names(numeric_data), numeric_data, character_data) %>%
    purrr::pmap(function(question, numcol, charcol){
      # If all charcol and numcol are identical to begin with, then treat as non-factor numeric
      if(identical(numcol[!is.na(numcol)], charcol[!is.na(charcol)])){
        return(list(tibble(), tibble()))
      }
      # Creating preliminary mapping of factors
      factors <- tibble(Question = question, NumVal = numcol, CharVal = charcol) %>%
        # Making unique mapping
        dplyr::distinct() %>%
        dplyr::filter(dplyr::if_all(c(NumVal, CharVal), ~!is.na(.))) %>%
        # Removing some standard qualtrics formatting
        dplyr::mutate(CharVal = CharReplace(CharVal))
      # Reporting error if mapping is non-unique
      if(length(unique(factors$CharVal)) != length(unique(factors$NumVal))){
        stop(paste0("Mismatch found in", question))
      }

      # Perform numeric coercion if desired, on by default
      if(.numcoerce){
        # First checking if character value is pseudo-numeric
        # Pseudo numeriic is satisfied if more than 1 case where charval is just a number
        # or satisfies regex below. Only works with factors of levels 0-99
        if(sum(stringr::str_detect(factors$CharVal,
                                   "^\\d{1,2}$|^\\d{1,2}(?=\\s)|(?<=\\s)\\d{1,2}$|\\d{1,2}(?=:\\s)"),
               na.rm = T) > 1){
          # Replacing all NumVal with number extracted from CharVal
          new_factors <- factors %>% dplyr::mutate(
            NumVal = stringr::str_extract(CharVal, "^\\d{1,2}$|^\\d{1,2}(?=\\s)|(?<=\\s)\\d{1,2}$|\\d{1,2}(?=:\\s)"))
          # If new numvals are unique and no NAs, then done
          if(!NA_character_ %in% new_factors$NumVal & !any(duplicated(new_factors$NumVal))){
            return(list(factors, new_factors))}
          # If not, then we will need to repair numbers
          # Taking all used numbers, excluding NAs
          used_numbers <- new_factors$NumVal %>%
            as.numeric() %>%
            unique() %>%
            .[!is.na(.)]
          # Creating vector of all numbers between minimal and maximal value
          used_numbers_seq <- min(used_numbers):max(used_numbers)
          # IF used_numbers_seq isn't long enough, add more to max until equivalent
          if(length(used_numbers_seq) != nrow(new_factors)){
            used_numbers_seq <- c(used_numbers_seq,
                                  seq(max(used_numbers) + 1,
                                      max(used_numbers) + nrow(new_factors) - length(used_numbers_seq))
            )
          }
          new_factors <- new_factors %>%
            # First removing any repeated numvals, keeping first instance
            dplyr::group_by(NumVal) %>%
            dplyr::mutate(NumVal = dplyr::if_else(row_number() > 1, NA_character_, NumVal))
          # Replacing NA NumVal with unused values, without replacement
          # Identifying index of NA values
          na_index <- which(is.na(new_factors$NumVal))
          # Vector of unused numbers of na_index length
          unused_numbers_seq <- used_numbers_seq[!used_numbers_seq %in% new_factors$NumVal] %>%
            .[seq_along(na_index)]
          # Performing replacment
          for(i in seq_along(na_index)){
            new_factors$NumVal[na_index[i]] <- unused_numbers_seq[i]
          }
          return(list(factors, new_factors))
        }
        # If not pseudo numeric, now check if dummy variable (treated as dummy if length 2)
        if(length(factors$CharVal) == 2){
          # If dummy variable, recode as a 0 or 1 with an ttempt made for intelligent ordering
          if(any(stringr::str_detect(factors$CharVal,
                                     stringr::regex("no|false|^f$|failure|not|never", ignore_case = T)))){
            # Taking first "no" index as the true no, assigning to a numval of 0
            no_index <- which(stringr::str_detect(factors$CharVal, stringr::regex("no|false|^f$|failure|not|never", ignore_case = T)))[1]
            new_factors <- factors %>%
              dplyr::mutate(NumVal = dplyr::if_else(row_number() == no_index, "0", "1"))
          } else{new_factors <- factors %>%
            dplyr::mutate(NumVal = dplyr::if_else(row_number() == 1, "0", "1"))}
          return(list(factors, new_factors))
        }
      }
      # No coercion performed
      return(list(tibble(), factors))
    })


  # Creating dataframe of all modified tibbles
  purrr::map_dfr(factors.list, function(factor_list){
    # If identical or unmodified, ignore
    if(nrow(factor_list[[1]]) == 0 | identical(factor_list[[1]], dplyr::ungroup(factor_list[[2]]))){
      return(tibble())} else{
        dplyr::full_join(factor_list[[1]] %>% rename("Original NumVal" = NumVal),
                         factor_list[[2]] %>% rename("New NumVal" = NumVal),
                         by = c("Question", "CharVal")) %>%
          dplyr::arrange(`Original NumVal`)}
  }) %>%
    # Printing all modified factors
    print()

  # Creating new factor database
  factordb <- purrr::map_dfr(factors.list, ~.[[2]]) %>%
    dplyr::mutate(
      Val = paste0(NumVal, ": ", CharVal),
      Dataset = data_name,
      Table_Name = .table_name,
      Schema_Name = .schema_name
    ) %>%
    dplyr::arrange(as.numeric(stringr::str_extract(Question, "\\d+")), as.numeric(NumVal))
  return(list("LabelDatabase" = labeldb, "FactorDatabase" = factordb))
}

