#' Generate formatted, ordered codebooks.
#'
#' Generate LaTeX formatted codebooks for use when knittiing to pdf or a list of
#' sumamry dataframes for all other use cases. Contains support for metadata
#' structures as described by \code{metaSQL}
#' This must be included in YAML header if specifying pdf output
#' \code{header-includes:}
#' \code{- `\usepackage{booktabs}`}
#' \code{- `\usepackage{longtable}`}
#' \code{- `\usepackage{float}`}
#' Need results='asis' in the chunk options for formatting
#'
#' @import magrittr
#' @param .freq_factors Optional vector containing all of the variables that
#' should have a frequency table displayed
#' @param .mean_factors Optional vector containing all of the variables that
#'  should have a means summary displayed
#' @param .data Tibble containing numeric data that you are attempting to
#'  process into a codebook
#' @param .labeldb Optional dataframe of labels.
#' @param .factordb Optional dataframe of factors.
#' @param .round.digits numeric containing how many decimal places to round to.
#'  Default = 2 See \code{kableExtra::column_spec()} width arg for more options
#' @param .output String specying output type. Currrently accepts pdf and list.
#' @returns A printout of the codebook in question order or list of dataframes
#' used to construct
#' @export

codebook_interleave <- function(.freq_factors = NULL,
                                .mean_factors = NULL,
                                .data,
                                .labeldb = NULL,
                                .factordb = NULL,
                                .round.digits = 2,
                                .output = "pdf"){

  freq_factors <- .freq_factors
  mean_factors <- .mean_factors
  data <- .data

  # Creating character data from factors if .factordb available
  if(!is.null(.factordb)){
    factordb <- dplyr::filter(.factordb, Question %in% names(data))
    if(nrow(.factordb) == 0){
      char_data <- data
      factordb <- .factordb
    } else{
      # Storing original factors
      factordb_legend <- factordb
      # Shortening strings to use as factors
      factordb <- dplyr::mutate(factordb, CharVal = str_repair_latex(
        stringr::str_trunc(Val, 50, "center")))
      #Converting to factor where possible
      char_data <- dplyr::mutate(
        data, dplyr::across(dplyr::all_of(factordb$Question), function(numval){
          factordb <- dplyr::filter(factordb, Question == dplyr::cur_column())
          factordb$Val[match(numval, factordb$NumVal)]
        }))
    }
  } else if(!is.null(freq_factors)){
    char_data <- data %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(freq_factors), as.character))
    factordb <- .factordb
  } else{
    char_data <- data
    factordb <- .factordb
  }

  if(!is.null(mean_factors)){
    # Make numeric any variables we are calculating the means of
    num_data <- data %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(mean_factors), as.numeric))
  } else{num_data <- data}

  # Defining functions used for frequencies and numeric summaries
  # Ideally could create this externally and load as a package

  # Frequencies
  freq <- function(x){
    freq_df <- unique(x) %>% sort(na.last = T) %>%
      purrr::map_dfr(function(value){
        if(!is.na(value)){
          col_1 <- sum(x == value, na.rm = T)
          col_2 <- col_1/sum(!is.na(x))
          col_3 <- col_1/length(x)
        } else{
          col_1 <- sum(is.na(x))
          col_2 <- NA_real_
          col_3 <- sum(is.na(x))/length(x)}
        dplyr::tibble(
          "Value" = as.character(value),
          "Freq" = col_1,
          "% Valid" = col_2,
          "% Total" = col_3)
      })
    freq_df <- freq_df %>%
      dplyr::add_row(
        "Value" = "Total",
        "Freq" = length(x),
        "% Valid" = sum(!is.na(x))/length(x), #JM: percent valid of total
        "% Total" = 1
      ) %>%
      dplyr::mutate(dplyr::across(3:4, ~round(.*100, .round.digits)))
    return(freq_df)
  }

  # Means and median
  descr <- function(x){
    n <- length(x)
    mean_df <- dplyr::tibble(
      " " = c("Min", "Max", "Mean", "Std.Dev", "N.Valid", "Pct.Valid"),
      "value" = c(round(min(x, na.rm = T),.round.digits),
                  round(max(x, na.rm = T),.round.digits),
                  round(mean(x, na.rm = T),.round.digits),
                  round(sd(x, na.rm = T),.round.digits),
                  round(sum(!is.na(x)),.round.digits),
                  round(sum(!is.na(x))/length(x) * 100,.round.digits)
      ))
    median_df <- dplyr::tibble(
      " " = c("Min", "Q1", "Median", "Q3", "Max", "IQR"),
      "value" = c(round(min(x, na.rm = T),.round.digits),
                  round(stats::quantile(x, .25, na.rm = T),.round.digits),
                  round(stats::quantile(x, .5, na.rm = T),.round.digits),
                  round(stats::quantile(x, .75, na.rm = T),.round.digits),
                  round(max(x, na.rm = T),.round.digits),
                  round(stats::quantile(x, .75, na.rm =  T) -
                          stats::quantile(x, .25, na.rm = T),
                       .round.digits)

      ))
    return(c(list(mean_df), list(median_df)))
  }

  # If "list" or "plain" listed as output, return list of dataframe output, named
  if(.output %in% c("plain", "list")){
    # Generating list
    codebook.list <- purrr::map(
      names(data) %>% set_names(., .), function(variable_name){
        if(!variable_name %in% c(mean_factors, freq_factors)){
          return(NULL)
        }
        # Initializing list to input codebook extracts into
        variable.list <- list()
        variable.list$variable_name <- variable_name
        # Storing verbose label if possible
        variable.list$variable_label <- dplyr::if_else(
          variable_name %in% .labeldb$variable_name,
          # Concatenate label if possible
          paste0(variable_name, ": ",
                 .labeldb$variable_label[.labeldb$variable_name == variable_name]),
          # Only variable name
          variable_name
        )
        # Storing frequencies if requested
        if(variable_name %in% freq_factors){
          variable.list$frequency <- freq(dplyr::pull(char_data, variable_name))
        }
        # Storing univariates if requested
        if(variable_name %in% mean_factors){
          univariates <- descr(dplyr::pull(num_data, variable_name))
          variable.list$mean <- univariates[[1]]
          variable.list$median <- univariates[[2]]
        }
        return(variable.list)
      })
    # Removing NULLs from list
    codebook.list <- codebook.list[which(!purrr::map_lgl(codebook.list, is.null))]
    return(codebook.list)
  }
  # Displaying formatted tables in order as they appear in data if output is "latex" or "pdf"
  else if(.output %in% c("LATEX", "LaTeX", "latex", "pdf", "PDF")){
    purrr::walk(names(data), function(variable_name){
      # IF not in names, skip to next question
      if(!variable_name %in% c(mean_factors, freq_factors)){
        return(NULL)
      }

      variable_print <- str_repair_latex(variable_name)

      # If selected for codebooking, print name and variable label if possible
      cat("\\begin{minipage}{\\textwidth}")
      cat("  \n")
      dplyr::if_else(variable_name %in% .labeldb$variable_name,
                     # Concatenate label if possible
                     paste0(str_repair_latex(variable_print), ": ",
                            str_repair_latex(
                              .labeldb$variable_label[
                                .labeldb$variable_name == variable_name])),
                     # Only variable name
                     str_repair(variable_print)
      ) %>% cat("\\subsubsection{", ., "}  \n")

      # First displaying frequencies if requested
      if(variable_name %in% freq_factors){
        kableExtra::kbl(freq(dplyr::pull(char_data, variable_name)),
                        booktabs = T, longtable = T, digits =.round.digits,
                        label = stringr::str_remove_all(
                          paste0(variable_print, "freq"), "[:graph:]")) %>%
          kableExtra::column_spec(1, width = "20em") %>%
          kableExtra::kable_styling(latex_options = "HOLD_position") %>%
          print()
      }

      # Displaying means if requested
      if(variable_name %in% mean_factors){
        numeric_summaries <- descr(dplyr::pull(num_data, variable_name)) %>%
          purrr::imap(~kableExtra::kable_styling(
            kableExtra::kbl(.,
                            booktabs = T, longtable = T, digits = .round.digits,
                            caption = dplyr::if_else(.y == 1, "Mean", "Median"),
                            label = stringr::str_remove_all(
                              dplyr::if_else(.y == 1,
                                             paste0(variable_print, "mean"),
                                             paste0(variable_print, "median")),
                              "[:graph:]")),
            latex_options = "HOLD_position"))
        cat("\\begin{figure}[H]  \n")
        cat("\\begin{minipage}[c]{0.5\\textwidth}  \n")
        cat("\\centering  \n")
        print(numeric_summaries[[1]])
        cat("  \n")
        cat("\\end{minipage}")
        cat("\\begin{minipage}[c]{0.5\\textwidth}  \n")
        cat("\\centering  \n")
        print(numeric_summaries[[2]])
        cat("  \n")
        cat("\\end{minipage}  \n")
        cat("\\end{figure}")
      }
      # Displaying legend if only means and factors exist
      # Displaying first 10 factors only
      if(variable_name %in% mean_factors & (!variable_name %in% freq_factors) &
         variable_name %in% factordb$Question){
        cat("\\textbf{Legend:}\\\\  \n")
        if(length(factordb$Val[factordb$Question == variable_name]) < 11){
          cat(factordb$Val[factordb$Question == variable_name],
              sep = "\\\\  \n")
        } else{
          cat(c(head(factordb$Val[factordb$Question == variable_name], 10),
                "Only first 10 factors displays, truncated"),
              sep = "\\\\  \n")
        }
      }
      # Adding newline
      cat("\\end{minipage}")
      cat("  \n\\vspace{5mm}  \n")
    })
  }
}


