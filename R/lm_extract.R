#' Generate dataframe of linear model extracts across a set of dependent variables
#'
#' @param .dataset    data frame name providing all necessary variables
#' @param .dependents   character vector of all column names for outcomes
#' @param .treatment    character vector of column name(s) for .treatment/independent variable
#' @param .covariates   character vector of column names for other predictors/explanatory variables
#' @param .asfactorcovs character vector of .covariates that should be treated as factors/categorical variables
#' @param .confintalpha decimal value of desired alpha level for confidence intervals, i.e. 0.90, 0.95, 0.99. Default is 0.95
#' @param .output_path optional, path to output location
#' @param .output_file optional, filename for output file
#'
#' @param .inc_trail   RTU parameter to include trail in output
#' @param .inc_sample  RTU parameter to include sample size row in output
#'
#'
#' @author Melvin Gutierrez
#' @return data frame with one row per dependent with desired regression statistics
#' @export
#'
#' @examples
#'# Functional Sample Call
#'# Copy outside of function and run after function
#'
#'# Convert COHORT to a factor for analysis
#'FinalBIF <- dplyr::mutate(bif_final, CohortFactor = factor(COHORT))
#'
#'# Create vector of dependent variables
#'## Using `setdiff()` to create vector of all variables except those listed
#' FinalBIF_dependents <- names(FinalBIF[setdiff(names(FinalBIF),
#'                                               c("SAMPLEID",
#'                                                 "RA_CODE",
#'                                                 "CohortFactor",
#'                                                 "COHORT",
#'                                                 "RA_DATE",
#'                                                 "AGE",
#'                                                 "DOB",
#'                                                 "blcurhrs_73plus"))])
#' FinalBIF_.treatment  <- "RA_CODE"
#' FinalBIF_.covariates <- "CohortFactor"
#'
#' extract <- lm_extract(FinalBIF,
#'                       FinalBIF_dependents,
#'                       FinalBIF_.treatment,
#'                       FinalBIF_.covariates,
#'                       .confintalpha = .9)
#'
#' extract
lm_extract <- function (.dataset,
                        .dependents,
                        .treatment,
                        .covariates = NA_character_,
                        .asfactorcovs = NA_character_,
                        .confintalpha = 0.95,
                        .inc_trail = TRUE,
                        .inc_sample = TRUE,
                        .output_path = NA_character_,
                        .output_file = NA_character_) {

  # Parameter checks
  if (!is.data.frame(.dataset))   {stop("ERROR: lm_extract expects a data frame as the first argument.")}
  if (!is.character(.dependents))   {stop("ERROR: lm_extract expects a character vector of dependents.")}
  if (!is.character(.covariates))   {stop("ERROR: lm_extract expects a character vector of .covariates.")}
  if (!is.character(.asfactorcovs)) {stop("ERROR: lm_extract expects a character vector of .covariates to be treated as factors.")}
  if (!is.character(.treatment))    {stop("ERROR: lm_extract expects a character vector of .treatment." )}
  if (length(attr(.dataset, "groups")) > 0) {stop("ERROR: data frame is grouped/rowwise. This may cause issues with calculations. Use ungroup() to remove.")}

  # Check for labels
  label_df <- data.frame()
  if (purrr::map(.dataset, attr, "label") %>% unlist %>% length() > 0) {
    label_df <-
      purrr::map(.dataset, attr, "label") %>%
      unlist() %>%
      as.data.frame() %>%
      tibble::rownames_to_column() %>%
      dplyr::rename("dependent" = "rowname", "label" = ".")
  }

  # curate to non missing for treatment variable
  .dataset_nomiss_treat <-
    .dataset %>%
    dplyr::filter(!is.na(!!sym(.treatment)))

  ##################################################################################
  # 5 number summary via mean_extract for ungrouped variables
  ##################################################################################

  UnadjStatsUnGrp <-
    mean_extract(
      .dataset = .dataset_nomiss_treat
      ,.inc_sample = FALSE
      ,.inc_trail = FALSE
    )

  ##################################################################################
  # trimming to the outcome and renaming variable names
  ##################################################################################

  UnadjStatsUnGrp <-
    UnadjStatsUnGrp %>%
      dplyr::filter(Variable %in% .dependents) %>%
      dplyr::select(dependent  = Variable
                   ,DepVarMean = Mean
                   ,DepVarSD   = Stddev
                   ,DepVarMin  = Min
                   ,DepVarMax  = Max
                   ,DepVarNMiss= NMiss
                   )

  ##################################################################################
  # 5 number summary via mean_extract for grouped variables
  ##################################################################################

  UnadjStatsGrp <-
    mean_extract(
       .dataset = .dataset_nomiss_treat
      ,.subgroup = .treatment
      ,.inc_sample = FALSE
      ,.inc_trail = FALSE
    ) %>%
    dplyr::filter(!is.na(!!rlang::sym(.treatment)))

  ##################################################################################
  # trimming to the outcome and renaming variable names
  ##################################################################################

  UnadjStatsGrp <-
    UnadjStatsGrp %>%
    dplyr::filter(Variable %in% .dependents) %>%
    dplyr::select(dependent   = Variable
                 ,DepVarN     = n_total
                 ,DepVarMean  = Mean
                 ,DepVarSD    = Stddev
                 ,!!rlang::sym(.treatment)
                 ) %>%
    tidyr::pivot_wider(id_cols = "dependent",
                       names_from = !!rlang::sym(.treatment),
                       values_from = c("DepVarN"
                                      ,"DepVarMean"
                                      ,"DepVarSD")
                       )

  ##################################################################################
  # Unadjusted Statistics combined
  ##################################################################################

  UnadjStats <-
    dplyr::full_join(UnadjStatsUnGrp, UnadjStatsGrp, by = "dependent") %>%
    janitor::clean_names()

  # Construct formulas
  ## Creating vector for unadjusted model and adjusted model
  ## Create list to hold each formula in advance (best practice); using as.formula to convert string to formula
  ## Create vector to hold values of dependents which have to be removed in advance (best practice)
  ## NOTE: Using seq_along to traverse each dependent; seq_along is recommended for loops (best practice)
  fms <- vector(mode="list", length=length(.dependents))
  unadjfms <- vector(mode="list", length=length(.dependents))
  fms_missing <- vector(mode="logical", length = length(.dependents))

  for (i in seq_along(.dependents)) {
    # Assign constructed formula to position in vector
    ## Note: if dependent doesn't have variation, the model will fail, note which elements would cause this issue
    if (length(table(.dataset[.dependents[i]])) > 1) {

      fms[[i]]      <- as.formula(paste(.dependents[i], paste(c(.treatment, .covariates), collapse = " + "), sep = " ~ "))
      unadjfms[[i]] <- as.formula(paste(.dependents[i], "~", .treatment))
    }
    else {
      print(paste("NOTE:",
                   .dependents[i],
                  "has no variation and will be removed",
                  "from the dependent list."))
      fms_missing[[i]] <- TRUE
      fms[[i]]         <- NA
      unadjfms[[i]]    <- NA
    }
  }
  # Remove formulas and dependents that have no variation
  if (length(which(fms_missing)) > 0) {
    fms        <- fms[-which(fms_missing)]
    unadjfms   <- unadjfms[-which(fms_missing)]
    .dependents <- .dependents[-which(fms_missing)]
  }
  else {
    fms        <- fms
    unadjfms   <- unadjfms
    .dependents <- .dependents
  }


  names(fms)      <- .dependents
  names(unadjfms) <- .dependents

  # Join on Unadj Statistics
  lm_unadj <-
    purrr::map2_dfr(unadjfms,
             .dependents,
             rowrecord,
             .dataset = .dataset,
             .treatment = .treatment,
             .asfactorcovars =  NA_character_,
             .confintalpha = .confintalpha,
             .output_path = .output_path,
             .output_file = .output_file)

  if (all(is.na(.covariates))) {

    lm_all <- lm_unadj

  } else {

    lm_all <-
      purrr::map2_dfr(fms,
               .dependents,
               rowrecord,
               .dataset = .dataset,
               .treatment = .treatment,
               .asfactorcovars = .asfactorcovs,
               .confintalpha = .confintalpha,
               .output_path = .output_path,
               .output_file = .output_file)

  }


  lm_unadj <-
    lm_unadj %>%
    janitor::clean_names() %>%
    dplyr::rename_with(~paste("unadj_",.), -all_of("dependent"))

  AllJoined <-
    lm_all %>%
    janitor::clean_names() %>%
    dplyr::full_join(lm_unadj, by = "dependent") %>%
    dplyr::full_join(UnadjStats, by = "dependent") %>%
    dplyr::relocate(formulaused, .after = last_col()) %>%
    dplyr::mutate(covarsused = toString(.covariates))


  # Adding Effect Sizes
  ## Effect Size Equations from: https://stats.idre.ucla.edu/other/mult-pkg/faq/general/effect-size-power/faqhow-is-effect-size-used-in-power-analysis/
  ## Create var with levels of .treatment variable
  treat_lvls <- tolower(names(table(.dataset[.treatment])))

  ## Overall Effect Size
  Overall_ES <-
    AllJoined %>%
    dplyr::mutate(EffectSize_Overall = impact/dep_var_sd)

  ## Pooled Effect Size
  Pooled_ES <-
    Overall_ES %>%
    dplyr::mutate(EffectSize_Pooled = impact/sqrt((((get(paste0("n_",treat_lvls[1]))-1)*get(paste0("dep_var_sd_",treat_lvls[1]))^2) +
                                              ((get(paste0("n_",treat_lvls[2]))-1)*get(paste0("dep_var_sd_",treat_lvls[2]))^2)) /
                                             (get(paste0("n_",treat_lvls[1])) + get(paste0("n_",treat_lvls[2])) - 2)))

  ## Control Effect Size
  Control_ES <-
    Pooled_ES %>%
    dplyr::mutate(EffectSize_Control = impact/get(paste0("dep_var_sd_",treat_lvls[1]))) %>%
    dplyr::relocate(starts_with("EffectSize"), .before = "r_squared")

  # clean names
  final_tib <- janitor::clean_names(Control_ES)

  # Add labels if they exist
  if (length(label_df) > 0) {
    final_tib <-
      final_tib %>%
      dplyr::left_join(label_df, by = "dependent") %>%
      dplyr::relocate("label", .after = "dependent")
  }

  # Adding the sample size function to the larger table
  if(.inc_sample == TRUE) {

    # pull out the sample counts
    sampleCounts <-
      sample_size(.dataset = .dataset,
                  .treatment = .treatment
      )

    # rename the sample count function first column to the first column of the Tib
    names(sampleCounts)[1] <- names(final_tib)[1]

    final_tib <-
      final_tib %>%
        dplyr::bind_rows(sampleCounts)

  }

  # Adding the trail function
  if(.inc_trail == TRUE) {

    ## Capturing the function parameters
    call <- match.call()

    final_tib <-
      final_tib %>%
      dplyr::bind_cols(trail(call))

  }

  # clean names
  final_tib <- janitor::clean_names(final_tib)

  # Outputting to an excel spreadsheet
  if (!is.na(.output_path) & !is.na(.output_file)) {

    create_excel(
      .output_path = .output_path,
      .output_file = .output_file,
      .x = final_tib
    )
  } #end of outputting to excel if statement

  # return extract tibble
  return(final_tib)

}


# Create function for running regression
## Function will take a single formula and create the output record
## Function will later be called on all formulas
rowrecord <- function(
    .formula,
    .dependent,
    .dataset,
    .treatment,
    .asfactorcovars = NA_character_,
    .confintalpha = 0.95,
    .output_path,
    .output_file) {


  # Run regression
  outlm <- lm(.formula, .dataset)

  # Apply tidy to outlm result
  ## Drop unnecessary rows (intercept/.covariates)
  outlmtidy <-
    broom::tidy(outlm) %>%
    dplyr::filter(stringr::str_detect(term, paste0("^", .treatment))) %>%
    dplyr::rename(dependent = term,
           impact = estimate,
           t_statistic = statistic) %>%
    dplyr::mutate(dependent = .dependent)

  ## Apply glance to lm result for r-squared
  outlmrsq <-
    broom::glance(outlm) %>%
    dplyr::select("r.squared", "adj.r.squared") %>%
    dplyr::mutate(dependent = .dependent)

  # Apply augment to lm to calculate sample sizes
  ## Calculating sample sizes
  SampleSizes <-
    broom::augment(outlm) %>%
    dplyr::count(.data[[.treatment]]) %>%
    tidyr::pivot_wider(names_from = 1, values_from = n, names_prefix = "N_") %>%
    dplyr::mutate(N_Total = sum(dplyr::c_across(everything()))) %>%
    dplyr::mutate(dependent = .dependent)

  # Apply emmeans to linear model
  # Change reference grid to treat only categorical variables as factors
  if (length(.asfactorcovars) > 0) {
    covkeep <- c(.treatment, .asfactorcovars)
  } else {
    covkeep <- .treatment
  }
  outlm_rg <- emmeans::ref_grid(outlm, cov.keep = covkeep)

  ## NOTE: In addition to providing the specs parameter (as with a single call), the data frame must also be specified
  ##       Otherwise, the emmeans functions do not know which data to use for the emmeans
  emm <- emmeans::emmeans(outlm_rg, data = .dataset, specs = .treatment)

  # Apply tidy to emmeans result
  ## Drop unnecessary columns
  emmtidy <-
    broom::tidy(emm) %>%
    dplyr::select(all_of(.treatment), "estimate", "std.error")

  # Apply confint to each linear model
  ## level based on argument provided
  ## default is .95

  if (.confintalpha == 0.90) {
    low_name <- "5 %"
    high_name <- "95 %"
  }
  else if (.confintalpha == 0.95) {
    low_name <- "2.5 %"
    high_name <- "97.5 %"
  }
  else if (.confintalpha == 0.99) {
    low_name <- "0.5 %"
    high_name <- "99.5 %"
  }

  ### estimated means CIs
  ### dropping unnecessary columns
  estCIs <-
    confint(emm, parm = "estimate", level = .confintalpha) %>%
    broom::tidy() %>%
    dplyr::select(-std.error, -df)

  ### estimated impact CIs
  impCIs <-
    confint(outlm, level = .confintalpha) %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(stringr::str_detect(rowname, .treatment)) %>%
    dplyr::rename(dependent = rowname,
           conf.low_impact = all_of(low_name),
           conf.high_impact = all_of(high_name)) %>%
    dplyr::mutate(dependent = .dependent)

  # Merge all emmtidys, estCIs, impCIs, outlmtidys
  outlmJoins <-
    dplyr::full_join(emmtidy, estCIs, by = c(.treatment, "estimate")) %>%
    tibble::add_column("dependent" = .dependent) %>%
    tidyr::pivot_wider(id_cols = "dependent", names_from = all_of(.treatment), values_from=-dependent) %>%
    dplyr::full_join(outlmtidy, by = "dependent") %>%
    dplyr::full_join(impCIs, by = "dependent") %>%
    dplyr::full_join(outlmrsq, by = "dependent") %>%
    dplyr::full_join(SampleSizes, by = "dependent")

  # Drop .treatment columns, add trail, and return final tibble
  final_outlm <-
    outlmJoins %>%
    dplyr::mutate(stars = dplyr::case_when(
        p.value < 0.01 ~ "***",
        p.value < 0.05 ~ "**",
        p.value < 0.10 ~ "*",
        TRUE ~ ""
    )) %>%
    dplyr::mutate(formulaused = toString(.formula)) %>%
    dplyr::select(dependent,
             dplyr::starts_with("estimate"),
               impact,
               p.value,
               stars,
               t_statistic,
             tidyselect::starts_with("std.error"),
             tidyselect::starts_with("conf"),
             tidyselect::starts_with("Effect"),
               r.squared,
               adj.r.squared,
             tidyselect::starts_with("n_"),
             tidyselect::starts_with("Dep"),
               formulaused
    )


return(final_outlm)

}
