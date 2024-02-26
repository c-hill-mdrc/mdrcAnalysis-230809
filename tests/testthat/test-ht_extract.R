#### All Row Variable Testing
testthat::test_that("HT Unit Testing",{

  # load ht_extract_sas output
  ht_extract_sas <-
    readxl::read_excel(
      path = testthat::test_path("sas_tests", "htsubgroupcomparison.xlsx"))

  # Convert COHORT to a factor for analysis
  finalBif <- bif_final %>%
    dplyr::mutate(CohortFactor = factor(COHORT))

  # Create vector of dependent variables
  finalBif_dependents <- names(finalBif[setdiff(names(finalBif),
                                                c("SAMPLEID", "RA_CODE", "CohortFactor",
                                                  "COHORT", "RA_DATE", "AGE", "DOB", "blcurhrs_73plus",
                                                  "blworkingyes", "blworkingno"))])
  finalBif_treatment  <- "RA_CODE"
  finalBif_covariates <- "CohortFactor"

  ht_extract_r <- ht_extract(.dataset = finalBif,
                             .subgroup = "COHORT",
                             .dependents = finalBif_dependents,
                             .treatment = "RA_CODE",
                             .covariates = NA_character_)

  ht_extract_r$subgrpval <- sub(".*\\D", "", ht_extract_r$UniqueId)

  ht_extract_sas$Dependent <- tolower(ht_extract_sas$Dependent)
  colnames(ht_extract_sas) <- tolower(colnames(ht_extract_sas))

  # filter ht extract without the full sample but just 1 and 2
  ht_extract_sas <- ht_extract_sas %>%
    dplyr::select_if(~!any(is.na(.))) %>%
    dplyr::filter(subgrpval != "Full")

  ht_extract_sas <- ht_extract_sas[ht_extract_sas$dependent %in% finalBif_dependents, ]

  ht_extract_r_matching <-
    ht_extract_r %>%
    dplyr::filter(dependent != "sample_size") %>%
    dplyr::select(
      subgrpval,
      dependent, estimate_p, estimate_c, impact, n_p, n_c, n_total, std_error,
      r_squared, unadj_impact, dep_var_mean, dep_var_sd
    ) %>%
    dplyr::arrange(subgrpval, dependent)

  ht_extract_sas_matching <-
    ht_extract_sas %>%
    dplyr::select(
      subgrpval,
      dependent, adjmean_p, adjmean_c,
      impact_p_c, nobs_p, nobs_c, nobs,
      stderr_p_c, rsquare, unadj_impact_p_c, depvarmean,
      depvarstddev
    ) %>%
    dplyr::rename(
      dependent = dependent,
      estimate_p = adjmean_p,
      estimate_c = adjmean_c,
      impact = impact_p_c,
      n_p = nobs_p,
      n_c = nobs_c,
      n_total = nobs,
      std_error = stderr_p_c,
      r_squared = rsquare,
      unadj_impact = unadj_impact_p_c,
      dep_var_mean = depvarmean,
      dep_var_sd = depvarstddev
    ) %>%
    dplyr::arrange(subgrpval, dependent)

  testthat::expect_equal(object = ht_extract_r_matching,
                         expected = ht_extract_sas_matching,
                         ignore_attr = TRUE,
                         tolerance = 0.01)

})
