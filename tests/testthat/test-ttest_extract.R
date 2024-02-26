#### All Row Variable Testing
testthat::test_that("T-Test Unit Testing No Subgroups",{

  # load ht_extract_sas output
  ttest_extract_sas <-
    readxl::read_excel(
      path = testthat::test_path("sas_tests", "ttestcomp.xlsx"),
      sheet = "NoBy")

  # select columns from ttest_extract_sas
  ttest_extract_sas <-
    ttest_extract_sas %>%
      dplyr::select(Variable, N_C, N_P, N_Total, Mean_C, Mean_P, Diff,
                    Mean_Total, EqualTValue, UnequalTValue, EqualProbt,
                    UnequalProbt, ProbF, TValue, ProbT)

  # Convert COHORT to a factor for analysis
  finalBif <- bif_final %>%
    dplyr::mutate(CohortFactor = factor(COHORT))

  # Pull t-test extract out
  ttest_extract_r <- ttest_extract(finalBif, .classvar = "RA_CODE")

  # Select columns from ttest_extract_r
  # ttest differences are aligned. What is going on with the sas version?
  ttest_extract_r <-
    ttest_extract_r %>%
      dplyr::select(Variable, N_C, N_P, N_Total, Mean_C, Mean_P, Diff,
                    Mean_Total, EqualTValue, UnequalTValue, EqualProbt,
                    UnequalProbt, ProbF, TValue, ProbT) %>%
      dplyr::filter(Variable %in% c("AGE", "bldiplomas_hs", "bldiplomas_ged", "bldiplomas_tec",
                                    "bldiplomas_as", "bldiplomas_4yr", "bldiplomas_md", "bldiplomas_non"))

  # comparing ttest for no group by object
  testthat::expect_equal(object = ttest_extract_sas,
                         expected = ttest_extract_r,
                         tolerance = 0.01,
                         ignore_attr = TRUE)

})

#### All Row Variable Testing
testthat::test_that("T-Test Unit Testing By Subgroup",{

  # load ht_extract_sas output
  ttest_extract_sas <-
    readxl::read_excel(
      path = testthat::test_path("sas_tests", "ttestcomp.xlsx"),
      sheet = "ByCohort")

  # select columns from ttest_extract_sas
  ttest_extract_sas <-
    ttest_extract_sas %>%
    dplyr::select(UniqueId, SubGrpVal, Variable, N_C, N_P,
                  N_Total, Mean_C, Mean_P, Diff, Mean_Total,
                  EqualTValue, UnequalTValue, EqualProbt, UnequalProbt,
                  ProbF, TValue, ProbT)

  # Convert COHORT to a factor for analysis
  finalBif <- bif_final %>%
    dplyr::mutate(CohortFactor = factor(COHORT))

  # Pull t-test extract out
  ttest_extract_r <- ttest_extract(finalBif, .classvar = "RA_CODE",
                                   .byx = "COHORT",
                                   .depvars = c("AGE", "bldiplomas_hs", "bldiplomas_ged",
                                                "bldiplomas_tec", "bldiplomas_as", "bldiplomas_4yr",
                                                "bldiplomas_md", "bldiplomas_non"),
                                   .depvar_labels = c("Age (years)", "High School Diploma",
                                                      "GED", "Technical School",
                                                      "Associate's Degree", "Bachelor's Degree",
                                                      "Master's Degree", "No Degree"))

  # Select columns from ttest_extract_r
  # ttest differences are aligned. What is going on with the sas version?
  ttest_extract_r <-
    ttest_extract_r %>%
    dplyr::select(
      UniqueID, COHORT, Variable, N_C, N_P,
      N_Total, Mean_C, Mean_P, Diff, Mean_Total,
      EqualTValue, UnequalTValue, EqualProbt, UnequalProbt,
      ProbF, TValue, ProbT) %>%
    dplyr::filter(Variable %in% c("AGE", "bldiplomas_hs", "bldiplomas_ged", "bldiplomas_tec",
                                  "bldiplomas_as", "bldiplomas_4yr", "bldiplomas_md", "bldiplomas_non"))


  ## Preparing each before subsetting
  ### Renaming the ttest_extract_r outputs
  ttest_extract_r <- ttest_extract_r %>%
    dplyr::rename("UniqueId" = "UniqueID",
                  "SubGrpVal" = "COHORT") %>%
    dplyr::arrange(UniqueId, SubGrpVal)

  ### Rearranging the ttest_extract_sas outputs
  ttest_extract_sas <- ttest_extract_sas %>%
    dplyr::arrange(UniqueId, SubGrpVal)

  # comparing ttest for no group by object
  testthat::expect_equal(object = ttest_extract_sas,
                         expected = ttest_extract_r,
                         tolerance = 0.01,
                         ignore_attr = TRUE)

})



