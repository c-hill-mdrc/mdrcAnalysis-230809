#### All Row Variable Testing
testthat::test_that("Mean Extract Unit Testing Without Subgroups",{

  # running mean_extract without subgroups
  meanFromFxWithoutSG <-
    mean_extract(.dataset = sim_data_robust_reg) %>%
    dplyr::filter(Variable != "sample_size") %>%
    dplyr::select(-filepath, -user, -datetime, -RVersion, -Packages, -FunctionParams)

  # Filtering Variables from the function
  filteredVars <- unique(meanFromFxWithoutSG$Variable)

  # arrange the variables
  meanFromFxWithoutSG <-
    meanFromFxWithoutSG %>%
    dplyr::arrange(Variable)

  # load mean extract without subgroups from sas
  meanFromSasWithoutSG <-
    readxl::read_excel(
      path = testthat::test_path("sas_tests/meansextractcomp.xlsx"),
      sheet = "full")

  # renaming and reshaping the output file
  meanFromSasWithoutSG <-
    meanFromSasWithoutSG %>%
      dplyr::select(Variable, N, MEAN, STDDEV, MIN, MAX, NMISS) %>%
      dplyr::rename(n_total = N
                   ,Mean = MEAN
                   ,Stddev = STDDEV
                   ,Min = MIN
                   ,Max = MAX
                   ,NMiss = NMISS)

  # arrange the variables
  meanFromSasWithoutSG <-
    meanFromSasWithoutSG %>%
      dplyr::arrange(Variable) %>%
      dplyr::filter(Variable %in% filteredVars)


  testthat::expect_equal(object = meanFromSasWithoutSG,
                         expected = meanFromFxWithoutSG,
                         tolerance = 0.01)

})


#### All Row Variable Testing
testthat::test_that("Mean Extract Unit Testing With Subgroups",{

  # running mean extract with subgroups
  meanFromFxWithSG <-
    mean_extract(
      .dataset = sim_data_robust_reg,
      .subgroup = "treatment") %>%
    dplyr::filter(Variable != "sample_size") %>%
    dplyr::select(-filepath, -user, -datetime, -RVersion, -Packages, -FunctionParams, -n_0, -n_1)

  # Filtering Variables from the function
  filteredVars <- unique(meanFromFxWithSG$Variable)

  # arrange the variables
  meanFromFxWithSG <-
    meanFromFxWithSG %>%
    dplyr::arrange(Variable, treatment)

  # load mean extract with subgroups from sas
  meanFromSasWithSG <-
    readxl::read_excel(
      path = testthat::test_path("sas_tests/meansextractcomp.xlsx"),
      sheet = "byTreatment")

  # renaming and reshaping the output file
  meanFromSasWithSG <-
    meanFromSasWithSG %>%
    dplyr::select(Variable, SubGrpVal, N, MEAN, STDDEV, MIN, MAX, NMISS) %>%
    dplyr::rename(treatment = SubGrpVal
                  ,n_total = N
                  ,Mean = MEAN
                  ,Stddev = STDDEV
                  ,Min = MIN
                  ,Max = MAX
                  ,NMiss = NMISS)

  # arrange the variables
  meanFromSasWithSG <-
    meanFromSasWithSG %>%
    dplyr::filter(Variable %in% filteredVars) %>%
    dplyr::arrange(Variable, treatment)

  testthat::expect_equal(object = meanFromSasWithSG,
                         expected = meanFromFxWithSG,
                         tolerance = 0.01)

})



