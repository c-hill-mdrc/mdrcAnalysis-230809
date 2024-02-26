#### All Row Variable Testing
testthat::test_that("RobustRegression Extract Unit Testing Without Cluster SE Adjustment",{

  # load robust regression without cluster standard error adjustment from sas
  sas_robustreg_withoutCluster <-
    readxl::read_excel(
      path = testthat::test_path("sas_tests", "robustreg_simdata.xlsx"),
      sheet = "Not Clustered")

  # applying the lower subgrp value function
  sas_robustreg_withoutCluster$Dependent <- tolower(sas_robustreg_withoutCluster$Dependent)

  sas_robustreg_withoutCluster <-
    sas_robustreg_withoutCluster %>%
      dplyr::select(UniqueId, SubGrpVal, Dependent, AdjMean_1, AdjMean_0, Impact_1_0,
                    ProbT_1_0, StdErr_1_0, NObs_1, NObs_0, NObs,
                    UnAdj_Mean_1, UnAdj_Mean_0, Unadj_Impact_1_0) %>%
      dplyr::arrange(UniqueId,
                     SubGrpVal)

  # running robust regression without clusters
  robustreg_withoutCluster <-
    robustreg_extract(sim_data_robust_reg,
                      .TmntVar = "treatment",
                      .DepVars = c("employed_01",
                                   "employed_02"),
                      .ByX = c("eduLevel"))

  # adjusting robust regression to match with sas robust regression outputs
  robustreg_withoutCluster$SubGrpVal <- as.character(robustreg_withoutCluster$SubGrpVal)

  robustreg_withoutCluster <-
    robustreg_withoutCluster %>%
      dplyr::select(UniqueId, SubGrpVal, Dependent, AdjMean_1, AdjMean_0, Impact_1_0,
                    ProbT_1_0, StdErr_1_0, NObs_1, NObs_0, NObs,
                    UnAdj_Mean_1, UnAdj_Mean_0, Unadj_Impact_1_0) %>%
      dplyr::arrange(UniqueId,
                     SubGrpVal)

  testthat::expect_equal(object = robustreg_withoutCluster,
                         expected = sas_robustreg_withoutCluster,
                         tolerance = 0.01)

})


#### All Row Variable Testing
testthat::test_that("RobustRegression Extract Unit Testing With Cluster SE Adjustment",{

  # load robust regression without cluster standard error adjustment from sas
  sas_robustreg_withCluster <-
    readxl::read_excel(
      path = testthat::test_path("sas_tests", "robustreg_simdata.xlsx"),
      sheet = "Clustered")

  # applying the lower subgrp value function
  sas_robustreg_withCluster$Dependent <- tolower(sas_robustreg_withCluster$Dependent)

  sas_robustreg_withCluster <-
    sas_robustreg_withCluster %>%
    dplyr::select(UniqueId, SubGrpVal, Dependent, AdjMean_1, AdjMean_0, Impact_1_0,
                  ProbT_1_0, StdErr_1_0, NObs_1, NObs_0, NObs,
                  UnAdj_Mean_1, UnAdj_Mean_0, Unadj_Impact_1_0) %>%
    dplyr::arrange(UniqueId,
                   SubGrpVal)

  # running robust regression without clusters
  robustreg_withCluster <-
    robustreg_extract(sim_data_robust_reg,
                      .TmntVar = "treatment",
                      .DepVars = c("employed_01",
                                   "employed_02"),
                      .ByX = c("eduLevel"),
                      .ClusterSE = "sites")

  # adjusting robust regression to match with sas robust regression outputs
  robustreg_withCluster$SubGrpVal <- as.character(robustreg_withCluster$SubGrpVal)

  robustreg_withCluster <-
    robustreg_withCluster %>%
    dplyr::select(UniqueId, SubGrpVal, Dependent, AdjMean_1, AdjMean_0, Impact_1_0,
                  ProbT_1_0, StdErr_1_0, NObs_1, NObs_0, NObs,
                  UnAdj_Mean_1, UnAdj_Mean_0, Unadj_Impact_1_0) %>%
    dplyr::arrange(UniqueId,
                   SubGrpVal)


  # number of rows of robust regression without cluster
  nrow(sas_robustreg_withCluster)

  # number of rows of robust regression with cluster
  nrow(robustreg_withCluster)

  testthat::expect_equal(object = robustreg_withCluster,
                         expected = sas_robustreg_withCluster,
                         tolerance = 0.01)

})

