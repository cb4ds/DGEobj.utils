context("DGEobj.utils - tests for runContrasts.R functions")


test_that('runContrasts.R: runContrasts()', {
    # Define the named contrasts from design matrix column names
    contrastList  <- list(Sham_vs_BDL     = "ReplicateGroupSham - ReplicateGroupBDL",
                          Sham_vs_EXT1024 = "ReplicateGroupSham  - ReplicateGroupBDL_EXT.1024",
                          Sham_vs_Nint    = "ReplicateGroupSham - ReplicateGroupBDL_Nint",
                          Sham_vs_Sora    = "ReplicateGroupSham - ReplicateGroupBDL_Sora")

    dgeObj_output <- runContrasts(dgeObj              = t_obj1,
                                  designMatrixName    = "ReplicateGroupDesign",
                                  contrastList        = contrastList,
                                  contrastSetName     = "ReplicateGroup_Contrasts")
    expect_s3_class(dgeObj_output, "DGEobj")

    contrastList  <- list(EXT1024_vs_Sham    = "ReplicateGroupBDL_EXT.1024 - ReplicateGroupSham",
                          BDL_vs_EXT1024     = "ReplicateGroupBDL  - ReplicateGroupBDL_EXT.1024",
                          EXT1024_vs_Nint    = "ReplicateGroupBDL_EXT.1024 - ReplicateGroupBDL_Nint",
                          EXT1024_vs_Sora    = "ReplicateGroupBDL_EXT.1024 - ReplicateGroupBDL_Sora")
    dgeObj_output <- runContrasts(dgeObj              = t_obj1,
                                  designMatrixName    = "ReplicateGroupDesign",
                                  contrastList        = contrastList,
                                  contrastSetName     = "ReplicateGroup_Contrasts",
                                  runTopTreat         = TRUE,
                                  qValue              = TRUE,
                                  IHW                 = TRUE,
                                  verbose             = TRUE)
    expect_s3_class(dgeObj_output, "DGEobj")

    # testing assert statements
    expect_error(runContrasts(dgeObj = NULL),
                 regexp = "dgeObj must be specified and should be of class 'DGEobj'.")
    expect_error(runContrasts(dgeObj = t_obj1),
                 regexp = "designMatrixName must be specified.")
    expect_error(runContrasts(dgeObj           = t_obj1,
                              designMatrixName = "ReplicateGroup",
                              contrastList     = "XYZ"),
                 regexp = "contrastList must specified and must be a named list.")
    expect_error(runContrasts(dgeObj              = t_obj1,
                              designMatrixName    = "ReplicateGroup",
                              contrastList        = contrastList,
                              foldChangeThreshold = -1),
                 regexp = "foldChangeThreshold must be greater than or equal to 0.")
    expect_error(runContrasts(dgeObj              = t_obj1,
                              designMatrixName    = "ReplicateGroup",
                              contrastList        = contrastList,
                              runTopTable         = FALSE,
                              runTopTreat         = FALSE),
                 regexp = "One of runTopTable or runTopTreat must be TRUE.")
    expect_error(runContrasts(dgeObj              = t_obj1,
                              designMatrixName    = "XYZ",
                              contrastList        = contrastList),
                 regexp = "The specified designMatrixName not found in dgeObj.")
})
