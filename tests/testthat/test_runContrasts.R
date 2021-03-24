context("DGEobj.utils - tests for runContrasts.R functions")


test_that('runContrasts.R: runContrasts()', {
    dgeObj <- t_obj1
    contrastList <- getType(dgeObj, "topTable")
    names(contrastList) <- colnames(dgeObj$ReplicateGroupDesign)[-1]

    dgeObj_output <- runContrasts(dgeObj              = dgeObj,
                                  designMatrixName    = "ReplicateGroupDesign",
                                  contrastList        = contrastList,
                                  contrastSetName     = "ReplicateGroup_Contrasts")
    expect_s3_class(dgeObj_output, "DGEobj")

    dgeObj_output <- runContrasts(dgeObj              = dgeObj,
                                  designMatrixName    = "ReplicateGroupDesign",
                                  contrastList        = contrastList,
                                  contrastSetName     = "ReplicateGroup_Contrasts",
                                  runTopTreat         = TRUE,
                                  qValue              = TRUE,
                                  IHW                 = TRUE)
    expect_s3_class(dgeObj_output, "DGEobj")

    # testing assert statements
    expect_error(runContrasts(dgeObj = NULL),
                 regexp = "dgeObj must be specified and should be of class 'DGEobj'.")
    expect_error(runContrasts(dgeObj = dgeObj),
                 regexp = "designMatrixName must be specified.")
    expect_error(runContrasts(dgeObj           = dgeObj,
                              designMatrixName = "ReplicateGroup",
                              contrastList     = "XYZ"),
                 regexp = "contrastList must specified and must be a named list.")
    expect_error(runContrasts(dgeObj              = dgeObj,
                              designMatrixName    = "ReplicateGroup",
                              contrastList        = contrastList,
                              foldChangeThreshold = -1),
                 regexp = "foldChangeThreshold must be greater than or equal to 0.")
    expect_error(runContrasts(dgeObj              = dgeObj,
                              designMatrixName    = "ReplicateGroup",
                              contrastList        = contrastList,
                              runTopTable         = FALSE,
                              runTopTreat         = FALSE),
                 regexp = "One of runTopTable or runTopTreat must be TRUE.")
    expect_error(runContrasts(dgeObj              = dgeObj,
                              designMatrixName    = "XYZ",
                              contrastList        = contrastList),
                 regexp = "The specified designMatrixName not found in dgeObj.")
})
