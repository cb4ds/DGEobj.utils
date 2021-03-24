context("DGEobj.utils - tests for runEdgeRNorm.R functions")


test_that('runEdgeRNorm: runEdgeRNorm()', {
    dgeobj <- t_obj1
    dgeobj$DGEList <- NULL
    runEdgeRNorm_one_test <- runEdgeRNorm(dgeobj, plotFile = FALSE)
    runEdgeRNorm_one_test_DGEList <- getType(runEdgeRNorm_one_test, "DGEList")

    expect_s3_class(runEdgeRNorm_one_test, "DGEobj")
    expect_true(is.list(runEdgeRNorm_one_test_DGEList))
    expect_equal(length(runEdgeRNorm_one_test$DGEList), 2)
    expect_equal(names(runEdgeRNorm_one_test$DGEList), c("counts", "samples"))

    runEdgeRNorm_two_test <- runEdgeRNorm(dgeobj, normMethod = "RLE", plotFile = FALSE)

    runEdgeRNorm_two_test_DGEList <- getType(runEdgeRNorm_two_test, "DGEList")

    expect_s3_class(runEdgeRNorm_two_test, "DGEobj")
    expect_true(is.list(runEdgeRNorm_two_test_DGEList))
    expect_equal(length(runEdgeRNorm_two_test$DGEList), 2)
    expect_equal(names(runEdgeRNorm_two_test$DGEList), c("counts", "samples"))

    runEdgeRNorm_three_test <- runEdgeRNorm(dgeobj, plotFile = FALSE)
    runEdgeRNorm_three_test_DGEList <- getType(runEdgeRNorm_three_test, "DGEList")

    expect_s3_class(runEdgeRNorm_three_test, "DGEobj")
    expect_true(is.list(runEdgeRNorm_three_test_DGEList))
    expect_equal(length(runEdgeRNorm_three_test$DGEList), 2)
    expect_equal(names(runEdgeRNorm_three_test$DGEList), c("counts", "samples"))

    expect_error(runEdgeRNorm(runEdgeRNorm_test),
                 regexp = "object 'runEdgeRNorm_test' not found")
})


test_that('runEdgeRNorm: incorrect usage', {
    expect_error(runEdgeRNorm(),
                 regexp = "argument \"dgeObj\" is missing, with no default")
})
