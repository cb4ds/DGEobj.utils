context("DGEobj.utils - tests for extractCol.R functions")

test_that('extractCol: extractCol()', {
    extractCol_contrastList <- getType(t_obj1, "topTable")[1:2]
    extractCol_one_test <- extractCol(extractCol_contrastList, colName = "P.Value")

    expect_true(is.data.frame(extractCol_one_test))
    expect_equal(nrow(extractCol_one_test), 959)
    expect_equal(ncol(extractCol_one_test), 2)
    expect_equal(names(extractCol_one_test), c("BDL_vs_Sham", "EXT1024_vs_BDL"))

    extractCol_two_test <- extractCol(extractCol_contrastList, colName = "P.Value", robust = FALSE)

    expect_true(is.data.frame(extractCol_two_test))
    expect_equal(nrow(extractCol_two_test), 959)
    expect_equal(ncol(extractCol_two_test), 2)
    expect_equal(names(extractCol_two_test), c("BDL_vs_Sham", "EXT1024_vs_BDL"))

    expect_error(extractCol(extractCol_three_test),
                 regexp = "object 'extractCol_three_test' not found")
    # Testing assert
    msg <- "robust must be a singular logical value. Assigning default value TRUE."
    expect_warning(extractCol(extractCol_contrastList, colName = "P.Value",
                              robust = NULL),
                   regexp = msg)
    expect_warning(extractCol(extractCol_contrastList, colName = "P.Value",
                              robust = "FALSE"),
                   regexp = msg)
    expect_warning(extractCol(extractCol_contrastList, colName = "P.Value",
                              robust = c(FALSE, FALSE)),
                   regexp = msg)
})


test_that('extractCol: incorrect usage', {
    expect_error(extractCol(),
                 regexp = "argument \"contrastList\" is missing, with no default")
})
