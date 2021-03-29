context("DGEobj.utils - tests for isoformFrac.R functions")


test_that("isoformFrac.R: isoformFrac()", {
    expect_s3_class(t_obj1, "DGEobj")
    expect_error(isoformFrac(t_obj1$DGEList),
                 regexp = "dgeObj must be of class 'DGEobj.")
    expect_error(isoformFrac(t_obj1),
                 regexp = "The levels attribute of dgeObj must be 'isoform'.")

    dgeObj <- addItem(dgeObj   = t_obj1,
                      item     = t_obj1$geneData,
                      itemName = "isoformData",
                      itemType = "meta")

    attr(dgeObj, "level") <- "isoform"
    iso_data <- isoformFrac(dgeObj)
    expect_equal(dim(iso_data), c(959, 48))
})
