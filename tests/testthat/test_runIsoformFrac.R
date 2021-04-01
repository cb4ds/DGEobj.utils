context("DGEobj.utils - tests for isoformFrac.R functions")


test_that("isoformFrac.R: isoformFrac()", {
    expect_s3_class(t_obj1, "DGEobj")
    msg <- "dgeObj must be of class 'DGEobj."
    expect_error(isoformFrac(),
                 regexp = msg)
    expect_error(isoformFrac(NULL),
                 regexp = msg)
    expect_error(isoformFrac(t_obj1$DGEList),
                 regexp = msg)
    expect_error(isoformFrac(t_obj1),
                 regexp = "The levels attribute of dgeObj must be 'isoform'.")

    dgeObj <- addItem(dgeObj   = t_obj1,
                      item     = t_obj1$geneData,
                      itemName = "isoformData",
                      itemType = "meta")

    attr(dgeObj, "level") <- "isoform"
    iso_data <- isoformFrac(dgeObj)
    expect_equal(dim(iso_data), c(959, 48))
    ## dataType
    msg <- "dataType must be only a singular value from 'fpkm', 'tpm'. Assigning default value 'fpkm'"
    expect_warning(isoformFrac(dgeObj   = dgeObj,
                               dataType = NULL),
                   regexp = msg)
    expect_warning(isoformFrac(dgeObj   = dgeObj,
                               dataType = 1),
                   regexp = msg)
    expect_warning(isoformFrac(dgeObj   = dgeObj,
                               dataType = "NULL"),
                   regexp = msg)
    expect_warning(isoformFrac(dgeObj   = dgeObj,
                               dataType = c("fpkm", "fpkm")),
                   regexp = msg)
    ## normalize
    msg <- "normalize must be only a singular value from 'TMM', 'RLE', 'upperquartile', 'none'. Assigning default value 'TMM'"
    expect_warning(isoformFrac(dgeObj   = dgeObj,
                               normalize = NULL),
                   regexp = msg)
    expect_warning(isoformFrac(dgeObj   = dgeObj,
                               normalize = 1),
                   regexp = msg)
    expect_warning(isoformFrac(dgeObj   = dgeObj,
                               normalize = "NULL"),
                   regexp = msg)
    expect_warning(isoformFrac(dgeObj   = dgeObj,
                               normalize = c("TMM", "TMM")),
                   regexp = msg)
})
