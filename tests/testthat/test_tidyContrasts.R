context("DGEobj.utils - tests for tidyContrasts.R functions")


test_that('tidyContrasts: tidyContrasts()', {
    tidyContrast_one_test <- tidyContrasts(t_obj1)

    expect_true(is.data.frame(tidyContrast_one_test))
    expect_equal(nrow(tidyContrast_one_test), 3836)
    expect_equal(ncol(tidyContrast_one_test), 15)

    tidyContrast_two_test <- tidyContrasts(t_obj1,
                                           rownameColumn = "rownameColumn")

    expect_true(is.data.frame(tidyContrast_two_test))
    expect_equal(nrow(tidyContrast_two_test), 3836)
    expect_equal(ncol(tidyContrast_two_test), 15)
    expect_true("rownameColumn" %in% names(tidyContrast_two_test))

    expect_warning(tidyContrasts(t_obj1,
                                 includeColumns = c("rownames", "logFC", "CI.L", "CI.R")),
                   regexp = "Some requested columns are not present in all dataframes.")

    expect_error(tidyContrasts(dgeObj_test),
                 regexp = "object 'dgeObj_test' not found")
})


test_that('tidyContrasts: incorrect usage', {
    expect_error(tidyContrasts(),
                 regexp = "argument \"DGEdata\" is missing, with no default")
})
