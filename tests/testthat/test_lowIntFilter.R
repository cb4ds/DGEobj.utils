context("DGEobj.utils - tests for lowIntFilter.R functions")


test_that('lowIntFilter: lowIntFilter()', {
    lowIntFilter_one_test <- lowIntFilter(t_obj1, countThreshold = 10)
    expect_s3_class(lowIntFilter_one_test, "DGEobj")
    expect_equal(nrow(lowIntFilter_one_test$counts), 959)

    #verbose is enabled
    lowIntFilter_one_test <- lowIntFilter(t_obj1, verbose = TRUE, zfpkmThreshold = -3.0)
    expect_s3_class(lowIntFilter_one_test, "DGEobj")
    expect_equal(nrow(lowIntFilter_one_test$counts), 950)

    # NULL gene length
    lowIntFilter_one_test <- lowIntFilter(t_obj1, tpmThreshold = 1, verbose = TRUE)
    expect_error(lowIntFilter(t_obj1,
                              tpmThreshold = 1,
                              verbose = TRUE,
                              geneLength = t_obj1$geneData$ExonLength[1:100]),
                 regexp = "geneLength must be specified and should be the same length as the number of rows in counts.")
    expect_s3_class(lowIntFilter_one_test, "DGEobj")
    expect_equal(nrow(lowIntFilter_one_test$counts), 959)
    # test TPM threshold
    lowIntFilter_one_test <- lowIntFilter(t_obj1,
                                          countThreshold = 10,
                                          verbose = TRUE)
    expect_s3_class(lowIntFilter_one_test, "DGEobj")
    expect_equal(nrow(lowIntFilter_one_test$counts), 959)

    lowIntFilter_two_test <- lowIntFilter(t_obj1, zfpkmThreshold = -3.0)

    expect_s3_class(lowIntFilter_two_test, "DGEobj")
    expect_equal(nrow(lowIntFilter_two_test$counts), 950)

    lowIntFilter_three_test <- lowIntFilter(t_obj1,
                                            fpkThreshold = 5,
                                            verbose = TRUE)

    expect_s3_class(lowIntFilter_three_test, "DGEobj")
    expect_equal(nrow(lowIntFilter_three_test$counts), 936)

    lowIntFilter_four_test <- lowIntFilter(t_obj1, countThreshold = 10, zfpkmThreshold = -3.0, sampleFraction = 0.75)

    expect_s3_class(lowIntFilter_four_test, "DGEobj")
    expect_equal(nrow(lowIntFilter_four_test$counts), 903)

    expect_error(lowIntFilter(lowIntFilter_five_test),
                 regexp = "object 'lowIntFilter_five_test' not found")
    expect_error(lowIntFilter(NULL),
                 regexp = "dgeObj must be of class 'DGEobj'.")
    expect_error(lowIntFilter(),
                 regexp = "dgeObj must be of class 'DGEobj'.")
})


test_that('lowIntFilter: incorrect usage', {
    expect_error(lowIntFilter(t_obj1, zfpkmThreshold = 3.0, tpmThreshold = 1),
                 regexp = "Must use zfpkmThreshold or tpmThreshold, but not both.")
})
