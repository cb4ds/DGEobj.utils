context("DGEobj.utils - tests for convertCounts.R functions")


test_that("convertCounts.R: convertCounts()", {
    # CPM
    ## Full Count Matrix "none" normalization
    count_matrix <- convertCounts(counts      = t_obj1$counts_orig,
                                  unit        = "CPM",
                                  log         = FALSE,
                                  normalize   = "none",
                                  prior.count = NULL)
    expect_true("matrix" %in% class(count_matrix))
    expect_identical(dim(t_obj1$counts_orig), dim(count_matrix))
    ## Full Count Matrix missing normalization
    count_matrix <- convertCounts(counts      = t_obj1$counts_orig,
                                  unit        = "CPM",
                                  normalize = NULL)
    expect_true("matrix" %in% class(count_matrix))
    expect_identical(dim(t_obj1$counts_orig), dim(count_matrix))
    ## Full Count Matrix "true" normalization, testing backward compatibility
    count_matrix <- convertCounts(counts      = t_obj1$counts_orig,
                                  unit        = "CPM",
                                  log         = FALSE,
                                  normalize   = TRUE,
                                  prior.count = NULL)
    expect_true("matrix" %in% class(count_matrix))
    expect_identical(dim(t_obj1$counts_orig), dim(count_matrix))
    expect_error(convertCounts(unit = "CPM"),
                 msg = "argument \"counts\" is missing, with no default")
    expect_error(convertCounts(counts = t_obj1$counts_orig),
                 msg = "argument \"unit\" is missing, with no default")
    ## Full Count Matrix "TMM" normalization
    count_matrix <- convertCounts(counts      = t_obj1$counts_orig,
                                  unit        = "CPM",
                                  log         = FALSE,
                                  normalize   = "TMM",
                                  prior.count = NULL)
    expect_true("matrix" %in% class(count_matrix))
    expect_identical(dim(t_obj1$counts_orig), dim(count_matrix))
    ## Full Count Matrix "RLE" normalization
    count_matrix <- convertCounts(counts      = t_obj1$counts_orig,
                                  unit        = "CPM",
                                  log         = TRUE,
                                  normalize   = "RLE",
                                  prior.count = NULL)
    expect_true("matrix" %in% class(count_matrix))
    expect_identical(dim(t_obj1$counts_orig), dim(count_matrix))
    ## Subset of the Count Matrix
    count_matrix <- convertCounts(counts = t_obj1$counts_orig[1:100, ], unit = "CPM")
    expect_true("matrix" %in% class(count_matrix))
    expect_identical(dim(t_obj1$counts[1:100,]), dim(count_matrix))

    # TPM
    ## Full Count Matrix normalize "none"
    genelength <- getItem(t_obj1, "geneData")$ExonLength
    ### FALSE log
    count_matrix <- convertCounts(counts      = t_obj1$counts,
                                  unit        = "TPM",
                                  geneLength  = genelength ,
                                  log         = FALSE,
                                  normalize   = "none",
                                  prior.count = NULL)
    expect_true("matrix" %in% class(count_matrix))
    expect_identical(dim(t_obj1$counts), dim(count_matrix))
    ### TRUE log
    count_matrix <- convertCounts(counts      = t_obj1$counts,
                                  unit        = "TPM",
                                  geneLength  = genelength ,
                                  log         = TRUE,
                                  normalize   = "none",
                                  prior.count = NULL)
    expect_true("matrix" %in% class(count_matrix))
    expect_identical(dim(t_obj1$counts), dim(count_matrix))
    ## testing genelength as a matrix
    count_matrix <- convertCounts(counts      = t_obj1$counts,
                                  unit        = "TPM",
                                  geneLength  = as.matrix(genelength) ,
                                  log         = FALSE,
                                  normalize   = "none",
                                  prior.count = NULL)
    expect_true("matrix" %in% class(count_matrix))
    expect_identical(dim(t_obj1$counts), dim(count_matrix))
    expect_error(convertCounts(counts = t_obj1$counts, unit = "TPM"), regexp = "geneLength is required for unit = FPK|FPKM|TPM")
    ## Full Count Matrix normalize "TMM"
    expect_warning(convertCounts(counts = t_obj1$counts,
                                 unit = "TPM",
                                 geneLength = genelength,
                                 normalize = "TMM"),
                   regexp = "TPM normalization overides TMM normalization!")
    ## Full Count Matrix normalize "upperquartile"
    expect_warning(convertCounts(counts = t_obj1$counts,
                                 unit = "TPM",
                                 geneLength = genelength,
                                 normalize = "upperquartile"),
                   regexp = "TPM normalization overides upperquartile normalization!")
    ## Full Count Matrix normalize "TMMWZP"
    expect_warning(convertCounts(counts = t_obj1$counts,
                                 unit = "TPM",
                                 geneLength = genelength,
                                 normalize = "TMMWZP"),
                   regexp = "TPM normalization overides TMMwzp normalization!")
    ## Full Count Matrix normalize "RLE"
    expect_warning(convertCounts(counts = t_obj1$counts,
                                 unit = "TPM",
                                 geneLength = genelength,
                                 normalize = "RLE"),
                   regexp = "TPM normalization overides RLE normalization!")
    expect_warning(convertCounts(counts = t_obj1$counts,
                                 unit = "TPM",
                                 geneLength = genelength,
                                 prior.count = 1,
                                 log = TRUE),
                   regexp = "Using a prior.count for logTPM calculations is not recommended and may produce unpredictable results!")
    ## Subset of the Count Matrix
    count_matrix <- convertCounts(counts = t_obj1$counts_orig[1:100, ],
                                  unit = "TPM",
                                  geneLength  = genelength[1:100])
    expect_true("matrix" %in% class(count_matrix))
    expect_identical(dim(t_obj1$counts[1:100,]), dim(count_matrix))

    # FPK
    ## Full Count Matrix normalize "none"
    count_matrix <- convertCounts(counts      = t_obj1$counts,
                                  unit        = "FPK",
                                  geneLength  = genelength ,
                                  log         = FALSE,
                                  normalize   = "none",
                                  prior.count = NULL)
    expect_true("matrix" %in% class(count_matrix))
    expect_identical(dim(t_obj1$counts), dim(count_matrix))
    ## Full Count Matrix normalize "none" with 0 genelength
    count_matrix <- convertCounts(counts      = t_obj1$counts,
                                  unit        = "FPK",
                                  geneLength  = integer(length(genelength)) ,
                                  log         = FALSE,
                                  normalize   = "none",
                                  prior.count = NULL)
    expect_true("matrix" %in% class(count_matrix))
    expect_identical(dim(t_obj1$counts), dim(count_matrix))
    ## Full Count Matrix normalize "TMM" with 0 genelength
    count_matrix <- convertCounts(counts      = t_obj1$counts,
                                  unit        = "FPK",
                                  geneLength  = integer(length(genelength)) ,
                                  log         = FALSE,
                                  normalize   = "TMM",
                                  prior.count = NULL)
    expect_true("matrix" %in% class(count_matrix))
    expect_identical(dim(t_obj1$counts), dim(count_matrix))
    ## Full Count Matrix normalize "TMM" log enabled
    count_matrix <- convertCounts(counts      = t_obj1$counts,
                                  unit        = "FPK",
                                  geneLength  = genelength,
                                  log         = TRUE,
                                  normalize   = "TMM",
                                  prior.count = NULL)
    expect_true("matrix" %in% class(count_matrix))
    expect_identical(dim(t_obj1$counts), dim(count_matrix))
    expect_error(convertCounts(counts = t_obj1$counts, unit = "FPK"), regexp = "geneLength is required for unit = FPK|FPKM|TPM")
    ## Full Count Matrix normalize "TMM" log disabled
    count_matrix <- convertCounts(counts = t_obj1$counts,
                                 unit = "FPK",
                                 geneLength = genelength,
                                 normalize = "TMM")
    expect_true("matrix" %in% class(count_matrix))
    expect_identical(dim(t_obj1$counts), dim(count_matrix))
    ## Full Count Matrix normalize "RLE"
    count_matrix <- convertCounts(counts = t_obj1$counts,
                                 unit = "FPK",
                                 geneLength = genelength,
                                 normalize = "RLE")
    expect_true("matrix" %in% class(count_matrix))
    expect_identical(dim(t_obj1$counts), dim(count_matrix))

    # FPKM
    ## Full Count Matrix normalize "none"
    count_matrix <- convertCounts(counts      = t_obj1$counts,
                                  unit        = "FPKM",
                                  geneLength  = genelength ,
                                  log         = FALSE,
                                  normalize   = "none",
                                  prior.count = NULL)
    expect_true("matrix" %in% class(count_matrix))
    expect_identical(dim(t_obj1$counts), dim(count_matrix))
    expect_error(convertCounts(counts = t_obj1$counts, unit = "FPKM"), regexp = "geneLength is required for unit = FPK|FPKM|TPM")
    ## Full Count Matrix normalize "TMM"
    count_matrix <- convertCounts(counts = t_obj1$counts,
                                  unit = "FPKM",
                                  geneLength = genelength,
                                  normalize = "TMM")
    expect_true("matrix" %in% class(count_matrix))
    expect_identical(dim(t_obj1$counts), dim(count_matrix))
    ## Full Count Matrix normalize "RLE"
    count_matrix <- convertCounts(counts = t_obj1$counts,
                                  unit = "FPKM",
                                  geneLength = genelength,
                                  normalize = "RLE")
    expect_true("matrix" %in% class(count_matrix))
    expect_identical(dim(t_obj1$counts), dim(count_matrix))

    # Subset of the Count Matrix
    count_matrix <- convertCounts(counts = t_obj1$counts_orig[1:100, ],
                                  unit = "FPK",
                                  geneLength  = genelength[1:100])
    expect_true("matrix" %in% class(count_matrix))
    expect_identical(dim(t_obj1$counts[1:100,]), dim(count_matrix))
    # Testing assert
    ## log
    msg <- "log must be a singular logical value. Assigning default values FALSE"
    expect_warning(convertCounts(counts      = t_obj1$counts_orig,
                                 unit        = "CPM",
                                 log         = NULL),
                   regexp = msg)
    expect_warning(convertCounts(counts      = t_obj1$counts_orig,
                                 unit        = "CPM",
                                 log         = "FALSE"),
                   regexp = msg)
    expect_warning(convertCounts(counts      = t_obj1$counts_orig,
                                 unit        = "CPM",
                                 log         = c(FALSE, FALSE)),
                   regexp = msg)

})

test_that("convertCounts.R: tpm.on.subset()", {
    geneLength <- NULL
    if (attr(t_obj1, "source") == "Omicsoft") {
        geneLength <- getItem(t_obj1, "geneData")$ExonLength
    } else if ("effectiveLength_orig" %in% names(t_obj1)) {
        geneLength <- rowMeans(getItem(t_obj1, "effectiveLength_orig"), na.rm = TRUE)
    }

    # testing level gene
    tpmObj <- tpm.on.subset(t_obj1)
    expect_true("matrix" %in% class(tpmObj))
    # testing level isoform
    isoform_dgeObj <- t_obj1
    isoform_dgeObj <- addItem(dgeObj   = isoform_dgeObj,
                              item     = isoform_dgeObj$geneData_orig,
                              itemName = "isoformData_orig",
                              itemType = "meta")
    attr(isoform_dgeObj, "level") <- "isoform"
    tpmObj <- tpm.on.subset(isoform_dgeObj)
    expect_true("matrix" %in% class(tpmObj))

    expect_error(tpm.on.subset("XYZ"), regexp = "dgeObj should be of class 'DGEobj'.")
    # testing level exon
    exon_dgeObj <- t_obj1
    attr(exon_dgeObj, "level") <- "exon"
    expect_error(tpm.on.subset(exon_dgeObj),
                 regexp = "The level of dgeObj should be of type 'isoform' or type 'gene'.")

    # testing for bad source attribute
    attr(isoform_dgeObj, "source") <- "XYZ"
    expect_error(tpm.on.subset(isoform_dgeObj),
                 regexp = "object 'geneLength' not found")
})

test_that("convertCounts.R: tpm.direct()", {
    genelength <- getItem(t_obj1, "geneData")$ExonLength
    tpmObj <- tpm.direct(t_obj1$counts, geneLength = genelength)
    expect_true("matrix" %in% class(tpmObj))

    # testing count as vector type
    tpmObj <- tpm.direct(counts = genelength, geneLength = genelength)
    expect_true("matrix" %in% class(tpmObj))

    # testing bad genelength parameter
    expect_error(tpm.direct(t_obj1$counts, geneLength = as.data.frame(genelength)),
                 regexp = "The dimensions of countsMatrix and geneLength should match.")

    # testing collapse parameter
    tpmObj <- tpm.direct(t_obj1$counts, geneLength = as.matrix(genelength), collapse = TRUE)
    expect_true("matrix" %in% class(tpmObj))
})
