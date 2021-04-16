#' Generate matrix of isoform fraction data
#'
#' Takes a DGEobj as input (transcript level data) and returns a matrix
#' containing isoform fraction data.
#'
#' Isoform fraction is calculated using length normalized data (FPKM or TPM).
#' Length normalized data is required because different isoforms have different
#' total exon lengths. If FPKM is specified, a normalization can be specified
#' via edgeR::calcNormFactors.
#'
#' Isoform fraction is calculated
#' as the isoform intensity divided by the summed gene intensity for all
#' isoforms of a given gene.
#'
#' TPM or FPKM are calculated directly from counts using all data in the DGEobj.
#' Performing low intensity filtering at the gene level before
#' running isoformFrac is recommended.
#'
#' @param dgeObj  An isoform level DGEobj created by function initDGEobj().
#'   Counts and isoformData must be present in the DGEobj (Required).
#'   isoformData$ExonLength must be present also.
#' @param dataType One of "fpkm" or "tpm" (Default = "fpkm")
#' @param normalize Default = "TMM" and invokes TMM normalization. Other allowed
#'   values are: "RLE", "upperquartile", "none". Invokes edgeR::calcNormFactors for
#'   normalization.  Only invoked when dataType = "fpkm".  This is because
#'   applying TPM overrides any prior column scaling.
#'
#' @return An isoform fraction dataframe
#'
#' @examples
#' \dontrun{
#'    myDGEobj <- isoformFrac(myDGEobj)
#' }
#'
#' @importFrom dplyr group_by mutate %>%
#' @importFrom assertthat assert_that
#' @importFrom tidyr gather spread
#' @importFrom DGEobj getItem addItem
#'
#' @export
isoformFrac <- function(dgeObj,
                        dataType = "fpkm",
                        normalize = "tmm") {

    assertthat::assert_that(!missing(dgeObj),
                            !is.null(dgeObj),
                            "DGEobj" %in% class(dgeObj),
                            msg = "dgeObj must be of class 'DGEobj.")
    assertthat::assert_that(attr(dgeObj, "level") == "isoform",
                            msg = "The levels attribute of dgeObj must be 'isoform'.")
    assertthat::assert_that(!is.null(dgeObj$isoformData$ExonLength),
                            msg = "An ExonLength column must be present in the isoformData table. ")

    if (any(is.null(dataType),
            !is.character(dataType),
            length(dataType) != 1,
            !tolower(dataType) %in% c("fpkm", "tpm"))) {
        warning("dataType must be only a singular value from 'fpkm', 'tpm'. Assigning default value 'fpkm'")
        dataType  <- "fpkm"
    }

    if (any(is.null(normalize),
            !is.character(normalize),
            length(normalize) != 1,
            !tolower(normalize) %in% c("tmm", "rle", "upperquartile", "none"))) {
        warning("normalize must be only a singular value from 'TMM', 'RLE', 'upperquartile', 'none'. Assigning default value 'TMM'")
        normalize  <- "tmm"
    }

    # Calculate sum of isoforms for each gene and sample
    counts      <- DGEobj::getItem(dgeObj, "counts")
    isoformData <- DGEobj::getItem(dgeObj, "isoformData")

    omicData <- switch(tolower(dataType),
                       "fpkm" = convertCounts(counts,
                                              unit = "fpkm",
                                              geneLength = isoformData$ExonLength,
                                              normalize = normalize),
                       "tpm" = convertCounts(counts,
                                             unit = "TPM",
                                             geneLength = isoformData$ExonLength,
                                             normalize = "none")
    ) %>% as.data.frame()

    omicData$GeneID <- isoformData$rgd_symbol
    omicData$TranscriptID <- rownames(omicData)

    # Calculate isoform fraction
    omictidy <- tidyr::gather(omicData,
                              key   = "sample",
                              value = "intensity",
                              -.data$GeneID, -.data$TranscriptID) %>%
        dplyr::group_by(.data$sample, .data$GeneID) %>%
        dplyr::mutate(geneTotal = sum(.data$intensity),
                      isofrac = .data$intensity / .data$geneTotal)

    # Drop uneeded columns
    omictidy$intensity <- NULL
    omictidy$geneTotal <- NULL

    # Now spread to an isoformPct matrix
    IsoformFrac <- tidyr::spread(omictidy, .data$sample, .data$isofrac) %>% as.data.frame
    # Set row names to transcript ID and remove ID columns
    rownames(IsoformFrac) <- IsoformFrac$TranscriptID
    IsoformFrac$GeneID <- NULL
    IsoformFrac$TranscriptID <- NULL

    return(IsoformFrac)
}
