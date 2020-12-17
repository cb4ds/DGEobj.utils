#' Create a tidy intensity dataframe
#'
#' This function takes a dataframe or matrix of intensity values with geneID row names
#' and sample name column names, e.g. from log2cpm.  It uses tidyr::gather to produce a
#' "tidy" dataframe.
#'
#' The rowIDcolname defines the name the resulting gene id column should have and will create the
#' column from the row names, which are typically ensemble gene IDs. If a gene symbol column is
#' already present in the intensity dataframe, it can be specified using rowIDcolname and will
#' be used instead of the rownames.
#'
#' @param intensityObj A dataframe or matrix with row names and col names. (Required)
#' @param rowIDcolname Column name for the rowID (Usually a geneID column).
#'   If it doesn't exist, it will be created from row names. (Required)
#' @param keyColname Defines the column name that will hold the sample names
#'   captured from the col names of the input dataframe.  (Required)
#' @param valueColname Defines the column name for the value column (intensity values). (Required)
#' @param group A grouping variable to define which columns go together (vector
#'   with length = ncol(intensityObj) of the input dataframe).  This is typically the name
#'   of a column from the design table. (Required)
#'
#' @return A tidy dataframe of intensity data.
#'
#' @examples
#' \dontrun{
#'   myTidyIntensity <- tidyIntensity(myLog2CPM,
#'                                    rowIDcolname = "GeneID",
#'                                    keyColname = "Sample"
#'                                    valueColname = "Log2CPM",
#'                                    group = dgeObj$design$ReplicateGroup)
#' }
#'
#' @importFrom assertthat assert_that
#' @importFrom dplyr left_join
#' @importFrom stats setNames
#' @importFrom tidyr gather
#' @importFrom tibble rownames_to_column tibble
#'
#' @export
tidyIntensity <- function(intensityObj,
                          rowIDcolname,
                          keyColname,
                          valueColname,
                          group) {

    assertthat::assert_that(any(c("data.frame", "matrix") %in% class(intensityObj)),
                            msg = "intensityObj must be of class 'data.frame' or 'matrix'.")
    assertthat::assert_that(!missing(rowIDcolname),
                            msg = "rowIDcolname must be specified.")
    assertthat::assert_that(!missing(keyColname),
                            msg = "keyColname must be specified.")
    assertthat::assert_that(!missing(valueColname),
                            msg = "valueColname must be specified.")
    assertthat::assert_that(!missing(group),
                            length(group) == ncol(intensityObj),
                            msg = "group must be specified and should be the same length as the number of columns in intensityObj.")

    if ("matrix" %in% class(intensityObj)) {
        intensityObj <- as.data.frame(intensityObj)
    }
    samples <- colnames(intensityObj)

    # Create a rownames column
    xcol <- ncol(intensityObj) + 1
    intensityObj <- tibble::rownames_to_column(intensityObj, var = rowIDcolname)
    intensityObj <- tidyr::gather(intensityObj, key = !!keyColname, value = !!valueColname, 2:xcol)

    # Join the group info
    groupinfo <- tibble::tibble(keycol = samples, group = group)
    intensityObj <- dplyr::left_join(intensityObj, groupinfo, by = stats::setNames(nm = keyColname, "keycol"))
    return(intensityObj)
}
