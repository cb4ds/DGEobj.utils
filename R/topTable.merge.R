#' Merge specified topTable df cols
#'
#' Take a list of topTable dataframes and consolidate output for specified
#' columns. Should work on any named list of dataframes where each member of the list
#' has the same columns.
#'
#' @param contrastList A named list of topTable data.frames which all have the same colnames and same row counts.
#' The dataframes in the list should have rownames (geneIDs).
#' @param colNames The list of column names of the data column to extract to a
#'   matrix (Default = c("logFC", "AveExpr", "P.Value", "adj.P.Val"))
#' @param digits Number of decimal places for specified columns. Should be same
#' length as colNames. (Default = c(2, 2, 4, 3)). If one value supplied, it is used
#' for all columns.
#'
#' @return A matrix containing the extracted columns.
#'
#' @examples
#' \dontrun{
#'     myContrastTable <- topTable.merge(topTablelist)
#' }
#'
#' @import magrittr
#' @importFrom stringr str_c
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom dplyr left_join
#'
#' @export
topTable.merge <- function(contrastList,
                           colNames = c("logFC",
                                        "AveExpr",
                                        "P.Value",
                                        "adj.P.Val"),
                           digits = c(2, 2, 4, 3)) {

    assertthat::assert_that(!missing(contrastList),
                            "list" %in% class(contrastList),
                            "data.frame" %in% class(contrastList[[1]]),
                            !is.null(names(contrastList)),
                            msg = "contrastList must be specified, be of class 'list' and be a named list specifically, and include items of class 'data.frame'.")
    assertthat::assert_that(length(digits) %in% c(1, length(colNames)),
                            msg = "digits must be either of length 1 or the same length as colNames.")

    if (length(digits) == 1) { # Expand the digits vector
        digits <- rep(digits, length(colNames))
    }

    # Add contrast suffix to each colname
    contrastNames <- names(contrastList)

    # Get the first set of columns
    dat <- extractCol(contrastList, colName = colNames[1], robust = TRUE) %>%
        as.data.frame
    colnames(dat) <- stringr::str_c(colNames[1], "_", colnames(dat))
    dat <- round(dat, digits[1])
    dat <- dat %>%
        tibble::rownames_to_column(var = "rowid")

    if (length(colNames)  > 1) {
        for (i in 1:length(colNames)) {
            dat2 <- extractCol(contrastList, colName = colNames[i], robust = TRUE) %>%
                as.data.frame
            # Add datatype as prefix on colname e.g. logFC_contrastname
            colnames(dat2) <- stringr::str_c(colNames[i], "_", colnames(dat2))
            dat2 <- round(dat2, digits[i])
            dat2 <- dat2 %>%
                tibble::rownames_to_column(var = "rowid")
            if (i == 1) {
                dat <- dat2
            } else {
                dat <- dat %>%
                    dplyr::left_join(dat2, by = "rowid")
            }
        }
    }

    dat <- dat %>%
        tibble::column_to_rownames(var = "rowid")
    return(dat)
}
