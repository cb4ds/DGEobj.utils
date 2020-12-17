#' Create tidy contrast list of topTable contrast dfs
#'
#' Takes a DGEobj or contrast list as input and merges them into one tidy dataframe.
#' A contrast list is simply a list of topTable contrast dataframes. For
#' example, DGEobj::getType(myDGEobj, "topTable") would retrieve a list of
#' contrasts from a DGEobj.  The contrast list must be a named list as the
#' contrast names are used during the merge operation.
#'
#' The input may or may not have rownames. If supplied rownameColumn does not exist as a colname in
#' the dataframes, it is created from the rownames. In tidy style, the output will have no rownames.
#'
#' The contrast names will be used as a new column in the tidy output format.
#'
#' @param DGEdata A DGEobj or named list of contrast dataframes. (Required)
#' @param rownameColumn Name of the rowname column. If a column by this
#'   name does not exist, it is created from the rownames property.
#' @param includeColumns A character vector of columns to include in the output.
#'   (Default = colnames of the first contrast)
#'
#' @return A dataframe with merged contrast data.
#'
#' @examples
#' \dontrun{
#'   # Get contrasts directly from a DGEobj
#'   myMergedTidyDF <- tidyContrasts(myDGEobj)
#'
#'   # Assemble a list of contrasts from two DGEobjs; just logFC and conf intervals
#'   myContrasts <- c(getType(DGEobj1, "topTable"), getType(DGEobj2, "topTable"))
#'   myMergedTidyDF <- tidyContrasts(myContrasts, includeColumns = c("logFC", "CI.R", "CI.L"))
#' }
#'
#' @importFrom assertthat assert_that
#' @importFrom dplyr select
#' @importFrom tibble rownames_to_column
#' @importFrom DGEobj getType
#'
#' @export
tidyContrasts <- function(DGEdata,
                          rownameColumn = "rownames",
                          includeColumns) {

    assertthat::assert_that(any(c("DGEobj", "list") %in% class(DGEdata)),
                            msg = "DGEdata must be of class 'DGEobj' or 'list'.")

    if ("DGEobj" %in% class(DGEdata)) {
        DGEdata <- DGEobj::getType(DGEdata, "topTable")
        assertthat::assert_that(!(length(DGEdata) == 0),
                                msg = "No topTable dataframes found in DGEdata. Please specify a DGEdata that contains topTable dataframes.")

    }

    # Make sure list contains only dataframes
    assertthat::assert_that(all(sapply(DGEdata, class) == "data.frame"),
                            msg = "DGEdata must only contain dataframes.")

    # Make sure each df has a name
    minNameLen <- min(sapply(names(DGEdata), nchar))
    assertthat::assert_that(!(minNameLen == 0),
                            msg = "All dataframes in DGEdata must be named (it must be a named list.)")

    # Set default columns
    if (missing(includeColumns)) {
        includeColumns <- colnames(DGEdata[[1]])
    }

    # Find the common set of columns present in all dataframes
    commonColumns <- colnames(DGEdata[[1]])
    for (i in 2:length(DGEdata))
        commonColumns <- intersect(commonColumns, colnames(DGEdata[[i]]))

    # Make sure user-requested columns are present
    if (!all(includeColumns %in% commonColumns)) {
        warning("Some requested columns are not present in all dataframes.")
    }
    commonColumns <- intersect(commonColumns, includeColumns)

    # Does the rownameColumn exist in df1?
    if (!rownameColumn %in% colnames(DGEdata[[1]])) {
        # Move rownames to rownameColumn
        DGEdata <- lapply(DGEdata, tibble::rownames_to_column, var = rownameColumn)
        commonColumns <- c(rownameColumn, commonColumns)
    }

    # Reduce all dataframes to the selected columns
    # This also insures same column order in each df
    DGEdata <- lapply(DGEdata, dplyr::select, commonColumns)

    # Add a contrast name column to each DF
    for (name in names(DGEdata)) {
        DGEdata[[name]]["Contrast"] <- name
    }

    # Now merge the dataframes vertically
    DGEdata <- do.call(rbind, DGEdata)

    return(DGEdata)
}
