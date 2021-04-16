#' Summarize a contrast list
#'
#' Takes a contrast list produced by runContrasts. Defaults are provided to specify columns to
#' summarize and thresholds for each column, though they can be adjusted. A fold change
#' threshold can optionally be specified. The function queries the topTable results and
#' returns a dataframe with the summary results, but only includes gene counts that meet
#' the specified conditions.
#'
#' Any specified column names that don't exist will be ignored. Normally the
#' defaults cover all the p-value and FDR related columns. However, a fcThreshold
#' can be added and the p-value/FDR thresholds can be modified using the fcThreshold
#' and sigThresholds arguments, respectively.
#'
#' @param contrastList A list of topTable dataframes.
#' @param columns Vector of column names to summarize from topTable dataframes.
#'   Default = c("P.Value", "adj.P.Val", "Qvalue", "qvalue.lfdr", "ihw.adj_pvalue")
#' @param sigThresholds Thresholds to use for each column specified in columns
#'   Must be same length at columns argument.
#'   Default = c(0.01, 0.1, 0.1, 0.1, 0.1)
#' @param fcThreshold Fold-change threshold (absolute value, not logged.)
#'
#' @return data.frame with one summary row per contrast.
#'
#' @importFrom assertthat assert_that
#'
#' @examples
#' \dontrun{
#'    # Get a contrast list from a dgeObj
#'    myContrastList <- DGEobj::getType(DGEobj, "topTable")
#'
#'    # All default thresholds, no fold change threshold
#'    mySigSummary <- summarizeSigCounts(myContrastList)
#'
#'    # All defaults with a foldchange threshold
#'    mySigSummary <- summarizeSigCounts(myContrastList, fcThreshold = 1.5)
#'
#'    # Change the p-value and fdr thresholds
#'    mySigSummary <- summarizeSigCounts(myContrastList, sigThresholds = c(0.05, 0.2, 0.2, 0.2, 0.2))
#' }
#'
#' @export
summarizeSigCounts <- function(contrastList,
                               columns = c("P.Value", "adj.P.Val", "Qvalue", "qvalue.lfdr", "ihw.adj_pvalue"),
                               sigThresholds = c(0.01, 0.1, 0.1, 0.1, 0.1),
                               fcThreshold = 0) {

    assertthat::assert_that(length(columns) == length(sigThresholds),
                            msg = "Supplied sigThresholds should be same length as supplied columns.")

    # Functions
    getSigCounts <- function(df,
                             columns = c("P.Value", "adj.P.Val", "Qvalue", "qvalue.lfdr", "ihw.adj_pvalue"),
                             thresholds = c(0.01, 0.1, 0.1, 0.1, 0.1),
                             fcThreshold = 0) {
        # Pass a topTable df, vector of field names, a threshold and optionally a fcThreshold)
        # Then get the sigcounts
        # Return a vector of named values
        CountSig <- function(df, column, threshold, fcThreshold = 0){
            # Supply a field and a threshold.
            # Return the number of samples <= threshold
            # If fcThreshold > 0, also filter on logFC field (convert the fcThreshold to
            # log2 and filter the abs of logFC)
            idx <- df[column] <= threshold
            if (fcThreshold > 0) {
                fcidx <- abs(df$logFC) >= log2(fcThreshold)
                idx <- idx & fcidx
            }
            return(sum(idx))
        }

        # Body getSigCounts
        counts <- list()

        for (i in 1:length(columns)) {
            counts[columns[i]] <- CountSig(df, columns[i], thresholds[i], fcThreshold)
        }
        return(unlist(counts))
    }

    # Reduce tableFields to only ones that exist
    columns <- columns[columns %in% colnames(contrastList[[1]])]
    sigThresholds <- sigThresholds[columns %in% colnames(contrastList[[1]])]

    # Collect the rows in a list
    myrows <- list()
    for (i in 1:length(contrastList)) {
        myrows[[i]] <- getSigCounts(contrastList[[i]], columns, sigThresholds, fcThreshold)
    }

    # Put rows into a matrix
    DF <- do.call("rbind", myrows)

    rownames(DF) <- names(contrastList)
    colnames(DF) <- columns

    return(DF)
}
