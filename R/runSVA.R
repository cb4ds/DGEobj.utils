#' Test for surrogate variables
#'
#' Takes an DGEobj from runVoom and tests for surrogate variables. Adds a new
#' design matrix to the DGEobj with the surrogate variable columns appended using cbind.
#' runVoom should then be run again with the new design matrix to complete the
#' analysis.
#'
#' @param dgeObj A DGEobj with normalized counts and a DesignMatrix.
#' @param designMatrixName The itemName of the design matrix in DGEobj.
#' @param n.sv  Optional; Use to override the default n.sv returned by num.sv
#'    for the number of SV to analyze.
#' @param method Method passed to num.sv. Supports "leek" or "be". (Default =
#'   "leek")
#'
#' @return dgeObj containing an updated design table, the svobj and a new design
#'   matrix.
#'
#' @examples
#' \dontrun{
#'    myDGEobj <- runSVA(myDGEobj)
#' }
#'
#' @importFrom sva sva num.sv
#' @import magrittr
#' @importFrom assertthat assert_that
#' @importFrom DGEobj getItem addItem
#' @importFrom stats model.matrix as.formula
#'
#' @export
runSVA <- function(dgeObj,
                   designMatrixName,
                   n.sv,
                   method = "leek") {

    assertthat::assert_that(!missing(dgeObj),
                            "DGEobj" %in% class(dgeObj),
                            with(dgeObj, exists("design")),
                            msg = "dgeObj must be specified, be of class 'DGEobj', and should have a 'design' attribute.")
    assertthat::assert_that(!missing(designMatrixName),
                            "character" %in% class(designMatrixName),
                            with(dgeObj, exists(designMatrixName)),
                            msg = "designMatrixName must be specified, should be of class 'character', and must exist as an attribute on the dgeObj.")
    assertthat::assert_that(tolower(method) %in% c("leek", "be"),
                            msg = "method must be one of 'leek' or 'be'.")

    method <- tolower(method)

    # Set up a NullFormula and DesignMatrix
    NullFormula = "~ 1"
    Design = DGEobj::getItem(dgeObj, "design")
    NullDesignMatrix = stats::model.matrix(as.formula(NullFormula), Design)

    log2cpm <- convertCounts(dgeObj$counts, unit = "cpm", log = TRUE, normalize = "tmm")
    designMatrix <- DGEobj::getItem(dgeObj, designMatrixName)
    if (missing(n.sv)) {
        n.sv <- sva::num.sv(log2cpm, designMatrix, method = method)
    } else { # can't have n.sv > number of residual degrees of freedom
        rdf <- ncol(dgeObj) - ncol(designMatrix)
        if (n.sv > rdf)
            n.sv <- rdf
    }
    svobj <- suppressWarnings(sva::sva(log2cpm, designMatrix, NullDesignMatrix, n.sv = n.sv))

    # Pull out the surrogate variables
    sv <- svobj$sv

    if (svobj$n.sv > 0) {
        # Give them a colname
        colnames(sv) <- paste("sv", 1:ncol(sv), sep = "")

        # Add the SVA colums to the DesignMatrix
        designMatrix_SVA <- cbind(designMatrix, sv)

        # Capture the function call
        FunArgs <- match.call()

        dgeObj <- addItem(dgeObj, svobj, paste(designMatrixName, "_svobj", sep = ""),
                          "svobj",
                          funArgs = FunArgs,
                          parent = designMatrixName)

        # Save the new designMatrix
        dgeObj <- addItem(dgeObj, designMatrix_SVA, paste(designMatrixName, "_sva", sep = ""),
                          "designMatrix",
                          funArgs = FunArgs,
                          parent = designMatrixName)
        # Add the SV columns to the Design table
        dgeObj$design <- cbind(dgeObj$design, sv)

    } else .tsmsg("No Surrogate Variables Found. DGEobj is unchanged.")

    return(dgeObj)
}
