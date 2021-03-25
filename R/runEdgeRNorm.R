#' Run edgeR TMM normalization on DGEobj
#'
#' Returns a DGEobj containing DGEList object representing the result of
#'  edgeR TMM normalization.
#'
#' @param dgeObj A DGEobj containing counts, design data, and gene annotation.
#' @param normMethod One of "TMM", "RLE", "upperquartile", or "none". (Default = "TMM")
#' @param includePlot Enable returning a bar plot of the norm.factors produced. (Default = FALSE)
#' @param plotLabels Sample labels for the plot. Length must equal the number of
#'   samples. (Default = NULL; sample number will be displayed)
#'
#' @return A DGEobj with a normalized DGEList added or a list containing the normalized DGEobj and a plot
#'
#' @examples
#' \dontrun{
#'    myDGEobj <- runEdgeRNorm(myDGEobj)
#' }
#'
#' @import magrittr ggplot2
#' @importFrom edgeR calcNormFactors DGEList
#' @importFrom DGEobj addItem getItem
#' @importFrom assertthat assert_that
#'
#' @export
runEdgeRNorm <- function(dgeObj,
                         normMethod = "TMM",
                         includePlot = FALSE,
                         plotLabels = NULL) {
    funArgs <- match.call()
    assertthat::assert_that(class(dgeObj) == "DGEobj",
                            msg = "dgeObj must be of class 'DGEobj'.")
    MyDGElist  <-  as.matrix(DGEobj::getItem(dgeObj, "counts")) %>%
        edgeR::DGEList() %>%
        edgeR::calcNormFactors(method = normMethod)


    # Capture the DGEList
    itemAttr <- list(normalization = normMethod)
    dgeObj   <- DGEobj::addItem(dgeObj,
                                item = MyDGElist,
                                itemName = "DGEList",
                                itemType = "DGEList",
                                funArgs = funArgs,
                                itemAttr = itemAttr,
                                parent = "counts")

    if (includePlot) { # Bar plot of norm.factors
        # Plot the Norm factors
        if (!is.null(plotLabels) && length(plotLabels == ncol(dgeObj))) {
            x = plotLabels
            angle = 45
        } else {
            x = 1:ncol(dgeObj)
            angle = 0
        }
        df <- data.frame(x = factor(x),
                         Norm.Factors = MyDGElist$samples$norm.factors)
        nfplot <- ggplot(df, aes(x = x, y = Norm.Factors)) +
            geom_bar(stat = "identity",
                     color = "dodgerblue4",
                     fill = "dodgerblue3",
                     width = 0.7) +
            geom_hline(yintercept = 1.0, color = "red") +
            xlab("Samples") +
            ylab("Norm Factors") +
            ggtitle("Normalization Factors") +
            theme_bw(12) +
            theme(axis.text.x = element_text(angle = angle, hjust = 1.0))
    }
    if (includePlot) {
        list(dgeObj = dgeObj, plot = nfplot)
    } else {
        dgeObj
    }
}
