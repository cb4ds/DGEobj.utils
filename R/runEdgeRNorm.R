#' Run edgeR TMM normalization on DGEobj
#'
#' Returns a DGEobj containing DGEList object representing the result of
#'  edgeR TMM normalization.
#'
#' @param dgeObj A DGEobj containing counts, design data, and gene annotation.
#' @param normMethod One of "TMM", "RLE", "upperquartile", or "none". (Default = "TMM")
#' @param includePlot Enable returning a "canvasXpress" or "ggplot" bar plot of the norm.factors produced (Default = FALSE).
#'  Possible values to pass:
#'  \itemize{
#'   \item \strong{FALSE or NULL}: Disable plot
#'   \item \strong{TRUE or "canvasXpress"}: returns "canvasXpress" bar plot of the norm.factors produced.
#'   \item \strong{"ggplot"}: returns "ggplot" bar plot of the norm.factors produced.
#' }
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
#' @importFrom canvasXpress canvasXpress
#'
#' @export
runEdgeRNorm <- function(dgeObj,
                         normMethod = "TMM",
                         includePlot = FALSE,
                         plotLabels = NULL) {
    funArgs <- match.call()
    assertthat::assert_that(class(dgeObj) == "DGEobj",
                            msg = "dgeObj must be of class 'DGEobj'.")
    assertthat::assert_that(!is.null(normMethod),
                            is.character(normMethod),
                            length(normMethod) == 1,
                            tolower(normMethod) %in% c("tmm", "rle", "upperquartile", "none"),
                            msg = "normMethod must be only one of the following values 'TMM', 'RLE', 'upperquartile', 'none'.")
    if (is.null(includePlot)) {
        plot_type <- "none"
    } else if (is.logical(includePlot) && length(includePlot) == 1) {
        plot_type <- ifelse(includePlot, "canvasxpress", "none")
    } else if (is.character(includePlot) && length(includePlot) == 1) {
        if (tolower(includePlot) %in% c("canvasxpress", "ggplot")) {
            plot_type <- tolower(includePlot)
        } else {
            warning("includePlot must be only one of the following values TRUE, FALSE, 'canvasXpress' or 'ggplot'.  Assigning default value FALSE.")
            plot_type <- "none"
        }
    } else {
        warning("includePlot must be only one of the following values TRUE, FALSE, 'canvasXpress' or 'ggplot'.  Assigning default value FALSE.")
        plot_type <- "none"
    }
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
    if (plot_type != "none") {
        if (!is.null(plotLabels) && length(plotLabels) == ncol(dgeObj)) {
            labels <-  plotLabels
            angle  <-  45
        } else {
            if (!is.null(plotLabels) && length(plotLabels) != ncol(dgeObj)) {
                warning(paste("plotLabels must be a character vectore and",
                              "its length must be equal to dgeobj number of columns.",
                              "Assiging default values from 1 to dgeobj columns number."))
            }
            labels <- 1:ncol(dgeObj)
            angle  <- ifelse(plot_type == "canvasxpress", 90, 0)
        }
        plot_data <- data.frame(row.names = factor(labels),
                         Norm.Factors = MyDGElist$samples$norm.factors)
    }

    if (plot_type == "canvasxpress") {
        plot <- canvasXpress::canvasXpress(data = as.data.frame(t(plot_data)),
                                           graphOrientation = "vertical",
                                           graphType = "Bar",
                                           showLegend = FALSE,
                                           smpLabelRotate = angle,
                                           smpTitle = "Samples",
                                           theme = "CanvasXpress",
                                           widthFactor = 1.5,
                                           title = "Normalization Factors",
                                           xAxisTitle = "Norm Factors",
                                           color    = "dodgerblue3",
                                           decorations = list(line = list(list(value = 1,
                                                                               width = 2,
                                                                               color = "rgb(255,0,0)"))),
                                           setMinX = 0)
        list(dgeObj = dgeObj, plot = plot)
    } else if (plot_type == "ggplot") {
        plot <- ggplot(plot_data, aes(x = labels, y = Norm.Factors)) +
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
        list(dgeObj = dgeObj, plot = plot)
    } else {
        dgeObj
    }
}
