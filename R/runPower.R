#' Run a power analysis on counts and design matrix
#'
#' Take a counts matrix and design matrix and return a power analysis using
#' the RNASeqPower package. The counts matrix should be pre-filtered to remove
#' non-expressed genes using an appropriate filtering criteria. The design matrix
#' should describe the major sources of variation so the procedure can dial
#' out those known effects for the power calculations.
#'
#' If includePlots = FALSE (the default) or NULL, the function will return a tall skinny dataframe
#' of power calculations for various requested combinations of N and significance
#' thresholds.
#'
#' If includePlots = TRUE, "canvasXpress" or "ggplot", a list is returned with an additional two
#' "canvasXpress" or ggplots (plots) to the dataframe.
#'
#' @param countsMatrix A counts matrix or dataframe of numeric data. (Required)
#' @param designMatrix A design matrix or dataframe of numeric data. (Required)
#' @param depth A set of depth to use in the calculations.  The default depths of
#'        c(10, 100, 1000) respectively represent a detection limit, below average
#'        expression, and median expression levels, expressed in readcount units.
#' @param N A set of N value to report power for. (Default = c(3, 6, 10, 20))
#' @param FDR FDR thresholds to filter for for FDR vs. Power graph. (Default = c(0.05, 0.1))
#' @param effectSize A set of fold change values to test. (Default = c(1.2, 1.5, 2))
#' @param includePlots controls adding tow plots to the returned dataframe (Default = FALSE).
#'        The two plots are; a ROC curve (FDR vs. Power) and a plot of N vs. Power.
#'        Possible values to pass:
#'        \itemize{
#'          \item \strong{FALSE or NULL}: Disable plots
#'          \item \strong{TRUE or "canvasXpress"}: returns "canvasXpress" plots.
#'          \item \strong{"ggplot"}: returns "ggplot" plots.}
#'
#' @return a dataframe of power calculations or a list of the dataframe and defined plots as defined by the "includePlots" argument.
#'
#' @examples
#' \dontrun{
#'    myPowerResults <- runPower(countsMatrix, designMatrix)
#' }
#'
#' @import magrittr ggplot2
#' @importFrom RNASeqPower rnapower
#' @importFrom edgeR estimateDisp DGEList calcNormFactors aveLogCPM
#' @importFrom dplyr filter
#' @importFrom stats approx
#' @importFrom assertthat assert_that
#'
#' @export
runPower <- function(countsMatrix,
                     designMatrix,
                     depth = c(10, 100, 1000),
                     N = c(3, 6, 10, 20),
                     FDR = c(0.05, 0.1),
                     effectSize = c(1.2, 1.5, 2),
                     includePlots = FALSE) {
    assertthat::assert_that(!missing(countsMatrix),
                            !is.null(countsMatrix),
                            class(countsMatrix)[[1]] %in% c("matrix","data.frame"),
                            msg = "countsMatrix must be specified and must be of class matrix or dataframe.")
    assertthat::assert_that(!missing(designMatrix),
                            !is.null(designMatrix),
                            class(designMatrix)[[1]] %in% c("matrix","data.frame"),
                            msg = "designMatrix must be specified and must be of class matrix or dataframe.")
    # Fit the BCV data and define the BCV for each depth requested.
    # Estimate dispersion
    dgelist <- countsMatrix %>%
        as.matrix() %>%
        edgeR::DGEList() %>%
        edgeR::calcNormFactors() %>%
        edgeR::estimateDisp(design = designMatrix, robust = TRUE)

    # Get a fitted CV values for each input value of depth
    # BCV is the sqrt of Dispersion
    GeoMeanLibSize  <- dgelist$samples$lib.size %>% log %>% mean %>% exp
    depth_avelogcpm <- edgeR::aveLogCPM(depth, GeoMeanLibSize)
    depthBCV        <- sqrt(approx(dgelist$AveLogCPM, dgelist$trended.dispersion,
                            xout = depth_avelogcpm, rule = 2, ties = mean)$y)

    n <- seq(min(N),max(N),1)   # For the N vs P plot
    alpha <- seq(0.05, 0.9, 0.05)  # Alpha is FDR levels

    # Initialize a dataframe for the results table
    pdat <- data.frame(depth = double(),
                       n = double(),
                       effect = double(),
                       alpha = double(),
                       powerVal = double(),
                       stringsAsFactors = FALSE)
    for (D in depth) {
        cv <- depthBCV[D == depth]
        for (Nf in n) {
            for (E in effectSize) {
                for (A in alpha) {
                    P    <- RNASeqPower::rnapower(depth = D, n = Nf, cv = cv, effect = E, alpha = A)
                    pdat <- rbind(pdat, c(depth = D, n = Nf, effect = E, alpha = A, powerVal = P))
                }
            }
        }
    }
    colnames(pdat) <- c("depth", "n", "effect", "alpha", "power")
    if (is.null(includePlots)) {
        plot_type <- "none"
    } else if (is.logical(includePlots) && length(includePlots) == 1) {
        plot_type <- ifelse(includePlots, "canvasxpress", "none")
    } else if (is.character(includePlots) && length(includePlots) == 1) {
        if (tolower(includePlots) %in% c("canvasxpress", "ggplot")) {
            plot_type <- tolower(includePlots)
        } else {
            warning("includePlots must be only one of the following values TRUE, FALSE, 'canvasXpress' or 'ggplot'.  Assigning default value FALSE.")
            plot_type <- "none"
        }
    } else {
        warning("includePlots must be only one of the following values TRUE, FALSE, 'canvasXpress' or 'ggplot'.  Assigning default value FALSE.")
        plot_type <- "none"
    }

    if (plot_type == "ggplot") {
        rocdat <- dplyr::filter(pdat, n %in% N)
        rocdat$depth <- as.factor(rocdat$depth)

        roc <- ggplot(rocdat, aes(x = alpha, y = powerVal, fill = depth, shape = depth, color = depth)) +
            geom_line(size = 1) +
            scale_x_continuous(breaks = seq(0, 1, 0.2)) +
            scale_y_continuous(breaks = seq(0, 1, 0.2)) +
            facet_grid(effect ~ n, labeller = label_both) +
            ggtitle("ROC curves") +
            xlab("\nFDR") +
            ylab("Power") +
            expand_limits(x = 0, y = 0) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
            theme_gray(18)

        # N vs Power
        # Filter to just a few FDR thresholds
        ndat <- dplyr::filter(pdat, alpha %in% FDR)

        ndat$depth <- as.factor(ndat$depth)
        ndat$FDR <- ndat$alpha

        NvP <- ggplot(ndat, aes(x = n, y = powerVal, fill = depth, shape = depth, color = depth)) +
            geom_line(size = 1) +
            scale_y_continuous(breaks = seq(0, 1, 0.2)) +
            facet_grid(FDR ~ effect, labeller = label_both) +
            ggtitle("N vs Power") +
            xlab("\nN") +
            ylab("Power") +
            expand_limits(x = 0, y = 0) +
            theme_gray()

        list(PowerData = pdat, ROC = roc, NvP = NvP)
    } else {
        pdat
    }
}
