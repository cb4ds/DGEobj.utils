#' Run a power analysis on counts and design matrix
#'
#' Take a counts matrix and design matrix and return a power analysis using
#' the RNASeqPower package. The counts matrix should be pre-filtered to remove
#' non-expressed genes using an appropriate filtering criteria. The design matrix
#' should describe the major sources of variation so the procedure can dial
#' out those known effects for the power calculations.
#'
#' If return = "dataframe," the function will return a tall skinny dataframe
#' of power calculations for various requested combinations of N and significance
#' thresholds.  If return = "plots" or "both", a list is returned with two
#' ggplots (plots) or the plots plus the dataframe (both).
#'
#' @param countsMatrix A counts matrix. (Required)
#' @param designMatrix A design matrix. (Required)
#' @param depth A set of depth to use in the calculations.  The default depths of
#'        c(10, 100, 1000) respectively represent a detection limit, below average
#'        expression, and median expression levels, expressed in readcount units.
#' @param N A set of N value to report power for. (Default = c(3, 6, 10, 20))
#' @param FDR FDR thresholds to filter for for FDR vs. Power graph. (Default = c(0.05, 0.1))
#' @param effectSize A set of fold change values to test. (Default = c(1.2, 1.5, 2))
#' @param return One of "dataframe", "plots", "both". (Default = "both").
#'        Two plots are generated; a ROC curve (FDR vs. Power) and a plot of N vs. Power.
#'
#' @return A list of result objects as defined by the "return" argument.
#'
#' @examples
#' \dontrun{
#'    myPowerResults <- runPower(countsMatrix, designMatrix)
#' }
#'
#' @import magrittr
#' @importFrom RNASeqPower rnapower
#' @importFrom edgeR estimateDisp DGEList calcNormFactors aveLogCPM
#' @importFrom dplyr filter
#' @importFrom stats approx
#'
#' @export
runPower <- function(countsMatrix,
                     designMatrix,
                     depth = c(10, 100, 1000),
                     N = c(3, 6, 10, 20),
                     FDR = c(0.05, 0.1),
                     effectSize = c(1.2, 1.5, 2),
                     return = "both") {

    # Fit the BCV data and define the BCV for each depth requested.
    # Estimate dispersion
    dgelist <- countsMatrix %>%
        as.matrix() %>%
        edgeR::DGEList() %>%
        edgeR::calcNormFactors() %>%
        edgeR::estimateDisp(design = designMatrix, robust = TRUE)

    # Get a fitted CV values for each input value of depth
    # BCV is the sqrt of Dispersion
    GeoMeanLibSize <- dgelist$samples$lib.size %>% log %>% mean %>% exp
    depth_avelogcpm <- edgeR::aveLogCPM(depth, GeoMeanLibSize)
    depthBCV <- sqrt(approx(dgelist$AveLogCPM, dgelist$trended.dispersion,
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
    cnames <- colnames(pdat)

    for (D in depth) {
        cv <- depthBCV[D == depth]
        for (Nf in n)
            for (E in effectSize)
                for (A in alpha) {
                    P <- RNASeqPower::rnapower(depth = D, n = Nf, cv = cv, effect = E, alpha = A)
                    pdat <- rbind(pdat, c(depth = D, n = Nf, effect = E, alpha = A, powerVal = P))
                }
    }
    colnames(pdat) <- cnames

    result <- list()
    if (tolower(return) %in% c("dataframe", "both")) {
        result$PowerData <- pdat
        colnames(result$PowerData) <- c("depth", "n", "effect", "alpha", "power")
    }

    if (tolower(return) %in% c("plot", "both")) {
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

        result$ROC <- roc

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

        result$NvP <- NvP
    }

    return(result)
}
