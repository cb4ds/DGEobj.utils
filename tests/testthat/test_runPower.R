context("DGEobj.utils - tests for runPower.R functions")


test_that("runPower.R: runPower()", {
    # data setup
    designMatrix <- model.matrix(~ 0 + ReplicateGroup, getItem(t_obj1, "design"))

    # no plots
    ## default value
    power_plot <- runPower(countsMatrix = t_obj1$counts, designMatrix = designMatrix)
    expect_s3_class(power_plot, "data.frame")
    ## NULL value
    power_plot <- runPower(countsMatrix = t_obj1$counts, designMatrix = designMatrix, includePlots = NULL)
    expect_s3_class(power_plot, "data.frame")
    ## FALSE value
    power_plot <- runPower(countsMatrix = t_obj1$counts, designMatrix = designMatrix, includePlots = FALSE)
    expect_s3_class(power_plot, "data.frame")
    #expect_equal(dim(power_plot), c())

    # with plots
    ## ggplot
    power_plot <- runPower(countsMatrix = t_obj1$counts, designMatrix = designMatrix, includePlots = "ggplot")
    expect_type(power_plot, "list")
    expect_s3_class(power_plot$ROC, c("gg", "ggplot"))
    expect_s3_class(power_plot$NvP, c("gg", "ggplot"))
    expect_s3_class(power_plot$PowerData, "data.frame")

    # Testing asserts
    ## countsMatrix
    ### Wrong data
    expect_error(runPower(counts = t_obj1),
                 regexp = "countsMatrix must be specified and must be of class matrix or dataframe.")
    ### missing data
    expect_error(runPower(),
                 regexp = "countsMatrix must be specified and must be of class matrix or dataframe.")
    ### NULL data
    expect_error(runPower(counts = NULL),
                 regexp = "countsMatrix must be specified and must be of class matrix or dataframe.")
    # designMatrix
    ### Wrong data
    expect_error(runPower(counts = t_obj1$counts, designMatrix = t_obj1),
                 regexp = "designMatrix must be specified and must be of class matrix or dataframe.")
    ### missing data
    expect_error(runPower(counts = t_obj1$counts),
                 regexp = "designMatrix must be specified and must be of class matrix or dataframe.")
    ### NULL data
    expect_error(runPower(counts = t_obj1$counts, designMatrix = NULL),
                 regexp = "designMatrix must be specified and must be of class matrix or dataframe.")
    ## includePlots
    msg <- "includePlots must be only one of the following values TRUE, FALSE, 'canvasXpress' or 'ggplot'.  Assigning default value FALSE."
    expect_warning(runPower(counts = t_obj1$counts, designMatrix = designMatrix, includePlots = c(TRUE, FALSE)),
                   regexp = msg)
    expect_warning(runPower(counts = t_obj1$counts, designMatrix = designMatrix, includePlots = "abc"),
                   regexp = msg)
    expect_warning(runPower(counts = t_obj1$counts, designMatrix = designMatrix, includePlots = c("canvasXpress", "ggplot")),
                   regexp = msg)
})
