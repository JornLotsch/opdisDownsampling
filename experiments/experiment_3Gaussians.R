#' Downsampling Analysis Function
#'
#' Performs downsampling analysis on provided data, generates plots,
#' and computes distribution comparison statistics.
#'
#' @param data_df Data frame with 'Data' and 'Cls' columns
#' @param nSamples Number of downsampling iterations
#' @param Size Target size for downsampled data
#' @param nTrials Number of trials for opdisDownsampling
#' @param CheckRemoved Logical; if TRUE, also optimize the removed part of the data for distribution equality with the original
#' @param path_functions Path to custom distribution comparison functions
#' @return List containing plots and statistical comparison results
downsample_analysis <- function(data_df, nSamples = 10, Size = 1500, nTrials = 1,
                                CheckRemoved = FALSE, OptimizeBetween = FALSE,
                                use_y_limits = FALSE, ylimits_list = NULL, MaxCores = getOption("mc.cores", 2L),
                                path_functions = "/home/joern/Aktuell/DownSamplingStructure/12RLibrary/opdisDownsampling/") {

  # Load required libraries
  library(ggplot2)
  library(ggthemes)
  library(reshape2)
  library(ggbeeswarm)
  library(ggh4x)
  library(cowplot)
  library(patchwork)

  # Source custom functions for distribution comparison
  source(paste0(path_functions, "R/", "dist_amrdd.R"))
  source(paste0(path_functions, "R/", "dist_compare.R"))
  source(paste0(path_functions, "R/", "dist_euclidean.R"))
  source(paste0(path_functions, "R/", "dist_kld.R"))
  source(paste0(path_functions, "R/", "density_smooth_hist.R"))
  source(paste0(path_functions, "R/", "utils.R"))
  source(paste0(path_functions, "experiments/", "pareto_density_estimation_ie.R"))

  # Perform downsampling iterations
  downsampled_results <- lapply(1:nSamples, function(i) {
    cat(paste0("  Iteration ", i, "/", nSamples, "\n"))

    result <- opdisDownsampling::opdisDownsampling(
      Data = data_df$Data,
      Cls = data_df$Cls,
      Size = Size,
      Seed = i + (i - 1) * 1000000,
      nTrials = nTrials,
      CheckRemoved = CheckRemoved,
      OptimizeBetween = OptimizeBetween,
      MaxCores = MaxCores
    )

    list(
      reduced_data = result$ReducedData$Data,
      removed_data = result$RemovedData$Data
    )
  })

  # Extract and organize downsampled data
  reduced_data_df <- do.call(cbind.data.frame, lapply(downsampled_results, "[[", "reduced_data"))
  removed_data_df <- do.call(cbind.data.frame, lapply(downsampled_results, "[[", "removed_data"))
  names(reduced_data_df) <- paste0("Sample_", 1:nSamples)
  names(removed_data_df) <- paste0("Sample_", 1:nSamples)

  # Helper function for PDE calculation
  calculate_pde <- function(data, pareto_radius = NULL) {
    pde_result <- pareto_density_estimation_ie(
      Data = data,
      paretoRadius = pareto_radius,
      kernels = seq(min(data_df$Data), max(data_df$Data), by = 0.001)
    )
    data.frame(x = pde_result$kernels, y = pde_result$paretoDensity)
  }

  # Calculate PDEs for original and downsampled data
  original_pde <- calculate_pde(data_df$Data)

  reduced_pdes <- apply(reduced_data_df, 2, calculate_pde)
  reduced_pdes <- mapply(function(df, name) {
    df$name <- name
    df
  }, reduced_pdes, names(reduced_pdes), SIMPLIFY = FALSE)

  removed_pdes <- apply(removed_data_df, 2, calculate_pde)
  removed_pdes <- mapply(function(df, name) {
    df$name <- name
    df
  }, removed_pdes, names(removed_pdes), SIMPLIFY = FALSE)

  # Prepare data for plotting
  reduced_plot_data <- rbind.data.frame(
    cbind.data.frame(original_pde, name = "Original", Category = "Original"),
    cbind.data.frame(do.call(rbind.data.frame, reduced_pdes), Category = "Samples")
  )

  removed_plot_data <- rbind.data.frame(
    cbind.data.frame(original_pde, name = "Original", Category = "Original"),
    cbind.data.frame(do.call(rbind.data.frame, removed_pdes), Category = "Samples")
  )

  # Create plots
  plot_theme <- theme_light() +
    theme(
      legend.position.inside = TRUE,
      legend.position = c(0.1, 0.8),
      legend.background = element_rect(fill = alpha("white", 0.5)),
      strip.background = element_rect(fill = "cornsilk"),
      strip.text = element_text(colour = "black")
    )

  p_reduced <- ggplot(reduced_plot_data, aes(x = x, y = y, color = Category, group = name)) +
    geom_line() +
    scale_color_colorblind() +
    plot_theme +
    labs(x = "Data", y = "PDE", title = "Reduced data: PDE comparison")

  p_removed <- ggplot(removed_plot_data, aes(x = x, y = y, color = Category, group = name)) +
    geom_line() +
    scale_color_colorblind() +
    plot_theme +
    labs(x = "Data", y = "PDE", title = "Removed data: PDE comparison")

  # Statistical comparison
  test_statistics <- c("ad", "kuiper", "cvm", "wass", "dts", "ks", "amrdd", "euc", "kld")

  statistical_results <- lapply(test_statistics, function(test_stat) {
    reduced_vs_orig <- apply(reduced_data_df, 2, function(x) {
      CompDistrib(vector1 = x, vector2 = data_df$Data, TestStat = test_stat)
    })

    removed_vs_orig <- apply(removed_data_df, 2, function(x) {
      CompDistrib(vector1 = x, vector2 = data_df$Data, TestStat = test_stat)
    })

    reduced_vs_removed <- mapply(CompDistrib, reduced_data_df, removed_data_df,
                                 MoreArgs = list(TestStat = test_stat))

    list(
      reduced_vs_orig = reduced_vs_orig,
      removed_vs_orig = removed_vs_orig,
      reduced_vs_removed = reduced_vs_removed
    )
  })
  names(statistical_results) <- test_statistics

  # Organize statistical results for plotting
  stats_df <- rbind.data.frame(
    cbind.data.frame(
      Comparison = "Reduced vs Original",
      Test = rownames(do.call(rbind, lapply(statistical_results, "[[", "reduced_vs_orig"))),
      do.call(rbind, lapply(statistical_results, "[[", "reduced_vs_orig"))
    ),
    cbind.data.frame(
      Comparison = "Removed vs Original",
      Test = rownames(do.call(rbind, lapply(statistical_results, "[[", "removed_vs_orig"))),
      do.call(rbind, lapply(statistical_results, "[[", "removed_vs_orig"))
    ),
    cbind.data.frame(
      Comparison = "Reduced vs Removed",
      Test = rownames(do.call(rbind, lapply(statistical_results, "[[", "reduced_vs_removed"))),
      do.call(rbind, lapply(statistical_results, "[[", "reduced_vs_removed"))
    )
  )

  stats_long <- suppressMessages(melt(stats_df))

  p_statistics <- ggplot(stats_long, aes(x = Comparison, y = value, color = Comparison)) +
    geom_beeswarm(method = "compactswarm", show.legend = FALSE) +
    facet_wrap(~Test, scales = "free", nrow = 1) +
    theme_light() +
    theme(
      legend.position = "bottom",
      legend.background = element_rect(fill = alpha("white", 0.5)),
      strip.background = element_rect(fill = "cornsilk"),
      strip.text = element_text(colour = "black"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)
    ) +
    scale_color_colorblind() +
    labs(y = "Test Statistic Value", title = "Distribution comparison statistics")

  n_facets <- length(unique(stats_long[["Test"]]))
  # Check if ylimits_list is a list of the right length (or NULL)
  use_y_limits <- is.list(ylimits_list) && length(ylimits_list) == n_facets && use_y_limits

  if (use_y_limits) {
    # Add individual scales per facet
    p_statistics <- p_statistics + facetted_pos_scales(y = ylimits_list)
  }

  # Return results
  return(list(
    plots = list(
      reduced_data_plot = p_reduced,
      removed_data_plot = p_removed,
      statistics_plot = p_statistics
    ),
    data = list(
      original_data = data_df,
      reduced_data = reduced_data_df,
      removed_data = removed_data_df,
      statistical_results = statistical_results
    )
  ))
}

# Generate synthetic 3-component GMM data (outside the function)
generate_gmm_data <- function() {
  N <- 3000
  weights <- c(.6, .1, .3)
  sds <- c(2, .001, .2)
  means <- c(0, 4, 6)

  set.seed(42)
  gmm_params <- rbind(weights * N, means, sds)
  gmm_params_list <- split(gmm_params, rep(1:ncol(gmm_params), each = nrow(gmm_params)))

  data_vector <- unlist(lapply(gmm_params_list, function(x) {
    do.call(rnorm, as.list(x))
  }))

  class_vector <- rep(1:3, weights * N)
  data.frame(Cls = class_vector, Data = data_vector)
}


# Experiment:

# y limits
my_limits <- list(
  scale_y_continuous(limits = c(0, 150000)),
  scale_y_continuous(limits = c(0, 1)),
  scale_y_continuous(limits = c(0, 12)),
  scale_y_continuous(limits = c(0, 80)),
  scale_y_continuous(limits = c(0, 0.65)),
  scale_y_continuous(limits = c(0, 18)),
  scale_y_continuous(limits = c(0, 0.15)),
  scale_y_continuous(limits = c(0, 0.25)),
  scale_y_continuous(limits = c(0, 0.55))
)

# Generate data
df_art_data_3GMM <- generate_gmm_data()

# Perform experiments with different number of trials
list_of_OptimizeBetween <- c(TRUE, FALSE)
list_of_CheckRemovedThreefold <- c(TRUE, FALSE)
list_of_sizes <- c(50, 1500, round(nrow(df_art_data_3GMM)) * 0.8)
list_of_nTrials <- c(1, 1000000)

experiment_3Gaussians_results <-
  lapply(list_of_OptimizeBetween, function(OptimizeBetween) {

    Res1 <- lapply(list_of_CheckRemovedThreefold, function(CheckRemovedThreefold) {
      Res2 <- lapply(list_of_sizes, function(Size) {
        experiment_3Gaussians_results_1 <- lapply(list_of_nTrials, function(nTrials) {
          results <- downsample_analysis(data_df = df_art_data_3GMM, nSamples = 10, Size = Size, nTrials = nTrials,
                                         CheckRemoved = CheckRemovedThreefold,
                                         OptimizeBetween = OptimizeBetween, use_y_limits = TRUE, ylimits_list = my_limits,
                                         MaxCores = parallel::detectCores() - 1)
          return(results)
        })
        names(experiment_3Gaussians_results_1) <- paste0("nTrials", list_of_nTrials)

        # Combine results plots
        p_experiment_3Gaussians_results <- cowplot::plot_grid(
          experiment_3Gaussians_results_1[[1]]$plots$reduced_data_plot,
          experiment_3Gaussians_results_1[[1]]$plots$removed_data_plot,
          experiment_3Gaussians_results_1[[1]]$plots$statistics_plot,
          experiment_3Gaussians_results_1[[2]]$plots$reduced_data_plot,
          experiment_3Gaussians_results_1[[2]]$plots$removed_data_plot,
          experiment_3Gaussians_results_1[[2]]$plots$statistics_plot,
          labels = "AUTO", nrow = 2, rel_widths = c(2, 2, 3)
        ) +
          plot_annotation(
            title = paste0("Data splits: Comparison of distributions with original data: ",
                           format(Size, scientific = FALSE), " points sampled from originally ", nrow(df_art_data_3GMM), "."),
            subtitle = ifelse(OptimizeBetween,
                              paste0("Top row: First split, Bottom row: Best of ",
                                     format(list_of_nTrials[2], scientific = FALSE), " splits; Training/testing versus validation optimized"),
                              ifelse(CheckRemovedThreefold,
                                     paste0("Top row: First split, Bottom row: Best of ",
                                            format(list_of_nTrials[2], scientific = FALSE), " splits; Threefold optimized"),
                                     paste0("Top row: First split, Bottom row: Best of ", format(list_of_nTrials[2], scientific = FALSE),
                                            " splits; Optimized for training/testing versus the original.")))
          ) &
          theme(
            plot.tag.position = c(0.5, 1),
            plot.tag = element_text(size = 14, face = "bold", vjust = 0)
          )

        print(p_experiment_3Gaussians_results)

        # Save combined plot
        ggsave(filename =
                 paste0("p_experiment_3Gaussians_results",
                        "_CheckRemovedThreefold", CheckRemovedThreefold,
                        "_Size", Size,
                        "_OptimizeBetween", OptimizeBetween,
                        ".svg"), plot = p_experiment_3Gaussians_results, width = 18, height = 12)

        return(list(
          results = experiment_3Gaussians_results_1,
          plot = p_experiment_3Gaussians_results
        ))
      })
      names(Res2) <- paste0("Size", list_of_sizes)
      return(Res2)
    })
    names(Res1) <- paste0("CheckRemovedThreefold", list_of_CheckRemovedThreefold)
    return(Res1)
  })

names(experiment_3Gaussians_results) <- paste0("OptimizeBetween", list_of_OptimizeBetween)

save.image(file = "experiment_3Gaussians.RData")
