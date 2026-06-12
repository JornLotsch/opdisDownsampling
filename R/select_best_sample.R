# ===== BEST TRIAL SELECTION FUNCTIONS =====

#' @title Select best trial for one matrix
#' @description Selects the best trial when optimizing for one matrix
#'   such as for the reduced data subset versus the original
#' @param AD_statMat Matrix of statistics comparing data
#' @return Integer. Index of the best trial
#' @keywords internal
select_best_trial_one_matrix <- function(AD_statMat) {
  if (is.null(AD_statMat) || nrow(AD_statMat) == 0) {
    stop("opdisDownsampling: No data available.")
  }

  max_vals <- apply(AD_statMat, 1, max, na.rm = TRUE)

  if (all(is.infinite(max_vals))) {
    warning("opdisDownsampling: All values are infinite in reduced vs removed comparison.",
      call. = FALSE
    )
    return(1)
  }

  which.min(max_vals)
}
