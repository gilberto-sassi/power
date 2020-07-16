#' Power and sample size to unbalanced ANOVA.
#'
#' \code{pwr_anova_unbalanced} computes the power and the sample size to
#' unbalanced ANOVA.
#'
#' @param group_means vector where each entry is the mean of the treatment
#' @param size_means vector where each entry is the sample size for each
#' treatment
#' @param mean population and general mean
#' @param sigma populational standard deviation
#' @param sig_level significance level (Type I error probability)
#' @keywords hypothesis testing, power, significance level,
#' unbalanced anova
#' @return \code{pwr_anova_unbalanced} returns a list with the following
#' components:
#' \describe{
#' \item{group_means}{vector where each entry is the mean of the treatment}
#' \item{size_means}{vector where each entry is the sample size for each
#' treatment}
#' \item{sigma}{population standard deviation}
#' \item{sig_level}{significance level}
#' \item{power_sampleSize}{A \code{tibble} with number of observations
#' \code{n}, \code{N} is the sample size and  power \code{pwr}}
#' }
#'
#' @usage pwr_anova_unbalanced(group_means, size_means, mean, sigma,
#' sig_level = 0.05)
#'
#' @examples
#' # Power
#' pwr_anova_unbalanced(group_means  = c(21, 14, 15, 25),
#' size_means = c(4, 5, 6, 4), mean = 21, sigma = sqrt(2), sig_level = 0.05)
pwr_anova_unbalanced <- function(group_means, size_means, mean, sigma,
                                 sig_level = 0.05) {
  # The user must give the effect size.
  if (missing(group_means) | missing(sigma) | missing(mean) |
      missing(size_means)) {
    stop("Means for each treatment,", sQuote("group_means"),
         ", population standard deviation,", sQuote("sigma"),
         ", populational and general mean, ", sQuote("mean"),
         ", and the number of observation in each treatment,",
         sQuote("size_means"),
         ", are required")
  }
  # the sample size must be greater or equal to 5
  if (!is.null(size_means)) {
    if (min(size_means) < 4) stop("Number of observations,",
                         sQuote("n"),
                         ", in each treatment must be at least 4.")
  }
  if (not_near(length(size_means), length(group_means))) {
    stop("size_means and group_means must have the same length.")
  }
  a <- length(size_means)
  ncp <- ((group_means - mean)^2 * size_means / (a * sigma^2)) %>%
    sum()
  1 - pf(qf(1 - sig_level, df1 = a - 1, df2 = sum(size_means) - a),
         df1 = a - 1, df2 = sum(size_means) - a, ncp = ncp)
}
