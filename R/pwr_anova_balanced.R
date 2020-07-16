#' Power and sample size to balanced ANOVA.
#'
#' \code{pwr_anova_balanced} computes the power and the sample size to ANOVA.
#'
#' @param group_means vector where each entry is the mean of the treatment
#' @param mean population and general mean
#' @param sigma populational standard deviation
#' @param n number of observations in each treatment.
#' @param pwr power of test \eqn{1 + \beta} (1 minus type II error probability)
#' @param sig_level significance level (Type I error probability)
#' @keywords hypothesis testing, power, significance level,
#' balanced anova, sample size
#' @return \code{pwr_anova_balanced} returns a list with the following
#' components:
#' \describe{
#' \item{group_means}{means for each treatment}
#' \item{sigma}{population standard deviation}
#' \item{sig_level}{significance level}
#' \item{power_sampleSize}{A \code{tibble} with number of observations
#' \code{n}, \code{N} is the sample size and  power \code{pwr}}
#' }
#'
#' @usage pwr_anova_balanced(group_means, mean, sigma, pwr = NULL, n = NULL,
#' sig_level = 0.05)
#'
#' @examples
#' # Power
#' pwr_anova_balanced(group_means  = c(15, 20, 25, 30, 35), mean = 15,
#' sigma = 5,  pwr = NULL, n = 10, sig_level = 0.05)
#' # Sample size
#' pwr_anova_balanced(group_means  = c(15, 20, 25, 30, 35), mean = 15,
#' sigma = 5,  pwr = 0.95, n = NULL, sig_level = 0.05)
pwr_anova_balanced <- function(group_means, mean, sigma, pwr = NULL, n = NULL,
                              sig_level = 0.05) {
  # The user gives the power ou the sample size. Just one option.
  if (not_near(sum(is.null(n), is.null(pwr)), 1)) {
    stop("Exactly one of n and pwr must be NULL")
  }
  # The user must give the effect size.
  if (missing(group_means) | missing(sigma) | missing(mean)) {
    stop("Means for each treatment,", sQuote("group_means"),
         ", population standard deviation,", sQuote("sigma"), ", and ",
         "populational and general mean, ", sQuote("mean"), ", are required")
  }
  # the sample size must be greater or equal to 5
  if (!is.null(n)) {
    if (min(n) < 5) stop("Number of observations,",
                         sQuote("n"),
                         ", in each treatment must be at least 5.")
  }
  # the pwr must belong to (0, 1)
  if (!is.null(pwr) & (!all(is.numeric(pwr)) | any(0 > pwr, pwr > 1))) {
    stop("Power, ",
         sQuote("pwr"),
         ", must be a real number belonging to (0,1).")
  }
  # Significance level must belong to (0, 1)
  if (is.null(sig_level) | !is.numeric(sig_level) |
      any(0 > sig_level, sig_level > 1)) {
    stop("Significance level, ",
         sQuote("sig_level"), ", must be a real number belonging to (0,1).")
  }
  es <- ((group_means - mean)^2 / sigma^2) %>%
    mean() %>%
    sqrt()
  a <- length(group_means)
  if (is.null(pwr)) {
    pwr <- n %>%
      map_dbl(function(n) {
        1 - pf(qf(1 - sig_level, df1 = a - 1, df2 = n * a - a),
             df1 = a - 1, df2 = n * a - a, ncp = n * es^2)
    })
  } else if (is.null(n)) {
    n <- pwr %>%
      map_int(function(pwr) {
        faux <- function(n) {
          n %>%
            map_dbl(function(n) {
              1 - pf(qf(1 - sig_level, df1 = a - 1, df2 = n * a - a),
                     df1 = a - 1, df2 = n * a - a, ncp = n * es^2)
            })
        }
        nlminb(5, faux, lower = 5, upper = Inf)$par %>%
          ceiling() %>%
          as.integer()
      })
  }
  list(power_sampleSize = tibble(n = as.integer(n), pwr, N = a * n),
       sig_level = sig_level,
       group_means = group_means,
       sigma = sigma)
}
