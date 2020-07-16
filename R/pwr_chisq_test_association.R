#' Power and sample size for chi squared test.
#'
#' \code{pwr_chisq_test_association} computes the power and the sample size
#' for the chi squared testing to check association between two qualitative
#' variables.
#'
#' @details  In this function, we use the effect size \code{es}
#' \eqn{=\sqrt{\frac{\chi_{obs}^2}{n}}}, where \eqn{n} is the sample size and
#' \eqn{\chi^2_{obs}} is the statistics test. \code{es} is a coefficient or
#' measure of association.
#'
#'The effect size \code{es} \eqn{=w=\sqrt{\frac{w^2}{n}}} is an strength measure
#'of association. Usually, we use
#'\describe{
#'\item{\code{w=0.1}}{to detect a weak association between two variables}
#'\item{\code{w=0.3}}{to detect a moderate association between two variables}
#'\item{\code{w=0.5}}{to detect a strong association between two variables}
#'}
#'
#' @usage pwr_chisq_test_association(es, nrow, ncol, n = NULL, pwr = NULL,
#' sig_level = 0.05)
#'
#' @param es effect size \code{w =} \eqn{\sqrt{\frac{\chi^2_{obs}}{n}}}
#' @param nrow number of categories of the first variable
#' @param ncol number of categories of the second variable
#' @param n number of observations (sample size)
#' @param pwr power of test \eqn{1 + \beta} (1 minus type II error probability)
#' @param sig_level significance level (Type I error probability)
#' @keywords hypothesis testing, power, significance level, chi-squared test,
#' sample size, association between two qualitative association
#' @return \code{pwr_chisq_test_association} returns a list with the following
#' components:
#' \describe{
#' \item{es}{effect size}
#' \item{nrow}{number of categories of the first variable}
#' \item{ncol}{number of categories of the second variable}
#' \item{sig_level}{significance level}
#' \item{power_sampleSize}{A \code{tibble} with sample size and \code{n}
#' is the number of pairs}
#' }
#' @examples
#' # Power
#' pwr_chisq_test_association(es = 0.5, nrow = 5, ncol = 5,
#' n = 300, pwr = NULL, sig_level = 0.05)
#' # Sample size
#' pwr_chisq_test_association(es = 0.5, nrow = 5, ncol = 5,
#' n = NULL, pwr = 0.99, sig_level = 0.05)
pwr_chisq_test_association <- function(es, nrow, ncol, n = NULL, pwr = NULL,
                                       sig_level = 0.05) {
  # The user gives the power ou the sample size. Just one option.
  if (sum(is.null(pwr), is.null(n)) %notin% 1) {
    stop("Exactly one of n and pwr must be NULL")
    }
  # The user must give the effect size.
  if (missing(es) | missing(nrow) | missing(ncol)) {
    stop("Effect size",
         " and number of categories for both variables",
         " are required.")
  }
  # the sample size must be greater or equal to 5 and should be equal
  if (!is.null(n)) {
    if (min(n) < 5) {
      stop("Number of observations must be at least 5.")
    }
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
  if (is.null(pwr)) {
    pwr <- n %>%
      map_dbl(function(n) {
        1 - pchisq(qchisq(1 - sig_level,
                          df = (nrow - 1) * (ncol - 1)),
                   df = (nrow - 1) * (ncol - 1),
                   ncp = es^2 * n)
      })
  } else if (is.null(n)) {
    n <- pwr %>%
      map_int(function(pwr) {
        faux <- function(n) n %>%
          map_dbl(~ (1 - pchisq(qchisq(1 - sig_level,
                                       df = (nrow - 1) * (ncol - 1)),
                                df = (nrow - 1) * (ncol - 1),
                                ncp = es^2  * .x) - pwr)^2)
        nlminb(5, faux, lower = 4, upper = Inf)$par %>%
          ceiling() %>%
          as.integer()
      })
  }
  if (is.null(n)) {
    list(power_sampleSize = tibble(n = n, pwr = pwr),
         sig_level = sig_level) %>%
      return()
  } else if (!is.null(n)) {
    list(power_sampleSize = tibble(n = as.integer(n), pwr = pwr),
         sig_level = sig_level) %>%
      return()
  } else {
    stop("Unexpected error. Check the arguments n1 and n2.")
  }
}
