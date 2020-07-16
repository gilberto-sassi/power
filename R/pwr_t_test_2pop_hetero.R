#' Power and sample size for t-test for two population.
#'
#' \code{pwr_t_test_2pop_hetero} computes the power and the sample size for
#' testing difference of means of two normal populations with different and
#' unknown variances.
#'
#' @param delta populational difference of means
#' @param delta0 difference of means under null hypothesis
#' @param sigma1 standard deviation for the first population
#' @param sigma2 standard deviation for the second population
#' @param n1 number of observations (sample size) for the first population
#' @param n2 number of observations (sample size) for the second population
#' @param pwr power of test \eqn{1 + \beta} (1 minus type II error probability)
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of "two.sided" (default), "greater" or "less"
#' @param sig_level significance level (Type I error probability)
#' @keywords hypothesis testing, power, significance level, t test,
#' two populations, two variables, difference of means, known standard and
#' different deviations, heteroscedastic t test, sample size
#' @return \code{pwr_t_test_2pop_hetero} returns a list with the following
#' components:
#' \describe{
#' \item{sigma1}{standard deviation for the first population}
#' \item{sigma2}{standard deviation for the second population}
#' \item{delta}{difference of means}
#' \item{delta0}{difference of means under null hypothesis}
#' \item{sig_level}{significance level}
#' \item{power_sampleSize}{A \code{tibble} with sample size, \code{n1} for the
#' first population and \code{n2} for the second population}
#' }
#'
#' @details Exactly one of the parameters samples sizes ('n1' and 'n2')
#' and 'pwr' must be passed as NULL, and that parameter is determined from
#' the other. Notice that the last one has non-NULL default so NULL must be
#' explicitly passed if you want to compute it.
#'
#' #' The parameters 'sigma1, 'sigma2' and 'delta' are required.
#' The effect size is computed internally.
#'
#' @usage pwr_t_test_2pop_hetero(sigma1, sigma2, delta, delta0 = 0, n1 = NULL,
#' n2 = NULL, pwr = NULL, alternative = "two.sided", sig_level = 0.05)
#'
#' @examples
#' # Power
#' pwr_t_test_2pop_hetero(sigma1 = 5, sigma2 = 10, delta = 5, delta0 = 0,
#' n1 = 10, n2 = 10, pwr = NULL, alternative = "two.sided", sig_level = 0.05)
#' # Sample size
#' pwr_t_test_2pop_hetero(sigma1 = 5, sigma2 = 10, delta = 5, delta0 = 0,
#' n1 = NULL, n2 = NULL, pwr = 0.99, alternative = "two.sided",
#' sig_level = 0.05)
pwr_t_test_2pop_hetero <- function(sigma1, sigma2, delta, delta0 = 0, n1 = NULL,
                                   n2 = NULL, pwr = NULL,
                                   alternative = "two.sided",
                                   sig_level = 0.05) {
  # The user gives the power ou the sample size. Just one option.
  if (sum(is.null(pwr), all(is.null(n1), is.null(n2))) %notin% 1) {
    stop("Exactly one of n1+n2 and pwr must be NULL")
  }
  # The user giver the sample for the two populations
  if ((is.null(n1) & !is.null(n2)) | (!is.null(n1) & is.null(n2))) {
    stop("The user should give n1 and n2.")
  }
  # The user must give the effect size.
  if (missing(delta) | missing(delta0) | missing(sigma1) | missing(sigma2)) {
    stop("Population standard deviation for the first standard deviation",
         "Population standard deviation for the second standard deviation",
         "Populational difference of means",
         " and difference of means under null hypothesis are required.")
  }
  # the sample size must be greater or equal to 5 and should be equal
  if (!is.null(n1) & !is.null(n2)) {
    if (min(n1, n2) < 5) {
      stop("Number of observations must be at least 5.")
    }
    if (not_near(length(n1), length(n2))) {
      stop("n1 and n2 must have the same length.")
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
  # Alternative should be in c("two.sided", "less", "greater")
  if (!(alternative %in% c("two.sided", "less", "greater"))) {
    stop("Alternative should be exactly one of options:",
         sQuote("two.sided"), ",",
         sQuote("less"), " and ",
         sQuote("greater"), ". Hint: check your spelling.")
  }
  if (is.null(pwr)) {
    es <- (delta - delta0) / sqrt((sigma1^2 + sigma2^2) / 2)
    pwr <- switch(alternative,
                  "two.sided" =
                    map2_dbl(n1, n2, ~
                               1 -
                               pt(qt(1 - sig_level / 2, df = .x + .y - 2),
                                  df = .x + .y - 2,
                                  ncp = abs(es) / sqrt(1 / .x + 1 / .y)) +
                               pt(qt(sig_level / 2, df = .x + .y - 2),
                                  df = .x + .y - 2,
                                  ncp = abs(es) / sqrt(1 / .x + 1 / .y))),
                  "less" =
                    map2_dbl(n1, n2, ~
                               pt(qt(sig_level, df = .x + .y - 2),
                                  df = .x + .y - 2,
                                  ncp = es / sqrt(1 / .x + 1 / .y))),
                  "greater" =
                    map2_dbl(n1, n2, ~
                               1 - pt(qt(1 - sig_level, df = .x + .y - 2),
                                      df = .x + .y - 2,
                                      ncp = es / sqrt(1 / .x + 1 / .y)))
    )
  } else if (is.null(n1) & is.null(n2)) {
    es <- (delta - delta0) / sqrt((sigma1^2 + sigma2^2) / 2)
    n <- switch(alternative,
                "two.sided" = {
                  pwr %>%
                    map_int(function(pwr) {
                      faux <- function(n) n %>%
                        map_dbl(~ (1 -
                                     pt(qt(1 - sig_level / 2, df = 2 * .x - 2),
                                        df = 2 * .x - 2,
                                        ncp = abs(es) / sqrt(2 / .x)) +
                                     pt(qt(sig_level / 2, df = 2 * .x - 2),
                                        df = 2 * .x - 2,
                                        ncp = abs(es) / sqrt(2 / .x)) -
                                     pwr)^2)
                      nlminb(5, faux, lower = 5, upper = Inf)$par %>%
                        ceiling() %>%
                        as.integer()
                    })
                },
                "less" = {
                  pwr %>%
                    map_int(function(pwr) {
                      faux <- function(n) n %>%
                        map_dbl(~ (pt(qt(sig_level, df = 2 * .x - 2),
                                      df = 2 * .x - 2,
                                      ncp = es / sqrt(2 / .x)) -
                                     pwr)^2)
                      nlminb(5, faux, lower = 5, upper = Inf)$par %>%
                        ceiling() %>%
                        as.integer()
                    })
                },
                "greater" = {
                  pwr %>%
                    map_int(function(pwr) {
                      faux <- function(n) n %>%
                        map_dbl(~ (1 - pt(qt(1 - sig_level, df = 2 * .x - 2),
                                          df = 2 * .x - 2,
                                          ncp = es / sqrt(2 / .x)) -
                                     pwr)^2)
                      nlminb(5, faux, lower = 5, upper = Inf)$par %>%
                        ceiling() %>%
                        as.integer()
                    })
                }
    )
  }
  if (is.null(n1) | is.null(n2)) {
    list(power_sampleSize = tibble(n1 = as.integer(n), n2 = as.integer(n),
                                   pwr = pwr),
         sig_level = sig_level,
         sigma1 = sigma1,
         sigma2 = sigma2,
         delta = delta,
         delta0 = delta0) %>%
      return()
  } else if (!is.null(n1) & !is.null(n2)) {
    list(power_sampleSize = tibble(n1 = n1, n2 = n2, pwr = pwr),
         sig_level = sig_level,
         sigma1 = sigma1,
         sigma2 = sigma2,
         delta = delta,
         delta0 = delta0) %>%
      return()
  } else {
    stop("Unexpected error. Check the arguments n1 and n2.")
  }
}
