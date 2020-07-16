#' Power and sample size for paired t-test for two population.
#'
#' \code{pwr_paired_t_test} computes the power and the sample size for testing
#' difference of means of two associates and normal populations .
#'
#' @param delta difference of means
#' @param delta0 difference of means under null hypothesis
#' @param sigma_d standard deviation of difference of each pair
#' @param n number of pairs
#' @param pwr power of test \eqn{1 + \beta} (1 minus type II error probability)
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of "two.sided" (default), "greater" or "less"
#' @param sig_level significance level (Type I error probability)
#' @keywords hypothesis testing, power, significance level, paired t test,
#' two dependent populations, difference of means, sample size
#' @return \code{pwr_paired_t_test} returns a list with the following
#' components:
#' \describe{
#' \item{sigma_d}{standard deviation of difference of each pair}
#' \item{delta}{difference of means}
#' \item{delta0}{dfference of means under null hypothesis}
#' \item{sig_level}{significance level}
#' \item{power_sampleSize}{A \code{tibble} with sample size and \code{n}
#' is the number of pairs}
#' }
#' @examples
#' # Power
#' pwr_paired_t_test(sigma_d = 5, delta = 5, delta0 = 0,
#' n1 = 10, n2 = 10, pwr = NULL, alternative = "two.sided", sig_level = 0.05)
#' # Sample size
#' pwr_paired_t_test(sigma_d = 5, delta = 5, delta0 = 0,
#' n1 = NULL, n2 = NULL, pwr = 0.99, alternative = "two.sided",
#' sig_level = 0.05)
pwr_paired_t_test <- function(sigma_d, delta, delta0 = 0, n = NULL, pwr = NULL,
                          alternative = "two.sided", sig_level = 0.05) {
  # The user gives the power ou the sample size. Just one option.
  if (sum(is.null(pwr), is.null(n)) %notin% 1) {
    stop("Exactly one of n and pwr must be NULL")
  }
  # The user must give the effect size.
  if (missing(delta) | missing(delta0) | missing(sigma_d)) {
    stop("Population standard deviation of difference of each pair",
         " Populational difference of means",
         " and difference of means under null hypothesis are required.")
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
  # Alternative should be in c("two.sided", "less", "greater")
  if (!(alternative %in% c("two.sided", "less", "greater"))) {
    stop("Alternative should be exactly one of options:",
         sQuote("two.sided"), ",",
         sQuote("less"), " and ",
         sQuote("greater"), ". Hint: check your spelling.")
  }
  if (is.null(pwr)) {
    es <- (delta - delta0) / sigma_d
    pwr <- switch(alternative,
                  "two.sided" =
                    map_dbl(n, ~
                              1 -
                              pt(qt(1 - sig_level / 2, df = .x - 1),
                                 df = .x - 1,
                                 ncp = abs(es) * sqrt(.x)) +
                              pt(qt(sig_level / 2, df = .x - 1),
                                 df = .x - 1,
                                 ncp = abs(es) * sqrt(.x))),
                  "less" =
                    map_dbl(n, ~
                              pt(qt(sig_level, df = .x - 1),
                                 df = .x - 1,
                                 ncp = es * sqrt(.x))),
                  "greater" =
                    map_dbl(n, ~
                              1 - pt(qt(1 - sig_level, df = .x - 1),
                                     df = .x - 1,
                                     ncp = es * sqrt(.x)))
    )
  } else if (is.null(n)) {
    es <- (delta - delta0) / sigma_d
    n <- switch(alternative,
                "two.sided" = {
                  pwr %>%
                    map_int(function(pwr) {
                      faux <- function(n) n %>%
                        map_dbl(~ (1 -
                                     pt(qt(1 - sig_level / 2, df = .x - 1),
                                        df = .x - 1,
                                        ncp = abs(es) * sqrt(.x)) +
                                     pt(qt(sig_level / 2, df = .x - 1),
                                        df = .x - 1,
                                        ncp = abs(es) * sqrt(.x)) -
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
                        map_dbl(~ (pt(qt(sig_level, df = .x - 1),
                                      df = .x - 1,
                                      ncp = es * sqrt(.x)) -
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
                        map_dbl(~ (1 - pt(qt(1 - sig_level, df = .x - 1),
                                          df = .x - 1,
                                          ncp = es * sqrt(.x)) -
                                     pwr)^2)
                      nlminb(5, faux, lower = 5, upper = Inf)$par %>%
                        ceiling() %>%
                        as.integer()
                    })
                }
    )
  }
  if (is.null(n)) {
    list(power_sampleSize = tibble(n = n, pwr = pwr),
         sig_level = sig_level,
         sigma_d = sigma_d,
         delta = delta,
         delta0 = delta0) %>%
      return()
  } else if (!is.null(n)) {
    list(power_sampleSize = tibble(n = as.integer(n), pwr = pwr),
         sig_level = sig_level,
         sigma_d = sigma_d,
         delta = delta,
         delta0 = delta0) %>%
      return()
  } else {
    stop("Unexpected error. Check the arguments n1 and n2.")
  }
}
