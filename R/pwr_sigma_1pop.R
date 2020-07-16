#' Power and sample size for chi-square test for variance.
#'
#' \code{pwr_sigma_1pop} computes the power and the sample size for testing
#' variance in a normal variable.
#'
#' @param sigma populational standard deviation
#' @param sigma0 standard deviation under null hypothesis
#' @param n number of observations (sample size)
#' @param pwr power of test \eqn{1 + \beta} (1 minus type II error probability)
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of "two.sided" (default), "greater" or "less"
#' @param sig_level significance level (Type I error probability)
#' @keywords hypothesis testing, power, significance level, chi-squared test
#' for variance, one population, one variable, sample size
#' @return \code{pwr_sigma_1pop} returns a list with the following components:
#' \describe{
#' \item{sigma}{populational standard deviation}
#' \item{sigma0}{standard deviation under null hypothesis}
#' \item{sig_level}{significance level}
#' \item{power_sampleSize}{A \code{tibble} with sample size \code{n} and
#' power \code{pwr}}
#' }
#'
#' @usage pwr_sigma_1pop(sigma, sigma0, n = NULL, pwr = NULL,
#' alternative = "two.sided", sig_level = 0.05)
#'
#' @details Exactly one of the parameters 'n' and 'pwr' must
#' be passed as NULL, and that parameter is determined from the other.
#' Notice that the last one has non-NULL default so NULL must be explicitly
#' passed if you want to compute it.
#'
#' This function computes internally the effect size, given the populational
#' standard deviation and the standard deviation under null hypothesis.
#' These parameters are required.
#'
#' @examples
#' # Power
#' pwr_sigma_1pop(sigma = 10, sigma0 = 20, n = 25, pwr = NULL,
#' alternative = "two.sided", sig_level = 0.05)
#' # Sample size
#' pwr_sigma_1pop(sigma = 10, sigma0 = 20, n = NULL, pwr = 0.95,
#' alternative = "two.sided", sig_level = 0.05)
pwr_sigma_1pop <- function(sigma, sigma0, n = NULL, pwr = NULL,
                            alternative = "two.sided", sig_level = 0.05) {
  # The user gives the power ou the sample size. Just one option.
  if (list(n, pwr) %>% map_lgl(~ is.null(.x)) %>% sum() %>% not_near(1)) {
    stop("Exactly one of n and pwr must be NULL")
  }
  # The user must give the effect size.
  if (missing(sigma0) | missing(sigma)) {
    stop("General standard deviation under null hypothesis",
         " or standard deviation are messing.")
  }
  # the sample size must be greater or equal to 5
  if (!is.null(n)) {
    if (min(n) < 5) stop("Number of observations,",
                         sQuote("n"),
                         ", in each group must be at least 5.")
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
    pwr <- switch(alternative,
                  "two.sided" = n %>%
                    map_dbl(~  1 -
                              pchisq(qchisq(1 - sig_level / 2,
                                            df = .x - 1) * sigma0^2 / sigma^2,
                                     df = .x - 1) +
                              pchisq(qchisq(sig_level / 2,
                                            df = .x - 1) * sigma0^2 / sigma^2,
                                     df = .x - 1)),
                  "less" = n %>%
                    map_dbl(~ pchisq(qchisq(sig_level,
                                            df = .x - 1) * sigma0^2 / sigma^2,
                                     df = .x - 1)),
                  "greater" = n %>%
                    map_dbl(~ 1 -
                              pchisq(qchisq(1 - sig_level,
                                            df = .x - 1) * sigma0^2 / sigma^2,
                                     df = .x - 1))
    )
  } else if (is.null(n)) {
    n <- switch(alternative,
                "two.sided" = {
                  pwr %>%
                    map_int(function(pwr) {
                      faux <- function(n) n %>%
                        map_dbl(~ (1 -
                                     pchisq(qchisq(1 - sig_level / 2,
                                                   df = .x - 1) *
                                              sigma0^2 / sigma^2,
                                            df = .x - 1) +
                                     pchisq(qchisq(sig_level / 2,
                                                   df = .x - 1) *
                                              sigma0^2 / sigma^2,
                                            df = .x - 1) - pwr)^2)
                      nlminb(5, faux, lower = 5, upper = Inf)$par %>%
                        ceiling() %>%
                        as.integer()
                    })
                },
                "less" = {
                  pwr %>%
                    map_int(function(pwr) {
                      faux <- function(n) n %>%
                        map_dbl(~ (pchisq(qchisq(sig_level,
                                                 df = .x - 1) *
                                            sigma0^2 / sigma^2,
                                          df = .x - 1) - pwr)^2)
                      nlminb(5, faux, lower = 5, upper = Inf)$par %>%
                        ceiling() %>%
                        as.integer()
                    })
                },
                "greater" = {
                  pwr %>%
                    map_int(function(pwr) {
                      faux <- function(n) n %>%
                        map_dbl(~ (1 -
                                     pchisq(qchisq(1 - sig_level,
                                                   df = .x - 1) *
                                              sigma0^2 / sigma^2,
                                            df = .x - 1) - pwr)^2)
                      nlminb(5, faux, lower = 5, upper = Inf)$par %>%
                        ceiling() %>%
                        as.integer()
                    })
                }
    )
  }
  list(power_sampleSize = tibble(n = as.integer(n), pwr),
       sig_level = sig_level,
       sigma0 = sigma0,
       sigma = sigma)
}
