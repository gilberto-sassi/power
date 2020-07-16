#' Power and sample size for testing proportion.
#'
#' \code{pwr_prop_1pop} computes the power and the sample size for testing
#' proportion.
#'
#' @param p populational proportion
#' @param p0 proportion under null hypothesis
#' @param n number of observations (sample size)
#' @param pwr power of test \eqn{1 + \beta} (1 minus type II error probability)
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of "two.sided" (default), "greater" or "less"
#' @param sig_level significance level (Type I error probability)
#' @keywords hypothesis testing, power, significance level, proportion test,
#' one population, one variable, sample size
#' @return \code{pwr_prop_1pop} returns a list with the following components:
#' \describe{
#' \item{p}{populational proportion}
#' \item{p0}{proportion under null hypothesis}
#' \item{sig_level}{significance level}
#' \item{power_sampleSize}{A \code{tibble} with sample size \code{n} and
#' power \code{pwr}}
#' }
#'
#' @usage pwr_prop_1pop(p, p0, n = NULL, pwr = NULL,
#' alternative = "two.sided", sig_level = 0.05)
#'
#' @details Exactly one of the parameters 'n' and 'pwr' must
#' be passed as NULL, and that parameter is determined from the other.
#' Notice that the last one has non-NULL default so NULL must be explicitly
#' passed if you want to compute it.
#'
#' This function computes internally the effect size, given the population
#' proportion and proportion under null hypothesis.
#' These parameters are required.
#'
#' @examples
#' # Power
#' pwr_prop_1pop(sigma = 10, sigma0 = 20, n = 25, pwr = NULL,
#' alternative = "two.sided", sig_level = 0.05)
#' # Sample size
#' pwr_prop_1pop(sigma = 10, sigma0 = 20, n = NULL, pwr = 0.95,
#' alternative = "two.sided", sig_level = 0.05)
pwr_prop_1pop <- function(p, p0, n = NULL, pwr = NULL,
                           alternative = "two.sided", sig_level = 0.05) {
  # The user gives the effect size ou the sample size. Just one option.
  if (list(n, pwr) %>% map_lgl(~ is.null(.x)) %>% sum() %>% not_near(1)) {
    stop("Exactly one of n and pwr must be NULL")
  }
  # The user must give the effect size.
  if (missing(p0) | missing(p)) {
    stop("Proportion under null hypothesis",
         " or population proportion are messing.")
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
                pnorm((-p + p0 + qnorm(1 - sig_level / 2) *
                         sqrt(p0 * (1 - p0) / .x)) / sqrt(p * (1 - p) / .x)) +
                pnorm((-p + p0 + qnorm(sig_level / 2) *
                         sqrt(p0 * (1 - p0) / .x)) / sqrt(p * (1 - p) / .x))),
    "less" = n %>%
      map_dbl(~ pnorm((-p + p0 + qnorm(sig_level) *
                         sqrt(p0 * (1 - p0) / .x)) / sqrt(p * (1 - p) / .x))),
    "greater" = n %>%
      map_dbl(~ 1 -
                pnorm((-p + p0 + qnorm(1 - sig_level) *
                         sqrt(p0 * (1 - p0) / .x)) / sqrt(p * (1 - p) / .x)))
    )
  } else if (is.null(n)) {
    n <- switch(alternative,
    "two.sided" = {
      pwr %>%
        map_int(function(pwr) {
          faux <- function(n) n %>%
            map_dbl(~ (1 -
                         pnorm((-p + p0 + qnorm(1 - sig_level / 2) *
                                  sqrt(p0 * (1 - p0) / .x)) /
                                 sqrt(p * (1 - p) / .x)) +
                         pnorm((-p + p0 + qnorm(sig_level / 2) *
                                  sqrt(p0 * (1 - p0) / .x)) /
                                 sqrt(p * (1 - p) / .x)) -
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
            map_dbl(~ (pnorm((-p + p0 + qnorm(sig_level) *
                                sqrt(p0 * (1 - p0) / .x)) /
                               sqrt(p * (1 - p) / .x)) -
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
            map_dbl(~ (1 -
                         pnorm((-p + p0 + qnorm(1 - sig_level) *
                                  sqrt(p0 * (1 - p0) / .x)) /
                                 sqrt(p * (1 - p) / .x)) -
                         pwr)^2)
          nlminb(5, faux, lower = 5, upper = Inf)$par %>%
            ceiling() %>%
            as.integer()
        })
    }
    )
  }
  list(power_sampleSize = tibble(n = as.integer(n), pwr),
       sig_level = sig_level,
       p = p,
       p0 = p0)
}
