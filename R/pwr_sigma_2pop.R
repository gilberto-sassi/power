#' Power and sample size for F test for variance.
#'
#' \code{pwr_sigma_2pop} computes the power and the sample size for testing
#' to compare variances of two normal population.
#'
#' @param sigma1 populational standard deviation for the first population
#' @param sigma2 populational standard deviation for the second population
#' @param n1 number of observations (sample size) for the first population
#' @param n2 number of observations (sample size) for the second population
#' @param pwr power of test \eqn{1 + \beta} (1 minus type II error probability)
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of "two.sided" (default), "greater" or "less"
#' @param sig_level significance level (Type I error probability)
#' @keywords hypothesis testing, power, significance level, F test
#' for variance, two population, two variables, sample size
#' @return \code{pwr_sigma_2pop} returns a list with the following components:
#' \describe{
#' \item{sigma1}{populational standard deviation}
#' \item{sigma2}{standard deviation under null hypothesis}
#' \item{sig_level}{significance level}
#' \item{power_sampleSize}{A \code{tibble} with sample size, \code{n1} for the
#' first population and \code{n2} for the second population}
#' }
#'
#' @usage pwr_sigma_2pop(sigma1, sigma2, n1 = NULL, n2 = NULL, pwr = NULL,
#' alternative = "two.sided", sig_level = 0.05)
#'
#' @details Exactly one of the parameters samples sizes ('n1' and 'n2')
#' and 'pwr' must be passed as NULL, and that parameter is determined from
#' the other. Notice that the last one has non-NULL default so NULL must be
#' explicitly passed if you want to compute it.
#'
#' The parameters 'sigma1 and 'sigma2' are required.
#' The effect size is computed internally.
#'
#' @examples
#' # Power
#' pwr_sigma_2pop(sigma1 = 10, sigma2 = 20, n1 = 25, n2 = 21, pwr = NULL,
#' alternative = "two.sided", sig_level = 0.05)
#' # Sample size
#' pwr_sigma_2pop(sigma = 10, sigma0 = 20, n = NULL, pwr = 0.95,
#' alternative = "two.sided", sig_level = 0.05)
pwr_sigma_2pop <- function(sigma1, sigma2, n1 = NULL, n2 = NULL, pwr = NULL,
                           alternative = "two.sided", sig_level = 0.05) {
  # The user gives the power ou the sample size. Just one option.
  if (sum(is.null(pwr), all(is.null(n1), is.null(n2))) %notin% 1) {
    stop("Exactly one of n1+n2 and pwr must be NULL")
  }
  # The user must give the effect size.
  if (missing(sigma1) | missing(sigma2)) {
    stop("Standard deviation for the populations are required.")
  }
  # The user giver the sample for the two populations
  if ((is.null(n1) & !is.null(n2)) | (!is.null(n1) & is.null(n2))) {
    stop("The user should give n1 and n2.")
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
    pwr <- switch(alternative,
                  "two.sided" = map2_dbl(n1, n2, ~
                 1 -
                   pf(qf(1  - sig_level / 2, df1 = .x - 1, df2 = .y - 1) *
                        sigma2^2 /  sigma1^2,
                      df1 = .x - 1, df2 = .y - 1) +
                   pf(qf(sig_level / 2, df1 = .x - 1, df2 = .y - 1) *
                        sigma2^2 /  sigma1^2,
                      df1 = .x - 1, df2 = .y - 1)),
                  "less" = map2_dbl(n1, n2, ~
                pf(qf(sig_level, df1 = .x - 1, df2 = .y - 1) *
                     sigma2^2 /  sigma1^2,
                   df1 = .x - 1, df2 = .y - 1)),
                  "greater" = map2_dbl(n1, n2, ~
                 1 -
                   pf(qf(1  - sig_level, df1 = .x - 1, df2 = .y - 1) *
                        sigma2^2 /  sigma1^2,
                      df1 = .x - 1, df2 = .y - 1))
    )
  } else if (is.null(n1) & is.null(n2)) {
    n <- switch(alternative,
                "two.sided" = {
                  pwr %>%
                    map_int(function(pwr) {
                      faux <- function(n) n %>%
                        map_dbl(~ (1 -
                   pf(qf(1  - sig_level / 2, df1 = .x - 1, df2 = .x - 1) *
                        sigma2^2 /  sigma1^2,
                      df1 = .x - 1, df2 = .x - 1) +
                   pf(qf(sig_level / 2, df1 = .x - 1, df2 = .x - 1) *
                        sigma2^2 /  sigma1^2,
                      df1 = .x - 1, df2 = .x - 1) -
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
                      map_dbl(~ (pf(qf(sig_level, df1 = .x - 1, df2 = .x - 1) *
                                      sigma2^2 /  sigma1^2,
                                    df1 = .x - 1, df2 = .x - 1) -
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
                 pf(qf(1  - sig_level, df1 = .x - 1, df2 = .x - 1) *
                      sigma2^2 /  sigma1^2,
                    df1 = .x - 1, df2 = .x - 1) -
                                     pwr)^2)
                      nlminb(5, faux, lower = 5, upper = Inf)$par %>%
                        ceiling() %>%
                        as.integer()
                    })
                }
    )
  }
  if (is.null(n1) & is.null(n2)) {
    list(power_sampleSize = tibble(n1 = as.integer(n), n2 = as.integer(n),
                                   pwr = pwr),
         sig_level = sig_level,
         sigma1 = sigma1,
         sigma2 = sigma2) %>%
      return()
  } else if (!is.null(n1) & !is.null(n2)) {
    list(power_sampleSize = tibble(n1 = n1, n2 = n2, pwr = pwr),
         sig_level = sig_level,
         sigma1 = sigma1,
         sigma2 = sigma2) %>%
      return()
  } else {
    stop("Unexpected error. Check the arguments n1 and n2.")
  }
}
