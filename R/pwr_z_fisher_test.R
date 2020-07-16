#' Power and sample size to test Pearson's correlation.
#'
#' \code{pwr_z_fisher_test} computes the power and the sample size for testing
#' Pearson's correlation using Fisher's Z transform.
#'
#' @details It can be proved that
#' \eqn{\frac{1}{2} \log\left(\frac{1 + r}{1 - r}\right)}
#' has asympotc normal distribution with mean
#' #' \eqn{\frac{1}{2} \log\left(\frac{1 + \rho}{1 - \rho}\right)}
#' and variance \eqn{\sqrt{\frac{1}{n - 3}}}, where \eqn{n} is the sample size,
#' \eqn{\rho} is the populational Pearson's correlation, and \eqn{r} is the
#' sample correlation.
#'
#' It's require to give the population Pearson's correlation and Pearson's
#' correlation under null hypothesis.
#'
#' @usage pwr_z_fisher_test(rho, rho0, n = NULL, pwr = NULL,
#' alternative = "two.sided", sig_level = 0.05)
#'
#' @param rho Pearson's correlation
#' @param rho0 Pearson's correlation under null hypothesis
#' @param n number of observations (sample size)
#' @param pwr power of test \eqn{1 + \beta} (1 minus type II error probability)
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of "two.sided" (default), "greater" or "less"
#' @param sig_level significance level (Type I error probability)
#' @keywords hypothesis testing, power, significance level,
#' Fisher's Z transform, correlation, sample size
#' @return \code{pwr_z_fisher_test} returns a list with the following
#' components:
#' \describe{
#' \item{rho}{Pearson's correlation}
#' \item{rho0}{Pearson's correlation under null hypothesis}
#' \item{sig_level}{significance level}
#' \item{power_sampleSize}{A \code{tibble} with sample size \code{n} and
#' power \code{pwr}}
#' }
#' @examples
#' # Power
#' pwr_z_fisher_test(rho = 0.8, rho0 = 0.7, n = 100, pwr = NULL,
#' alternative = "two.sided", sig_level = 0.05)
#' # Sample size
#' pwr_z_fisher_test(rho = 0.8, rho0 = 0.7, n = NULL, pwr = 0.99,
#' alternative = "two.sided", sig_level = 0.05)
pwr_z_fisher_test <- function(rho, rho0, n = NULL, pwr = NULL,
                        alternative = "two.sided", sig_level = 0.05) {
  # The user gives the power ou the sample size. Just one option.
  if (sum(is.null(n), is.null(pwr)) %notin% 1) {
    stop("Exactly one of n and pwr must be NULL")
  }
  # The user must give the effect size.
  if (missing(rho) | missing(rho0)) {
    stop("Pearson's correlation and",
         " Pearson's correlation under null hypothesis.")
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
                    map_dbl(function(n) {
                      m <- (0.5 * log((1 + rho) / (1 - rho)) -
                              0.5 * log((1 + rho0) / (1 - rho0))) /
                        sqrt(1 / (n - 3))
                      1 - pnorm(qnorm(1 - sig_level / 2) - m) +
                        pnorm(qnorm(sig_level / 2) - m)
                    }),
                  "less" = n %>%
                    map_dbl(function(n) {
                      m <- (0.5 * log((1 + rho) / (1 - rho)) -
                              0.5 * log((1 + rho0) / (1 - rho0))) /
                        sqrt(1 / (n - 3))
                      pnorm(qnorm(sig_level) - m)
                    }),
                  "greater" = n %>%
                    map_dbl(function(n) {
                      m <- (0.5 * log((1 + rho) / (1 - rho)) -
                              0.5 * log((1 + rho0) / (1 - rho0))) /
                        sqrt(1 / (n - 3))
                      1 - pnorm(qnorm(1 - sig_level) - m)
                    })
    )
  } else if (is.null(n)) {
    n <- switch(alternative,
                "two.sided" = {
                  pwr %>%
                    map_int(function(pwr) {
                      faux <- function(n) n %>%
                        map_dbl(function(n) {
                          m <- (0.5 * log((1 + rho) / (1 - rho)) -
                                   0.5 * log((1 + rho0) / (1 - rho0))) /
                             sqrt(1 / (n - 3))
                           (1 - pnorm(qnorm(1 - sig_level / 2) - m) +
                             pnorm(qnorm(sig_level / 2) - m) - pwr)^2
                        })
                      nlminb(5, faux, lower = 5, upper = Inf)$par %>%
                        ceiling() %>%
                        as.integer()
                    })
                },
                "less" = {
                  pwr %>%
                    map_int(function(pwr) {
                      faux <- function(n) n %>%
                        map_dbl(function(n) {
                          m <- (0.5 * log((1 + rho) / (1 - rho)) -
                                  0.5 * log((1 + rho0) / (1 - rho0))) /
                            sqrt(1 / (n - 3))
                          (pnorm(qnorm(sig_level) - m) - pwr)^2
                        })
                      nlminb(5, faux, lower = 5, upper = Inf)$par %>%
                        ceiling() %>%
                        as.integer()
                    })
                },
                "greater" = {
                  pwr %>%
                    map_int(function(pwr) {
                      faux <- function(n) n %>%
                        map_dbl(function(n) {
                          m <- (0.5 * log((1 + rho) / (1 - rho)) -
                                  0.5 * log((1 + rho0) / (1 - rho0))) /
                            sqrt(1 / (n - 3))
                          (1 - pnorm(qnorm(1 - sig_level) - m) - pwr)^2
                        })
                      nlminb(5, faux, lower = 5, upper = Inf)$par %>%
                        ceiling() %>%
                        as.integer()
                    })
                }
    )
  }
  list(power_sampleSize = tibble(n = as.integer(n), pwr),
       sig_level = sig_level,
       rho = rho,
       rho0 = rho0)
}
