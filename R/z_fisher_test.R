#' Fisher's Z test for Pearson's correlation.
#'
#' \code{z_fisher_test} performs test on Pearson's correlation using Fisher's Z
#' transform.
#'
#' @details It can be proved that
#' \eqn{\frac{1}{2}\ln\left(\frac{1 + r}{1- r} \right)}
#' has normal distribution with mean
#' \eqn{\frac{1}{2}\ln\left(\frac{1 + \rho}{1- \rho} \right)}
#' and variance \eqn{\sigma^2 = \frac{1}{n - 3}}.
#'
#'@usage z_fisher_test(x, y, rho0 = 0, alternative = "two.sided",
#'conf_level = 0.95)
#'
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param rho0 Pearson's correlation under null hypothesis.
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of "two.sided" (default), "greater" or "less".
#' @param conf_level confidence level of the interval
#'
#' @return \code{z_fisher_test} return a list with the following components:
#' \describe{
#' \item{statistic}{statistic of test}
#' \item{p_value}{p-value}
#' \item{alternative}{a character string specifying the alternative hypothesis}
#' \item{lower_ci}{lower bound in the confidence interval}
#' \item{upper_ci}{upper bound in the confidence interval}
#' \item{conf_level}{confidence level of the interval}
#' }
#'
#' @examples
#' # x and y are not correlated
#' Sigma <- matrix(c(1, 0, 0, 1), ncol = 2)
#' m <- MASS::mvrnorm(500, mu = c(0, 1), Sigma = Sigma)
#' z_fisher_test(m[, 1], m[, 2])
#' # x and y are correlated
#' Sigma <- matrix(c(1, 0.85, 0.85, 1), ncol = 2)
#' m <- MASS::mvrnorm(500, mu = c(0, 1), Sigma = Sigma)
#' z_fisher_test(m[, 1], m[, 2])
z_fisher_test <- function(x, y, rho0 = 0, alternative = "two.sided",
                          conf_level = 0.95) {
  # x and y are required
  if (missing(x) | missing(y)) {
    stop("x and y are required.")
  }
  # x and y must be numeric vector
  if (!rlang::is_bare_numeric(x) | !rlang::is_bare_numeric(y)) {
    stop("x and y must be numeric vector.")
  }
  # x and y must have the same length
  if (not_near(length(x), length(y))) {
    stop("x and y must have the same length.")
  }
  if (length(x) < 5) {
    stop("x and y must have at least 5 entries.")
  }
  # Alternative should be in c("two.sided", "less", "greater")
  if (alternative %notin% c("two.sided", "less", "greater")) {
    stop("Alternative should be exactly one of options:",
         sQuote("two.sided"), ",",
         sQuote("less"), " and ",
         sQuote("greater"), ". Hint: check your spelling.")
  }
  r <- cor(x, y)
  z0  <- (0.5 * log((1 + r) / (1 - r)) - 0.5 * log((1 + rho0) / (1 - rho0))) /
    sqrt(1 / (length(x) - 3))
  p_value <- base::switch(alternative,
                          "two.sided" = 2 * (1 - pnorm(abs(z0))),
                          "less" = pnorm(z0),
                          "greater" = 1 - pnorm(z0))
  xi <- 0.5 * log((1 + r) / (1 - r))
  z_q <- qnorm(1 - conf_level / 2)
  sd <- sqrt(1 / (length(x) - 3))
  lower_ci <- (exp(2 * (-z_q * sd + xi)) - 1) /
    (exp(2 * (-z_q * sd + xi)) + 1)
  upper_ci <- (exp(2 * (z_q * sd + xi)) - 1) /
    (exp(2 * (z_q * sd + xi)) + 1)
  tibble(statistic = z0, p_value = p_value, alternative = alternative,
         lower_ci = lower_ci, upper_ci = upper_ci, conf_level = conf_level,
         cor = r)
}
