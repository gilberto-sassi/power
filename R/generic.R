#' Value not matching.
#'
#' \code{%notin%} is a more intuitive interface as a binary operator, which
#' returns a logical vector indicating if there is a match or not for its
#' left operand.
#'
#'
#' @inherit base::`%in%`
`%notin%` <- function(x, table){
  !(x %in% table)
}

#' Difference between two numbers.
#'
#' This is a safe way of comparing if two vectors of floating point numbers.
#'
#' @keywords comparison, equality
#' @inherit dplyr::near
#' @examples
#' not_near(sqrt(2), sqrt(3))
not_near <- function(x, y, tol = 1e-8) !dplyr::near(x, y, tol)

#' Convert numeric to character
#'
#' @description \code{cortar} divides the range of x into intervals and codes
#' the values in x according to which interval they fall. The leftmost
#' interval corresponds to level one, the next leftmost to level two and so on.
#'
#' @details \code{cortar} is similar to \code{base::cut}, but we have modify to use
#' with the package \code{bookdown}
#'
#' @usage cortar(x)
#'
#' @param x  numeric vector which is to be converted to a character vector by
#' cutting.
#' @param breaks either a numeric vector of two or more unique cut points or a
#' single number (greater than or equal to 2) giving the number of intervals
#' into which x is to be cut.
#' @param open_left logical, \code{TRUE} indicates to include in the interval the
#' left limit
#' @param open_right logica, \code{TRUE} indicates to include in the interval the
#' right limit
#' @param dig_lab nteger which is used when labels are not given. It determines
#' the number of digits used in formatting the break numbers.
#'
#' @return a character vector with the new codification
#'
#' @examples
#' x <- rnorm(1000)
#' cortar(x, breaks = -10:10)
cortar <- function(x, breaks, open_left = T, open_right = F, dig_lab = 3){

  cut(x, breaks = breaks, include.lowest = !open_left,
      right = !open_right, dig.lab = dig_lab) %>%
    map_chr(function(u){
      valor <- u
      if (str_detect(valor, '^\\[') & str_detect(valor, '\\)$')) {
        valor <- str_replace_all(valor, ',', ' |=== ')
      } else if (str_detect(valor, '^\\[') & str_detect(valor, '\\]$')){
        valor <- str_replace_all(valor, ',', ' |===| ')
      } else if (str_detect(valor, '^\\(') & str_detect(valor, '\\]$')) {
        valor <- str_replace_all(valor, ',', ' ===| ')
      }
      str_replace_all(valor, '(\\[|\\]|\\(|\\))', '')
    })
}

#' Distribution table
#'
#' \code{table_distribution} computes the distribution table for a character,
#' integer or discrete vector
#'
#' @param df a tibble
#' @param variable a scalar character
#'
#' @return a tibble that is distribution table
#'
#' @usage table_distribution(df, variable)
#'
#' @examples
#' n <- 1000
#'df <- tibble(nome = sample(c('A', 'B', 'C'), size = n, replace = T,
#'                           prob = c(0.6, 0.2, 0.2)),
#'            numero = rnorm(n))
#'table_distribution(df, 'nome')
table_distribution <- function(df, variable) {

  suppressWarnings({
    tab <- df %>%
      group_by_(variable) %>%
      summarise(frequencia = n()) %>%
      mutate(frequencia_relativa = frequencia / sum(frequencia),
             porcentagem = frequencia_relativa * 100) %>%
      arrange(desc(porcentagem))
  })

  tab <- tab %>%
    add_case(frequencia = sum(tab$frequencia),
             frequencia_relativa = sum(tab$frequencia_relativa),
             porcentagem = sum(tab$porcentagem))

  v_log <- is.na(tab[[1]])
  tab[[1]][v_log] <- 'Total'

  colnames(tab) <- c(colnames(tab)[1], 'Frequência', 'Frequencia relativa', 'Porcentagem')

  tab
}

#' Concatenating strings
#'
#'
#'\code{`%colar%`} concatenates vectors afeter converting to character.
#'
#'@param vetor1 character vector
#'@param vetor2 character vector
#'
#'@usage vetor1 %colar% vetor2
#'
#' @examples
#' "string a" %colar% "string b"
`%colar%` <- function(vetor1, vetor2) paste0(vetor1, vetor2)

#' Bowley's coefficient
#'
#' \code{bowley_coeff} computes the Bowley's Coefficient.
#'
#' @details Bowley's coefficient assess the symmetry of the data.
#'
#' @param x a numerical vector
#' @param na.rm a logical value indicating whether NA values should be stripped
#' before the computation proceeds.
#'
#'@usage bowley_coeff(x)
#'
#'@examples
#'bowley_coeff(runif(1000))
bowley_coeff <- function(x, na.rm = F){
  quartis <- numeric(3)
  quartis[1] <- quantile(x, probs = 1 / 4, na.rm = na.rm)
  quartis[2] <- quantile(x, probs = 2 / 4, na.rm = na.rm)
  quartis[3] <- quantile(x, probs = 3 / 4, na.rm = na.rm)

  (quartis[3] - 2 * quartis[2] + quartis[1]) / (quartis[3] - quartis[1])
}

