#' @name poisson_test_pv
#'
#' @title
#' Poisson Test
#'
#' @description
#' `poisson_test_pv()` performs an exact or approximate Poisson test about the
#' rate parameter of a Poisson distribution. In contrast to
#' [`stats::poisson.test()`], it is vectorised, only calculates p-values and
#' offers a normal approximation of their computation. Furthermore, it is
#' capable of returning the discrete p-value supports, i.e. all observable
#' p-values under a null hypothesis. Multiple tests can be evaluated
#' simultaneously. In two-sided tests, several procedures of obtaining the
#' respective p-values are implemented.
#'
#' `r lifecycle::badge('deprecated')`\cr
#' **Note**: Please use `poisson_test_pv()`! The older `poisson.test.pv()` is
#' deprecated in order to migrate to snake case. It will be removed in a future
#' version.
#'
#' @param x        integer vector giving the number of events.
#' @param lambda   non-negative numerical vector of hypothesised rate(s).
#'
#' @template param
#' @templateVar alternative TRUE
#' @templateVar ts_method TRUE
#' @templateVar exact TRUE
#' @templateVar correct TRUE
#' @templateVar simple_output TRUE
#'
#' @details
#' The parameters `x`, `lambda` and `alternative` are vectorised. They are
#' replicated automatically to have the same lengths. This allows multiple null
#' hypotheses to be tested simultaneously.
#'
#' Since the Poisson distribution itself has an infinite support, so do the
#' p-values of exact Poisson tests. Thus supports only contain p-values that
#' are not rounded off to 0.
#'
#' @template details_two_sided
#'
#' @template return
#'
#' @seealso
#' [`stats::poisson.test()`], [`binom_test_pv()`]
#'
#' @references
#' Blaker, H. (2000) Confidence curves and improved exact confidence intervals
#'   for discrete distributions. *Canadian Journal of Statistics*,
#'   **28**(4), pp. 783-798. \doi{10.2307/3315916}
#'
#' Hirji, K. F. (2006). *Exact analysis of discrete data*. New York: Chapman
#'   and Hall/CRC. pp. 55-83. \doi{10.1201/9781420036190}
#'
#' @examples
#' # Constructing
#' k <- c(4, 2, 2, 14, 6, 9, 4, 0, 1)
#' lambda <- c(3, 2, 1)
#'
#' # Computation of exact two-sided p-values ("blaker") and their supports
#' results_ex  <- poisson_test_pv(k, lambda, ts_method = "blaker")
#' raw_pvalues <- results_ex$get_pvalues()
#' pCDFlist    <- results_ex$get_pvalue_supports()
#'
#' # Computation of normal-approximated one-sided p-values ("less") and their supports
#' results_ap  <- poisson_test_pv(k, lambda, "less", exact = FALSE)
#' raw_pvalues <- results_ap$get_pvalues()
#' pCDFlist    <- results_ap$get_pvalue_supports()
#'
#' @importFrom stats pnorm qnorm
#' @importFrom checkmate assert_integerish qassert
#' @export
poisson_test_pv <- function(
  x,
  lambda = 1,
  alternative = "two.sided",
  ts_method = "minlike",
  exact = TRUE,
  correct = TRUE,
  simple_output = FALSE
) {
  # plausibility checks of input parameters
  len_x <- length(x)
  len_l <- length(lambda)
  len_a <- length(alternative)
  len_g <- max(len_x, len_l, len_a)

  qassert(x, "A+")
  assert_integerish(x, lower = 0)
  x <- round(x)
  if(len_x < len_g) x <- rep_len(x, len_g)

  qassert(x = lambda, "N+[0,)")
  if(len_l < len_g) lambda <- rep_len(lambda, len_g)

  qassert(exact, "B1")
  if(!exact) qassert(correct, "B1")

  ts_method <- match.arg(
    ts_method,
    c("minlike", "blaker", "absdist", "central")
  )

  for(i in seq_len(len_a)){
    alternative[i] <- match.arg(
      alternative[i],
      c("two.sided", "less", "greater")
    )
    if(exact && alternative[i] == "two.sided")
      alternative[i] <- ts_method
  }
  if(len_a < len_g) alternative <- rep_len(alternative, len_g)

  qassert(simple_output, "B1")

  # find unique parameter sets
  params   <- unique(data.frame(lambda, alternative))
  lambda_u <- params$lambda
  alt_u    <- params$alternative
  len_u    <- length(lambda_u)

  # prepare output
  res <- numeric(len_g)
  if(!simple_output) {
    supports <- vector("list", len_u)
    indices  <- vector("list", len_u)
  }

  # begin computations
  for(i in 1:len_u) {
    idx <- which(lambda_u[i] == lambda & alt_u[i] == alternative)
    N <- max(x[idx])

    if(exact) {
      # generate all probabilities under current lambda
      d <- generate_poisson_probs(lambda_u[i])
      # difference between number of probabilities and largest desired x-value
      len_diff <- N - length(d) + 1
      # add 0-probabilities, if necessary
      if(len_diff > 0) d <- c(d, rep(0, len_diff))
      # compute p-value support
      pv_supp <- support_exact(
        alternative = alt_u[i],
        probs = d,
        expectations = abs(seq_along(d) - 1 - lambda_u[i])
      )
    } else {
      # find quantile with (approx.) smallest probability > 0
      q <- -floor(qnorm(2^-1074, -lambda_u[i], sqrt(lambda_u[i])))
      # greatest observation
      M <- max(q, N)
      # compute p-value support
      if(lambda_u[i] == 0){
        pv_supp <- switch(
          EXPR      = alt_u[i],
          less      = rep(1, M + 1),
          greater   = c(1, rep(0, M)),
          two.sided = c(1, rep(0, M))
        )
      } else {
        # possible observations (minus expectation)
        z <- 0:M - lambda_u[i]
        # standard deviation
        std <- sqrt(lambda_u[i])
        # compute p-value support
        pv_supp <- switch(
          EXPR      = alt_u[i],
          less      = rev(c(1, pnorm(rev(z)[-1] + correct * 0.5, 0, std))),
          greater   = c(1, pnorm(z[-1] - correct * 0.5, 0, std, FALSE)),
          two.sided = pmin(1, 2 * pnorm(-abs(z) + correct * 0.5, 0, std))
        )
      }
    }

    res[idx] <- pv_supp[x[idx] + 1]
    if(!simple_output) {
      supports[[i]] <- sort(unique(pv_supp))
      indices[[i]]  <- idx
    }
  }

  out <- if(!simple_output) {
    dnames <- sapply(match.call(), deparse1)
    DiscreteTestResults$new(
      test_name = ifelse(
        exact,
        "Exact Poisson test",
        paste0(
          "Normal-approximated Poisson test",
          ifelse(correct, " with continuity correction", "")
        )
      ),
      inputs = list(
        observations = data.frame(
          `number of events` = x,
          check.names = FALSE
        ),
        nullvalues = data.frame(
          `event rate` = lambda,
          check.names = FALSE
        ),
        parameters = data.frame(alternative = alternative)
      ),
      p_values = res,
      pvalue_supports = supports,
      support_indices = indices,
      data_name = paste(dnames["x"], "and", dnames["lambda"])
    )
  } else res

  return(out)
}

#' @rdname poisson_test_pv
#' @export
#' @importFrom lifecycle deprecate_soft
poisson.test.pv <- function(
    x,
    lambda = 1,
    alternative = "two.sided",
    ts.method = "minlike",
    exact = TRUE,
    correct = TRUE,
    simple.output = FALSE
) {
  deprecate_soft("0.2", "poisson.test.pv()", "poisson_test_pv()")

  poisson_test_pv(x, lambda, alternative, ts.method, exact, correct, simple.output)
}
