createNAs <- function (x, pctNA = 0.1) {
  n <- nrow(x)
  p <- ncol(x)
  NAloc <- rep(FALSE, n * p)
  NAloc[sample.int(n * p, floor(n * p * pctNA))] <- TRUE
  x[matrix(NAloc, nrow = n, ncol = p)] <- NA
  return(x)
}

mttest <- createNAs(tibble::rownames_to_column(mtcars), .05)
