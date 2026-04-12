#' Title
#'
#' @inheritParams infoxtr::kocmi
#' @param regulator Integer vector of column indices for the regulator variable.
#' @param monte
#' @param contain_null (optional) Logical. If `TRUE`, the test statistic is
#'   computed using knockoffs generated under the null model. In this case
#'   the difference is defined as \eqn{I(Y; X_{null} | Z) - I(Y; X_{knockoff} | Z)}.
#'   If `FALSE`, the original conditional mutual information \eqn{I(Y; X | Z)} is
#'   used instead and compared against the knockoff estimates \eqn{I(Y; X_{knockoff} | Z)}.
#'
#' @returns A named numeric vector.
#' @export
#'
#' @examples
kocmi_single = \(data, target, regulator, conds, type = c("cont", "disc"),
                 monte = 40, nboots = 1e4, k = 3, threads = 1, seed = 42,
                 base = exp(1), method = "equal", contain_null = TRUE) {
  type = match.arg(type)
  mat = infoxtr:::.convert2mat(data, type)
  knockoff = construct_knockoff(mat, regulator, conds, monte, seed)
  null_knockoff = NULL
  if (contain_null) {
    null_knockoff = construct_knockoff(mat, regulator, c(target,conds), monte, seed)
  }
  return(infoxtr:::RcppKOCMI(
            mat, target, regulator, conds, knockoff, null_knockoff, type,
            nboots, k, 0, threads, seed, base, method, contain_null))
}
