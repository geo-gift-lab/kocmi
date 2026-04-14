#' Generate Ghost Knockoff Samples
#'
#' @inheritParams construct_ghostknockoff
#'
#' @details
#' The function first estimates the empirical correlation matrix of the input
#' variables and then applies the Ghost Knockoff construction using the
#' `"asdp"` method. A set of Gaussian perturbations is generated and combined
#' with the original data through a linear transformation to produce knockoff
#' variables.
#'
#' The output is a three-dimensional array where each slice along the third
#' dimension corresponds to one Monte Carlo knockoff realization.
#'
#' @returns A numeric array of dimension
#'   `(n_samples, n_variables, monte)` containing knockoff samples.
#' @export
#'
#' @examples
#' set.seed(42)
#' z = stats::rnorm(100)
#' x = sin(z) + 0.01 * stats::rnorm(100)
#' y = sin(x) + exp(z) + 0.01 * stats::rnorm(100)
#' data = cbind(x, y, z)
#' kocmi::x_ghostknockoff(data)
#'
x_ghostknockoff = \(x, monte = 50, seed = 42) {
  set.seed(seed)

  cor.G = stats::cor(x)
  n.sample = nrow(x)
  n.G = ncol(x)
  x.knockoff = array(NA, dim = c(n.sample, n.G, monte),
                     dimnames = list(rownames(x), colnames(x), 1:monte))

  knockoff = GhostKnockoff.prelim(cor.G, M = monte, method = "asdp")
  P.each = knockoff$P.each
  V.left = as.matrix(knockoff$V.left)
  permute.index = knockoff$permute.index
  permute.V.index = rep(permute.index, monte)

  Normal_50Studies = as.matrix(V.left %*% matrix(stats::rnorm(ncol(V.left) * n.sample), ncol(V.left), n.sample))  # 50 -> sample size
  Normal_50Studies[permute.V.index == 1, ] = matrix(
    stats::rnorm(sum(permute.index) * monte * n.sample), sum(permute.index) * monte, n.sample)

  for (i in seq_len(n.sample)) {
    Normal_k = matrix(Normal_50Studies[, i], nrow = n.G)
    x_ik = as.vector(P.each %*% t(x[i, , drop = F])) + Normal_k
    for (j in seq_len(monte)) {
      x.knockoff[i,,j] = x_ik[,j]
    }
  }

  return(x.knockoff)
}

#' Construct Target-wise Ghost Knockoff Samples
#'
#' @inheritParams infoxtr::surd
#' @param x A numeric matrix of observations. Rows correspond to samples
#'   and columns correspond to variables.
#' @param monte (optional) Number of Monte Carlo knockoff samples to generate.
#'   Default is 50.
#' @param seed (optional) Random seed for reproducibility of knockoff generation.
#'
#' @details
#' The function applies ghost knockoff construction to the joint set of
#' target and agent variables. From each Monte Carlo realization, only
#' the knockoff copy of the target variable is retained.
#'
#' @returns A numeric matrix.
#' @export
#'
#' @examples
#' set.seed(42)
#' z = stats::rnorm(100)
#' x = sin(z) + 0.01 * stats::rnorm(100)
#' y = sin(x) + exp(z) + 0.01 * stats::rnorm(100)
#' data = cbind(x, y, z)
#' kocmi::construct_ghostknockoff(data, 1, 2:3)
#'
construct_ghostknockoff = \(x, target, agent, monte = 50, seed = 42) {
  agent = sort(unique(agent))
  if(target %in% agent) stop("`target` and `agent` must be disjoint.")
  k.array = x_ghostknockoff(x[,c(target,agent),drop = FALSE], monte, seed)
  return(k.array[,1,])
}
