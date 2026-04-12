kocmi_single = \(data, target, agent, conds, type = c("cont", "disc"),
                 nboots = 1e4, k = 3, threads = 1, seed = 42,
                 base = exp(1), method = "equal", contain_null = TRUE) {
  type = match.arg(type)
  mat = infoxtr:::.convert2mat(data, type)
  knockoff = .convert2mat(knockoff, contain_type = FALSE)
  null_knockoff = .convert2mat(null_knockoff, contain_type = FALSE)
  return(infoxtr:::RcppKOCMI(
            mat, target, agent, conds, knockoff, null_knockoff, type,
            nboots, k, 0, threads, seed, base, method, contain_null))
}
