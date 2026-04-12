kocmi_single = \(data, target, agent, conds, type = c("cont", "disc"),
                 monte = 40, nboots = 1e4, k = 3, threads = 1, seed = 42,
                 base = exp(1), method = "equal", contain_null = TRUE) {
  type = match.arg(type)
  mat = infoxtr:::.convert2mat(data, type)
  knockoff = construct_knockoff(mat, agent, conds, monte, seed)
  null_knockoff = NULL
  if (contain_null) null_knockoff = construct_knockoff(mat, agent, c(target,conds), monte, seed)
  return(infoxtr:::RcppKOCMI(
            mat, target, agent, conds, knockoff, null_knockoff, type,
            nboots, k, 0, threads, seed, base, method, contain_null))
}
