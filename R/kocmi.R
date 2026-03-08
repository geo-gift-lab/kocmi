KOCMI.net = \(expr, k=3, M, tf=colnames(expr), pcc = 0){
  expr = as.data.frame(expr)
  expr = expr[,tf]
  weightNet = transMatrix(expr)
  if (pcc == 0) {
    df.ko0 = x_knockoff(expr,M)
    a = apply(weightNet,1,function(x){KOCMI(expr,x[1],x[2],k,M,df.ko0)})
  } else {
    sigma = stats::cor(expr)
    a = apply(weightNet,1,function(x){KOCMI(expr,x[1],x[2],k,M,sigma=sigma,pcc=pcc)})
  }
  p.value = lapply(a,\(x) x$pvalue)
  t1 = lapply(a,\(x) x$t1)
  weightNet$pvalue = as.numeric(p.value)
  weightNet$p_adj = as.numeric(stats::p.adjust(p.value,method = 'BH'))
  weightNet$t1 = as.numeric(t1)
  weightNet$cs = abs(as.numeric(t1))
  return(weightNet)
}

transMatrix = \(expr){
  tf = colnames(expr)
  n = length(tf)
  re = matrix(0, n*(n-1), 4, dimnames = list(c(), c("regulator", "target", "cs", "pvalue")))
  re = as.data.frame(re)
  re$regulator = rep(tf, each = (n-1))
  target = c()
  for(i in seq_len(n)){
    target = c(target, tf[-i])
  }
  re$target = target
  return(re)
}

KOCMI = \(data, cause, effect, k=3, M = 50, df.ko0 = NULL, sigma = NULL, pcc = NULL, seed = 42){
  if(!is.null(pcc) & !is.null(sigma)){
    sigma = sigma[, colnames(sigma) != effect]
    pccnode = names(which(abs(sigma[cause,]) >= pcc))
    pccnode = c(pccnode, effect)
    if (length(pccnode) <= 2) {
      pccnode = names(sort(abs(sigma[cause,]), decreasing = T)[1:2])
      pccnode = c(pccnode, effect)
    }
    condition = setdiff(pccnode, c(cause, effect))

    data = data[, pccnode]
  } else {
    condition = setdiff(colnames(data), c(cause, effect))
  }

  df.cause = data[, cause]
  df.effect = data[, effect]
  df.condition = data[, condition]

  if(is.null(df.ko0)){df.ko0 = x_knockoff(data, M, seed = seed)}

  cmi0 = KNN.CMI(data, cause, effect, k)

  cmi0.knockoff = apply(df.ko0, 3, \(x){KNN.CMI(cbind(cause = x[,cause],effect = df.effect, df.condition), "cause", "effect", k)})

  df.knockoff = x_knockoff(cbind(df.cause, df.condition), M)

  cmi.knockoff = apply(df.knockoff, 3, \(x){KNN.CMI(cbind(cause = x[,1], effect = df.effect, df.condition), "cause", "effect", k)})

  cmi.knockoff = round(cmi.knockoff, digits = 8)
  cmi0.knockoff = round(cmi0.knockoff, digits = 8)

  D = cmi0.knockoff-cmi.knockoff
  p = permutation_test_mean(D)$p_value
  t1 = mean(D) / stats::sd(D)
  #t3 = mean(cmi0-cmi.knockoff) / stats::sd(cmi0-cmi.knockoff)
  #t2 = (mean(cmi0.knockoff) - mean(cmi.knockoff)) / sqrt(var(cmi0.knockoff)/M + var(cmi0.knockoff)/M)

  return(list(pvalue = p,t1 = t1,cmi0 = cmi0,
              cmi0.knockoff = cmi0.knockoff,
              cmi.knockoff = cmi.knockoff))
}

KNN.CMI = \(data, cause = "x", effect = "y", k = 3){
  condition = setdiff(colnames(data), c(cause, effect))
  df_xz = data[,c(cause, condition)]
  df_yz = data[,c(effect, condition)]
  df_z = as.matrix(data[,condition])

  D_all = as.matrix(stats::dist(data, method = 'maximum'))
  D_z = as.matrix(stats::dist(df_z, method = 'maximum'))
  D_xz = as.matrix(stats::dist(df_xz, method = 'maximum'))
  D_yz = as.matrix(stats::dist(df_yz, method = 'maximum'))

  N = nrow(data)
  n_xz = c(); n_yz = c(); n_z = c()

  for(i in seq_len(N)){
    epsilon = sort(D_all[i,])[k+1]
    n_xz = c(n_xz, length(which(D_xz[i,] < epsilon)))
    n_yz = c(n_yz, length(which(D_yz[i,] < epsilon)))
    n_z = c(n_z, length(which(D_z[i,] < epsilon)))
  }

  c = digamma(k) - mean(digamma(n_xz)) - mean(digamma(n_yz)) + mean(digamma(n_z))
  return(c)
}

permutation_test_mean = \(x, n_perm = 10000) {
  observed_stat = abs(mean(x))
  permuted_stats = numeric(n_perm)

  for (i in seq_len(n_perm)) {
    signs = sample(c(-1, 1), size = length(x), replace = TRUE)
    permuted_x = x * signs
    permuted_stats[i] = abs(mean(permuted_x))
  }

  p_value = mean(permuted_stats >= observed_stat)

  return(list(
    observed_statistic = observed_stat,
    p_value = p_value
  ))
}
