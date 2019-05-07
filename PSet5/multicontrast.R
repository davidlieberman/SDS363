multicontrast <- function(contrast, data, grouping){
  
  # Groups
  groups <- as.vector(as.matrix(unique(grouping)))
  
  # Multivariate Means for Each Variable in Each Group
  M <- apply(data, 2, function(y) tapply(y, grouping, mean))
  M <- M[match(groups, row.names(M)),]
  
  # Counts for Each Group
  N <- table(grouping)
  N <- N[match(groups, row.names(N))]
  
  # Calculate Weighted Sum of Squared Weights
  SW <- sum(contrast^2 / N)
  
  # Calculate SSCP Between (Hypothesis Matrix)
  C_hat <- colSums(contrast*M)
  SSCP_between <- (cbind(C_hat) %*% rbind(C_hat)) / SW
  
  # Calculate SSCP Within (Error Matrix)
  SSCP_within_each_group <- list()
  for (i in seq_along(groups)){
    X <- subset(data, grouping == groups[i])
    deviations <- matrix(0, nrow(X), ncol(X))
    for (j in 1:ncol(X)){
      deviations[,j] <- X[,j] - M[i,j]
    }
    SSCP_within_each_group[[i]] <- t(deviations) %*% deviations
  }
  SSCP_within <- Reduce("+", SSCP_within_each_group)
  
  # Calculate Wilks' Lambda
  lambda <- det(SSCP_within) / det(SSCP_between + SSCP_within)
  
  # Calculate Degrees of Freedom
  p <- ncol(data) # no.variables
  m <- nrow(data) - length(groups)  # no.observations - no.groups
  
  df1 <- p
  df2 <- m-p+1
  
  # Calculate approx. F
  F.stat <- (1 - lambda) / lambda * df2 / df1
  
  # Calculate p-value from F distribution
  p.value <- 1 - pf(F.stat, df1, df2)
  
  # Output
  out <- c(lambda, F.stat, df1, df2, p.value)
  names(out) <- c("Wilks", "approx.F", "df1", "df2", "p.value")
  return(out) }
