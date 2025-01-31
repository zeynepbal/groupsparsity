library(arrangements)

logdet <- function(Sigma) {
  return(2 * sum(log(diag(chol(Sigma)))))
}

bayesian_group_sparse_linear_classification_variational_train <- function(X, y, membership, parameters) {    
  set.seed(parameters$seed)
  
  log2pi <- log(2 * pi)
  N <- dim(X)[1]
  D <- dim(X)[2]
  G <- max(membership)

  # one-of-K-encoding
  binary_encoding <- matrix(0, length(membership), max(membership))
  binary_encoding[cbind(1:length(membership), membership)] <- 1 

  binary_indices <- matrix(0, 2^G, G)
  combination_table <- list()
  for(m in 0:(2^G - 1)) {
    m2 <- as.numeric(intToBits(m))[1:G] + 1
    binary_indices[m + 1,] <- m2
    selected <- which(m2 == 2)
    if(length(selected) > 1) {
      combination_table[[m + 1]] <- combinations(selected, 2)
    } 
  }
  
  selection_table <- array(0, c(G, G, 2^(G - 2)))
  for(g1 in 1:(G - 1)) {
    for(g2 in (g1 + 1):G) {
      selection_table[g1, g2,] <- which(binary_indices[,g1] == 2 & binary_indices[,g2] == 2)
    }
  }
  
  lambda <- list(alpha = matrix(parameters$lambda_alpha + 0.5, D, 1), beta = matrix(parameters$lambda_beta, D, 1))
  gamma <- list(alpha = parameters$gamma_alpha + 0.5, beta = parameters$gamma_beta)
  bw <- list(mu = matrix(rnorm(D + 1), D + 1, 1), sigma = diag(1, D + 1, D + 1))
  f <- list(mu = (abs(matrix(rnorm(N), N, 1)) + parameters$margin) * sign(y), sigma = matrix(1, N, 1))
  
  theta_mat <- matrix(1, G, G)
  diag(theta_mat) <- parameters$p / (1 - parameters$p)
  theta_mat <- log(theta_mat)
  
  P_initial <- rep(0, 2^G)
  for(m in 0:(2^G - 1)) {
    m2 <- binary_indices[m + 1,]
    selected <- which(m2 == 2)
    if(length(selected) > 0) { #first-order
      P_initial[m + 1] <- P_initial[m + 1] + sum(diag(theta_mat)[selected])
    }
    if(length(selected) > 1) { #second-order
      P_initial[m + 1] <- P_initial[m + 1] + sum(theta_mat[combination_table[[m + 1]]])
    }
  }
  
  P_initial <- exp(P_initial - max(P_initial))
  P_initial <- P_initial / sum(P_initial)
  
  zeta <- list(theta = theta_mat)
  
  XX <- crossprod(X, X)
  
  lower <- matrix(-1e40, N, 1)
  lower[which(y > 0)] <- +parameters$margin
  upper <- matrix(+1e40, N, 1)
  upper[which(y < 0)] <- -parameters$margin
  
  btimesbT.mu <- bw$mu[1]^2 + bw$sigma[1, 1]
  wtimeswT.mu <- tcrossprod(bw$mu[2:(D + 1)], bw$mu[2:(D + 1)]) + bw$sigma[2:(D + 1), 2:(D + 1)]
  wtimesb.mu <- bw$mu[2:(D + 1)] * bw$mu[1] + bw$sigma[2:(D + 1), 1]

  bounds <- c()
  for(iter in 1:parameters$iteration) {
    P <- rep(0, 2^G)
    for(m in 0:(2^G - 1)) {
      m2 <- binary_indices[m + 1,]
      selected <- which(m2 == 2)
      if(length(selected) > 0) { #first-order
        P[m + 1] <- P[m + 1] + sum(diag(zeta$theta)[selected])
      }
      if(length(selected) > 1) { #second-order
        P[m + 1] <- P[m + 1] + sum(zeta$theta[combination_table[[m + 1]]])
      }
    }
    
    P <- exp(P - max(P))
    P <- P / sum(P)
    exp_zeta <- matrix(0, G, G)
    for(g in 1:G) {
      exp_zeta[g, g] <- sum(P[which(binary_indices[,g] == 2)])
    }
    for(g1 in 1:(G - 1)) {
      for(g2 in (g1 + 1):G) {
        exp_zeta[g1, g2] <- sum(P[selection_table[g1, g2,]])
        exp_zeta[g2, g1] <- exp_zeta[g1, g2]
      }
    }

    # update lambda
    lambda$beta <- 1 / (1/parameters$lambda_beta + 0.5 * diag(wtimeswT.mu))

    # update gamma
    gamma$beta <- 1 / (1 / parameters$gamma_beta + 0.5 * btimesbT.mu)

    # update b and w
    bw$sigma <- chol2inv(chol(rbind(cbind(N + gamma$alpha * gamma$beta, t((binary_encoding %*% diag(exp_zeta)) * colSums(X))), cbind((binary_encoding %*% diag(exp_zeta)) * colSums(X), diag(as.vector(lambda$alpha * lambda$beta), D, D) + ((binary_encoding %*% exp_zeta %*% t(binary_encoding)) * XX)))))
    bw$mu <- bw$sigma %*% (rbind(matrix(1, 1, N), t(X) * matrix((binary_encoding %*% diag(exp_zeta)), D, N)) %*% f$mu)

    btimesbT.mu <- bw$mu[1]^2 + bw$sigma[1, 1]
    wtimeswT.mu <- tcrossprod(bw$mu[2:(D + 1)], bw$mu[2:(D + 1)]) + bw$sigma[2:(D + 1), 2:(D + 1)]
    wtimesb.mu <- bw$mu[2:(D + 1)] * bw$mu[1] + bw$sigma[2:(D + 1), 1]
    
    A_mat <- wtimeswT.mu * XX
    A_mat_tilde <- matrix(0, G, G)
    for(a1 in 1:G) {
      for(a2 in 1:G) {
        A_mat_tilde[a1, a2] <- sum((membership == a1) %*% t(membership == a2) * A_mat)
      }
    }

    diag(zeta$theta) <- diag(theta_mat) - diag(A_mat_tilde) / 2 + (t(binary_encoding) %*% (bw$mu[2:(D + 1)] * colSums(X * matrix(f$mu, N, D)))) - t(binary_encoding) %*% (wtimesb.mu * colSums(X))
    theta_mn <- (-A_mat_tilde + theta_mat)
    zeta$theta[upper.tri(zeta$theta)] <- theta_mn[upper.tri(theta_mn)]
    zeta$theta[lower.tri(zeta$theta)] <- theta_mn[lower.tri(theta_mn)]
    
    # update f
    output <- cbind(matrix(1, N, 1), X) %*% (bw$mu * c(1, (binary_encoding %*% diag(exp_zeta))))
    alpha_norm <- lower - output
    beta_norm <- upper - output
    normalization <- pnorm(beta_norm) - pnorm(alpha_norm)
    normalization[which(normalization == 0)] <- 1
    f$mu <- output + (dnorm(alpha_norm) - dnorm(beta_norm)) / normalization
    f$sigma <- 1 + (alpha_norm * dnorm(alpha_norm) - beta_norm * dnorm(beta_norm)) / normalization - (dnorm(alpha_norm) - dnorm(beta_norm))^2 / normalization^2
  
    if(parameters$progress == 1){
      # lower bound
      lb <- 0
      
      # p(f | w, b, zeta, X)
      lb <- lb - 0.5 * (crossprod(f$mu, f$mu) + sum(f$sigma)) + (t(bw$mu[2:(D + 1)]) %*% ((binary_encoding %*% diag(exp_zeta)) * colSums(X * matrix(f$mu, N, D)))) + sum(bw$mu[1] * f$mu) - 0.5 * sum(diag(((binary_encoding %*% exp_zeta %*% t(binary_encoding)) * XX) %*% wtimeswT.mu)) -(t(wtimesb.mu) %*% ((binary_encoding %*% diag(exp_zeta)) * colSums(X))) - 0.5 * N * btimesbT.mu - 0.5 * N * log2pi

      # p(lambda)
      lb <- lb + sum((parameters$lambda_alpha - 1) * (digamma(lambda$alpha) + log(lambda$beta)) - lambda$alpha * lambda$beta / parameters$lambda_beta - lgamma(parameters$lambda_alpha) - parameters$lambda_alpha * log(parameters$lambda_beta))
      
      # p(w | lambda)
      lb <- lb - 0.5 * sum(as.vector(lambda$alpha * lambda$beta) * diag(wtimeswT.mu)) - 0.5 * (D * log2pi - sum(digamma(lambda$alpha) + log(lambda$beta)))
  
      # p(gamma)
      lb <- lb + (parameters$gamma_alpha - 1) * (digamma(gamma$alpha) + log(gamma$beta)) - gamma$alpha * gamma$beta / parameters$gamma_beta - lgamma(parameters$gamma_alpha) - parameters$gamma_alpha * log(parameters$gamma_beta)
      
      # p(b | gamma)
      lb <- lb - 0.5 * gamma$alpha * gamma$beta * btimesbT.mu - 0.5 * (log2pi - (digamma(gamma$alpha) + log(gamma$beta)))
      
      # p(zeta)
      lb <- lb + sum(P * log(P_initial) , na.rm = TRUE)

      # q(lambda)
      lb <- lb + sum(lambda$alpha + log(lambda$beta) + lgamma(lambda$alpha) + (1 - lambda$alpha) * digamma(lambda$alpha))  
      
      # q(gamma)
      lb <- lb + gamma$alpha + log(gamma$beta) + lgamma(gamma$alpha) + (1 - gamma$alpha) * digamma(gamma$alpha)
      
      # q(b, w)
      lb <- lb + 0.5 * ((D + 1) * (log2pi + 1) + logdet(bw$sigma))
      
      # q(f)
      lb <- lb + 0.5 * sum(log2pi + f$sigma) + 0.5 + sum(log(normalization)) + sum((alpha_norm * dnorm(alpha_norm) - beta_norm * dnorm(beta_norm)) / 2*normalization)
      
      # q(zeta)
      lb <- lb - sum(P * log(P) , na.rm = TRUE)
      
      bounds[iter] <- lb
    }
  }  
  obj <- list(lambda = lambda, gamma = gamma, bw = bw, zeta = list(pi = diag(exp_zeta)), parameters = parameters, bounds = bounds)
  return(obj)
}
