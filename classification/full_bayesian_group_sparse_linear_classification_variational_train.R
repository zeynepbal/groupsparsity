library(arrangements)

logdet <- function(Sigma) {
  return(2 * sum(log(diag(chol(Sigma)))))
}

full_bayesian_group_sparse_linear_classification_variational_train <- function(X, y, membership, parameters) {    
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
  eta <- list(alpha = matrix(parameters$eta_alpha, G, 1), beta = matrix(parameters$eta_beta, G, 1))
  Theta <- list(alpha = matrix(parameters$theta_alpha, G, G), beta = matrix(parameters$theta_beta, G, G))
  bw <- list(mu = matrix(rnorm(D + 1), D + 1, 1), sigma = diag(1, D + 1, D + 1))
  f <- list(mu = (abs(matrix(rnorm(N), N, 1)) + parameters$margin) * sign(y), sigma = matrix(1, N, 1))
  
  theta_mat <- Theta$alpha * Theta$beta
  diag(theta_mat) <- eta$alpha * eta$beta
  theta_mat <- log(theta_mat)
  
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

    # update eta
    eta$alpha <- parameters$eta_alpha + diag(exp_zeta)
    
    # update Theta
    Theta$alpha <- parameters$theta_alpha + exp_zeta
    
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

    diag(zeta$theta) <- (digamma(eta$alpha) + log(eta$beta)) - diag(A_mat_tilde) / 2 + (t(binary_encoding) %*% (bw$mu[2:(D + 1)] * colSums(X * matrix(f$mu, N, D)))) - t(binary_encoding) %*% (wtimesb.mu * colSums(X))
    theta_mn <- (-A_mat_tilde + (digamma(Theta$alpha) + log(Theta$beta)))
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
      
      # p(zeta | eta, Theta)
      theta_zeta <- exp_zeta * (digamma(Theta$alpha) + log(Theta$beta))
      lb <- lb +(P[1] + sum(diag(exp_zeta) * (digamma(eta$alpha) + log(eta$beta))) + sum(theta_zeta[upper.tri(theta_zeta)]))
      
      # p(eta)
      lb <- lb + sum((parameters$eta_alpha - 1) * (digamma(eta$alpha) + log(eta$beta)) - eta$alpha * eta$beta / parameters$eta_beta - lgamma(parameters$eta_alpha) - parameters$eta_alpha * log(parameters$eta_beta))
      
      # p(Theta)
      p_theta <- (parameters$theta_alpha - 1) * (digamma(Theta$alpha) + log(Theta$beta)) - Theta$alpha * Theta$beta / parameters$theta_beta - lgamma(parameters$theta_alpha) - parameters$theta_alpha * log(parameters$theta_beta)
      lb <- lb + sum(p_theta[upper.tri(p_theta)])
      
      # q(lambda)
      lb <- lb + sum(lambda$alpha + log(lambda$beta) + lgamma(lambda$alpha) + (1 - lambda$alpha) * digamma(lambda$alpha))  
      
      # q(gamma)
      lb <- lb + gamma$alpha + log(gamma$beta) + lgamma(gamma$alpha) + (1 - gamma$alpha) * digamma(gamma$alpha)
      
      # q(b, w)
      lb <- lb + 0.5 * ((D + 1) * (log2pi + 1) + logdet(bw$sigma))
      
      # q(f)
      lb <- lb + 0.5 * sum(log2pi + f$sigma) + 0.5 + sum(log(normalization)) + sum((alpha_norm * dnorm(alpha_norm) - beta_norm * dnorm(beta_norm)) / 2*normalization)
      
      # q(zeta)
      lb <- lb - sum(P * log(P), na.rm = TRUE)
      
      # q(eta)
      lb <- lb + sum(eta$alpha + log(eta$beta) + lgamma(eta$alpha) + (1 - eta$alpha) * digamma(eta$alpha))  
      
      # q(Theta)
      q_theta <- Theta$alpha + log(Theta$beta) + lgamma(Theta$alpha) + (1 - Theta$alpha) * digamma(Theta$alpha)
      lb <- lb + sum(q_theta[upper.tri(q_theta)])  
      
      bounds[iter] <- lb
    }
  }  
  obj <- list(lambda = lambda, gamma = gamma, bw = bw, zeta = list(pi = diag(exp_zeta)), exp_zeta = exp_zeta, eta = eta, Theta = Theta, parameters = parameters, bounds = bounds)
  return(obj)
}
