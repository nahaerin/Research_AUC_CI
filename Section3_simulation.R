###############################
############ packages #########
# install.packages("pROC")
# install.packages("VGAM")
# install.packages("clinfun")
# install.packages("mvtnorm")
# install.packages("expm")
###############################


################################
########## functions ###########
################################


############ Hanley & McNeil

Q1 <- function(A) {  A/(2-A)  }
Q2 <- function(A) {  2*(A^2) / (1+A)  }

H_Mc <- function(data, y){
  nd <- sum(y); n <- length(y); nn <- n-nd
  glm_f <- glm(y ~ data, family = binomial(link = "logit"))
  z <- summary(glm_f)$coefficients[3,4]
  # Linear score of full model
  full_model <- glm.fit(cbind(1,data), y, family = binomial())$linear.predictors
  
  # Linear score of reduced model (p=1 -> coef = 1)
  red_model <- glm.fit(cbind(1,data[,1]), y, family = binomial())$linear.predictors
  
  # AUC using Mann-whitney
  A1 <- pROC::auc(y, as.vector(full_model), levels = c(0,1), direction = "<")
  A2 <- pROC::auc(y, as.vector(red_model), levels = c(0,1), direction = "<")  
  delta <- A1 - A2
  
  # Calculate sd(delta)
  r1 <- VGAM::kendall.tau(full_model[which(y==1)],red_model[which(y==1)])
  r2 <- VGAM::kendall.tau(full_model[which(y==0)],red_model[which(y==0)])
  r <- (r1+r2)/2
  Var_A1 <- ( A1*(1-A1) + (nd - 1)*(Q1(A1)-A1^2) + (nn-1)*(Q2(A1) - A1^2) ) / (nd*nn)
  Var_A2 <- ( A2*(1-A2) + (nd - 1)*(Q1(A2)-A2^2) + (nn-1)*(Q2(A2) - A2^2) ) / (nd*nn)
  sd_delta <- sqrt(Var_A1 + Var_A2 - 2*r*sqrt(Var_A1)*sqrt(Var_A2))
  
  # 95% Confidence interval
  CI_L <- delta - 1.96 * sd_delta
  CI_U <- delta + 1.96 * sd_delta
  
  return(list(A_f = as.vector(A1), A_r = as.vector(A2), delta = delta, CI_L=CI_L, CI_U=CI_U, z=z))
}


############ DeLong
DeLong <- function(data, y) {
  nd <- sum(y); n <- length(y); nn <- n-nd
  glm_f <- glm(y ~ data, family = binomial(link = "logit"))
  z <- summary(glm_f)$coefficients[3,4]
  
  # Linear score of full model
  full_model <- glm.fit(cbind(1,data), y, family = binomial())$linear.predictors
  
  # Linear score of reduced model (p=1 -> coef = 1)
  red_model <- glm.fit(cbind(1,data[,1]), y, family = binomial())$linear.predictors
  
  # Delong's variance and mann-whitney auc
  DELONG <- clinfun::roc.area.test(cbind(full_model,red_model),y)
  L <- t(c(1,-1))
  delta <- DELONG$area[1]-DELONG$area[2]
  
  # 95% Confidence interval
  CI_L <- L %*% DELONG$area - 1.96 * sqrt(L %*% DELONG$var %*% t(L))
  CI_U <- L %*% DELONG$area + 1.96 * sqrt(L %*% DELONG$var %*% t(L))
  
  return(list(A_f = as.vector(DELONG$area[1]), A_r = as.vector(DELONG$area[2]), delta = delta, CI_L=CI_L, CI_U=CI_U, z=z))
}



############ Bandos
Bandos <- function(data, y) {
  nd <- sum(y); n <- length(y); nn <- n-nd
  glm_f <- glm(y ~ data, family = binomial(link = "logit"))
  z <- summary(glm_f)$coefficients[3,4]
  
  # Linear score of full model
  full_model <- glm.fit(cbind(1,data), y, family = binomial())$linear.predictors
  
  # Linear score of reduced model (p=1 -> coef = 1)
  red_model <- glm.fit(cbind(1,data[,1]), y, family = binomial())$linear.predictors
  
  # Mann-whitney AUC
  A1 <- pROC::auc(y, as.vector(full_model), levels = c(0,1), direction = "<")
  A2 <- pROC::auc(y, as.vector(red_model), levels = c(0,1), direction = "<")  
  delta <- A1 - A2
  
  # Caculate sd(delta)
  psi_full <- psi_red <- matrix(NA, nrow = nd, ncol = nn)
  full_model_1 <- full_model[which(y==1)]; full_model_0 <- full_model[which(y==0)]
  red_model_1 <- red_model[which(y==1)]; red_model_0 <- red_model[which(y==0)]
  for(v in 1:nd) {
    for(w in 1:nn){
      psi_full[v,w] <-  as.numeric(full_model_1[v]>full_model_0[w]) + 0.5*as.numeric(full_model_1[v]==full_model_0[w])
      psi_red[v,w] <-  as.numeric(red_model_1[v]>red_model_0[w]) + 0.5*as.numeric(red_model_1[v]==red_model_0[w])
    }
  }
  W <- psi_full - psi_red
  W.i <- apply(W, 1, mean)
  W.j <- apply(W, 2, mean)
  W.. <- delta 
  last_term <- 0
  for(s in 1:nd) {
    for(t in 1:nn) {
      last_term <- last_term + (W[s,t]-W.i[s]-W.j[t]+W..)^2
    }
  }
  sd_delta <- sqrt(sum((W.i-W..)^2)/nd^2 + sum((W.j-W..)^2)/nn^2 + last_term/(nd^2 * nn^2)) 
  
  # 95% Confidence interval  
  CI_L <- delta - 1.96 * sd_delta
  CI_U <- delta + 1.96 * sd_delta
  
  return(list(A_f = as.vector(A1), A_r = as.vector(A2), delta = delta, CI_L=CI_L, CI_U=CI_U, z=z))
}


############ Li
Li <- function(data, y) {
  nd <- sum(y); n <- length(y); nn <- n-nd
  glm_f <- glm(y ~ data, family = binomial(link = "logit"))
  z <- summary(glm_f)$coefficients[3,4]
  
  # Linear score of full model
  full_model <- glm.fit(cbind(1,data), y, family = binomial())$linear.predictors
  
  # Linear score of reduced model (p=1 -> coef = 1)
  red_model <- glm.fit(cbind(1,data[,1]), y, family = binomial())$linear.predictors
  
  # Mann-whitney AUC
  A1 <- pROC::auc(y, as.vector(full_model), levels = c(0,1), direction = "<")
  A2 <- pROC::auc(y, as.vector(red_model), levels = c(0,1), direction = "<")  
  delta <- A1 - A2
  
  full_model_1 <- full_model[which(y==1)]; red_model_1 <- red_model[which(y==1)]
  full_model_0 <- full_model[which(y==0)]; red_model_0 <- red_model[which(y==0)]
  
  # Monte Carlo Algorithm
  x_11_bar <- mean(full_model_1); x_21_bar <- mean(red_model_1);
  x_10_bar <- mean(full_model_0); x_20_bar <- mean(red_model_0)
  d <- c( x_11_bar-x_10_bar, x_21_bar-x_20_bar )
  ssx_1.1 <- sum( (full_model_1-x_11_bar)^2 )
  ssx_2.1 <- sum( (red_model_1-x_21_bar)^2 )
  ssx_12.1 <- sum( (full_model_1-x_11_bar)*(red_model_1-x_21_bar) )
  ssx_1.0 <- sum( (full_model_0-x_10_bar)^2 )
  ssx_2.0 <- sum( (red_model_0-x_20_bar)^2 )
  ssx_12.0 <- sum( (full_model_0-x_10_bar)*(red_model_0-x_20_bar) )
  b <- ssx_12.1 / ssx_1.1
  g <- ssx_12.0 / ssx_1.0
  ssx_given.1 <- ssx_2.1 - ssx_12.1^2/ssx_1.1
  ssx_given.0 <- ssx_2.0 - ssx_12.0^2/ssx_1.0
  R_delta_MC <- c()
  
  for(K in 1:3000) {
    # Step 1
    U1 <- rchisq(1, nd-1) ; U_given <- rchisq(1, nd-2)
    V1 <- rchisq(1, nn-1) ; V_given <- rchisq(1, nn-2)
    ZB <- rnorm(1); ZG <- rnorm(1)
    
    # Step 2
    R_sg_1 <- ssx_1.1 / U1
    R_sg_given <- ssx_given.1 / U_given
    R_beta <- b - ZB*sqrt( ssx_given.1/(U_given*ssx_1.1) )
    R_lmb_1 <- ssx_1.0 / V1
    R_lmb_given <- ssx_given.0 / V_given
    R_gamma <- g - ZG*sqrt( ssx_given.0/(V_given*ssx_1.0) )
    
    # Step 3
    R_sg_2 <- R_beta^2*R_sg_1 + R_sg_given
    R_sg_12 <- R_beta*R_sg_1
    R_lmb_2 <- R_gamma^2*R_lmb_1 + R_lmb_given
    R_lmb_12 <- R_gamma*R_lmb_1
    
    # Step 4
    ZD <- mvtnorm::rmvnorm(1,rep(0,2),diag(1,2))
    R_SIG <- matrix(c(R_sg_1, R_sg_12, R_sg_12, R_sg_2), nrow = 2, byrow = T)
    R_LMB <- matrix(c(R_lmb_1, R_lmb_12, R_lmb_12, R_lmb_2), nrow = 2, byrow = T)
    R_tau <- d - as.vector(expm::sqrtm(R_SIG/nd + R_LMB/nn) %*% t(ZD))
    
    R_A1 <- pnorm(R_tau[1] / sqrt(R_SIG[1,1]+R_LMB[1,1]))
    R_A2 <- pnorm(R_tau[2] / sqrt(R_SIG[2,2]+R_LMB[2,2]))
    
    R_delta <- R_A1 - R_A2
    
    R_delta_MC[K] <- R_delta  
  }
  
  # Step 5  95% CI
  CI_Li <- quantile(R_delta_MC,c(0.025,0.975))
  CI_L <- CI_Li[1]
  CI_U <- CI_Li[2]
  
  return(list(A_f = as.vector(A1), A_r = as.vector(A2), delta = delta, CI_L=as.numeric(CI_L), CI_U=as.numeric(CI_U), z=z))
}


############ Heller
Heller <- function(data,y) {
  nd <- sum(y); n <- length(y); nn <- n-nd
  glm_f <- glm(y ~ data, family = binomial(link = "logit"))
  z <- summary(glm_f)$coefficients[3,4]
  
  full_mat <- data
  red_mat <- data[,-ncol(data)]
  full_mat <- as.matrix(full_mat); red_mat <- as.matrix(red_mat)
  
  MRC <- clinfun::deltaAUC(y, red_mat, full_mat[,ncol(data)])
  
  # Linear score of full model
  Linear_full <- full_mat %*% MRC$par.full
  
  # Linear score of reduced model (p=1 -> coef = 1)
  Linear_red <- red_mat
  
  # Smoothing AUC
  A1 <- MRC$results[2,1] 
  A2 <- MRC$results[2,2] 
  delta <- A1 - A2
  
  # bandwidth for smoothed AUC, variance etc.
  hf <- sd(Linear_full) / n^(1/3)
  
  # Caculate sd(delta)
  var1 <- var2 <- 0
  
  for(i in 1:n) {
    for(j in 1:n) {
      for(k in 1:n) {
        if(y[i]==1 & y[j]==0 & y[k]==0 & k != j) {
          var1 <- var1 + (pnorm((Linear_full[i]-Linear_full[j])/hf) - pnorm((Linear_red[i]-Linear_red[j])/hf) - delta) * (pnorm((Linear_full[i]-Linear_full[k])/hf) - pnorm((Linear_red[i]-Linear_red[k])/hf) - delta)
        }
        if(y[i]==1 & y[j]==0 & y[k]==1 & k != j) {
          var2 <- var2 + (pnorm((Linear_full[i]-Linear_full[j])/hf) - pnorm((Linear_red[i]-Linear_red[j])/hf) - delta) * (pnorm((Linear_full[k]-Linear_full[j])/hf) - pnorm((Linear_red[k]-Linear_red[j])/hf) - delta)
        }
      }
    }
  }
  
  var1 <- var1 / (nn*nd*(nn-1))
  var2 <- var2 / (nn*nd*(nd-1))
  var_delta <- var1/nd + var2/nn   # variance of delta
  var_tau <- var_delta/(4*delta)   # variance of tau
  
  # 95% CI using delta
  CI_DIFF <- c(delta - 1.96*sqrt(var_delta), delta + 1.96*sqrt(var_delta))
  
  # 95% CI using tau
  if(delta < 0) {
    tau <- 0 
    CI_DIFFvst <- c(0,0)
  } else {
    tau <- sqrt(delta)
    lower <- tau - 1.96*sqrt(var_tau)
    ifelse(lower < 0, lower <- 0, lower <- lower)
    upper <- tau + 1.96*sqrt(var_tau)
    CI_DIFFvst <- c((lower)^2, (upper)^2)
  }
  return(list(A_f = as.vector(A1), A_r = as.vector(A2), delta = delta, 
              CI_L=CI_DIFF[1], CI_U=CI_DIFF[2], 
              CIvst_L = CI_DIFFvst[1], CIvst_U = CI_DIFFvst[2], z=z)) 
}


################################
########## Simulation ##########
################################

###### mud = 0.4 ###### 2020.09.16

mud <- 0.4
rho <- seq(-0.8,0.8,0.1)
iter <- 1000
HM <- DL <- BA <- LI <- array(NA, dim = c(iter,3))
HE <- array(NA, dim = c(iter,5))
RESULT <- array(NA, dim = c(length(rho), 15))
colnames(RESULT) <- c("delta1", "HM_CI_L", "HM_CI_U", "DL_CI_L", "DL_CI_U", "BA_CI_L", "BA_CI_U",
                      "delta2", "LI_CI_L", "LI_CI_U",
                      "delta3", "HE_CI_L", "HE_CI_U", "HE_CIvst_L", "HE_CIvst_U")

y <- rep(c(1,0), c(30,30))
nd <- nn <- 30; n <- 60
k <- 1
for(j in 1:length(rho)) {
  i <- 1
  while(i <= iter) {
    set.seed(k)
    
    nondefault <- mvtnorm::rmvnorm(nn, c(0,0), matrix(c(1,rho[j],rho[j],1), nrow = 2))
    default <- mvtnorm::rmvnorm(nd, c(0.5, mud), matrix(c(1,rho[j],rho[j],1), nrow = 2))
    data <- rbind(default, nondefault)
    
    glm_f <- glm(y ~ data, family = binomial(link = "logit")) 
    z <- summary(glm_f)$coefficients[3,4]  ## p-value
    
    if(z < 0.05) {
      A <- H_Mc(data, y)
      B <- DeLong(data, y)
      C <- Bandos(data, y)
      D <- Li(data, y)
      E <- Heller(data,y)
      
      HM[i,] <- c(A$delta, A$CI_L, A$CI_U)
      DL[i,] <- c(B$delta, B$CI_L, B$CI_U)
      BA[i,] <- c(C$delta, C$CI_L, C$CI_U)
      LI[i,] <- c(D$delta, D$CI_L, D$CI_U)
      HE[i,] <- c(E$delta, E$CI_L, E$CI_U, E$CIvst_L, E$CIvst_U)
      
      cat("rho = ", rho[j], ", iteration = ", i, "\n")
      
      i <- i+1
    } else {  i <- i  }
    
    k <- k + 1    
    
  }
  RESULT[j,] <- c(apply(HM, 2, mean), apply(DL, 2, mean)[-1], apply(BA, 2, mean)[-1], 
                  apply(LI, 2, mean), apply(HE, 2, mean))
}

RESULT <- cbind.data.frame(cbind(mud, rho, RESULT))

write.csv(RESULT, "C:/Users/nahaerin/Desktop/git/Research_AUC_CI/mud_0.4.csv")


###### mud = 0.45 ###### 2020.09.21

mud <- 0.45
HM <- DL <- BA <- LI <- array(NA, dim = c(iter,3))
HE <- array(NA, dim = c(iter,5))
RESULT <- array(NA, dim = c(length(rho), 15))
colnames(RESULT) <- c("delta1", "HM_CI_L", "HM_CI_U", "DL_CI_L", "DL_CI_U", "BA_CI_L", "BA_CI_U",
                      "delta2", "LI_CI_L", "LI_CI_U",
                      "delta3", "HE_CI_L", "HE_CI_U", "HE_CIvst_L", "HE_CIvst_U")

k <- 1
for(j in 1:length(rho)) {
  i <- 1
  while(i <= iter) {
    set.seed(k)
    
    nondefault <- mvtnorm::rmvnorm(nn, c(0,0), matrix(c(1,rho[j],rho[j],1), nrow = 2))
    default <- mvtnorm::rmvnorm(nd, c(0.5, mud), matrix(c(1,rho[j],rho[j],1), nrow = 2))
    data <- rbind(default, nondefault)
    
    glm_f <- glm(y ~ data, family = binomial(link = "logit")) 
    z <- summary(glm_f)$coefficients[3,4]  ## p-value
    
    if(z < 0.05) {
      A <- H_Mc(data, y)
      B <- DeLong(data, y)
      C <- Bandos(data, y)
      D <- Li(data, y)
      E <- Heller(data,y)
      
      HM[i,] <- c(A$delta, A$CI_L, A$CI_U)
      DL[i,] <- c(B$delta, B$CI_L, B$CI_U)
      BA[i,] <- c(C$delta, C$CI_L, C$CI_U)
      LI[i,] <- c(D$delta, D$CI_L, D$CI_U)
      HE[i,] <- c(E$delta, E$CI_L, E$CI_U, E$CIvst_L, E$CIvst_U)
      
      cat("rho = ", rho[j], ", iteration = ", i, "\n")
      
      i <- i+1
    } else {  i <- i  }
    
    k <- k + 1    
    
  }
  RESULT[j,] <- c(apply(HM, 2, mean), apply(DL, 2, mean)[-1], apply(BA, 2, mean)[-1], 
                  apply(LI, 2, mean), apply(HE, 2, mean))
}

RESULT <- cbind.data.frame(cbind(mud, rho, RESULT))

write.csv(RESULT, "C:/Users/nahaerin/Desktop/git/Research_AUC_CI/mud_0.45.csv")


###### mud = 0.5 ###### 2020.09.28

mud <- 0.5
HM <- DL <- BA <- LI <- array(NA, dim = c(iter,3))
HE <- array(NA, dim = c(iter,5))
RESULT <- array(NA, dim = c(length(rho), 15))
colnames(RESULT) <- c("delta1", "HM_CI_L", "HM_CI_U", "DL_CI_L", "DL_CI_U", "BA_CI_L", "BA_CI_U",
                      "delta2", "LI_CI_L", "LI_CI_U",
                      "delta3", "HE_CI_L", "HE_CI_U", "HE_CIvst_L", "HE_CIvst_U")

k <- 1
for(j in 1:length(rho)) {
  i <- 1
  while(i <= iter) {
    set.seed(k)
    
    nondefault <- mvtnorm::rmvnorm(nn, c(0,0), matrix(c(1,rho[j],rho[j],1), nrow = 2))
    default <- mvtnorm::rmvnorm(nd, c(0.5, mud), matrix(c(1,rho[j],rho[j],1), nrow = 2))
    data <- rbind(default, nondefault)
    
    glm_f <- glm(y ~ data, family = binomial(link = "logit")) 
    z <- summary(glm_f)$coefficients[3,4]  ## p-value
    
    if(z < 0.05) {
      A <- H_Mc(data, y)
      B <- DeLong(data, y)
      C <- Bandos(data, y)
      D <- Li(data, y)
      E <- Heller(data,y)
      
      HM[i,] <- c(A$delta, A$CI_L, A$CI_U)
      DL[i,] <- c(B$delta, B$CI_L, B$CI_U)
      BA[i,] <- c(C$delta, C$CI_L, C$CI_U)
      LI[i,] <- c(D$delta, D$CI_L, D$CI_U)
      HE[i,] <- c(E$delta, E$CI_L, E$CI_U, E$CIvst_L, E$CIvst_U)
      
      cat("rho = ", rho[j], ", iteration = ", i, "\n")
      
      i <- i+1
    } else {  i <- i  }
    
    k <- k + 1    
    
  }
  RESULT[j,] <- c(apply(HM, 2, mean), apply(DL, 2, mean)[-1], apply(BA, 2, mean)[-1], 
                  apply(LI, 2, mean), apply(HE, 2, mean))
}

RESULT <- cbind.data.frame(cbind(mud, rho, RESULT))

write.csv(RESULT, "C:/Users/nahaerin/Desktop/git/Research_AUC_CI/mud_0.5.csv")


###### mud = 0.55 ###### 2020.10.01

mud <- 0.55
HM <- DL <- BA <- LI <- array(NA, dim = c(iter,3))
HE <- array(NA, dim = c(iter,5))
RESULT <- array(NA, dim = c(length(rho), 15))
colnames(RESULT) <- c("delta1", "HM_CI_L", "HM_CI_U", "DL_CI_L", "DL_CI_U", "BA_CI_L", "BA_CI_U",
                      "delta2", "LI_CI_L", "LI_CI_U",
                      "delta3", "HE_CI_L", "HE_CI_U", "HE_CIvst_L", "HE_CIvst_U")

k <- 1
for(j in 1:length(rho)) {
  i <- 1
  while(i <= iter) {
    set.seed(k)
    
    nondefault <- mvtnorm::rmvnorm(nn, c(0,0), matrix(c(1,rho[j],rho[j],1), nrow = 2))
    default <- mvtnorm::rmvnorm(nd, c(0.5, mud), matrix(c(1,rho[j],rho[j],1), nrow = 2))
    data <- rbind(default, nondefault)
    
    glm_f <- glm(y ~ data, family = binomial(link = "logit")) 
    z <- summary(glm_f)$coefficients[3,4]  ## p-value
    
    if(z < 0.05) {
      A <- H_Mc(data, y)
      B <- DeLong(data, y)
      C <- Bandos(data, y)
      D <- Li(data, y)
      E <- Heller(data,y)
      
      HM[i,] <- c(A$delta, A$CI_L, A$CI_U)
      DL[i,] <- c(B$delta, B$CI_L, B$CI_U)
      BA[i,] <- c(C$delta, C$CI_L, C$CI_U)
      LI[i,] <- c(D$delta, D$CI_L, D$CI_U)
      HE[i,] <- c(E$delta, E$CI_L, E$CI_U, E$CIvst_L, E$CIvst_U)
      
      cat("rho = ", rho[j], ", iteration = ", i, "\n")
      
      i <- i+1
    } else {  i <- i  }
    
    k <- k + 1    
    
  }
  RESULT[j,] <- c(apply(HM, 2, mean), apply(DL, 2, mean)[-1], apply(BA, 2, mean)[-1], 
                  apply(LI, 2, mean), apply(HE, 2, mean))
}

RESULT <- cbind.data.frame(cbind(mud, rho, RESULT))

write.csv(RESULT, "C:/Users/nahaerin/Desktop/git/Research_AUC_CI/mud_0.55.csv")

###### mud = 0.6 ###### 2020.10.04

mud <- 0.6
HM <- DL <- BA <- LI <- array(NA, dim = c(iter,3))
HE <- array(NA, dim = c(iter,5))
RESULT <- array(NA, dim = c(length(rho), 15))
colnames(RESULT) <- c("delta1", "HM_CI_L", "HM_CI_U", "DL_CI_L", "DL_CI_U", "BA_CI_L", "BA_CI_U",
                      "delta2", "LI_CI_L", "LI_CI_U",
                      "delta3", "HE_CI_L", "HE_CI_U", "HE_CIvst_L", "HE_CIvst_U")

k <- 1
for(j in 1:length(rho)) {
  i <- 1
  while(i <= iter) {
    set.seed(k)
    
    nondefault <- mvtnorm::rmvnorm(nn, c(0,0), matrix(c(1,rho[j],rho[j],1), nrow = 2))
    default <- mvtnorm::rmvnorm(nd, c(0.5, mud), matrix(c(1,rho[j],rho[j],1), nrow = 2))
    data <- rbind(default, nondefault)
    
    glm_f <- glm(y ~ data, family = binomial(link = "logit")) 
    z <- summary(glm_f)$coefficients[3,4]  ## p-value
    
    if(z < 0.05) {
      A <- H_Mc(data, y)
      B <- DeLong(data, y)
      C <- Bandos(data, y)
      D <- Li(data, y)
      E <- Heller(data,y)
      
      HM[i,] <- c(A$delta, A$CI_L, A$CI_U)
      DL[i,] <- c(B$delta, B$CI_L, B$CI_U)
      BA[i,] <- c(C$delta, C$CI_L, C$CI_U)
      LI[i,] <- c(D$delta, D$CI_L, D$CI_U)
      HE[i,] <- c(E$delta, E$CI_L, E$CI_U, E$CIvst_L, E$CIvst_U)
      
      cat("rho = ", rho[j], ", iteration = ", i, "\n")
      
      i <- i+1
    } else {  i <- i  }
    
    k <- k + 1    
    
  }
  RESULT[j,] <- c(apply(HM, 2, mean), apply(DL, 2, mean)[-1], apply(BA, 2, mean)[-1], 
                  apply(LI, 2, mean), apply(HE, 2, mean))
}

RESULT <- cbind.data.frame(cbind(mud, rho, RESULT))

write.csv(RESULT, "C:/Users/nahaerin/Desktop/git/Research_AUC_CI/mud_0.6.csv")

###### mud = 0.65 ###### 2020.10.08

mud <- 0.65
HM <- DL <- BA <- LI <- array(NA, dim = c(iter,3))
HE <- array(NA, dim = c(iter,5))
RESULT <- array(NA, dim = c(length(rho), 15))
colnames(RESULT) <- c("delta1", "HM_CI_L", "HM_CI_U", "DL_CI_L", "DL_CI_U", "BA_CI_L", "BA_CI_U",
                      "delta2", "LI_CI_L", "LI_CI_U",
                      "delta3", "HE_CI_L", "HE_CI_U", "HE_CIvst_L", "HE_CIvst_U")

k <- 1
for(j in 1:length(rho)) {
  i <- 1
  while(i <= iter) {
    set.seed(k)
    
    nondefault <- mvtnorm::rmvnorm(nn, c(0,0), matrix(c(1,rho[j],rho[j],1), nrow = 2))
    default <- mvtnorm::rmvnorm(nd, c(0.5, mud), matrix(c(1,rho[j],rho[j],1), nrow = 2))
    data <- rbind(default, nondefault)
    
    glm_f <- glm(y ~ data, family = binomial(link = "logit")) 
    z <- summary(glm_f)$coefficients[3,4]  ## p-value
    
    if(z < 0.05) {
      A <- H_Mc(data, y)
      B <- DeLong(data, y)
      C <- Bandos(data, y)
      D <- Li(data, y)
      E <- Heller(data,y)
      
      HM[i,] <- c(A$delta, A$CI_L, A$CI_U)
      DL[i,] <- c(B$delta, B$CI_L, B$CI_U)
      BA[i,] <- c(C$delta, C$CI_L, C$CI_U)
      LI[i,] <- c(D$delta, D$CI_L, D$CI_U)
      HE[i,] <- c(E$delta, E$CI_L, E$CI_U, E$CIvst_L, E$CIvst_U)
      
      cat("rho = ", rho[j], ", iteration = ", i, "\n")
      
      i <- i+1
    } else {  i <- i  }
    
    k <- k + 1    
    
  }
  RESULT[j,] <- c(apply(HM, 2, mean), apply(DL, 2, mean)[-1], apply(BA, 2, mean)[-1], 
                  apply(LI, 2, mean), apply(HE, 2, mean))
}

RESULT <- cbind.data.frame(cbind(mud, rho, RESULT))

write.csv(RESULT, "C:/Users/nahaerin/Desktop/git/Research_AUC_CI/mud_0.65.csv")

###### LRT ######
rm(list=ls())

m <- seq(0.4,0.65,0.05)
r <- seq(-0.8,0.8,0.1)
mud <- rep(m, rep(length(r),length(m)))
rho <- rep(r,length(m))
y <- rep(c(1,0), c(30,30))
iter <- 1000
loglik_f <- loglik_r <- array(NA, dim = c(length(mud), iter))

for(j in 1:length(mud)){
  for(i in 1:iter){
    nondefault <- mvtnorm::rmvnorm(30, c(0,0), matrix(c(1,rho[j],rho[j],1),nrow = 2))
    default <- mvtnorm::rmvnorm(30, c(0.5,mud[j]), matrix(c(1,rho[j],rho[j],1), nrow =2))
    data <- rbind(default, nondefault)
    glm_f <- glm(y ~ data, family = binomial(link = "logit"))
    glm_r <- glm(y ~ data[,1], family = binomial(link = "logit"))
    test <- lmtest::lrtest(glm_r, glm_f)
    loglik_f[j,i] <- test$LogLik[2]    # df=3
    loglik_r[j,i] <- test$LogLik[1]    # df=2
  }
  cat("iter = ", j, "\n")  
}

lrt_result <- cbind(mud, 
                    rho,
                    apply(loglik_r, 1, mean),
                    apply(loglik_f, 1, mean),
                    rep(NA, length(mud)),
                    rep(NA, length(mud)),
                    rep(NA, length(mud)),
                    rep(NA, length(mud)))

colnames(lrt_result) <- c("mud", "rho", "l_r", "l_f", "-2l_r","-2l_f", "-2(l_r-l_f)", "pval")

for(i in 1:length(mud)){
  lrt_result[i,5] <- -2*lrt_result[i,3]
  lrt_result[i,6] <- -2*lrt_result[i,4]
  lrt_result[i,7] <- lrt_result[i,5]-lrt_result[i,6]
  lrt_result[i,8] <- 1-lrt_pchisq(result[i,7],1)
}

write.csv(lrt_result, "C:\\Users\\nahaerin\\Desktop\\git\\Research_AUC_CI\\LRT.csv")


##### same result
#nondefault <- mvtnorm::rmvnorm(30, c(0,0), matrix(c(1,rho[j],rho[j],1),nrow = 2))
#default <- mvtnorm::rmvnorm(30, c(0.5,mud[j]), matrix(c(1,rho[j],rho[j],1), nrow =2))
#data <- rbind(default, nondefault)
#glm_f <- glm(y ~ data, family = binomial(link = "logit"))
#glm_r <- glm(y ~ data[,1], family = binomial(link = "logit"))
#f <- cbind(1,data) %*% as.vector(glm_f$coefficients)
#r <- cbind(1,data[,1]) %*% as.vector(glm_r$coefficients)
#pf <- exp(f)/(1+exp(f))
#pr <- exp(r)/(1+exp(r))
#Lf <- prod(pf^y * (1-pf)^(1-y))
#Lr <- prod(pr^y * (1-pr)^(1-y))
#lf <- log(Lf)
#lr <- log(Lr)
#test <- lmtest::lrtest(glm_r, glm_f)
#loglik_f[j,i] <- test$LogLik[2]    # df=3
#loglik_r[j,i] <- test$LogLik[1]    # df=2
#diff[j,i] <- loglik_f[j,i] - loglik_r[j,i]
#
#lrt_result <- cbind(mud, 
#                rho,
#                apply(loglik_f, 1, mean),
#                apply(loglik_r, 1, mean),
#                apply(diff, 1, mean),
#                rep(NA, length(mud)))
#
#colnames(lrt_result) <- c("mud", "rho", "loglik_f", "loglik_r", "l_f-l_r", "pval")
#
#for(i in 1:length(mud)){
#  lrt_result[i,6] <- 1-pchisq(2*lrt_result[i,5],1)
#  
#}


###################### combine ########################

LRT <- read.csv("C:\\Users\\nahaerin\\Desktop\\git\\Research_AUC_CI\\LRT.csv")
mud4 <- read.csv("C:\\Users\\nahaerin\\Desktop\\git\\Research_AUC_CI\\mud_0.4.csv")
mud45 <- read.csv("C:\\Users\\nahaerin\\Desktop\\git\\Research_AUC_CI\\mud_0.45.csv")
mud5 <- read.csv("C:\\Users\\nahaerin\\Desktop\\git\\Research_AUC_CI\\mud_0.5.csv")
mud55 <- read.csv("C:\\Users\\nahaerin\\Desktop\\git\\Research_AUC_CI\\mud_0.55.csv")
mud6 <- read.csv("C:\\Users\\nahaerin\\Desktop\\git\\Research_AUC_CI\\mud_0.6.csv")
mud65 <- read.csv("C:\\Users\\nahaerin\\Desktop\\git\\Research_AUC_CI\\mud_0.65.csv")


FINAL <- cbind(LRT$pval,rbind(mud4[,-1], mud45[,-1], mud5[,-1], mud55[,-1], mud6[,-1], mud65[,-1]))

View(round(FINAL,4))
HM <- DL <- BA <- LI <- HE1 <- HE2 <- rep("X", nrow(FINAL))
FINAL2 <- data.frame(FINAL$mud, FINAL$rho, FINAL$`LRT$pval`, FINAL$delta1, HM, DL, BA, FINAL$delta2, LI, FINAL$delta3, HE1, HE2)
FINAL2

FINAL2$HM[which(FINAL$HM_CI_L < 0)] <- "O"
FINAL2$DL[which(FINAL$DL_CI_L < 0)] <- "O"
FINAL2$BA[which(FINAL$BA_CI_L < 0)] <- "O"
FINAL2$LI[which(FINAL$LI_CI_L < 0)] <- "O"
FINAL2$HE1[which(FINAL$HE_CI_L < 0)] <- "O"  
FINAL2$HE2[which(FINAL$HE_CIvst_L < 0)] <- "O"

table(FINAL2$HM[FINAL2$FINAL..LRT.pval. >= 0.05])
table(FINAL2$DL[FINAL2$FINAL..LRT.pval. >= 0.05])
table(FINAL2$BA[FINAL2$FINAL..LRT.pval. >= 0.05])
table(FINAL2$LI[FINAL2$FINAL..LRT.pval. >= 0.05])
table(FINAL2$HE1[FINAL2$FINAL..LRT.pval. >= 0.05])
table(FINAL2$HE2[FINAL2$FINAL..LRT.pval. >= 0.05])


table(FINAL2$HM[FINAL2$FINAL..LRT.pval. < 0.05])
table(FINAL2$DL[FINAL2$FINAL..LRT.pval. < 0.05])
table(FINAL2$BA[FINAL2$FINAL..LRT.pval. < 0.05])
table(FINAL2$LI[FINAL2$FINAL..LRT.pval. < 0.05])
table(FINAL2$HE1[FINAL2$FINAL..LRT.pval. < 0.05])
table(FINAL2$HE2[FINAL2$FINAL..LRT.pval. < 0.05])


write.csv(FINAL2, "C:\\Users\\nahaerin\\Desktop\\git\\Research_AUC_CI\\total.csv")
