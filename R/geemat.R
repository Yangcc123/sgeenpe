
#' Title Calculate the  correlation matrix for clustered failure time
#'
#' @param formula a formula expression, of the form \code{response ~ predictors}. The \code{response} is a \code{Surv} object with right censoring.
#' @param data a data frame in which to interpret the variables named in the \code{formula}.
#' @param beta_Lambda A dataframe contains regression parameters and cumulative risk functions when using an independent working matrix
#'
#' @return correlation matrix for clustered failure time
#' @export
#'
#' @import MASS
#' @importFrom ncvreg ncvfit
geemat= function(formula, data,beta_Lambda){
  library(MASS)
  library(ncvreg)

  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required but not installed.")
  }

  # 构建模型框架（model frame）
  mf <- stats::model.frame(formula, data)

  # 提取设计矩阵 X（协变量）
  Covariates<- stats::model.matrix(attr(mf, "terms"), mf)[,-1]

  # 提取响应变量 Y（Surv 对象）
  Y <- stats::model.extract(mf, "response")

  # 检查 Y 是否是 Surv 对象
  if (!inherits(Y, "Surv")) {
    stop("The response variable must be a 'Surv' object.")
  }
  time <- Y[, 1]
  cens <- Y[, 2]


  cluster2= data
  t2<- time
  c1<- cens
  xxx = matrix(Covariates,nrow(Covariates),ncol(Covariates))

  t2<-cluster2$time
  c1<-cluster2$cens
  id<-cluster2$id

  K=length(unique(id))
  n=length(id)/K


  betainit=matrix(beta_Lambda$df_beta$betainit,ncol=1)

  new_lambda = beta_Lambda$df_Lambda
  new_lambda_vec = unlist(new_lambda[1,])
  indepence_y = c1/new_lambda_vec

  X1=xxx

  Lambda = new_lambda_vec
  newY1=c1/Lambda

  W1=diag(Lambda)
  mu=exp(X1%*%betainit)
  S1=newY1-mu

  D1=diag(mu[id==1])%*%diag(rep(1,n))%*%(X1[id==1,])
  for(i in 2:K)
  { D1=rbind(D1,diag(mu[id==i])%*%diag(rep(1,n))%*%(X1[id==i,]))
  }


  b = mu
  B =sqrt(ginv(diag(as.vector(b))))
  B_diag=diag(B)

  sample_correlative_matrix = function(K,n,t2){
    mat_time = matrix(0,nrow =K,ncol = n )
    for (i in 1:n) {
      for (j in 1:K) {
        #cluster2["id"==1][1]
        mat_time[j,i] = t2[id==j][i]
      }
    }
    return(cor(mat_time,method = "kendall"))
  }



  cor_matrix = sample_correlative_matrix(K,n,(indepence_y))

  cor_matrix_inverse = ginv(cor_matrix)

  eigen_vectors_multiply = function(cor_matrix,m,n){
    #Solve the eigenvectors of the correlation matrix and construct the m related basis matrices
    eigen_result <- eigen(cor_matrix)


    eigen_values <- eigen_result$values
    eigen_vectors <- eigen_result$vectors


    sorted_indices <- order(eigen_values, decreasing = TRUE)
    sorted_eigen_values <- eigen_values[sorted_indices]
    sorted_eigen_vectors <- eigen_vectors[, sorted_indices]

    selected_eigen_vectors <- sorted_eigen_vectors[, 1:m]


    result_matrices <- lapply(1:m, function(i) selected_eigen_vectors[, i] %*% t(selected_eigen_vectors[, i]))

    mid = c(list(diag(n)),result_matrices)
    return(mid)
  }

  Y_variabel = function(K,D1,B_diag,cor_matrix_inverse,Lambda,S1){

    Y =numeric(0)
    for (i in 1:K){
      tem_y = (t(matrix((D1[id==i,]),ncol=ncol(xxx)))) %*% (diag(B_diag[id==i])) %*% cor_matrix_inverse %*%  (diag(B_diag[id==i])) %*% diag(Lambda[id==i]) %*% S1[id==i]
      Y = c(Y,tem_y)
    }
    return(matrix(Y,ncol = 1))
  }

  #----------------------------------------------#
  X_variable = function(m,K,D1,B_diag,mid,Lambda,S1)
  {
    m=m+1
    X  = matrix(0,ncol(xxx)*K,m)
    for (j in 1:m){
      x1_col = numeric(0)
      for(i in 1:K){
        tem_x = (t(matrix((D1[id==i,]),ncol=ncol(xxx)))) %*% (diag(B_diag[id==i])) %*% mid[[j]] %*% (diag(B_diag[id==i])) %*% diag(Lambda[id==i]) %*% S1[id==i]
        x1_col = c(x1_col,tem_x)
      }
      X[,j] = x1_col
    }
    return(X)
  }

  matrix_estimated = function(n,m,beta_eatimate,mid)
  {
    corr_matrix_estimated = matrix(0,ncol=n,nrow = n)
    m = m+1
    for (i in 1:m){
      corr_matrix_estimated = corr_matrix_estimated + beta_eatimate[i]*mid[[i]]
    }
    return(corr_matrix_estimated)
  }
  ######################################################

  BIC_lamada = function(K,beta_eatimate,cor_matrix_inverse,corr_matrix_estimated)

  {
    non_zero_count = length(which(unlist(list(beta_eatimate)) != 0))

    BIC = K*(norm(as.vector(cor_matrix_inverse-corr_matrix_estimated), type = "2"))/norm(as.vector(cor_matrix_inverse+corr_matrix_estimated), type = "2")+log(K)*non_zero_count
    return(list(BIC,corr_matrix_estimated))
  }

  result_k_bda = data.frame()
  for (m in 2:n){
    mid = eigen_vectors_multiply(cor_matrix =cor_matrix,m = m,n = n)
    Y_1 =Y_variabel(K = K,D1 = D1,B_diag = B_diag,cor_matrix_inverse = cor_matrix_inverse,Lambda = Lambda,S1 = S1)
    X_1 = X_variable(m=m,K=K,D1 = D1,B_diag = B_diag,mid = mid,Lambda = Lambda,S1 = S1)
    for (bda in seq(0,1,0.002)){

      fit <- ncvfit(X_1, Y_1, penalty="SCAD",lambda=bda,penalty.factor = c(0, rep(0.5, ncol(X_1)-1))) #c(1,1,1)c(0,rep(1,m))
      corr_matrix_estimated = matrix_estimated(n=n,m=m,beta_eatimate = fit$beta,mid = mid)
      bic_bda = BIC_lamada(K = K ,beta_eatimate = fit$beta,cor_matrix_inverse = cor_matrix_inverse,corr_matrix_estimated = corr_matrix_estimated)

      pair = c(m,bda,bic_bda[[1]])
      result_k_bda = rbind(result_k_bda, pair)
    }
  }
  # The minimum value of the third column
  min_value <- min(result_k_bda[,3])
  row_index <- which(result_k_bda[,3] == min_value)[1]
  m = result_k_bda[row_index,1]
  mid = eigen_vectors_multiply(cor_matrix =cor_matrix,m = m,n = n)
  X_1 = X_variable(m=m,K=K,D1 = D1,B_diag = B_diag,mid = mid,Lambda = Lambda,S1 = S1)
  Y_1 =Y_variabel(K = K,D1 = D1,B_diag = B_diag,cor_matrix_inverse = cor_matrix_inverse,Lambda = Lambda,S1 = S1)
  fit <- ncvfit(X_1, Y_1, penalty="SCAD",lambda=result_k_bda[row_index,2],penalty.factor = c(0, rep(0.5, ncol(X_1)-1)))
  matrix_inverse_estimated = matrix_estimated(n = n,m=m,beta_eatimate = fit$beta,mid = mid)
  corr_mat= ginv(matrix_inverse_estimated)
  return(corr_mat)
}
