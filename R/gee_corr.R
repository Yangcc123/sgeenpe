#' @title Calculate the inverse of the correlation matrix
#'
#' @param cluster2 ensure that the first three columns are the ID, deletion status, and expiration time of the cluster failure time, followed by the covariates.
#' @param dataframe_indepence In an independent situation, the cumulative hazard function and parameter beta
#'
#' @return inverse of the correlation matrix
#' @export
#' @import MASS
#' @importFrom ncvreg ncvfit
#' @examples gee_corr(cluster2,dataframe_indepence)
gee_corr= function(cluster2,dataframe_indepence){
  library(MASS)
  library(ncvreg)
  t2<-cluster2$time
  c1<-cluster2$cens
  #xxx<-matrix(c(cluster2$Gender,cluster2$Decayed ,cluster2$bleed_ave,cluster2$Smoking, cluster2$sAge),dim(cluster2)[1],5)
  #loadmid = as.matrix(cluster2[,c(4:ncol(cluster2))])
  loadmid = as.matrix(cluster2[, -c(1:3)])
  xxx = matrix(loadmid,nrow(loadmid),ncol(loadmid))


  id<-cluster2$id

  K=length(unique(id))
  n=length(id)/K

  #这里的betainint为列向量
  betainit=matrix(dataframe_indepence$betainit,ncol=1)#设定的beta值，需要在先使用在独立结构下的beta系数

  new_lambda = dataframe_indepence[1,-(1:8)]
  new_lambda_vec = unlist(new_lambda[1,])
  indepence_y = c1/new_lambda_vec

  X1=xxx

  Lambda = new_lambda_vec  #失效时间400*1
  newY1=c1/Lambda #c1为删失指标400*1， Lambda为失效时间400*1，对应的Ki

  W1=diag(Lambda)
  mu=exp(X1%*%betainit)#对应相关的均值
  S1=newY1-mu#就是ki-u_new(Xi)但是这里计算的所有的差值



  #第一个变量值---------------------u_new关于参数beta求导的值
  D1=diag(mu[id==1])%*%diag(rep(1,n))%*%(X1[id==1,])
  #计算的是u_new关于参数beta求导的值
  #rbind是按列进行拼接，最终获得400*1
  for(i in 2:K)
  { D1=rbind(D1,diag(mu[id==i])%*%diag(rep(1,n))%*%(X1[id==i,]))
  }
  #---------------------------------------------------------------#D1

  #第二个变量值--------------------对角矩阵b的负二之一逆
  b = mu #直接就是均值矩阵，是多维的情况下，将方差矩阵分解出来的散布参数和矩阵矩阵
  #b是一个400*1的列矩阵，首先转换成向量，转化成对角矩阵，之后求逆和开根号，因为对角矩阵的缘故，可以整体进行操作
  B =sqrt(solve(diag(as.vector(b))))
  #对角矩阵，使用diag函数，可以直接获得对角元素，会变成一个向量
  B_diag=diag(B)
  #diag(B_diag[id==1])#使用布尔索引，索引出来第一个类中的两个元素
  #---------------------------------------------------------------#B_diag


  ##计算样本相关性矩阵的逆
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



  cor_matrix = sample_correlative_matrix(K,n,(indepence_y-mu))
  # 求解相关性矩阵的逆
  cor_matrix_inverse = ginv(cor_matrix)

  eigen_vectors_multiply = function(cor_matrix,m,n){
    #求解相关性矩阵的特征向量并构造相关的m个基矩阵
    # 计算特征值和特征向量
    # n为群内的个体数，下面需要添加维度为n的单位矩阵在result_matrices，方便后续的运算
    eigen_result <- eigen(cor_matrix)

    # 提取特征值和特征向量
    eigen_values <- eigen_result$values
    eigen_vectors <- eigen_result$vectors

    # 根据特征值的大小排序特征向量
    sorted_indices <- order(eigen_values, decreasing = TRUE)
    sorted_eigen_values <- eigen_values[sorted_indices]
    sorted_eigen_vectors <- eigen_vectors[, sorted_indices]

    # 选择前k个特征向量
    selected_eigen_vectors <- sorted_eigen_vectors[, 1:m]

    # 计算相乘得到的矩阵，%*%矩阵相乘，t是转置
    result_matrices <- lapply(1:m, function(i) selected_eigen_vectors[, i] %*% t(selected_eigen_vectors[, i]))

    # 输出结果
    #for (i in 1:k) {
    #print(result_matrices[[i]])
    #}
    mid = c(list(diag(n)),result_matrices)#这里的n表示群内个体数
    return(mid)
  }

  Y_variabel = function(K,D1,B_diag,cor_matrix_inverse,Lambda,S1){
    # 第 一个小k为特征及矩阵的个数
    # K为群的个数
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
    X  = matrix(0,ncol(xxx)*K,m)############################################3个协变量，这里需要更改
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
    for (i in 1:m){#m+1,基矩阵个数的m，多一个
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
      #X1 是协变变量矩阵，X_1 是基矩阵的值
      fit <- ncvfit(X_1, Y_1, penalty="SCAD",lambda=bda,penalty.factor = c(0, rep(0.5, ncol(X_1)-1))) #c(1,1,1)c(0,rep(1,m))
      corr_matrix_estimated = matrix_estimated(n=n,m=m,beta_eatimate = fit$beta,mid = mid)
      bic_bda = BIC_lamada(K = K ,beta_eatimate = fit$beta,cor_matrix_inverse = cor_matrix_inverse,corr_matrix_estimated = corr_matrix_estimated)
      #print(bic_bda[[2]])
      pair = c(m,bda,bic_bda[[1]])
      result_k_bda = rbind(result_k_bda, pair)
    }
  }
  # 找到第三列的最小值
  min_value <- min(result_k_bda[,3])
  row_index <- which(result_k_bda[,3] == min_value)[1]
  m = result_k_bda[row_index,1]
  mid = eigen_vectors_multiply(cor_matrix =cor_matrix,m = m,n = n)
  X_1 = X_variable(m=m,K=K,D1 = D1,B_diag = B_diag,mid = mid,Lambda = Lambda,S1 = S1)
  Y_1 =Y_variabel(K = K,D1 = D1,B_diag = B_diag,cor_matrix_inverse = cor_matrix_inverse,Lambda = Lambda,S1 = S1)
  fit <- ncvfit(X_1, Y_1, penalty="SCAD",lambda=result_k_bda[row_index,2],penalty.factor = c(0, rep(0.5, ncol(X_1)-1)))
  matrix_inverse_estimated = matrix_estimated(n = n,m=m,beta_eatimate = fit$beta,mid = mid)
  return(matrix_inverse_estimated)
}
