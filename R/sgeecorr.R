
#' Title Calculate the parameters of the clustered failure data
#'
#' @param formula a formula expression, of the form \code{response ~ predictors}. The \code{response} is a \code{Surv} object with right censoring.
#' @param data a data frame in which to interpret the variables named in the \code{formula}.
#' @param corrstr a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{def}.
#' @param corr_matix When corrstr is \code{independence}, it does not need to be specified. When \code{corrstr} is \code{def}, it needs to be specified.When corrstr is independence, it does not need to be specified. When corrstr is def, it needs to be specified.
#' @param eps tolerance for convergence. The default is \code{eps = 1e-6}
#' @param itermax specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{itermax} iterations and
#' the estimates will be based on the last iteration. The default \code{itermax = 100}.
#'
#' @return the coefficient of the covariate, robust variance
#' @export
#'
#' @import Matrix
#' @importFrom MASS ginv
#' @importFrom survival Surv
sgeecorr= function(formula, data,corrstr = "independence" ,corr_matix=NULL,eps = 1e-06,itermax =100){

  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required but not installed.")
  }

  mf <- stats::model.frame(formula, data)

  Covariates<- stats::model.matrix(attr(mf, "terms"), mf)[,-1]
  beta_name <- colnames(Covariates)
  Y <- stats::model.extract(mf, "response")

  if (!inherits(Y, "Surv")) {
    stop("The response variable must be a 'Surv' object.")
  }
  time <- Y[, 1]
  cens <- Y[, 2]


  cluster2= data
  t2<- time
  c1<- cens
  xxx = matrix(Covariates,nrow(Covariates),ncol(Covariates))

  id<-cluster2$id

  K=length(unique(id))
  n=length(id)/K
  if(corrstr =="independence"){
    matrix_inverse_estimated = diag(n)
  }
  else if(corrstr =="def"){
    matrix_inverse_estimated = ginv(corr_matix)
  }
  cens=c1

  t11=sort(t2)
  c11=c1[order(t2)]
  x111=xxx[order(t2),]

  tt1=unique(t11[c11==1])
  kk=length(table(t11[c11==1]))
  dd=as.matrix(table(t11[c11==1]))

  gg1=rep(1,length(c1))
  g11=gg1[order(t2)]

  X1=xxx
  Y1=matrix(cluster2$cens,ncol=1)


  betainit = matrix(rep(0,ncol(X1)),ncol=1)
  beta2 = matrix(rep(0,ncol(X1)),ncol=1)

  gSSS1=rep(0,kk)
  KK1=1

  repeat{

    gSS=rep(0,kk)
    gSS1=rep(1,kk)

    gSS[1]=dd[1]/(sum(g11[min((1:(K*n))[t11==tt1[1]]):(K*n)]*exp(x111[(min((1:(K*n))[t11==tt1[1]]):(K*n)),]%*%betainit)))

    for (i in 1:(kk-1))
    {
      gSS[i+1]=gSS[i]+dd[i+1]/
        (sum(g11[min((1:(K*n))[t11==tt1[i+1]]):(K*n)]*exp(x111[min((1:(K*n))[t11==tt1[i+1]]):(K*n),]%*%betainit)))
    }

    gSS1=exp(-gSS)

    gSS2=rep(0,K*n)
    gSS3=rep(0,K*n)

    for(i in 1:(K*n))
    {  kk1=1

    if(t2[i]<tt1[1])
    {
      gSS2[i]=1
      gSS3[i]=0.00000001
    }
    else {
      if(t2[i]>=tt1[kk])
      {
        gSS2[i]=0
        gSS3[i]=gSS[kk]
      }
      else {
        repeat{
          if(t2[i]>=tt1[kk1]) kk1=kk1+1
          else break
        }
        {      gSS2[i]=(gSS1[kk1-1])^(exp(xxx[i,]%*%betainit))
          gSS3[i]=gSS[kk1-1]
        }
      }
    }
    }



    Lambda<-gSS3

    W1=diag(Lambda)

    SK1=1
    beta1 =  matrix(rep(0,ncol(X1)),ncol=1)



    repeat{
      mu=exp(X1%*%betainit)
      newY1=c1/Lambda

      res=as.vector((newY1-mu)/sqrt(mu))
      rres=0

      pphi=(sum(res^2)/(K*n-dim(X1)[2]))

      res=matrix(res,ncol=K)
      res=t(res)
      for(i in 1:K)
      {
        for(j in 1:(n-1))
          rres=rres+res[i,j]*sum(res[i,(j+1):n])
      }
      rho=(pphi^(-1))*rres/(K*n*(n-1)/2-dim(X1)[2])

      SK=1

      repeat{
        D1=diag(mu[id==1,])%*%diag(rep(1,n))%*%(X1[id==1,])
        for(i in 2:K)
        { D1=rbind(D1,diag(mu[id==i,])%*%diag(rep(1,n))%*%(X1[id==i,])) }

        S1=newY1-mu

        R1 = ginv(matrix_inverse_estimated)

        V1=sqrt(diag(mu[id==1,]))%*%R1%*%sqrt(diag(mu[id==1,]))*pphi
        for(i in 2:K)
        { V1=bdiag(V1,sqrt(diag(mu[id==i,]))%*%R1%*%sqrt(diag(mu[id==i,]))*pphi) }

        V1=as.matrix(V1)
        Z1=D1%*%betainit+S1

        geebeta=ginv(t(D1)%*%ginv(V1)%*%W1%*%D1)%*%t(D1)%*%ginv(V1)%*%W1%*%Z1

        if(any(abs(geebeta-betainit)>eps) && (SK<=500))
        {
          betainit<-geebeta
          mu=exp(X1%*%betainit)
          SK=SK+1
        }
        else break

      }



      if(any(abs(betainit-beta1)>eps) && (SK1<30))

      {  beta1=betainit
      SK1=SK1+1
      }

      else break
    }


    if ((any(abs(betainit-beta2)>eps) || any(abs(gSS1-gSSS1)>eps) ) && KK1<=itermax)
    {
      beta2<-betainit
      gSSS1<-gSS1
      if(KK1%%5==0){
        print(c(KK1, betainit))
      }
      KK1<-KK1+1
    }


    else  break
  }

  print(c(KK1, betainit))





  ############################
  ##    variance estimate   ##
  ############################

  ###
  #  for beta variance
  ###

  betascale=1
  be=betainit
  gS=c(gSS[1],gSS[2:kk]-gSS[1:(kk-1)])

  gg1=rep(1,K*n)
  xxxx=xxx

  #######################################
  Q1 = ginv(matrix_inverse_estimated)
  IQ1 = matrix_inverse_estimated

  #####################################
  z2=cbind(xxx)

  c2=c1

  mu2=exp(z2%*%betainit)

  B2=matrix(0,n,n)
  ABC1=rep(0,K)
  VA1=matrix(0,dim(z2)[2],dim(z2)[2])

  for(v in 1:(dim(z2)[2]))
  {
    for(w in 1:(dim(z2)[2]))
    {
      for(i in 1:K)
      {
        z22=z2[id==i,]

        A2=t(z22[,v])

        c22=c2[id==i]
        Lam22=Lambda[id==i]
        # t22=t2[id==i]

        mu22=mu2[id==i,]

        BB1=(mu22^(1/2))%*%((t(mu22))^(-1/2))*IQ1

        for(s in 1:n)
        {
          for(l in 1:n)
          {
            B2[s,l]=(1/2)*(z22[s,w]-z22[l,w])*BB1[s,l]
          }
        }

        C2=(c22/Lam22)-mu22
        D2=BB1
        E2=z22[,w]*mu22
        G2=diag(Lam22)
        ABC1[i]=A2%*%(B2%*%G2%*%C2-D2%*%G2%*%E2)

      }
      VA1[v,w]=sum(ABC1)*(betascale^(-1))

      ABC1=rep(0,K)
    }

  }

  sdm=VA1

  BBC=matrix(0,kk,dim(xxxx)[2])

  for(s in 1:(kk))
  {
    BCm=gS[s]*exp(xxxx[(c1==1)&(t2==tt1[s]),]%*%be)

    for(j in 1:dim(xxxx)[2])
    {
      BBC[s,j]=sum(exp(xxxx[(c1==1)&(t2==tt1[s]),]%*%be)*(exp(-BCm)+BCm*exp(-BCm)-1)/((1-exp(-BCm))^2)*xxxx[(c1==1)&(t2==tt1[s]),j])+
        sum(gg1[t2>=tt1[s]]*exp(xxxx[t2>=tt1[s],]%*%be)*xxxx[t2>=tt1[s],j])
    }
  }

  CCC=rep(0,(kk))

  for(s in 1:(kk))
  {
    CCm=gS[s]*exp(xxxx[(c1==1)&(t2==tt1[s]),]%*%be)
    CCC[s]=sum(exp(2*(xxxx[(c1==1)&(t2==tt1[s]),]%*%be)-CCm)/(1-exp(-CCm))^2)
  }




  BC=matrix(0,dim(xxxx)[2],kk)

  for(r in 1:dim(xxxx)[2])
  {

    for(s in 1:(kk))
    {
      elem=0
      for(i in 1:K)
      {
        mu22=mu2[id==i,]
        xxx1=xxx[id==i,r]
        t21=t2[id==i]

        for(j in 1:n)
        {
          if(t21[j]>=tt1[s])
            elem=elem+sum(xxx1*((mu22)^(1/2))*((mu22[j])^(-1/2))*IQ1[,j])*mu22[j]*(betascale^(-1))
        }
      }
      BC[r,s]=elem
    }
  }

  M22=-sdm
  M23=BC
  M32=BBC
  M33=diag(CCC)

  M=rbind(cbind(M22,M23),cbind(M32,M33))

  ###### first derivative vector (fdv) #####

  fdm=0

  for(i in 1:K)
  {
    xxx1=xxxx[id==i,]

    gg11=gg1[id==i]
    c111=c1[id==i]
    t21=t2[id==i]

    g11=Lambda[id==i]

    z22=z2[id==i,]
    mu22=mu2[id==i,]
    mu22m=diag(mu22)

    G2=diag(g11)
    c22=c2[id==i]
    C2=(c22/g11)-mu22

    fdv=t(mu22m%*%z22)%*%ginv(sqrt(mu22m)%*%Q1%*%sqrt(mu22m)*(betascale))%*%G2%*%C2

    d11=rep(0,(kk))

    for(s in 1:(kk))
    {
      d11[s]=sum(exp(xxx1[(c111==1)&(t21==tt1[s]),]%*%be)/(1-exp(-gS[s]*exp(xxx1[(c111==1)&(t21==tt1[s]),]%*%be))))-
        sum(gg11[t21>=tt1[s]]*exp(xxx1[t21>=tt1[s],]%*%be))
    }


    fdm=fdm+t(t(c(fdv,d11)))%*%t(c(fdv,d11))

  }

  vcmR=ginv(M)%*%fdm%*%t(ginv(M))

  V1=diag(vcmR)[1:dim(xxx)[2]]
  V2=(diag(ginv(M)))[1:dim(xxx)[2]]



  sandv=V1
  naivv=V2

  niuCI_lower = betainit-1.96*sqrt(sandv)
  niuCI_upper = betainit+1.96*sqrt(sandv)
  p_value = 2*(1-pnorm(abs(betainit/sqrt(sandv))))
  df_Lambda = data.frame(t(Lambda))
  df_beta = data.frame(betainit,p_value,sandv,rho,pphi,KK1,naivv,niuCI_lower,niuCI_upper)
  row.names(df_beta) = beta_name
  print(df_beta)
  if(corrstr =="independence"){
    beta_Lambda = list(df_Lambda = df_Lambda,df_beta = df_beta)
    return(beta_Lambda)
  }
  else{
    return(df_beta)
  }
}
