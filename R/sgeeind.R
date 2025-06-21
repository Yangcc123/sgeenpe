#' @title Calculate the parameters of the correlated failure data
#' @description
#' Calculate the parameters of the correlated failure data
#' in the case of the inverse of the independent working matrix or the inverse of the given working matrix.
#'
#' @param cluster2 ensure that the first three columns are the ID, deletion status, and expiration time of the cluster failure time, followed by the covariates.
#' @param corr_inv Inverse of the corresponding working matrix
#' @param return_ind Return the cumulative hazard function in an independent situation.
#'
#' @return the coefficient of the covariate, robust variance
#' @export
#' @import Matrix
#' @importFrom MASS ginv
#'
#' @examples sgeeind(data,corr_inv,TRUE)
sgeeind = function(cluster2, corr_inv,return_ind =TRUE){
  #library(Matrix)
  #library(MASS)
  t2<-cluster2$time
  c1<-cluster2$cens
  #xxx<-matrix(c(cluster2$Gender,cluster2$Decayed ,cluster2$bleed_ave,cluster2$Smoking, cluster2$sAge),dim(cluster2)[1],5)
  #loadmid = as.matrix(cluster2[,c(4:ncol(cluster2))])
  loadmid = as.matrix(cluster2[, -c(1:3)])
  xxx = matrix(loadmid,nrow(loadmid),ncol(loadmid))

  id<-cluster2$id

  K=length(unique(id))
  n=length(id)/K
  matrix_inverse_estimated = corr_inv
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

  ###################


  # betainit=matrix(c(-0.425,0.341,-0.846),ncol=1)

  betainit=matrix(c(0,0,0,0),ncol=1)
  #betainit = matrix(coxphresult$coef,ncol=1)

  beta2=matrix(c(0,0,0,0),ncol=1)
  #beta2 = matrix(coxphresult$coef,ncol=1)
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
    beta1=matrix(c(0,0,0,0),ncol=1)
    #beta1 = matrix(coxphresult$coef,ncol=1)


    repeat{
      mu=exp(X1%*%betainit)
      newY1=c1/Lambda

      res=as.vector((newY1-mu)/sqrt(mu))
      rres=0

      pphi=(sum(res^2)/(K*n-dim(X1)[2]))
      #  pphi=1

      res=matrix(res,ncol=K)
      res=t(res)
      for(i in 1:K)
      {
        for(j in 1:(n-1))
          rres=rres+res[i,j]*sum(res[i,(j+1):n])
      }

      rho=(pphi^(-1))*rres/(K*n*(n-1)/2-dim(X1)[2])
      #  rho=0

      SK=1

      repeat{
        D1=diag(mu[id==1,])%*%diag(rep(1,n))%*%(X1[id==1,])
        for(i in 2:K)
        { D1=rbind(D1,diag(mu[id==i,])%*%diag(rep(1,n))%*%(X1[id==i,])) }

        S1=newY1-mu

        #R1=matrix(rho,n,n)
        #diag(R1)=1
        #####################################，使用估计的相关结构
        #R1 = diag(n)
        R1 = ginv(matrix_inverse_estimated)
        #############################################
        V1=sqrt(diag(mu[id==1,]))%*%R1%*%sqrt(diag(mu[id==1,]))*pphi
        for(i in 2:K)
        { V1=bdiag(V1,sqrt(diag(mu[id==i,]))%*%R1%*%sqrt(diag(mu[id==i,]))*pphi) }

        V1=as.matrix(V1)
        Z1=D1%*%betainit+S1

        geebeta=solve(t(D1)%*%solve(V1)%*%W1%*%D1)%*%t(D1)%*%solve(V1)%*%W1%*%Z1

        if(any(abs(geebeta-betainit)>1e-6) && (SK<=500))
        {
          betainit<-geebeta
          mu=exp(X1%*%betainit)
          SK=SK+1
        }
        else break

      }



      if(any(abs(betainit-beta1)>0.000001) && (SK1<30))

      {  beta1=betainit
      # mu=exp(X1%*%betainit)
      SK1=SK1+1
      }

      else break
    }


    if (any(abs(betainit-beta2)>0.000001) || any(abs(gSS1-gSSS1)>0.000001) )
    {
      beta2<-betainit
      gSSS1<-gSS1
      # Lambda<-gSS3
      print(c(KK1, betainit))
      KK1<-KK1+1
    }


    else  break
  }







  ############################
  ##    variance estimate   ##
  ############################

  ###
  #  for beta variance
  ###

  betacorr=rho
  # betascale=pphi
  betascale=1
  # betacorr=0

  be=betainit
  gS=c(gSS[1],gSS[2:kk]-gSS[1:(kk-1)])

  gg1=rep(1,K*n)
  xxxx=xxx

  #  weight1=Lambda

  #Q1=matrix(betacorr,n,n)
  #diag(Q1)=1
  #IQ1=solve(Q1)
  #######################################,使用估计的出来的相关性结构和逆
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

  #  BBC=rep(0,(kk))

  #  for(s in 1:(kk))
  #   {
  #      BCm=gS[s]*exp(be*xxxx[(c1==1)&(t2==tt1[s])])
  #      BBC[s]=sum(exp(be*xxxx[(c1==1)&(t2==tt1[s])])*(exp(-BCm)+BCm*exp(-BCm)-1)/((1-exp(-BCm))^2)*xxxx[(c1==1)&(t2==tt1[s])])+
  #             sum(gg1[t2>=tt1[s]]*exp(be*xxxx[t2>=tt1[s]])*xxxx[t2>=tt1[s]])
  #   }

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

  #  BC=rep(0,kk)

  #  for(s in 1:(kk))
  #   {
  #      elem=0
  #      for(i in 1:K)
  #        {
  #           mu22=mu2[id==i,]
  #           xxx1=xxx[id==i]
  #           t21=t2[id==i]

  #           for(j in 1:n)
  #            {
  #               if(t21[j]>=tt1[s])
  #               elem=elem+sum(xxx1*((mu22)^(1/2))*((mu22[j])^(-1/2))*IQ1[,j])*mu22[j]*(betascale^(-1))
  #            }
  #        }
  #      BC[s]=elem
  #   }


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

    fdv=t(mu22m%*%z22)%*%solve(sqrt(mu22m)%*%Q1%*%sqrt(mu22m)*(betascale))%*%G2%*%C2

    d11=rep(0,(kk))

    for(s in 1:(kk))
    {
      d11[s]=sum(exp(xxx1[(c111==1)&(t21==tt1[s]),]%*%be)/(1-exp(-gS[s]*exp(xxx1[(c111==1)&(t21==tt1[s]),]%*%be))))-
        sum(gg11[t21>=tt1[s]]*exp(xxx1[t21>=tt1[s],]%*%be))
    }


    fdm=fdm+t(t(c(fdv,d11)))%*%t(c(fdv,d11))

  }

  vcmR=solve(M)%*%fdm%*%t(solve(M))

  V1=diag(vcmR)[1:dim(xxx)[2]]
  V2=(diag(solve(M)))[1:dim(xxx)[2]]



  sandv=V1
  naivv=V2

  niuCI_lower = betainit-1.96*sqrt(sandv)
  niuCI_upper = betainit+1.96*sqrt(sandv)
  p_value = 2*(1-pnorm(abs(betainit/sqrt(sandv))))

  if(return_ind){
    df1 = data.frame(t(Lambda))
    df2 = data.frame(betainit,rho,pphi,KK1,sandv,naivv,niuCI_lower,niuCI_upper)
    return(dataframe_indepence = cbind(df2, df1))
  }
  else{
    return(dataframe_QU = data.frame(betainit,p_value,sandv,rho,pphi,KK1,naivv,niuCI_lower,niuCI_upper))
  }
}
