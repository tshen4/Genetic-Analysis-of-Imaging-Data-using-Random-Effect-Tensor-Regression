
herit_non<-function(ROI1=c(1,2),ROI2=c(3,4),SNP_data,Y,freq_band=1){
  

  n=dim(Y)[1];dim.set=c(length(ROI1),length(ROI2)); D=length(dim.set)
  G=dim(SNP_data)[2];noise=1
  
  
  #get the genetic variants Z:n-by-p (fixed)
  Z=SNP_data
  Z=scale(Z,scale=FALSE)
  
  #get the value of y(fixed)
  Y=dat_array[,ROI1,ROI2,freq_band]
  y1=unfold(as.tensor(Y),row_idx=1,col_idx = seq(2,D+1));y1=y1@data
  y1=scale(y1,scale=FALSE)
  #b.arr[,,iter]=b1
  #y.arr[,,iter]=y1
  
  y1.vec=as.matrix(as.vector(y1));
  
  Z11=diag(prod(dim.set))%x%Z
  
  #the vectorized b(gamma) and epsilon
  #gamma1=as.matrix(as.vector(b1))
  #e1=as.matrix(as.vector(epsilon1))
  
  ##initial values
  E.set=list()
  for(i in 1:1) E.set[[i]]=diag(prod(dim.set))*noise^(length(dim.set))
  
  E=diag(n)
  for(i in 1:1) E=E.set[[i]]%x%E
  
  Sigma.set=list()
  for(i in 1:1) Sigma.set[[i]]=diag(prod(dim.set))
  
  Sigma=diag(G)
  for(i in 1:1) Sigma=Sigma.set[[i]]%x%Sigma
  
  ##the true value of h2
  h2.true=tr(Z%*%t(Z)/n)*tr(Sigma)/(tr(Z%*%t(Z)/n)*tr(Sigma)+tr(E))
  
  V=Z11%*%Sigma%*%t(Z11)+E
  
  Sigma.old=matrix(0,prod(dim.set)*G,prod(dim.set)*G)
  E.old=matrix(0,prod(dim.set)*n,prod(dim.set)*n)
  count=0
  
  while(abs((max(Sigma-Sigma.old))>0.01|abs(max(E-E.old))>0.01)&count<=50){
    count=count+1;
    Sigma.old=Sigma
    E.old=E
    
    ##EM algorithm to update variance components
    E.gamma1=Sigma%*%t(Z11)%*%solve(V)%*%y1.vec
    Var.gamma1=Sigma-Sigma%*%t(Z11)%*%solve(V)%*%Z11%*%Sigma
    E.e=E%*%solve(V)%*%y1.vec
    Var.e=E-E%*%solve(V)%*%E
    A=E.gamma1%*%t(E.gamma1)+Var.gamma1
    B=E.e%*%t(E.e)+Var.e
    
    #store the scalar of total variance 
    var.scalar.sig1=1;var.scalar.sig2=1;
    var.scalar.e1=1;var.scalar.e2=1;
    #update sigma1,2
    L=diag(G)
    #if(D>=2){for(i in 1:(D-1)) L=solve(Sigma.set[[i]])%x%L}
    
    
    for(i in 1:prod(dim.set)){
      for(j in i:prod(dim.set)){
        Sigma.set[[1]][i,j]=tr(L%*%A[((i-1)*dim(L)[1]+1):(i*dim(L)[1]),((j-1)*dim(L)[1]+1):(j*dim(L)[1])])/G
      }
    }
    for(i in 1:prod(dim.set)){
      for(j in i:prod(dim.set)){
        Sigma.set[[1]][j,i]=Sigma.set[[1]][i,j]
      }
    }
    
    
    #Sigma2[1,1]=1
    #var.scalar.sig2=var.scalar.sig2*det(Sigma2)
    #Sigma2=Sigma2/det(Sigma2)
    #Sigma2=Sigma2*(matrix(2,K,K)-diag(K))
    
    
    #Sigma1[1,1]=1
    #var.scalar.sig1=var.scalar.sig1*det(Sigma1)
    #Sigma1=Sigma1/det(Sigma1)
    #Sigma1=Sigma1*(matrix(2,K,K)-diag(K))
    
    #update E's
    L=diag(n)
    #if(D>=2){for(i in 1:(D-1)) L=solve(E.set[[i]])%x%L}
    
    
    for(i in 1:prod(dim.set)){
      for(j in i:prod(dim.set)){
        E.set[[1]][i,j]=tr(L%*%B[((i-1)*dim(L)[1]+1):(i*dim(L)[1]),((j-1)*dim(L)[1]+1):(j*dim(L)[1])])/n
      }
    }
    
    for(i in 1:prod(dim.set)){
      for(j in i:prod(dim.set)){
        E.set[[1]][j,i]=E.set[[1]][i,j]
      }
    }
    #E2[1,1]=1
    #var.scalar.e2=var.scalar.e2*det(E2)
    #E2=E2/det(E2)
    #E2=E2*(matrix(2,K,K)-diag(K))
    
    
    #E1[1,1]=1
    #var.scalar.e1=var.scalar.e1*det(E1)
    #E1=E1/det(E1)
    #E1=E1*(matrix(2,K,K)-diag(K))
    
    ##reconstruct Sigma and E
    Sigma=diag(G)
    for(i in 1:1) Sigma=Sigma.set[[i]]%x%Sigma
    
    E=diag(n)
    for(i in 1:1) E=E.set[[i]]%x%E
    
    V=Z11%*%Sigma%*%t(Z11)+E
    
    #print(Sigma.set[[3]])
    #print(Sigma2%x%Sigma1)
    #print(det(V))
  }
  #loglike[iter]=-0.5*log(det(Sigma))-0.5*log(det(E))-0.5*t(gamma1)%*%solve(Sigma)%*%gamma1-0.5*t(e1)%*%solve(E)%*%e1
  
  #calculate heritability
  h2.non=tr(Z%*%t(Z)/n)*tr(Sigma)/(tr(Z%*%t(Z)/n)*tr(Sigma)+tr(E))
  return(h2.non)
}
