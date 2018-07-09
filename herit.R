
herit<-function(ROI1=c(1,2),ROI2=c(3,4),SNP_data,Y,freq_band=1){
  
  
  n=dim(Y)[1];dim.set=c(length(ROI1),length(ROI2)); D=length(dim.set)
  G=dim(SNP_data)[2];noise=1
  
  #get the genetic variants Z:n-by-p (fixed)
  Z=SNP_data
  Z=scale(Z,scale=FALSE)
  
  #get the value of response(y)
  Y=dat_array[,ROI1,ROI2,freq_band]
  y1=unfold(as.tensor(Y),row_idx=1,col_idx = seq(2,D+1));y1=y1@data
  y1=scale(y1,scale=FALSE)
  y1.vec=as.matrix(as.vector(y1));
  
  Z11=diag(prod(dim.set))%x%Z
  
  #the vectorized b(gamma) and epsilon
  #gamma1=as.matrix(as.vector(b1))
  #e1=as.matrix(as.vector(epsilon1))
  
  ##initial values
  E.set=list()
  for(i in 1:D) E.set[[i]]=diag(dim.set[i])*noise
  
  E=diag(n)
  for(i in 1:D) E=E.set[[i]]%x%E
  
  Sigma.set=list()
  for(i in 1:D) Sigma.set[[i]]=diag(dim.set[i])
  
  Sigma=diag(G)
  for(i in 1:D) Sigma=Sigma.set[[i]]%x%Sigma
  
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
    if(D>=2){for(i in 1:(D-1)) L=solve(Sigma.set[[i]])%x%L}
    
    
    for(i in 1:dim.set[D]){
      for(j in i:dim.set[D]){
        Sigma.set[[D]][i,j]=dim.set[D]*tr(L%*%A[((i-1)*dim(L)[1]+1):(i*dim(L)[1]),((j-1)*dim(L)[1]+1):(j*dim(L)[1])])/G/prod(dim.set)
      }
    }
    
    for(i in 1:dim.set[D]){
      for(j in 1:dim.set[D]){
        if(i>j)
          Sigma.set[[D]][i,j]=Sigma.set[[D]][j,i]
      }
    }
    #Sigma2[1,1]=1
    #var.scalar.sig2=var.scalar.sig2*det(Sigma2)
    #Sigma2=Sigma2/det(Sigma2)
    #Sigma2=Sigma2*(matrix(2,K,K)-diag(K))
    
    Sigma.temp=Sigma.set
    if(D>=2){
      for(p in 1:(D-1)){
        H=solve(Sigma.set[[p+1]])
        if(p<(D-1)){
          for(temp in (p+2):D) H=solve(Sigma.set[[temp]])%x%H
        }
        L=diag(G)
        if(p>1){
          for(temp in 2:p) L=solve(Sigma.set[[temp-1]])%x%L
        }
        
        for(i in 1:dim.set[p]){
          for(j in 1:dim.set[p]){
            Sigma.temp[[p]][i,j]=0
            for(r in 1:dim(H)[1]){
              for(s in 1:dim(H)[1]){
                A.temp=A[((r-1)*dim(L)[1]*dim.set[p]+1):(r*dim(L)[1]*dim.set[p]),((s-1)*dim(L)[1]*dim.set[p]+1):(s*dim(L)[1]*dim.set[p])]
                Sigma.temp[[p]][i,j]=Sigma.temp[[p]][i,j]+dim.set[p]*tr(L%*%A.temp[((i-1)*dim(L)[1]+1):(i*dim(L)[1]),((j-1)*dim(L)[1]+1):(j*dim(L)[1])])/prod(dim.set)/G*H[r,s]
              }
            }
          }
        }
      }
    }
    Sigma.set=Sigma.temp
    
    for(p in 1:D){
      for(i in 1:dim.set[p]){
        for(j in 1:dim.set[p]){
          if(i>j)
            Sigma.set[[p]][i,j]=Sigma.set[[p]][j,i]
        }
      }
    }
    #Sigma1[1,1]=1
    #var.scalar.sig1=var.scalar.sig1*det(Sigma1)
    #Sigma1=Sigma1/det(Sigma1)
    #Sigma1=Sigma1*(matrix(2,K,K)-diag(K))
    
    #update E's
    L=diag(n)
    if(D>=2){for(i in 1:(D-1)) L=solve(E.set[[i]])%x%L}
    
    
    for(i in 1:dim.set[D]){
      for(j in i:dim.set[D]){
        E.set[[D]][i,j]=dim.set[D]*tr(L%*%B[((i-1)*dim(L)[1]+1):(i*dim(L)[1]),((j-1)*dim(L)[1]+1):(j*dim(L)[1])])/n/prod(dim.set)
      }
    }
    
    for(i in 1:dim.set[D]){
      for(j in 1:dim.set[D]){
        if(i>j)
          E.set[[D]][i,j]=E.set[[D]][j,i]
      }
    }
    #E2[1,1]=1
    #var.scalar.e2=var.scalar.e2*det(E2)
    #E2=E2/det(E2)
    #E2=E2*(matrix(2,K,K)-diag(K))
    
    E.temp=E.set
    if(D>=2){
      for(p in 1:(D-1)){
        H=solve(E.set[[p+1]])
        if(p<(D-1)){
          for(temp in (p+2):D) H=solve(E.set[[temp]])%x%H
        }
        L=diag(n)
        if(p>1){
          for(temp in 2:p) L=solve(E.set[[temp-1]])%x%L
        }
        
        for(i in 1:dim.set[p]){
          for(j in 1:dim.set[p]){
            E.temp[[p]][i,j]=0
            for(r in 1:dim(H)[1]){
              for(s in 1:dim(H)[1]){
                B.temp=B[((r-1)*dim(L)[1]*dim.set[p]+1):(r*dim(L)[1]*dim.set[p]),((s-1)*dim(L)[1]*dim.set[p]+1):(s*dim(L)[1]*dim.set[p])]
                E.temp[[p]][i,j]=E.temp[[p]][i,j]+dim.set[p]*tr(L%*%B.temp[((i-1)*dim(L)[1]+1):(i*dim(L)[1]),((j-1)*dim(L)[1]+1):(j*dim(L)[1])])/prod(dim.set)/n*H[r,s]
              }
            }
          }
        }
        
      }
    }
    E.set=E.temp
    
    for(p in 1:D){
      for(i in 1:dim.set[p]){
        for(j in 1:dim.set[p]){
          if(i>j)
            E.set[[p]][i,j]=E.set[[p]][j,i]
        }
      }
    }
    #E1[1,1]=1
    #var.scalar.e1=var.scalar.e1*det(E1)
    #E1=E1/det(E1)
    #E1=E1*(matrix(2,K,K)-diag(K))
    
    ##reconstruct Sigma and E
    Sigma=diag(G)
    for(i in 1:D) Sigma=Sigma.set[[i]]%x%Sigma
    
    E=diag(n)
    for(i in 1:D) E=E.set[[i]]%x%E
    
    V=Z11%*%Sigma%*%t(Z11)+E
    #print(det(V))
  }
  #loglike[iter]=-0.5*log(det(Sigma))-0.5*log(det(E))-0.5*t(gamma1)%*%solve(Sigma)%*%gamma1-0.5*t(e1)%*%solve(E)%*%e1
  
  #calculate heritability
  h2=tr(Z%*%t(Z)/n)*tr(Sigma)/(tr(Z%*%t(Z)/n)*tr(Sigma)+tr(E))
  return(h2)
}


