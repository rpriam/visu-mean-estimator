
getQ <- function(Ys,C0_2,C1_2,C2_2,C01,C02,C12,lb) {
  
  Q = matrix(c(
    1       ,   0                  , 0                   , 0          , lb*C1_2        ,  lb*C2_2            ,  lb*C01       , lb*C02  ,  lb*C12 ,
    0       ,   lb*C0_2            , lb*C01              , lb*C02    ,    0           ,   0                 ,  0              , 0      ,  0 ,
    0       ,   lb*C01            , lb*C1_2             , lb*C12     ,    0           ,   0                 ,  0              , 0      ,  0 ,
    0       ,   lb*C02           , lb*C12              , lb*C2_2    ,    0           ,   0                 ,  0              , 0      ,  0 ,
    lb*C1_2  ,   0                  , 0                   , 0         ,  lb^2*C1_2^2     , lb^2*C1_2*C2_2      ,  lb^2*C01*C1_2  , lb^2*C02*C1_2  ,  lb^2*C1_2*C12 ,
    lb*C2_2 ,   0                  , 0                   , 0          , lb^2*C1_2*C2_2 , lb^2*C2_2^2           ,  lb^2*C01*C2_2  , lb^2*C02*C2_2  ,  lb^2*C12*C2_2 ,
    lb*C01   ,   0                  , 0                   , 0          , lb^2*C01*C1_2  ,  lb^2*C01*C2_2      ,  lb^2*C01^2      , lb^2*C01*C02 ,  lb^2*C01*C12 ,
    lb*C02  ,   0                  , 0                   , 0          ,  lb^2*C02*C1_2 ,  lb^2*C02*C2_2      ,  lb^2*C01*C02   ,lb^2*C02^2  ,  lb^2*C02*C12 ,
    lb*C12  ,   0                  , 0                   , 0          , lb^2*C1_2*C12  ,  lb^2*C12*C2_2      ,  lb^2*C01*C12   , lb^2*C02*C12 ,  lb^2*C12^2),
    nrow=9,ncol=9,byrow = TRUE);    
  return(Q);
}

Y_RPR03 <- function(a,b,Y_bar,X1_bar,X2_bar,Cy_2,Cx1_2,Cx2_2,Cx1x2,lb) {
  A = Y_bar^2 + Y_bar^2*lb*(Cy_2+(a^2+2*b)*Cx2_2+4*a*Cyx2);
  B = lb*Cx1_2;
  C = 2*Y_bar*lb*(a*Cx1x2+0.5*Cyx1);
  D0 = +Y_bar^2 + Y_bar^2*lb*(b*Cx2_2+a*Cyx2);
  D1 = +a*Y_bar*lb*Cx1x2;
  E = Y_bar^2;
  
  k0_ = (B*D0-C*D1)/(A*B-C^2);
  k1_ = (A*D1-C*D0)/(A*B-C^2);
  
  MSE_ = k0_^2*A+k1_^2*B+2*k0_*k1_*C-2*k0_*D0-2*k1_*D1+E;
  BIAS_ = k0_*Y_bar*( (k0_-1)/k0_ + b*lb*Cx2_2+a*lb*Cyx2 ) + k1_*a*lb*Cx1x2;
  return(list(MSE=MSE_,BIAS=BIAS_,k0=k0_, k1=k1_));
}

Y_RPR03_k01 <- function(a,b,Y_bar,X1_bar,X2_bar,Cy_2,Cx1_2,Cx2_2,Cx1x2,lb) {
  A = Y_bar^2 + Y_bar^2*lb*(Cy_2+(a^2+2*b)*Cx2_2+4*a*Cyx2);
  B = lb*Cx1_2;
  C = 2*Y_bar*lb*(a*Cx1x2+0.5*Cyx1);
  D0 = +Y_bar^2 + Y_bar^2*lb*(b*Cx2_2+a*Cyx2);
  D1 = +a*Y_bar*lb*Cx1x2;
  E = Y_bar^2;
  
  k0_ = 1;
  k1_ = (D1-C)/B;
  
  MSE_ = k0_^2*A+k1_^2*B+2*k0_*k1_*C-2*k0_*D0-2*k1_*D1+E;
  BIAS_ = k0_*Y_bar*( (k0_-1)/k0_ + b*lb*Cx2_2+a*lb*Cyx2 ) + k1_*a*lb*Cx1x2;
  return(list(MSE=MSE_,BIAS=BIAS_,k0=k0_, k1=k1_));
}

Y_MSE_MIN_SUMf1f2 <- 
  function(a1,b1,a2,b2,Y_bar,X1_bar,X2_bar,Cy_2,Cx1_2,Cx2_2,Cx1x2,lb) { # k1, k2 no constrained
    A = a1^2*Cx1_2;
    B = a2^2*Cx2_2;
    C = a1*a2*Cx1x2;
    D1 = -a1*Cyx1;
    D2 = -a2*Cyx2;
    E = Cy_2;
    
    k1 = (B*D1-C*D2)/(A*B-C^2);
    k2 = (A*D2-C*D1)/(A*B-C^2);
    
    MSE_ = 
      Y_bar^2 * lb * (k1^2*A+k2^2*B+2*k1*k2*C-2*k1*D1-2*k2*D2+E );
    MSE2_ = 
      Y_bar^2 * lb * ( Cy_2 - (B*D1^2-2*C*D1*D2+A*D2^2)/(A*B-C^2) );
    MSE3_ = Y_bar^2 * lb * 
      ( Cy_2 - (Cx2_2*Cyx1^2-2*Cx1x2*Cyx1*Cyx2+Cx1_2*Cyx2^2)/(Cx1_2*Cx2_2-Cx1x2^2) );
    BIAS_ = Y_bar * lb * 
      ( k1* (a1*Cyx1 + b1*Cx1_2 ) + k2* (a2*Cyx2 + b2*Cx2_2 ) ) ;
    return(list(MSE=MSE_,MSE2_=MSE2_,MSE3_=MSE3_,BIAS=BIAS_,alpha1=k1,alpha2=k2));
  }

Y_MSE_MIN_SUMf1f2_usualGEN <- 
  function(a1,b1,a2,b2,Ys,C0_2,C1_2,C2_2,C01,C02,C12,lb) {
    A  = Ys^2 + Ys^2*lb*(C0_2+(a1^2+2*b1)*C1_2+4*a1*C01)
    B  = Ys^2 + Ys^2*lb*(C0_2+(a2^2+2*b2)*C2_2+4*a2*C02)
    C  = Ys^2 + Ys^2*lb*(C0_2+b1*C1_2+b2*C2_2+2*a1*C01+2*a2*C02+a1*a2*C12);
    D0 = Ys^2 + Ys^2*lb*(b1*C1_2+a1*C01)
    D1 = Ys^2 + Ys^2*lb*(b2*C2_2+a2*C02)
    E  = Ys^2;
    
    k0_ = (B*D0-C*D1)/(A*B-C^2);
    k1_ = (A*D1-C*D0)/(A*B-C^2);
    
    MSE_ = k0_^2*A+k1_^2*B+2*k0_*k1_*C-2*k0_*D0-2*k1_*D1+E;
    BIAS_ = -1;
    return(list(MSE=MSE_,BIAS=BIAS_,k0=k0_, k1=k1_));
    
  }

Y_MSE_MIN_SUMf1f2_extendolkin <- 
  function(a1,b1,a2,b2,Y_bar,X1_bar,X2_bar,Cy_2,Cx1_2,Cx2_2,Cx1x2,lb) {
    A = a1^2*Y_bar^2*lb*Cx1_2;
    B = a2^2*Y_bar^2*lb*Cx2_2;
    C = a1*a2*Y_bar^2*lb*Cx1x2;
    D1 = -a1*Y_bar^2*lb*Cyx1;
    D2 = -a2*Y_bar^2*lb*Cyx2;
    E = Y_bar^2*lb*Cy_2;
    
    U = a2*Cyx2-a1*Cyx1+a2^2*Cx2_2-a1*a2*Cx1x2;
    V = a1^2*Cx1_2+a2^2*Cx2_2-2*a1*a2*Cx1x2;
    
    k1 = U/V;
    k2 = 1-k1;
    
    MSE_ = ( k1^2*A+k2^2*B+2*k1*k2*C-2*k1*D1-2*k2*D2+E );
    BIAS_ = Y_bar * lb * ( k1* (a1*Cyx1 + b1*Cx1_2 ) + 
                             k2* (a2*Cyx2 + b2*Cx2_2 ) ) ;
    return(list(MSE=MSE_,BIAS=BIAS_,alpha1=k1,alpha2=k2));
  }

plotpdf <- function(xs,ys,xlab_,abval_,namefile_) {
  pdf(namefile_);
  par(lwd=3)
  plot(xs,ys,type='l',xlab=xlab_, ylab= " ",
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  abline(v=abval_,col="red");
  dev.off();
}

Y_RPR_RC <- function(a,b,Ys,C0_2,C1_2,C2_2,C01,C02,C12,lb) {
  A  = Ys^2 + Ys^2*lb*(C0_2+(a^2+2*b)*C2_2+4*a*C02);
  B  = lb*C1_2;
  C  = 2*Ys*lb*(a*C12+0.5*C01);
  D0 = Ys^2 + Ys^2*lb*(b*C2_2+a*C02);
  D1 = a*Ys*lb*C12;
  E  = Ys^2;
  
  k0_ = (B*D0-C*D1)/(A*B-C^2);
  k1_ = (A*D1-C*D0)/(A*B-C^2);
  
  MSE_ = k0_^2*A+k1_^2*B+2*k0_*k1_*C-2*k0_*D0-2*k1_*D1+E;
  BIAS_ = k0_*Y_bar*( (k0_-1)/k0_ + b*lb*Cx2_2+a*lb*Cyx2 ) + k1_*a*lb*Cx1x2;
  return(list(MSE=MSE_,BIAS=BIAS_,k0=k0_, k1=k1_));
}

Y_RPR_RC_L <- function(a,b,Ys,C0_2,C1_2,C2_2,C01,C02,C12,lb) {
  A  = Ys^2 + Ys^2*lb*(C0_2+(a^2+2*b)*C2_2+4*a*C02) + Ys^2*lb^2*(a*C02+b*C2_2)^2;
  B  = lb*C1_2 + a^2*lb^2*C12^2;
  C  = 2*Ys*lb*(a*C12+0.5*C01) + a*b*Ys*lb^2*C2_2*C12 + a^2*Ys*lb^2*C02*C12;
  D0 = Ys^2 + Ys^2*lb*(b*C2_2+a*C02);
  D1 = a*Ys*lb*C12;
  E  = Ys^2;
  
  k0_ = (B*D0-C*D1)/(A*B-C^2);
  k1_ = (A*D1-C*D0)/(A*B-C^2);
  
  MSE_ = k0_^2*A+k1_^2*B+2*k0_*k1_*C-2*k0_*D0-2*k1_*D1+E;
  BIAS_ = -1;
  return(list(MSE=MSE_,BIAS=BIAS_,k0=k0_, k1=k1_));
}

Y_RPR_RC4_L <- function(a1,b1,a2,b2,Ys,C0_2,C1_2,C2_2,C01,C02,C12,lb) {
  A = Ys^2 * ( 1 + lb*C0_2 );
  B = lb*C1_2 + lb^2*(a2^2*C12^2);
  C = lb*C2_2 + lb^2*(a1^2*C12^2);
  D = Ys*lb*(a2*C12+C01);
  E = Ys*lb*(a1*C12+C02);
  F = lb*(C12+a1*a2*lb*C12^2);
  G0 = Ys^2;
  G1 = Ys*lb*(a2*C12);
  G2 = Ys*lb*(a1*C12);
  H = Ys^2;
  
  M = matrix(c(A,D,E,D,B,F,E,F,C),ncol = 3,nrow=3);
  k012 = solve(M)%*%cbind(c(G0,G1,G2));
  k0=k012[1]
  k1=k012[2]
  k2=k012[3]
  MSE_  = A*k0^2+B*k1^2+C*k2^2+2*D*k0*k1+2*E*k0*k2+2*F*k1*k2-2*G0*k0-2*G1*k1-2*G2*k2+H;
  BIAS_ = -1;
  return(list(MSE=MSE_,BIAS=BIAS_,k0=k0,k1=k1,k2=k2));
}

MSE_2vars_3param_GEN <- function(NN,Ys,C0_2,C1_2,C2_2,C01,C02,C12,lb) {
  Q = getQ(Ys,C0_2,C1_2,C2_2,C01,C02,C12,lb);
  
  MMgen = t(NN[,2:4])%*%Q%*%(NN[,2:4])
  Ggen  = -t(NN[,1])%*%Q%*%(NN[,2:4])
  
  k012_gen = solve(MMgen)%*%t(Ggen);
  k0_gen=k012_gen[1]
  k1_gen=k012_gen[2]
  k2_gen=k012_gen[3]
  MSE_gen_est  = c(t(NN%*%rbind(1,k012_gen))%*%Q%*%(NN%*%rbind(1,k012_gen)));
  
  PP=NN%*%as.matrix(c(1,k0_gen,k1_gen,k2_gen))
  T =PP[1];
  U0=PP[2];
  U1=PP[3];
  U2=PP[4];
  V11=PP[5];
  V22=PP[6];
  V01=PP[7];
  V02=PP[8];
  V12=PP[9];
  MSEgen_ = lb*U0^2*C0_2+2*lb*(T*V11+0.5*U1^2)*C1_2+2*lb*(T*V22+0.5*U2^2)*C2_2+T^2+
    2*lb*(T*V12+U1*U2)*C12+2*lb*(T*V01+U1*U0)*C01+2*lb*(T*V02+U0*U2)*C02+
    lb^2*(V11*C1_2+V22*C2_2+V01*C01+V02*C02+V12*C12)^2;
  
  x0=as.matrix(NN[,1])
  xK=NN[,-1]
  MSEgen_formula  = c(t(x0)%*%Q%*%x0-t(x0)%*%Q%*%xK %*%solve(t(xK)%*%Q%*%xK)%*% t(xK)%*%Q%*%x0);
  
  Kopt = k012_gen;
  MSEgen_formula2 = t(x0)%*%Q%*%x0 + 2*t(x0)%*%Q%*%xK%*%Kopt + t(Kopt )%*%(t(xK)%*%Q%*%xK)%*%Kopt;
  
  return(list(k0_gen=k0_gen,k1_gen=k1_gen,k2_gen=k2_gen,MSE1=MSEgen_,MSE2=MSE_gen_est));
}

Y_RPR_RC4_L_GEN <- function(a1,b1,a2,b2,Ys,C0_2,C1_2,C2_2,C01,C02,C12,lb) {
  NN = matrix(c(-Ys,Ys,0,0,
                0,Ys,0,0,
                0,0,1,0,
                0,0,0,1,
                0,0,0,0,
                0,0,0,0,
                0,0,0,0,
                0,0,0,0,
                0,0,a2,a1),nrow=9,ncol=4,byrow = TRUE)
  
  MSE_est  = MSE_2vars_3param_GEN(NN,Ys,C0_2,C1_2,C2_2,C01,C02,C12,lb);
  k0_gen = MSE_est$k0_gen
  k1_gen = MSE_est$k1_gen
  k2_gen = MSE_est$k2_gen
  MSEgen1=MSE_est$MSE1;
  MSEgen2=MSE_est$MSE2;
  BIAS_ = -1;
  return(list(MSE=MSEgen1,BIAS=BIAS_,k0=k0_gen,k1=k1_gen,k2=k2_gen));
}

MSE_2vars_2param_GEN <- function(NN,Ys,C0_2,C1_2,C2_2,C01,C02,C12,lb) {
  Q = getQ(Ys,C0_2,C1_2,C2_2,C01,C02,C12,lb);
  
  MMgen = t(NN[,2:3])%*%Q%*%(NN[,2:3])
  Ggen  = -t(NN[,1])%*%Q%*%(NN[,2:3])
  
  k012_gen = solve(MMgen)%*%t(Ggen);
  k0_gen=k012_gen[1]
  k1_gen=k012_gen[2]
  MSE_gen_est  = c(t(NN%*%rbind(1,k012_gen))%*%Q%*%(NN%*%rbind(1,k012_gen)));
  
  PP=NN%*%as.matrix(c(1,k0_gen,k1_gen))
  T =PP[1];
  U0=PP[2];
  U1=PP[3];
  U2=PP[4];
  V11=PP[5];
  V22=PP[6];
  V01=PP[7];
  V02=PP[8];
  V12=PP[9];
  MSEgen_ = lb*U0^2*C0_2+2*lb*(T*V11+0.5*U1^2)*C1_2+2*lb*(T*V22+0.5*U2^2)*C2_2+T^2+
    2*lb*(T*V12+U1*U2)*C12+2*lb*(T*V01+U1*U0)*C01+2*lb*(T*V02+U0*U2)*C02+
    lb^2*(V11*C1_2+V22*C2_2+V01*C01+V02*C02+V12*C12)^2;
  
  return(list(k0_gen=k0_gen,k1_gen=k1_gen,MSE1=MSEgen_,MSE2=MSE_gen_est));
}

Y_RPR_RC_L_GEN <- function(a,b,Ys,C0_2,C1_2,C2_2,C01,C02,C12,lb) {
  NN = matrix(c(-Ys,Ys,0,
                0,Ys,0,
                0,0,1,
                0,a*Ys,0,
                0,0,0,
                0,b*Ys,0,
                0,0,0,
                0,a*Ys,0,
                0,0,a),nrow=9,ncol=3,byrow = TRUE)
  
  MSE_est  = MSE_2vars_2param_GEN(NN,Ys,C0_2,C1_2,C2_2,C01,C02,C12,lb);
  k0_gen = MSE_est$k0_gen
  k1_gen = MSE_est$k1_gen
  k2_gen = MSE_est$k2_gen
  MSEgen1=MSE_est$MSE1;
  MSEgen2=MSE_est$MSE2;
  BIAS_ = -1;
  return(list(MSE=MSEgen1,BIAS=BIAS_,k0=k0_gen,k1=k1_gen,k2=k2_gen));
}

MSESALL2 <- function (Y_PRL_FUN, as1, bs1, as2, bs2, name_) {
  T=length(as);
  mses1_L_=numeric(T);
  for (t in 1:T) {
    resu_L_ = Y_PRL_FUN(as1[t],bs1[t],as2[t],bs2[t],Y_bar,Cy_2,Cx1_2,Cx2_2,Cyx1,Cyx2,Cx1x2,lb);
    mses1_L_[t] = resu_L_$MSE;
  }
  tmin = which.min(mses1_L_);
  MSE_min1_L_ = mses1_L_[tmin]
  cat("---------------------------------");
  cat(name_,"\n");
  cat("V_YO/MSE_min1_L_=",V_YO/MSE_min1_L_,"\n");
  aopt_L_=as[tmin];
  cat("aopt_L_=",as[tmin],"  ( b=",b,")\n");
  cat("---------------------------------");
  rl = list(mses1_L_=mses1_L_,as=as,bs=bs,aopt_L_=aopt_L_,name=name_);
  return(rl);
}

MSESALL <- function (Y_PRL_FUN, as, bs, name_) {
  T=length(as);
  mses1_L_=numeric(T);
  for (t in 1:T) {
    resu_L_ = Y_PRL_FUN(as[t],bs[t],Y_bar,Cy_2,Cx1_2,Cx2_2,Cyx1,Cyx2,Cx1x2,lb);
    mses1_L_[t] = resu_L_$MSE;
  }
  tmin = which.min(mses1_L_);
  MSE_min1_L_ = mses1_L_[tmin]
  cat("---------------------------------");
  cat(name_,"\n");
  cat("V_YO/MSE_min1_L_=",V_YO/MSE_min1_L_,"\n");
  aopt_L_=as[tmin];
  cat("aopt_L_=",as[tmin],"  ( b=",b,")\n");
  cat("---------------------------------");
  rl = list(mses1_L_=mses1_L_,as=as,bs=bs,aopt_L_=aopt_L_,name=name_);
  return(rl);
}
