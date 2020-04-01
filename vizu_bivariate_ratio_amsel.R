rm(list = ls())
library(latex2exp)
#
directory_="~/OUTPUTRATIO/";
source("~/CODES_csmt2019/vizu_lib_ratio.R");

nb_l  = 50; # number of values for a, a1, a2 for drawing the curves
            # to be reduced eventually in the case of simulations

datas     = matrix(0,nrow=6,ncol=11);
## ----------------------------------------------------------------
dd_ = 0;
for (data_ in c("D1","D2","D3","D4","D6","D8")) { 
  names_res = NULL;
  dd_ = dd_+1;
  
  cat(("###################################################"))
  
  as_left  = seq(-1.2,-0.05, length.out = nb_l);
  as_right = seq(0.05,+1.2,  length.out = nb_l);

  as = c(as_left, as_right);
  as00=as;
  
  source("~/CODES_csmt2019/getdata_art.R"); 
  
  datas[dd_,]=c(N,n,Y_bar,X1_bar,X2_bar,Cy_2,Cx1_2,Cx2_2,rhoyx1,rhoyx2,rhox1x2)
  
  ## ----------------------------------------------------------------
  
  alpha=0;
  a=-1/2;
  #b=3/8-alpha/4;
  b=0.5;
  b00=b;
  
  V_YO        = lb*Y_bar^2*Cy_2
  
  ## ----------------------------------------------------------------
  
  resu2_=Y_RPR03(a,b,Y_bar,X1_bar,X2_bar,Cy_2,Cx1_2,Cx2_2,Cx1x2,lb);
  BIAS_CRR_ab2 = resu2_$BIAS;
  MSE_CRR_ab2 = resu2_$MSE;
  cat("V_YO/MSE_CRR_ab2=",V_YO/MSE_CRR_ab2,"\n");
  cat("BIAS_CRR_ab2/Y=",BIAS_CRR_ab2/Y_bar,"\n");
  cat("BIAS_CRR_ab2^2/MSE_CRR_ab2=",BIAS_CRR_ab2^2/MSE_CRR_ab2,"\n");
  # 
  resuRC_L_GEN=Y_RPR_RC_L_GEN(a,b,Y_bar,Cy_2,Cx1_2,Cx2_2,Cyx1,Cyx2,Cx1x2,lb);
  MSE_CRR_ab2_L_GEN = resuRC_L_GEN$MSE;
  cat("V_YO/MSE_CRR_ab2_L_GEN=",V_YO/MSE_CRR_ab2_L_GEN," (a=",a," b=",b,"\n");
  # 
  resu4_L_GEN = Y_RPR_RC4_L_GEN (a,b,a,b,Y_bar,Cy_2,Cx1_2,Cx2_2,Cyx1,Cyx2,Cx1x2,lb);
  MSE_yRc4_L_GEN = resu4_L_GEN$MSE;
  cat("V_YO/MSE_yRc4_L_GEN=",V_YO/MSE_yRc4_L_GEN," (a=",a," b=",b,"\n");
  
  ## ----------------------------------------------------------------

  bs=rep(b,2*nb_l)
  resu_test = MSESALL(Y_RPR_RC, as, bs, "yRc");
  mses1 =   resu_test$mses1_L_;
  
  bs=rep(b,2*nb_l)
  resu_test = MSESALL(Y_RPR_RC_L, as, bs, "yRc_L");
  mses1_L =   resu_test$mses1_L_;
  
  gs_left = seq( (-1-(min(as_left)))/(min(as_left))*X2_bar,(-1-(max(as_left)))/(max(as_left))*X2_bar,length.out = 5000);
  gs_right = seq( (-1-(min(as_right)))/(min(as_right))*X2_bar,(-1-(max(as_right)))/(max(as_right))*X2_bar,length.out = 5000);
  gs = c(gs_left, gs_right);
  asg = -X2_bar/(X2_bar+gs);
  bsg = ( (X2_bar)/(X2_bar+gs) )^2;
  resu_test = MSESALL(Y_RPR_RC_L, asg, bsg, "yRc_L_gamma");
  mses3_L = resu_test$mses1_L_;
  
  ## ----------------------------------------------------------------
  
  bs=rep(b,2*nb_l)
  resu_test = MSESALL2(Y_RPR_RC4_L_GEN, as, bs, as, bs, "yRc4_GEN");
  mses_Rc4_L_GEN =   resu_test$mses1_L_;
  
  ## ----------------------------------------------------------------
  
  MSE_yreg = lb*Sy_2*(1-rhoyx1^2-rhoyx2^2+2*rhoyx1*rhoyx2*rhox1x2);
  cat("V_YO/MSE_yreg=",V_YO/MSE_yreg,"\n");
  
  MSE_ydif = lb*Sy_2*(1-(rhoyx1^2+rhoyx2^2-2*rhoyx1*rhoyx2*rhox1x2)/(1-rhox1x2^2));
  cat("V_YO/MSE_ydif=",V_YO/MSE_ydif,"\n");   
  
  ## ----------------------------------------------------------------
  
  L = 1-rhox1x2^2-(rhoyx1^2+rhoyx2^2-2*rhoyx1*rhoyx2*rhox1x2);
  A = 1-rhox1x2^2
  dt = (1-f)/n
  
  MSE_yRE = dt*Sy_2*L/(A+dt*Cy_2*L);
  cat("V_YO/MSE_yRE=",V_YO/MSE_yRE,"\n");   
  
  MSE_ypr = 0.25*dt*Y_bar^2*(4*Cy_2*L-dt*A*(Cx1_2+Cx2_2)^2)/(A+dt*Cy_2*L+dt*A*(Cx1_2+Cx2_2));
  cat("V_YO/MSE_ypr=",V_YO/MSE_ypr,"\n");   
  
  ## ----------------------------------------------------------------
  
  resu_add12_ = 
    Y_MSE_MIN_SUMf1f2(-0.5,1,-0.5,1,Y_bar,X1_bar,X2_bar,Cy_2,Cx1_2,Cx2_2,Cx1x2,lb)
  MSE_add12 =resu_add12_$MSE;
  BIAS_add12 =resu_add12_$BIAS;
  cat("V_YO/MSE_add12=",V_YO/MSE_add12,"\n");
  
  resu_extolk_ = 
    Y_MSE_MIN_SUMf1f2_extendolkin(-0.5,1,-0.5,1,
                                  Y_bar,X1_bar,X2_bar,Cy_2,Cx1_2,Cx2_2,Cx1x2,lb);
  MSE_extolk   = resu_extolk_$MSE;
  BIAS_extolk  = resu_extolk_$BIAS;
  cat("V_YO/MSE_extolk=",V_YO/MSE_extolk,"\n");
  
  ## ----------------------------------------------------------------
  
  adeb=0;afin=nb_l*2;
  if (data_=="D5") {
    ymin_=-0.25;
    ymax_=8;
   } else {
     ymin_=V_YO/max(mses1,mses1_L,mses_Rc4_L_GEN); #
     ymax_=V_YO/min(mses1,mses1_L,mses_Rc4_L_GEN); #
   }
  namefile_=paste(directory_,"PRE_Rc_RcL__",data_,".png",sep="");
  png(namefile_);
  par(lwd=3)
  plot(as[adeb:afin],V_YO/mses1_L[adeb:afin],type='l',xlab="a", ylab= " ",
       cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,lty=0,ylim = c(ymin_, ymax_),col="white")
  points(as[adeb:afin],V_YO/mses1_L[adeb:afin],type='l',xlab="a", ylab= " ",
         cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,lty=1,col="black")
  points(as[adeb:afin],V_YO/mses1[adeb:afin],type='l',xlab="a", ylab= " ",
         cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,lty=3,col="black")
  points(asg[adeb:afin],V_YO/mses3_L[adeb:afin],type='l',xlab="a", ylab= " ",
         cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,lty=4,col="black")
  points(as[adeb:afin],V_YO/mses_Rc4_L_GEN[adeb:afin],type='l',xlab="a", ylab= " ",
         cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,lty=5,col="black")
  abline(v=-0.5,col="red");
  abline(h = V_YO/MSE_add12,col="blue",lty=8);
  abline(h = V_YO/MSE_yRE,col="blue",lty=8);
  
  if (data_=="D6")
    legend("topright",bg="white",legend=c(TeX('$\\bar{y}_{1}$'),TeX('$\\bar{y}_{3}$'),TeX('$\\bar{y}_{5}$'),
                                TeX('$\\bar{y}_{5;L}$'),TeX('$\\bar{y}_{6;L}$'),TeX('$\\bar{y}_{10;L}$')), 
         col=c("blue","blue","black","black","black","black"),lty=c(8,8,3,1,4,5), 
         ncol=1, cex=1.8);
  dev.off();

  ## ----------------------------------------------------------------
}
