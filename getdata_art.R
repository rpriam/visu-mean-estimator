if (data_=="D1") {
  load("~/CODES_csmt2019/D1art.RData")
} else if (data_=="D2") {
  load("~/CODES_csmt2019/D2art.RData")
} else if (data_=="D3") {
  load("~/CODES_csmt2019/D3art.RData")
} else if (data_=="D4") {
  load("~/CODES_csmt2019/D4art.RData")
} else if (data_=="D5") {
  load("~/CODES_csmt2019/D5art.RData")
} else if (data_=="D6") {
  load("~/CODES_csmt2019/D6art.RData")
} else if (data_=="D7") {
  load("~/CODES_csmt2019/D7art.RData")
} else if (data_=="D8") {
  load("~/CODES_csmt2019/D8art.RData")
}
##
Cy_2  = Sy_2/Y_bar^2
Cx1_2 = Sx1_2 / X1_bar^2
Cx2_2 = Sx2_2 / X2_bar^2
Cyx1  = Sx1y/Y_bar/X1_bar;
Cyx2  = Sx2y/Y_bar/X2_bar;
Cx1x2 = Sx1x2/X1_bar/X2_bar;
rhoyx1 = Sx1y / sqrt(Sx1_2) / sqrt(Sy_2)
rhoyx2 = Sx2y / sqrt(Sx2_2) / sqrt(Sy_2)
rhox1x2 = Sx1x2 / sqrt(Sx1_2) / sqrt(Sx2_2)
f  = n/N;
lb = (N-n)/n/N;
cat("data=",data_,"\n");
cat("N =",N,"|","n =",n,"|","f =",f,"|","Y_bar =",Y_bar,"|","X1_bar =",
    X1_bar,"|","X2_bar =",X2_bar,"\n");
cat("Cy_2 =",Cy_2,"|","Cx1_2 =",Cx1_2,"|","Cx2_2 =",Cx2_2,"Cx1x2 =",Cx1x2,"\n");
cat("rhoyx1 =",rhoyx1,"|","rhoyx2 =",rhoyx2,"|","rhox1x2 =",rhox1x2)

