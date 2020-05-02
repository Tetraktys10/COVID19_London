function dydt = seir_eq_v3(t,y,beta1,beta2,beta3,c,k,p,gammaa,gammas,gammah,mus,muh1,muh2,muh3,psi1,psi2,psi3,alpha,q,xi,d1,d2,m1,m2,m3,d1s,d2s,ps)  %(S E Iu Ia Is Ih Ru Rh D Nc)
    
    if t<4
        d=d1;
        ds=d1s;
    else
        d=d2;
        ds=d2s;
    end

   if t<9
       beta=beta1;
       a1=1;
   elseif t>=9&&t<15
      beta=beta2*beta1;%
      a1=1;
   elseif t>=15&&t<61
       beta=beta3*beta1;%
       a1=1;
   elseif t>=61&&t<720
       beta=beta3*beta1;%
       %beta=beta1;
       %beta=beta2*beta1;
       a1=1;
   else
       beta=beta3*beta1;%
       a1=1;
   end

   
   
   if t<9
       f=1.343e-05*exp(-0.09117*t) -1.168e-05*exp(-0.2641*t);
       g=7.056e-06*exp(-0.1528*t) -6.188e-06*exp(-0.7754*t);
   elseif t>=9&&t<61
       f=0;
       g=0;
   else
       f=1.343e-05*exp(-0.09117*t) -1.168e-05*exp(-0.2641*t);%0;%
       g=7.056e-06*exp(-0.1528*t) -6.188e-06*exp(-0.7754*t);%0;%
   end
    
    %0-14
    dydt1 = - c(1,1)*beta*y(1)*(y(3)+alpha*y(4)+a1*xi*y(5))-c(2,1)*beta*y(1)*(y(13)+0.92*934000*g+alpha*(y(14)+0.92*934000*f)+a1*xi*y(15))-c(3,1)*beta*y(1)*(y(23)+0.08*934000*g+alpha*(y(24)+0.08*934000*f)+a1*xi*y(25));%
    dydt2 = c(1,1)*beta*y(1)*(y(3)+alpha*y(4)+a1*xi*y(5))+c(2,1)*beta*y(1)*(y(13)+0.92*934000*g+alpha*(y(14)+0.92*934000*f)+a1*xi*y(15))+c(3,1)*beta*y(1)*(y(23)+alpha*y(24)+a1*xi*y(25))-k*y(2);%
    dydt3 = p*k*y(2)-q*y(3);%Iu
    dydt4 = (1-p)*k*y(2)-gammaa*y(4);%Ia
    dydt5 = q*y(3)-(psi1+gammas+mus)*y(5);%Is
    dydt6 = psi1*y(5)-(gammah+muh1)*y(6)-m1*y(6);%Ih
    dydt7 = gammaa*y(4)+gammas*y(5);%Ru
    dydt8 = gammah*y(6);%Rh
    dydt9 = mus*y(5)+muh1*y(6)+m1*y(6);%D -y(9)
    dydt10 = d*y(5)+psi1*y(5);%NC
    
    %15-59
    dydt11 = - c(1,2)*beta*y(11)*(y(3)+alpha*y(4)+a1*xi*y(5))-c(2,2)*beta*y(11)*(y(13)+0.92*934000*g+alpha*(y(14)+0.92*934000*f)+a1*xi*y(15))-c(3,2)*beta*y(11)*(y(23)+0.08*934000*g+alpha*(y(24)+0.08*934000*f)+a1*xi*y(25));%
    dydt12 = c(1,2)*beta*y(11)*(y(3)+alpha*y(4)+a1*xi*y(5))+c(2,2)*beta*y(11)*(y(13)+0.92*934000*g+alpha*(y(14)+0.92*934000*f)+a1*xi*y(15))+c(3,2)*beta*y(11)*(y(23)+0.08*934000*g+alpha*(y(24)+0.08*934000*f)+a1*xi*y(25))-k*y(12);%
    dydt13 = p*k*y(12)-q*y(13);%Iu
    dydt14 = (1-p)*k*y(12)-gammaa*y(14);%Ia
    dydt15 = q*y(13)-(psi2+gammas+mus)*y(15);%Is
    dydt16 = psi2*y(15)-(gammah+muh2)*y(16)-m2*y(16);%Ih
    dydt17 = gammaa*y(14)+gammas*y(15);%Ru
    dydt18 = gammah*y(16);%Rh
    dydt19 = mus*y(15)+muh2*y(16)+m2*y(16);%D -y(19)
    dydt20 = d*y(15)+psi2*y(15);%NC
    
    %60+
    dydt21 = - c(1,3)*beta*y(21)*(y(3)+alpha*y(4)+a1*xi*y(5))-c(2,3)*beta*y(21)*(y(13)+0.92*934000*g+alpha*(y(14)+0.92*934000*f)+a1*xi*y(15))-c(3,3)*beta*y(21)*(y(23)+0.08*934000*g+alpha*(y(24)+0.08*934000*f)+a1*xi*y(25));%
    dydt22 = c(1,3)*beta*y(21)*(y(3)+alpha*y(4)+a1*xi*y(5))+c(2,3)*beta*y(21)*(y(13)+0.92*934000*g+alpha*(y(14)+0.92*934000*f)+a1*xi*y(15))+c(3,3)*beta*y(21)*(y(23)+0.08*934000*g+alpha*(y(24)+0.08*934000*f)+a1*xi*y(25))-k*y(22);%
    dydt23 = ps*k*y(22)-q*y(23);%Iu
    dydt24 = (1-ps)*k*y(22)-gammaa*y(24);%Ia
    dydt25 = q*y(23)-(psi3+gammas+mus)*y(25);%Is
    dydt26 = psi3*y(25)-(gammah+muh3+m3)*y(26);%Ih 
    dydt27 = gammaa*y(24)+gammas*y(25);%Ru
    dydt28 = gammah*y(26);%Rh
    dydt29 = mus*y(25)+muh3*y(26)+m3*y(26);%D -y(29)
    dydt30 = ds*y(25)+psi3*y(25);%NC
    
    
    dydt = [dydt1; dydt2; dydt3; dydt4; dydt5; dydt6; dydt7; dydt8; dydt9; dydt10; dydt11; dydt12; dydt13; dydt14; dydt15; dydt16; dydt17; dydt18; dydt19; dydt20; dydt21; dydt22; dydt23; dydt24; dydt25; dydt26; dydt27; dydt28; dydt29; dydt30];
    
end

