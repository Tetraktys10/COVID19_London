function dydt = seir_eq_v5(t,y,beta1,beta2,beta3,c,k,p,gammaa,gammas,gammah,mus,muh1,muh2,muh3,psi1,psi2,psi3,alpha,q,xi,d1,d2,m1,m2,m3)  %(S E Iu Ia Is Ih Ru Rh D Nc)
    
    if t<4
        d=d1;
    else
        d=d2;
    end

   if t<9
       beta=beta1;
       test=0;
       a1=1; %facemasks
       a2=1; %face covering
       b1=0; %percentage of early symptomatic who go in isolation
   elseif t>=9&&t<15
      beta=beta2*beta1;%
      test=0;
      a1=1;
      a2=1;
      b1=0;
   elseif t>=15&&t<61
       beta=beta3*beta1;%
       %beta=beta1;
       test=0;
       a1=1;
       a2=1;
       b1=0;
   elseif t>=61&&t<720
%      beta=beta1;
       beta=beta3*beta1;
      test=0;
       a1=1; 
       a2=1;
      b1=0;


%       test=(1/(7/2))*0.8;  %%3.5  2.3  1.75  
%       a1=0.7;%1;%
%       a2=0.9;%1;%
%       b1=1*1/7;
   else
       %beta=beta1;
       beta=beta3*beta1;
       test=0;
       a1=1;
       a2=1;
       b1=0;
   end
   
   if t<9
       f=1.343e-05*exp(-0.09117*t) -1.168e-05*exp(-0.2641*t);
       g=7.056e-06*exp(-0.1528*t) -6.188e-06*exp(-0.7754*t);
   elseif t>=9&&t<53
       f=0;
       g=0;
   else
       f=0;
       g=0;
%        f=1.343e-05*exp(-0.09117*t) -1.168e-05*exp(-0.2641*t);
%        g=7.056e-06*exp(-0.1528*t) -6.188e-06*exp(-0.7754*t);
   end
    
    %0-14
    dydt1 = - c(1,1)*beta*y(1)*(a2*y(3)+alpha*a2*y(4)+a1*xi*y(5)+a1*xi*alpha*y(31))-c(2,1)*beta*y(1)*(a2*y(13)+0.92*934000*g+alpha*a2*(y(14)+0.92*934000*f)+a1*xi*y(15)+a1*xi*alpha*y(32))-c(3,1)*beta*y(1)*(a2*y(23)+0.08*934000*g+alpha*a2*(y(24)+0.08*934000*f)+a1*xi*y(25)+a1*xi*alpha*y(33));%
    dydt2 = c(1,1)*beta*y(1)*(a2*y(3)+alpha*a2*y(4)+a1*xi*y(5)+a1*xi*alpha*y(31))+c(2,1)*beta*y(1)*(a2*y(13)+0.92*934000*g+alpha*a2*(y(14)+0.92*934000*f)+a1*xi*y(15)+a1*xi*alpha*y(32))+c(3,1)*beta*y(1)*(a2*y(23)+0.08*934000*g+alpha*a2*(y(24)+0.08*934000*f)+a1*xi*y(25)+a1*xi*alpha*y(33))-k*y(2)-b1*y(2);%
    dydt3 = p*k*y(2)-q*y(3);%Iu
    dydt4 = (1-p)*k*y(2)-gammaa*y(4)-test*y(4);%Ia
    dydt5 = q*y(3)+p*k/2*y(34)-(psi1+gammas+mus)*y(5);%Is
    dydt6 = psi1*y(5)-(gammah+muh1)*y(6)-m1*y(6);%Ih
    dydt7 = gammaa*y(4)+gammas*y(5)+1/3*y(31);%Ru
    dydt8 = gammah*y(6);%Rh
    dydt9 = mus*y(5)+muh1*y(6)+m1*y(6);%D -y(9)
    dydt10 = d*y(5)+y(6)+y(8);%+y(9)
    
    %15-59
    dydt11 = - c(1,2)*beta*y(11)*(a2*y(3)+alpha*a2*y(4)+a1*xi*y(5)+a1*xi*alpha*y(31))-c(2,2)*beta*y(11)*(a2*y(13)+0.92*934000*g+alpha*a2*(y(14)+0.92*934000*f)+a1*xi*y(15)+a1*xi*alpha*y(32))-c(3,2)*beta*y(11)*(a2*y(23)+0.08*934000*g+alpha*a2*(y(24)+0.08*934000*f)+a1*xi*y(25)+a1*xi*alpha*y(33));%
    dydt12 = c(1,2)*beta*y(11)*(a2*y(3)+alpha*a2*y(4)+a1*xi*y(5)+a1*xi*alpha*y(31))+c(2,2)*beta*y(11)*(a2*y(13)+0.92*934000*g+alpha*a2*(y(14)+0.92*934000*f)+a1*xi*y(15)+a1*xi*alpha*y(32))+c(3,2)*beta*y(11)*(a2*y(23)+0.08*934000*g+alpha*a2*(y(24)+0.08*934000*f)+a1*xi*y(25)+a1*xi*alpha*y(33))-k*y(12)-b1*y(12);%
    dydt13 = p*k*y(12)-q*y(13);%Iu
    dydt14 = (1-p)*k*y(12)-gammaa*y(14)-test*y(14);%Ia
    dydt15 = q*y(13)+p*k/2*y(35)-(psi2+gammas+mus)*y(15);%Is
    dydt16 = psi2*y(15)-(gammah+muh2)*y(16)-m2*y(16);%Ih
    dydt17 = gammaa*y(14)+gammas*y(15)+1/3*y(32);%Ru
    dydt18 = gammah*y(16);%Rh
    dydt19 = mus*y(15)+muh2*y(16)+m2*y(16);%D -y(19)
    dydt20 = d*y(15)+y(16)+y(18);%+y(19)
    
    %60+
    dydt21 = - c(1,3)*beta*y(21)*(a2*y(3)+alpha*a2*y(4)+a1*xi*y(5)+a1*xi*alpha*y(31))-c(2,3)*beta*y(21)*(a2*y(13)+0.92*934000*g+alpha*a2*(y(14)+0.92*934000*f)+a1*xi*y(15)+a1*xi*alpha*y(32))-c(3,3)*beta*y(21)*(a2*y(23)+0.08*934000*g+alpha*a2*(y(24)+0.08*934000*f)+a1*xi*y(25)+a1*xi*alpha*y(33));%
    dydt22 = c(1,3)*beta*y(21)*(a2*y(3)+alpha*a2*y(4)+a1*xi*y(5)+a1*xi*alpha*y(31))+c(2,3)*beta*y(21)*(a2*y(13)+0.92*934000*g+alpha*a2*(y(14)+0.92*934000*f)+a1*xi*y(15)+a1*xi*alpha*y(32))+c(3,3)*beta*y(21)*(a2*y(23)+0.08*934000*g+alpha*a2*(y(24)+0.08*934000*f)+a1*xi*y(25)+a1*xi*alpha*y(33))-k*y(22)-b1*y(22);%
    dydt23 = p*k*y(22)-q*y(23);%Iu
    dydt24 = (1-p)*k*y(22)-gammaa*y(24)-test*y(24);%Ia 3.5
    dydt25 = q*y(23)+p*k/2*y(36)-(psi3+gammas+mus)*y(25);%Is
    dydt26 = psi3*y(25)-(gammah+muh3+m3)*y(26);%Ih
    dydt27 = gammaa*y(24)+gammas*y(25)+1/3*y(33);%Ru
    dydt28 = gammah*y(26);%Rh
    dydt29 = mus*y(25)+muh3*y(26)+m3*y(26);%D -y(29)
    dydt30 = d*y(25)+y(26)+y(28);%+y(29)
    
    %asymptomatic isolated
    dydt31= test*y(4)+(1-p)*k/2*y(34)-1/3*y(31);
    dydt32=test*y(14)+(1-p)*k/2*y(35)-1/3*y(32);
    dydt33=test*y(24)+(1-p)*k/2*y(36)-1/3*y(33);
    
    %contacts traced
    dydt34=b1*y(2)-k/2*y(34); 
    dydt35=b1*y(12)-k/2*y(35);
    dydt36=b1*y(22)-k/2*y(36);
    
    
    dydt = [dydt1; dydt2; dydt3; dydt4; dydt5; dydt6; dydt7; dydt8; dydt9; dydt10; dydt11; dydt12; dydt13; dydt14; dydt15; dydt16; dydt17; dydt18; dydt19; dydt20; dydt21; dydt22; dydt23; dydt24; dydt25; dydt26; dydt27; dydt28; dydt29; dydt30; dydt31; dydt32; dydt33; dydt34; dydt35; dydt36];
    
end

