%% Next Generation Matrix method   E,Ia,Is,Ii,Ih,Iai Ei

syms beta S E Iu Ia Is Ih S1 E1 Iu1 Ia1 Is1 Ih1 S2 E2 Iu2 Ia2 Is2 Ih2 S3 E3 Iu3 Ia3 Is3 Ih3 Iai1 Iai2 Iai3 Ei1 Ei2 Ei3 c11 c12 c13 c21 c22 c23 c31 c32 c33 alpha xi k p ps gammaa q psi1 psi2 psi3 gammas mus gammah muh1 muh2 muh3 m1 m2 m3 f g N1 N2 N3 a1 a2 b1 test % define x and y as symbolic https://uk.mathworks.com/help/symbolic/syms.html

eqns = [- c11*beta*S1*(Ia1+alpha*Iu1+xi*Is1)-c21*beta*S1*(Ia2+0.92*934000*g+alpha*(Iu2+0.92*934000*f)+xi*Is2)-c31*beta*S1*(Ia3+0.08*934000*g+alpha*(Iu3+0.08*934000*f)+xi*Is3)==0, c11*beta*S1*(Ia1+alpha*Iu1+xi*Is1)+c21*beta*S1*(Ia2+0.92*934000*g+alpha*(Iu2+0.92*934000*f)+xi*Is2)+c31*beta*S1*(Ia3+alpha*Iu3+xi*Is3)-k*E1==0, p*k*E1-q*Ia1==0, (1-p)*k*E1-gammaa*Iu1==0, q*Ia1-(psi1+gammas+mus)*Is1==0, psi1*Is1-(gammah+muh1)*Ih1-m1*Ih1==0, - c12*beta*S2*(Ia1+alpha*Iu1+xi*Is1)-c22*beta*S2*(Ia2+0.92*934000*g+alpha*(Iu2+0.92*934000*f)+xi*Is2)-c32*beta*S2*(Ia3+0.08*934000*g+alpha*(Iu3+0.08*934000*f)+xi*Is3)==0, c12*beta*S2*(Ia1+alpha*Iu1+xi*Is1)+c22*beta*S2*(Ia2+0.92*934000*g+alpha*(Iu2+0.92*934000*f)+xi*Is2)+c32*beta*S2*(Ia3+0.08*934000*g+alpha*(Iu3+0.08*934000*f)+xi*Is3)-k*E2==0, p*k*E2-q*Ia2==0, (1-p)*k*E2-gammaa*Iu2==0, q*Ia2-(psi2+gammas+mus)*Is2==0, psi2*Is2-(gammah+muh2)*Ih2-m2*Ih2==0, - c13*beta*S3*(Ia1+alpha*Iu1+xi*Is1)-c23*beta*S3*(Ia2+0.92*934000*g+alpha*(Iu2+0.92*934000*f)+xi*Is2)-c33*beta*S3*(Ia3+0.08*934000*g+alpha*(Iu3+0.08*934000*f)+xi*Is3)==0, c13*beta*S3*(Ia1+alpha*Iu1+xi*Is1)+c23*beta*S3*(Ia2+0.92*934000*g+alpha*(Iu2+0.92*934000*f)+xi*Is2)+c33*beta*S3*(Ia3+0.08*934000*g+alpha*(Iu3+0.08*934000*f)+xi*Is3)-k*E3==0, ps*k*E3-q*Ia3==0, (1-ps)*k*E3-gammaa*Iu3==0, q*Ia3-(psi3+gammas+mus)*Is3==0, psi3*Is3-(gammah+muh3+m3)*Ih3==0];
Sol=solve(eqns, [S1 E1 Ia1 Iu1 Is1 Ih1 S2 E2 Ia2 Iu2 Is2 Ih2 S3 E3 Ia3 Iu3 Is3 Ih3 Iai1 Iai2 Iai3 Ei1 Ei2 Ei3]); %y(3)=Iu, y(4)=Ia, y(5)=Is
%Sol.S=N;

%%
%matrix of derivatives of infected equations respect to infected
%compartments %E Ia Iu Is Ih Iai Ei

%c11*beta*S1*(Ia1+alpha*Iu1+xi*Is1)-c21*beta*S1*(Ia2+0.92*934000*g+alpha*(Iu2+0.92*934000*f)+xi*Is2)-c31*beta*S1*(Ia3+0.08*934000*g+alpha*(Iu3+0.08*934000*f)+xi*Is3)



V=[-k-b1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    (1-p)*k -gammaa-test 0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0 0 0 0
    p*k 0 -q 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0 0 0 0
    0 0 q -(psi1+gammas+mus) 0 0 0 0 0 0  0 0 0 0 0 0 0 0 0 0 0
    0 0 0 psi1 -(gammah+muh1+m1) 0 0 0 0 0  0 0 0 0 0 0 0 0 0 0 0
     0 0 0 0 0 -k-b1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
     0 0 0 0 0 (1-p)*k -gammaa-test 0 0 0 0 0 0 0 0 0 0 0 0 0 0
     0 0 0 0 0 p*k 0 -q 0 0 0 0 0 0 0 0 0 0 0 0 0
     0 0 0 0 0 0 0 q -(psi2+gammas+mus) 0 0 0 0 0 0 0 0 0 0 0 0
     0 0 0 0 0 0 0 0 psi2 -(gammah+muh2+m2) 0 0 0 0 0 0 0 0 0 0 0
     0 0 0 0 0  0 0 0 0 0 -k-b1 0 0 0 0 0 0 0 0 0 0
     0 0 0 0 0  0 0 0 0 0 (1-ps)*k -gammaa-test 0 0 0 0 0 0 0 0 0
     0 0 0 0 0  0 0 0 0 0 ps*k 0 -q 0 0 0 0 0 0 0 0
     0 0 0 0 0  0 0 0 0 0 0 0 q -(psi3+gammas+mus) 0 0 0 0 0 0 0
     0 0 0 0 0  0 0 0 0 0 0 0 0 psi3 -(gammah+muh3+m3) 0 0 0 0 0 0
     0 test 0 0 0  0 0 0 0 0  0 0 0 0 0  -1/3 0 0 (1-p)*k/2 0 0
     0 0 0 0 0  0 test 0 0 0  0 0 0 0 0  0 -1/3 0 0 (1-p)*k/2 0
     0 0 0 0 0  0 0 0 0 0  0 test 0 0 0 0 0 -1/3 0 0 (1-p)*k/2
     b1 0 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0 -k/2 0 0
     0 0 0 0 0  b1 0 0 0 0  0 0 0 0 0  0 0 0  0 -k/2 0
     0 0 0 0 0  0 0 0 0 0  b1 0 0 0 0  0 0 0 0 0 -k/2]; %matrix of disease progression components

%evaluate F and V in the disease free equilibrium x0=(N,0,0,0,0,0)


F=[0 c11*beta*N1*alpha*a2 c11*beta*N1*a2 c11*beta*N1*xi*a1 0 0 c21*beta*N1*alpha*a2 c21*beta*N1*a2 c21*beta*N1*xi*a1 0 0 c31*beta*N1*alpha*a2 c31*beta*N1*a2 c31*beta*N1*xi*a1 0 c11*beta*N1*alpha*a1*xi c21*beta*N1*alpha*a1*xi c31*beta*N1*alpha*a1*xi 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 c12*beta*N2*alpha*a2 c12*beta*N2*a2 c12*beta*N2*xi*a1 0 0 c22*beta*N2*alpha*a2 c22*beta*N2*a2 c22*beta*N2*xi*a1 0 0 c32*beta*N2*alpha*a2 c32*beta*N2*a2 c32*beta*N2*xi*a1 0 c12*beta*N2*alpha*a1*xi c22*beta*N2*alpha*a1*xi c32*beta*N2*alpha*a1*xi 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 c13*beta*N3*alpha*a2 c13*beta*N3*a2 c13*beta*N3*xi*a1 0 0 c23*beta*N3*alpha*a2 c23*beta*N3*a2 c23*beta*N3*xi*a1 0 0 c33*beta*N3*alpha*a2 c33*beta*N3*a2 c33*beta*N3*xi*a1 0 c13*beta*N3*alpha*a1*xi c23*beta*N3*alpha*a1*xi c33*beta*N3*alpha*a1*xi 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]; %matrix of disease transmission components


F=[0 c11*beta*N1*alpha*a2 c11*beta*N1*a2 c11*beta*N1*xi*a1 0 0 c21*beta*N1*alpha*a2 c21*beta*N1*a2 c21*beta*N1*xi*a1 0 0 c31*beta*N1*alpha*a2 c31*beta*N1*a2 c31*beta*N1*xi*a1 0 c11*beta*N1*alpha*xi*a1 1.26*beta*N1*alpha*xi*a1 c31*beta*N1*alpha*xi*a1 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 c12*beta*N2*alpha*a2 c12*beta*N2*a2 c12*beta*N2*xi*a1 0 0 c22*beta*N2*alpha*a2 c22*beta*N2*a2 c22*beta*N2*xi*a1 0 0 c32*beta*N2*alpha*a2 c32*beta*N2*a2 c32*beta*N2*xi*a1 0 c12*beta*N2*alpha*xi*a1 c22*beta*N2*alpha*xi*a1 c32*beta*N2*alpha*xi*a1 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 c13*beta*N3*alpha*a2 c13*beta*N3*a2 c13*beta*N3*xi*a1 0 0 c23*beta*N3*alpha*a2 c23*beta*N3*a2 c23*beta*N3*xi*a1 0 0 c33*beta*N3*alpha*a2 c33*beta*N3*a2 c33*beta*N3*xi*a1 0 c13*beta*N3*alpha*xi*a1 c23*beta*N3*alpha*xi*a1 c33*beta*N3*alpha*xi*a1 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]; %matrix of disease transmission components
%%

beta1 = 1.9e-07; 
beta2 = 0.9539*beta1;  
beta3 = 0.315*beta1; 

NextGen = F/V;
NextGen = subs(NextGen,[beta S1 E1 Iu1 Ia1 Is1 Ih1 S2 E2 Iu2 Ia2 Is2 Ih2 S3 E3 Iu3 Ia3 Is3 Ih3 Iai1 Iai2 Iai3 Ei1 Ei2 Ei3 c11 c12 c13 c21 c22 c23 c31 c32 c33 alpha xi k p ps gammaa q psi1 psi2 psi3 gammas mus gammah muh1 muh2 muh3 m1 m2 m3 N1 N2 N3 a1 a2 b1 test],[beta1 1.7478e+06 235 23.16 10 25 0 5.8138e+06 781.3 77 31 64 14 1.4852e+06 200 20 8 10 10 0 0 0 0 0 0 2.2500 0.3806 0.2537 1.2599 0.9800 0.9635 0.2080 0.2469 0.7100 0.5 0.04047 1/4.6 0.66 0.66 1/6.5 1/1.5 0.0017 0.046 0.227 1/6.5 0 1/10.4 0.05*0.5 0.054*0.5 0.472*0.5 0 0 0.4 1.7478e+06+235+23.16+10+25  5.8138e+06+781.3+77+31+64+14  1.4852e+06+200+20+8+10+10 0.2 0.2 0 0]); %1.7478e+06+235+23.16+10+25  5.8138e+06+781.3+77+31+64+14  1.4852e+06+200+20+8+10+10   1739700 5800400 1481700
%NextGen = subs(NextGen,[beta S1 E1 Iu1 Ia1 Is1 Ih1 S2 E2 Iu2 Ia2 Is2 Ih2 S3 E3 Iu3 Ia3 Is3 Ih3 Iai1 Iai2 Iai3 Ei1 Ei2 Ei3 c11 c12 c13 c21 c22 c23 c31 c32 c33 alpha xi k p ps gammaa q psi1 psi2 psi3 gammas mus gammah muh1 muh2 muh3 m1 m2 m3 N1 N2 N3 a1 a2 b1 test],[beta1 1.7478e+06 235 23.16 10 25 0 5.8138e+06 781.3 77 31 64 14 1.4852e+06 200 20 8 10 10 0 0 0 0 0 0 2.2500 0.3806 0.0799 1.2599 0.9800 0.3035 0.0655 0.0778 0.2236 0.5 0.04047 1/4.6 0.66 0.66 1/6.5 1/1.5 0.0017 0.046 0.227 1/6.5 0 1/10.4 0.05*0.5 0.054*0.5 0.472*0.5 0 0 0.4 1.7478e+06+235+23.16+10+25  5.8138e+06+781.3+77+31+64+14  1.4852e+06+200+20+8+10+10 1 1 0 0]); %1.7478e+06+235+23.16+10+25  5.8138e+06+781.3+77+31+64+14  1.4852e+06+200+20+8+10+10   1739700 5800400 1481700


e = eig(NextGen);

beta1 = 1.9e-07; %1.047e-07;
beta2 = 0.9539;  
beta3 = 0.3803; %3803
