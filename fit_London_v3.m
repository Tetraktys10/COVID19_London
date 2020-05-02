function SOL=fit_London_v3(m1,m2,m3,b,x)%m1,m2,


y0= [1.7478e+06    235    23.16    10    25    0         0         0         0    11    5.8138e+06    781.3    77    31    64    14    0     0         0    34    1.4852e+06    200    20    8    10    10         0         0         0    15];
tspan=[0 35];

%define the parameters



beta1 = 1.9e-07;
beta2 = 0.9539;  
beta3 =  0.315;%0.3311;%  (0.299, 0.362) %0.315;
d2 =     0.055;%  (0.055, 0.055)
xi=0.04668; %(0,0.12)
d1 = 0.01539;
d1s=d1;
d2s=d2;

c = [2.2500    0.3806    0.2537;   1.2599    0.9800    0.9635;   0.2080    0.2469    0.7100];       %       %contact matrix
alpha = 0.5;        %reduction of infectiousness from asymptomatic people 50%
k = 1/4.6;          %progression from exposed to infectious initial stage, 4.6 days
q = 1/1.5;          %time where the symptomatic case remains undetected and where the asymptomatic case remains mixed with the symptomatic before moving to its own compartment
p = 0.66;           %probability of being symptomatic 
ps=0.66;
gammaa= 1/6.5;      %recovery rate, infectious period = 1/6.5 where 6.5 days is the mean infectiousness time
gammas= 1/6.5;      %recovery rate, infectious period = 1/6.5 where 6.5 days is the mean infectiousness time
gammah = 1/10.4;    %recovery rate, infectious period = 1/10.4 where 10.4 is the mean hospitalization period
psi1 = 0.0017;      %percentage of 0-14 cases that require hospitalization 0.17%
psi2 = 0.044;       %percentage of 15-59 cases that require hospitalization 4.4%
psi3 = 0.227;       %percentage of 60+ cases that require hospitalization 22.7%
mus = 0;            %mortality rate of symptomatic people
muh1 = 0.05*0.5;%mortality rate of hospitalised people: 5% of hospitalised require critical care, of which 50% will die
muh2 = 0.054*0.5;%mortality rate of hospitalised people: 5.4% of hospitalised require critical care, of which 50% will die
muh3 = b*0.472*0.5;%mortality rate of hospitalised people: 47.2% of hospitalised require critical care, of which 50% will die
% m1 = 0;
% m2 = 0;
% m3 = 0.29;

[t,y] = ode45(@(t,y) seir_eq_v3(t,y,beta1,beta2,beta3,c,k,p,gammaa,gammas,gammah,mus,muh1,muh2,muh3,psi1,psi2,psi3,alpha,q,xi,d1,d2,m1,m2,m3,d1s,d2s,ps),tspan,y0);



%SOL = zeros(size(years, 1), size(years, 2));
 SOL = zeros(size(x, 1), size(x, 2));
% %Fill the output matrix
    for r=1:size(SOL,1)
        for c=1:size(SOL,2)
           SOL(r,c) = interp1q(t,y(:,9)+y(:,19)+y(:,29),x(r, c));  %deaths       y(:,10)+y(:,20)+y(:,30) %notifications
        end
    end

end