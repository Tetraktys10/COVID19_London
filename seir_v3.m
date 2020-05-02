%COVID-19 London Model with ages

clear all

% All assumptions and parameters come from Imperial report "Impact of non-pharmaceutical interventions (NPIs) to reduce COVID19 mortality and healthcare demand"
% https://www.imperial.ac.uk/media/imperial-college/medicine/sph/ide/gida-fellowships/Imperial-College-COVID19-NPI-modelling-16-03-2020.pdf

%Initialise differential equations

%main model
y0= [1.7478e+06    235    23.16    10    25    0         0         0         0    11    5.8138e+06    781.3    77    31    64    14    0     0         0    34    1.4852e+06    200    20    8    10    10         0         0         0    15];

%scenarios
%y0= [1.7478e+06    235    23.16    10    25    0         0         0         0    11    5.8138e+06    781.3    77    31    64    14    0     0         0    34    1.4852e+06    200    20    8    10    10         0         0         0    15 0 0 0 0 0 0];

tspan=[0 60];

%define the parameters

%From calibration
beta1 = 1.9e-07; 
beta2 = 0.9539;  
beta3 =  0.315;  %(0.299, 0.362) 
d1 = 0.01539;
d2 = 0.055; %(0.055, 0.055)
xi=0.04668; %(0,0.12)

d1s=d1;
d2s=d2;


c = [2.2500    0.3806    0.2537;   1.2599    0.9800    0.9635;   0.2080    0.2469    0.7100];    %contact matrix from Mossong et al. adjusted to 2019 London population
alpha = 0.5;        %reduction of infectiousness from asymptomatic people 50%
k = 1/4.6;          %progression from exposed to infectious initial stage, 4.6 days
q = 1/1.5;          %time where the symptomatic case remains undetected and where the asymptomatic case remains mixed with the symptomatic before moving to its own compartment
p = 0.66;           %probability of being symptomatic 
ps=0.66;            %probability of being symptomatic in the 60+ age group
gammaa= 1/6.5;      %recovery rate, infectious period = 1/6.5 where 6.5 days is the mean infectiousness time
gammas= 1/6.5;      %recovery rate, infectious period = 1/6.5 where 6.5 days is the mean infectiousness time
gammah = 1/10.4;    %recovery rate, infectious period = 1/10.4 where 10.4 is the mean hospitalization period
psi1 = 0.0017;      %percentage of 0-14 cases that require hospitalization 0.17%
psi2 = 0.044;       %percentage of 15-59 cases that require hospitalization 4.4%
psi3 = 0.227;       %percentage of 60+ cases that require hospitalization 22.7%
mus = 0;            %mortality rate of symptomatic people
muh1 = 0.05*0.5;    %mortality rate of hospitalised people in ICU (0-14): 5% of hospitalised require critical care, of which 50% will die
muh2 = 0.054*0.5;   %mortality rate of hospitalised people in ICU (15-59): 5.4% of hospitalised require critical care, of which 50% will die
muh3 = 0.9*0.472*0.5;%mortality rate of hospitalised people in ICU (60+): 47.2% of hospitalised require critical care, of which 50% will die
m1 = 0;             %mortality rate of hospitalised people not in ICU (0-14): assumption
m2 = 0;             %mortality rate of hospitalised people not in ICU (15-59): assumption
m3 = 0.29;          %mortality rate of hospitalised people not in ICU (60+): against NHS data on COVID deaths in London

%main model
[t,y] = ode45(@(t,y) seir_eq_v3(t,y,beta1,beta2,beta3,c,k,p,gammaa,gammas,gammah,mus,muh1,muh2,muh3,psi1,psi2,psi3,alpha,q,xi,d1,d2,m1,m2,m3,d1s,d2s,ps),tspan,y0);

%scenarios
%[t,y] = ode45(@(t,y) seir_eq_v5(t,y,beta1,beta2,beta3,c,k,p,gammaa,gammas,gammah,mus,muh1,muh2,muh3,psi1,psi2,psi3,alpha,q,xi,d1,d2,m1,m2,m3),tspan,y0);

%total_notified_cases = 0.6*y(:,5)+y(:,6)+0.6*y(:,13)+y(:,14)+0.6*y(:,21)+y(:,22);   %50% of infections may not be identified as cases

%data from PHE on cumulative notified cases in London 9 March 2020 to 27th March 2020 (before the change in their recording system)
%https://coronavirus.data.gov.uk/#regions
data=[60 91	104	136	167	313	407	480	621	953	1221 1588 1965 2189	2433 2872 3247 3916	4634 5299 5957 6521 7121 8341 9291 10247 10764	11978 12636 13378 14355 15217 16011]; 
days=0:1:length(data)-1;

%data on daily deaths from NHS
%https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-daily-deaths/
deaths=[1	0	4	2	6	13	9	16	19	24	19	41	46	45	34	55	87	118	122	120	176	166	112	163	152	172	210	168	172	189	199	173	149	156	132	142];

for i=2:length(deaths)
    deaths_cumul(1)=deaths(1);
    deaths_cumul(i)=deaths_cumul(i-1)+deaths(i);
end
days2=0:1:length(deaths)-1;



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOTS

%%
figure(05) %(S E Iu Ia Is Ih Ru Rh D Nc)
plot(t,y(:,9)+y(:,19)+y(:,29),'-', 'LineWidth', 2.5, 'color', 'm') %'*-','MarkerIndices',1:10:length(y)
hold on
%plot(days2,deaths_cumul, 'o')%
set(gca, 'FontSize', 16);
xlabel('Days')
ylabel('Populations')
title(['Cumulative deaths'] )
% set(gca,'XTick',[0 50 100 160 365 547]);
xticks([0 50 100 160 365 547])
xticklabels({'9/03','26/04','17/06','16/08','9/03/21','7/09/21',''})
%legend('Deaths')
%%
figure(06) %(S E Iu Ia Is Ih Ru Rh D Nc)
plot(t,y(:,4)+y(:,5)+y(:,6)+y(:,14)+y(:,15)+y(:,16)+y(:,24)+y(:,25)+y(:,26),'-', 'LineWidth', 2.5, 'color', 'r')%,'*-','MarkerIndices',1:10:length(y),
hold on
%plot(days,data, 'o')
set(gca, 'FontSize', 16);
xlabel('Days')
ylabel('Populations')
title(['Infections'] )
xticks([0 50 100 160 365 547])
xticklabels({'9/03','26/04','17/06','16/08','9/03/21','7/09/21',''})
%legend('Infections (all types)')

%%
% figure(05) %(S E Iu Ia Is Ih Ru Rh D Nc)
% plot(t,y(:,9)+y(:,19)+y(:,29),'-*','MarkerIndices',1:10:length(y), 'LineWidth', 0.5, 'color', 'm') %'-.','MarkerIndices',1:10:length(y)
% hold on
% %plot(days2,deaths, 'o')
% set(gca, 'FontSize', 16);
% xlabel('Time')
% ylabel('Populations')
% title(['Cumulative deaths'] )
% xticks([0 50 100 160 365 547])
% xticklabels({'9/03','26/04','17/06','16/08','9/03/21','7/09/21',''})
% %legend('Days')
% %%
% figure(06) %(S E Iu Ia Is Ih Ru Rh D Nc)
% plot(t,y(:,4)+y(:,5)+y(:,6)+y(:,14)+y(:,15)+y(:,16)+y(:,24)+y(:,25)+y(:,26),'-*','MarkerIndices',1:10:length(y), 'LineWidth', 0.5, 'color', 'r')%
% hold on
% set(gca, 'FontSize', 16);
% xlabel('Days')
% ylabel('Populations')
% title(['Infections'] )
% xticks([0 50 100 160 365 547])
% xticklabels({'9/03','26/04','17/06','16/08','9/03/21','7/09/21',''})
% %legend('Infections (all types)')



%%% (y(:,1)+y(:,11)+y(:,21))/(y(:,1)+y(:,2)+y(:,3)+y(:,4)+y(:,5)+y(:,6)+y(:,7)+y(:,8)+y(:,9)+y(:,11)+y(:,12)+y(:,13)+y(:,14)+y(:,15)+y(:,16)+y(:,17)+y(:,18)+y(:,19)+y(:,21)+y(:,22)+y(:,23)+y(:,24)+y(:,25)+y(:,26)+y(:,27)+y(:,28)+y(:,29))
