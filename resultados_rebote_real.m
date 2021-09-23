clc
clear all 
close all

B_r=10^6;
theta_r=10^-2;
theta_v=10;
theta_d=100;

load('curva_RH_def_100');  %Cargo los vectores T, R_def, alpha_def, P_def, M1n_def y M2n_def
%%
% Si M1_i y beta_i son conocidos, conoczco todo de la onda incidente
M1_i=13;
beta_i=45;
M1n=M1_i*sind(beta_i);
for i=1:length(M1n_def)
    if M1n_def(i)>M1n
        index=i;   
        break;           
    end
end
T_i=T(index);
R_i=R_def(index);
alpha_2=alpha_def(index);
P_i=P_def(index);
M1n_i=M1n_def(index);
M2n_i=M2n_def(index);
theta=atand((R_i-1)/(tand(beta_i)+R_i/tand(beta_i)));
M2_i=M2n_i/sind(beta_i-theta);
c2_c1=M1n_i/(R_i*M2n_i);    




%%
%Onda reflejada

syms alpha_3 T_r

%Fórmulas sin sustituir:
% (23) P_r=T_r*R_r*(1+alpha_3)/(1+alpha_2);
% (25) R_r=tand(beta_r)/tand(beta_r-theta);
% (26) T_r=(6+alpha_2-R_r^-1*(1+alpha_2)+2*(1-alpha_2)*theta_v/T_i/(exp(theta_v/T_i)-1)-2*(1-alpha_3)*theta_v/T_i/(exp(theta_v/(T_i*T_r))-1)-2*(alpha_3-alpha_2)*theta_d/T_i)/(6+2*alpha_3-R_r*(1+alpha_3));
% (27) alpha_3^2/(1-alpha_3)-B_r*exp(-theta_d/(T_i*T_r))*sqrt(T_i*T_r)/(R_i*R_r)*(1-exp(-theta_v/(T_i*T_r)));
% (29) ((1+alpha_3)*R_r*T_r-(1+alpha_2))/(7*R_i/(5*P_i)*c2_c1^2*(1+alpha_2)*(1-R_r^-1)*M2_i^2)=(tand(beta_r))^2/(1+(tand(beta_r))^2);
% (30) tand(beta_r)=((R_r-1)/(2*tand(theta))*(1-sqrt(1-4*R_r*(tand(theta))^2/(R_r-1)^2)));

%De (27), despejo R_r->R_r=f(alpha_3,T_r)
% R_r=((1-alpha_3)/alpha_3^2*B_r*exp(-theta_d/(T_i*T_r))*sqrt(T_i*T_r)/R_i*(1-exp(-theta_v/(T_i*T_r)))); 

%Sustituyo R_r en (26)->(26)'=f(alpha_3,T_r)
% fun1=@(alpha_3,T_r) T_r-(6+alpha_2-((1-alpha_3)/alpha_3^2*B_r*exp(-theta_d/(T_i*T_r))*sqrt(T_i*T_r)/R_i*(1-exp(-theta_v/(T_i*T_r))))^-1*(1+alpha_2)+2*(1-alpha_2)*theta_v/T_i/(exp(theta_v/T_i)-1)-2*(1-alpha_3)*theta_v/T_i/(exp(theta_v/(T_i*T_r))-1)-2*(alpha_3-alpha_2)*theta_d/T_i)/(6+2*alpha_3-((1-alpha_3)/alpha_3^2*B_r*exp(-theta_d/(T_i*T_r))*sqrt(T_i*T_r)/R_i*(1-exp(-theta_v/(T_i*T_r))))*(1+alpha_3));

%Combino (29) y (30)->ecuacion (31)
% ((1+alpha_3)*R_r*T_r-(1+alpha_2))/(7*R_i/(5*P_i)*c2_c1^2*(1+alpha_2)*(1-R_r^-1)*M2_i^2)=(((R_r-1)/(2*tand(theta))*(1-sqrt(1-4*R_r*(tand(theta))^2/(R_r-1)^2))))^2/(1+(((R_r-1)/(2*tand(theta))*(1-sqrt(1-4*R_r*(tand(theta))^2/(R_r-1)^2))))^2);
%Y sustituyo R_r en (31)->(31')=f(alpha_3,T_r)
% fun2=@(alpha_3,T_r) ((1+alpha_3)*((1-alpha_3)/alpha_3^2*B_r*exp(-theta_d/(T_i*T_r))*sqrt(T_i*T_r)/R_i*(1-exp(-theta_v/(T_i*T_r))))*T_r-(1+alpha_2))/(7*R_i/(5*P_i)*c2_c1^2*(1+alpha_2)*(1-((1-alpha_3)/alpha_3^2*B_r*exp(-theta_d/(T_i*T_r))*sqrt(T_i*T_r)/R_i*(1-exp(-theta_v/(T_i*T_r))))^-1)*M2_i^2)-(((((1-alpha_3)/alpha_3^2*B_r*exp(-theta_d/(T_i*T_r))*sqrt(T_i*T_r)/R_i*(1-exp(-theta_v/(T_i*T_r))))-1)/(2*tand(theta))*(1-sqrt(1-4*((1-alpha_3)/alpha_3^2*B_r*exp(-theta_d/(T_i*T_r))*sqrt(T_i*T_r)/R_i*(1-exp(-theta_v/(T_i*T_r))))*(tand(theta))^2/(((1-alpha_3)/alpha_3^2*B_r*exp(-theta_d/(T_i*T_r))*sqrt(T_i*T_r)/R_i*(1-exp(-theta_v/(T_i*T_r))))-1)^2))))^2/(1+(((((1-alpha_3)/alpha_3^2*B_r*exp(-theta_d/(T_i*T_r))*sqrt(T_i*T_r)/R_i*(1-exp(-theta_v/(T_i*T_r))))-1)/(2*tand(theta))*(1-sqrt(1-4*((1-alpha_3)/alpha_3^2*B_r*exp(-theta_d/(T_i*T_r))*sqrt(T_i*T_r)/R_i*(1-exp(-theta_v/(T_i*T_r))))*(tand(theta))^2/(((1-alpha_3)/alpha_3^2*B_r*exp(-theta_d/(T_i*T_r))*sqrt(T_i*T_r)/R_i*(1-exp(-theta_v/(T_i*T_r))))-1)^2))))^2);

%% A partir de aqui no avanzo
%Calculo alpha_3 y T_r


fun = @(x)[x(2)-(6+alpha_2-((1-x(1))/x(1)^2*B_r*exp(-theta_d/(T_i*x(2)))*...
    sqrt(T_i*x(2))/R_i*(1-exp(-theta_v/(T_i*x(2)))))^-1*(1+alpha_2)+2*(1-alpha_2)*...
    theta_v/T_i/(exp(theta_v/T_i)-1)-2*(1-x(1))*theta_v/T_i/(exp(theta_v/(T_i*x(2)))-1)...
    -2*(x(1)-alpha_2)*theta_d/T_i)/(6+2*x(1)-((1-x(1))...
    /x(1)^2*B_r*exp(-theta_d/(T_i*x(2)))...
    *sqrt(T_i*x(2))/R_i*(1-exp(-theta_v/(T_i*x(2)))))*(1+x(1))),... %2a eq
    ((1+x(1))*((1-x(1))/x(1)^2*B_r*exp(-theta_d/(T_i*x(2)))*sqrt(T_i*x(2))...
    /R_i*(1-exp(-theta_v/(T_i*x(2)))))*x(2)-(1+alpha_2))/(7*R_i/(5*P_i)*c2_c1^2*...
    (1+alpha_2)*(1-((1-x(1))/x(1)^2*B_r*exp(-theta_d/(T_i*x(2)))*sqrt(T_i*x(2))...
    /R_i*(1-exp(-theta_v/(T_i*x(2)))))^-1)*M2_i^2)-(((((1-x(1))/x(1)^2*...
    B_r*exp(-theta_d/(T_i*x(2)))*sqrt(T_i*x(2))/R_i*...
    (1-exp(-theta_v/(T_i*x(2)))))-1)/(2*tand(theta))*(1-sqrt(1-4*((1-x(1))...
    /x(1)^2*B_r*exp(-theta_d/(T_i*x(2)))*sqrt(T_i*x(2))/R_i*...
    (1-exp(-theta_v/(T_i*x(2)))))*(tand(theta))^2/(((1-x(1))/x(1)^2*...
    B_r*exp(-theta_d/(T_i*x(2)))*sqrt(T_i*x(2))/R_i*(1-exp(-theta_v/(T_i*x(2)))))...
    -1)^2))))^2/(1+(((((1-x(1))/x(1)^2*B_r*exp(-theta_d/(T_i*x(2)))*...
    sqrt(T_i*x(2))/R_i*(1-exp(-theta_v/(T_i*x(2)))))-1)/(2*tand(theta))*...
    (1-sqrt(1-4*((1-x(1))/x(1)^2*B_r*exp(-theta_d/(T_i*x(2)))*...
    sqrt(T_i*x(2))/R_i*(1-exp(-theta_v/(T_i*x(2)))))*(tand(theta))^2/(((1-x(1))...
    /x(1)^2*B_r*exp(-theta_d/(T_i*x(2)))*sqrt(T_i*x(2))/R_i*...
    (1-exp(-theta_v/(T_i*x(2)))))-1)^2))))^2)];

[solucion] = fsolve(fun,[0.2,3.3]);
%solucion(1)=alpha_3

%% Una vez que consiga alpha_3 y T_r, puedo calcular todo
%Calculo R_r

R_r=((1-solucion(1))/solucion(1)^2*B_r*exp(-theta_d/(T_i*solucion(2)))*sqrt(T_i*solucion(2))/R_i*(1-exp(-theta_v/(T_i*solucion(2))))); 

%Calculo beta_r

beta_r_sol=atand((R_r-1)/(2*tand(theta))*(1-sqrt(1-4*R_r*(tand(theta))^2/(R_r-1)^2)));

%Calculo theta

theta_sol=atand((R_r-1)/(R_r/tand(beta_r_sol)+tand(beta_r_sol)));

%Calculo P_r

P_r=(1+solucion(1))/(1-alpha_2)*solucion(2)*R_r;


