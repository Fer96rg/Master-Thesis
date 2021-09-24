clc
clear all 
% close all

% DATOS:

B_r=10^6;
theta_r=10^-2;
theta_v=10;
theta_d=100;

%INPUTS: T y beta_i
%Onda incidente: calculo TODOS los parámetros
load('curva_RH_def_100');  %Cargo los vectores T, R_def, alpha_def, P_def, M1n_def y M2n_def
% Para cierto salto de temperaturas T(m)

m=553;  %Para M1n_i(m)=4.91 -> M1_i=3

% m=320;   %Para M1n_i(m)=1.5  -> M1_i=3
% m=1623;  %Para M1n_i(m)=3    -> M1_i=6
% m=3867;  %Para M1n_i(m)=5    -> M1_i=10

T_i=T(m);
R_i=R_def(m);
alpha_2=alpha_def(m);
P_i=P_def(m);
M1n_i=M1n_def(m);
M2n_i=M2n_def(m);
%Para beta_i determinado

beta_i=37.62;

M1_i=M1n_i./sind(beta_i);
theta=atand((R_i-1)./(tand(beta_i)+R_i./tand(beta_i)));
M2_i=M2n_i./sind(beta_i-theta);

% Onda reflejada
syms x beta_r alpha_3
c2_c1=M1n_i./(R_i.*M2n_i);

syms alpha_3 beta_r

% Cálculo de los parámetros adimensionales por dos formas (en comentarios la segunda)

R_r=((tand(theta).*tand(beta_r)+1)./(1-tand(theta)./tand(beta_r)));
T_r=((1+alpha_2+M2_i.^2.*(sind(beta_r)).^2.*7.*R_i./(5.*P_i).*c2_c1.^2.*(1+alpha_2).*(1-R_r.^-1))./(R_r.*(1+alpha_3)));
fun(1)=@(alpha_3,beta_r) alpha_3.^2./(1-alpha_3)-B_r.*exp(-theta_d./(T_i.*T_r)).*sqrt(T_i.*T_r)./(R_i.*R_r).*(1-exp(-theta_v./(T_i.*T_r)));
fun(2)=@(alpha_3,beta_r) T_r-((1+alpha_2).*(1-R_r.^-1)+alpha_2+5+(2.*theta_v.*(1-alpha_2))./(T_i.*(exp(theta_v./T_i)-1))-2.*theta_d.*(alpha_3-alpha_2)./T_i)./((1+alpha_3).*(1+R_r)+alpha_3+5+(2.*theta_v.*(1-alpha_3))./(T_i.*T_r.*(exp(theta_v./(T_i.*T_r))-1)));

% R_r=((tand(theta)*tand(beta_r)+1)/(1-tand(theta)/tand(beta_r)));
% T_r=((1+alpha_2+M2_i^2*(sind(beta_r))^2*7*R_i/(5*P_i)*c2_c1^2*(1+alpha_2)*(1-R_r^-1))/(R_r*(1+alpha_3)));
% fun(1)=@(alpha_3,beta_r) alpha_3^2/(1-alpha_3)-B_r*exp(-theta_d/(T_i*T_r))*sqrt(T_i*T_r)/(R_i*R_r)*(1-exp(-theta_v/(T_i*T_r)));
% fun(2)=@(alpha_3,beta_r) T_r-((1+alpha_2)*(1-R_r^-1)+alpha_2+5+(2*theta_v*(1-alpha_2))/(T_i*(exp(theta_v/T_i)-1))-2*theta_d*(alpha_3-alpha_2)/T_i)/((1+alpha_3)*(1+R_r)+alpha_3+5+(2*theta_v*(1-alpha_3))/(T_i*T_r*(exp(theta_v/(T_i*T_r))-1)));

% % Sustituyo R_r en T_r 
T_r=((1+alpha_2+M2_i.^2.*(sind(beta_r)).^2.*7.*R_i./(5.*P_i).*c2_c1.^2.*(1+alpha_2).*(1-((tand(theta).*tand(beta_r)+1)./(1-tand(theta)./tand(beta_r))).^-1))./(((tand(theta).*tand(beta_r)+1)./(1-tand(theta)./tand(beta_r))).*(1+alpha_3)));

% T_r=((1+alpha_2+M2_i^2*(sind(beta_r))^2*7*R_i/(5*P_i)*c2_c1^2*(1+alpha_2)*(1-((tand(theta)*tand(beta_r)+1)/(1-tand(theta)/tand(beta_r)))^-1))/(((tand(theta)*tand(beta_r)+1)/(1-tand(theta)/tand(beta_r)))*(1+alpha_3)));

% % Sustituyo R_r y T_r en fun1 y fun2
fun1=@(alpha_3,beta_r) alpha_3.^2./(1-alpha_3)-B_r.*exp(-theta_d./(T_i.*((1+alpha_2+M2_i.^2.*(sind(beta_r)).^2.*7.*R_i./(5.*P_i).*c2_c1.^2.*(1+alpha_2).*(1-((tand(theta).*tand(beta_r)+1)./(1-tand(theta)./tand(beta_r))).^-1))./(((tand(theta).*tand(beta_r)+1)./(1-tand(theta)./tand(beta_r))).*(1+alpha_3))))).*sqrt(T_i.*((1+alpha_2+M2_i.^2.*(sind(beta_r)).^2.*7.*R_i./(5.*P_i).*c2_c1.^2.*(1+alpha_2).*(1-((tand(theta).*tand(beta_r)+1)./(1-tand(theta)./tand(beta_r))).^-1))./(((tand(theta).*tand(beta_r)+1)./(1-tand(theta)./tand(beta_r))).*(1+alpha_3))))./(R_i.*((tand(theta).*tand(beta_r)+1)./(1-tand(theta)./tand(beta_r)))).*(1-exp(-theta_v./(T_i.*((1+alpha_2+M2_i.^2.*(sind(beta_r)).^2.*7.*R_i./(5.*P_i).*c2_c1.^2.*(1+alpha_2).*(1-((tand(theta).*tand(beta_r)+1)./(1-tand(theta)./tand(beta_r))).^-1))./(((tand(theta).*tand(beta_r)+1)./(1-tand(theta)./tand(beta_r))).*(1+alpha_3))))));
fun2=@(alpha_3,beta_r) ((1+alpha_2+M2_i.^2.*(sind(beta_r)).^2.*7.*R_i./(5.*P_i).*c2_c1.^2.*(1+alpha_2).*(1-((tand(theta).*tand(beta_r)+1)./(1-tand(theta)./tand(beta_r))).^-1))./(((tand(theta).*tand(beta_r)+1)./(1-tand(theta)./tand(beta_r))).*(1+alpha_3)))-((1+alpha_2).*(1-((tand(theta).*tand(beta_r)+1)./(1-tand(theta)./tand(beta_r))).^-1)+alpha_2+5+(2.*theta_v.*(1-alpha_2))./(T_i.*(exp(theta_v./T_i)-1))-2.*theta_d.*(alpha_3-alpha_2)./T_i)./((1+alpha_3).*(1+((tand(theta).*tand(beta_r)+1)./(1-tand(theta)./tand(beta_r))))+alpha_3+5+(2.*theta_v.*(1-alpha_3))./(T_i.*((1+alpha_2+M2_i.^2.*(sind(beta_r)).^2.*7.*R_i./(5.*P_i).*c2_c1.^2.*(1+alpha_2).*(1-((tand(theta).*tand(beta_r)+1)./(1-tand(theta)./tand(beta_r))).^-1))./(((tand(theta).*tand(beta_r)+1)./(1-tand(theta)./tand(beta_r))).*(1+alpha_3))).*(exp(theta_v./(T_i.*((1+alpha_2+M2_i.^2.*(sind(beta_r)).^2.*7.*R_i./(5.*P_i).*c2_c1.^2.*(1+alpha_2).*(1-((tand(theta).*tand(beta_r)+1)./(1-tand(theta)./tand(beta_r))).^-1))./(((tand(theta).*tand(beta_r)+1)./(1-tand(theta)./tand(beta_r))).*(1+alpha_3)))))-1)));
% 
% fun1=@(alpha_3,beta_r) alpha_3^2/(1-alpha_3)-B_r*exp(-theta_d/(T_i*((1+alpha_2+M2_i^2*(sind(beta_r))^2*7*R_i/(5*P_i)*c2_c1^2*(1+alpha_2)*(1-((tand(theta)*tand(beta_r)+1)/(1-tand(theta)/tand(beta_r)))^-1))/(((tand(theta)*tand(beta_r)+1)/(1-tand(theta)/tand(beta_r)))*(1+alpha_3)))))*sqrt(T_i*((1+alpha_2+M2_i^2*(sind(beta_r))^2*7*R_i/(5*P_i)*c2_c1^2*(1+alpha_2)*(1-((tand(theta)*tand(beta_r)+1)/(1-tand(theta)/tand(beta_r)))^-1))/(((tand(theta)*tand(beta_r)+1)/(1-tand(theta)/tand(beta_r)))*(1+alpha_3))))/(R_i*((tand(theta)*tand(beta_r)+1)/(1-tand(theta)/tand(beta_r))))*(1-exp(-theta_v/(T_i*((1+alpha_2+M2_i^2*(sind(beta_r))^2*7*R_i/(5*P_i)*c2_c1^2*(1+alpha_2)*(1-((tand(theta)*tand(beta_r)+1)/(1-tand(theta)/tand(beta_r)))^-1))/(((tand(theta)*tand(beta_r)+1)/(1-tand(theta)/tand(beta_r)))*(1+alpha_3))))));
% fun2=@(alpha_3,beta_r) ((1+alpha_2+M2_i^2*(sind(beta_r))^2*7*R_i/(5*P_i)*c2_c1^2*(1+alpha_2)*(1-((tand(theta)*tand(beta_r)+1)/(1-tand(theta)/tand(beta_r)))^-1))/(((tand(theta)*tand(beta_r)+1)/(1-tand(theta)/tand(beta_r)))*(1+alpha_3)))-((1+alpha_2)*(1-((tand(theta)*tand(beta_r)+1)/(1-tand(theta)/tand(beta_r)))^-1)+alpha_2+5+(2*theta_v*(1-alpha_2))/(T_i*(exp(theta_v/T_i)-1))-2*theta_d*(alpha_3-alpha_2)/T_i)/((1+alpha_3)*(1+((tand(theta)*tand(beta_r)+1)/(1-tand(theta)/tand(beta_r))))+alpha_3+5+(2*theta_v*(1-alpha_3))/(T_i*((1+alpha_2+M2_i^2*(sind(beta_r))^2*7*R_i/(5*P_i)*c2_c1^2*(1+alpha_2)*(1-((tand(theta)*tand(beta_r)+1)/(1-tand(theta)/tand(beta_r)))^-1))/(((tand(theta)*tand(beta_r)+1)/(1-tand(theta)/tand(beta_r)))*(1+alpha_3)))*(exp(theta_v/(T_i*((1+alpha_2+M2_i^2*(sind(beta_r))^2*7*R_i/(5*P_i)*c2_c1^2*(1+alpha_2)*(1-((tand(theta)*tand(beta_r)+1)/(1-tand(theta)/tand(beta_r)))^-1))/(((tand(theta)*tand(beta_r)+1)/(1-tand(theta)/tand(beta_r)))*(1+alpha_3)))))-1)));
% 
% % Supón que conoces Beta_r, y usas una ecuación para calcular Tr.
% % Usas la otra ecuación que no has usado para obtener un nuevo Beta_r.

beta_r_inicial=30.34;
alpha_3_sol=vpasolve(fun1(alpha_3,beta_r_inicial),alpha_3,alpha_2);
beta_r_sol=vpasolve(fun2(alpha_3_sol,beta_r),beta_r,beta_r_inicial);

%Método alternativo de resolución (por defecto en comentarios)
%x(1)=alpha_3
%x(2)=beta_r
% 
% %En esta función están las dos ecuaciones:
% F =@(x) [x(1)^2/(1-x(1))-B_r*exp(-theta_d/(T_i*((1+alpha_2+M2_i^2*(sind(x(2)))^2*7*R_i/(5*P_i)*c2_c1^2*(1+alpha_2)*(1-(tand(x(2))*(tand(x(2))*tand(theta)-1)/(tand(x(2))-tand(theta)))^-1))/((tand(x(2))*(tand(x(2))*tand(theta)-1)/(tand(x(2))-tand(theta)))*(1+x(1))))))*sqrt(T_i*((1+alpha_2+M2_i^2*(sind(x(2)))^2*7*R_i/(5*P_i)*c2_c1^2*(1+alpha_2)*(1-(tand(x(2)).*(tand(x(2)).*tand(theta)-1)./(tand(x(2))-tand(theta)))^-1))/((tand(x(2)).*(tand(x(2)).*tand(theta)-1)./(tand(x(2))-tand(theta)))*(1+x(1)))))/(R_i*(tand(x(2)).*(tand(x(2))*tand(theta)-1)/(tand(x(2))-tand(theta))))*(1-exp(-theta_v/(T_i*((1+alpha_2+M2_i^2*(sind(x(2)))^2*7*R_i/(5*P_i)*c2_c1^2*(1+alpha_2)*(1-(tand(x(2))*(tand(x(2))*tand(theta)-1)/(tand(x(2))-tand(theta)))^-1))/((tand(x(2))*(tand(x(2))*tand(theta)-1)/(tand(x(2))-tand(theta)))*(1+x(1)))))));
% ((1+alpha_2+M2_i^2*(sind(x(2)))^2*7*R_i/(5*P_i)*c2_c1^2*(1+alpha_2)*(1-(tand(x(2))*(tand(x(2))*tand(theta)-1)/(tand(x(2))-tand(theta)))^-1))/((tand(x(2))*(tand(x(2))*tand(theta)-1)/(tand(x(2))-tand(theta)))*(1+x(1))))-((1+alpha_2)*(1-(tand(x(2))*(tand(x(2))*tand(theta)-1)/(tand(x(2))-tand(theta)))^-1)+alpha_2+5+(2*theta_v*(1-alpha_2))/(T_i*(exp(-theta_v/T_i)-1))-2*theta_d*(x(1)-alpha_2)/T_i)/((1+x(1))*(1+(tand(x(2))*(tand(x(2))*tand(theta)-1)/(tand(x(2))-tand(theta))))+x(1)+5+(2*theta_v*(1-x(1)))/(T_i*((1+alpha_2+M2_i^2*(sind(x(2)))^2*7*R_i/(5*P_i)*c2_c1^2*(1+alpha_2)*(1-(tand(x(2))*(tand(x(2))*tand(theta)-1)/(tand(x(2))-tand(theta)))^-1))/((tand(x(2))*(tand(x(2))*tand(theta)-1)/(tand(x(2))-tand(theta)))*(1+x(1))))*(exp(-theta_v/(T_i*((1+alpha_2+M2_i^2*(sind(x(2)))^2*7*R_i/(5*P_i)*c2_c1^2*(1+alpha_2)*(1-(tand(x(2))*(tand(x(2))*tand(theta)-1)/(tand(x(2))-tand(theta)))^-1))/((tand(x(2))*(tand(x(2))*tand(theta)-1)/(tand(x(2))-tand(theta)))*(1+x(1))))))-1)))];
% 
% % fun=@sistema;
% x0=[alpha_2,beta_i];
% %Para de iterar en 400, por lo que quiero aumentar las iteraciones hasta
% %1500: el problema es que no me deja cambiarlo
% options = optimoptions(@fsolve,'MaxIterations',1500);
% x = fsolve(F,x0);

