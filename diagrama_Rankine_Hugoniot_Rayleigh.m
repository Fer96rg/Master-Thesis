clear all
clc
close all

%Diagrama RH y Linea de Rayleigh

gamma=1.4;
Gamma=(1.4+1)/(1.4-1);

rel_rho_inv=0:0.001:8;
rel_p1=(6-rel_rho_inv)./(6.*rel_rho_inv-1);
rel_p1(rel_p1<0)=NaN;
% rel_p2=1+gamma*M1^2*(1-rel_rho_inv);

a=1+gamma*2^2*(1-rel_rho_inv);
a(a<=1)=NaN;
b=1+gamma*1^2*(1-rel_rho_inv);
c=1+gamma*0.7^2*(1-rel_rho_inv);
c(c>=1)=NaN;


figure()
plot(rel_rho_inv,rel_p1,'r');
hold on
plot(rel_rho_inv,a,'b');
hold on
plot(rel_rho_inv,b,'k');
hold on
plot(rel_rho_inv,c,'g');
hold on
plot([Gamma^-1 Gamma^-1],[0 7],'k--');
hold on
plot([1 1],[0 7],'k--');
hold on
plot([0 7],[1 1],'k--');
axis([0 7 0 7]);
legend('Curva de Hugoniot','Linea de Rayleigh para $M_{1}>1$','Linea de Rayleigh para $M_{1}=1$','Linea de Rayleigh para $M_{1}<1$','interpreter','latex')
xlabel('$\frac{\rho_1}{\rho_2}$','interpreter','latex','FontSize',17);
ylabel('$\frac{p_2}{p_1}$','interpreter','latex','FontSize',17,'Rotation',0);
text(5.95,0.3,'$\Gamma$','interpreter','latex');
text(0.2,0.3,'$\Gamma^-1$','interpreter','latex');
text(0.5,4.5,'Punto 1');
text(1.1,1.4,'Punto 2');
text(2,0.65,'Punto 3');


% '$\rho_{2}/\rho_{1}$','$p_{2}/p_{1}$','$T_{2}/T_{1}$','interpreter','latex'

