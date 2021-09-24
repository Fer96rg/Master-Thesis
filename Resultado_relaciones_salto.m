%% Gráficas de las relaciones de salto
clc
clear all
close all

M1=1:0.5:10;
gamma=1.4;
for i=1:length(M1)
    rel_rho(i)=(gamma+1)*M1(i)^2/((gamma-1)*M1(i)^2+2);
    rel_p(i)=(2*gamma*M1(i)^2-(gamma-1))/(gamma+1);
    rel_T(i)=(2*gamma*M1(i)^2-(gamma-1))*((gamma-1)*M1(i)^2+2)/((gamma+1)^2*M1(i)^2);
    M2(i)=sqrt((2+(gamma-1)*M1(i)^2)/(2*gamma*M1(i)^2-gamma+1));
    rel_v(i)=1/rel_rho(i);
end

figure()
subplot(2,2,1);
plot(M1,rel_rho,'r');
hold on
plot([1 10],[(gamma+1)/(gamma-1) (gamma+1)/(gamma-1)],'--k');
xlabel('$M_1$','interpreter','latex');
ylh=ylabel('${\frac{\rho_{2}}{\rho_{1}}}$','interpreter','latex','Rotation',0,'FontSize',15);
ylh.Position(1) = ylh.Position(1) + 0.5;
text(2,5.3, '$\frac{\gamma+1}{\gamma-1}$','interpreter','latex','FontSize',15);
axis([1 10 1 6.5]);
grid on
% xlabel('${\mu}$','interpreter','latex', 'FontWeight','bold')


subplot(2,2,2);
semilogy(M1,rel_p,'b');
xlabel('${M_1}$','interpreter','latex');
ylabel('p_{2}/p_{1}','Rotation',0);
y2h=ylabel('${\frac{p_{2}}{p_{1}}}$','interpreter','latex','Rotation',0,'FontSize',15);
y2h.Position(1) = y2h.Position(1)+ 0.8;
axis([1 10 1 200]);
grid on

subplot(2,2,3);
semilogy(M1,rel_T,'g');
xlabel('${M_1}$','interpreter','latex');
y3h=ylabel('${\frac{T_{2}}{T_{1}}}$','interpreter','latex','Rotation',0,'FontSize',15);
y3h.Position(1) = y3h.Position(1)+ 1.2;
axis([1 10 0 20]);
grid on


subplot(2,2,4);
plot(M1,rel_v,'y');
hold on
plot([1 10],[(gamma-1)/(gamma+1) (gamma-1)/(gamma+1)],'--k');
xlabel('${M_1}$','interpreter','latex');
y4h=ylabel('${\frac{u_{2}}{u_{1}}}$','interpreter','latex','Rotation',0,'FontSize',15);
y4h.Position(1) = y4h.Position(1)+0.8;
text(6,0.3, '$\frac{\gamma-1}{\gamma+1}$','interpreter','latex','FontSize',15);
axis([1 10 0 1]);
grid on


figure()
subplot(1,2,1)
plot(M1,rel_rho,'r',M1,rel_p,'b',M1,rel_T,'g');
hold on
plot([1 10],[(gamma+1)/(gamma-1) (gamma+1)/(gamma-1)],'--k');
xlabel('${M_1}$','interpreter','latex');
axis([1 10 0 10]);
legend('$\rho_{2}/\rho_{1}$','$p_{2}/p_{1}$','$T_{2}/T_{1}$','interpreter','latex')
text(3,7, '$\frac{\gamma+1}{\gamma-1}$','interpreter','latex','FontSize',15);
grid on


subplot(1,2,2)
plot(M1,M2);
hold on
plot([1 10],[sqrt((gamma-1)/(2*gamma)) sqrt((gamma-1)/(2*gamma))],'--k');
xlabel('${M_1}$','interpreter','latex');
ylabel('${M_2}$','interpreter','latex');
axis([1 10 0 1]);
text(4,0.25, '$\sqrt{\frac{\gamma-1}{2\gamma}}$','interpreter','latex','FontSize',15);
grid on

