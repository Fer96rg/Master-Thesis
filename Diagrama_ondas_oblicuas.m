clc
clear all
% close all

% A continuación se representa el diagrama de ondas de choque oblicuas. Este
% diagrama se rige por una ecuación que relaciona:
%     -El ángulo de deflexión, theta.
%     -El ángulo de inclinación de la onda, beta.
%     -El Mach del fluido que incide sobre la onda, M1.

%Datos y definición de parámetros:
gamma=1.4;   % Gas diatómico ideal (aire compuesto por N2 y O2)
M1=[1.1 1.2 1.3 1.4 1.6 1.8 2 2.25 2.5 3 4 6 8 15 1000];    % M1=1000 equivale a Infinito
beta=NaN(length(M1));

for j=1:length(M1)
    rango=length(asind(1/M1(j)):1:90);
    beta(1:rango,j)=asind(1/M1(j)):1:90;
end
beta(beta==0)=NaN;
for j=1:length(M1)
    for i=1:90 
       theta_x(i,j)=atand((M1(j)^2*(sind(beta(i,j)))^2-1)*2*cotd(beta(i,j))/(gamma*M1(j)^2+M1(j)^2*cosd(2*beta(i,j))+2));
       M2(i,j)=sqrt((2+(gamma-1)*M1(j)^2)/(2*gamma*M1(j)^2*sind(beta(i,j))^2-gamma+1)+(2*M1(j)^2*(cosd(beta(i,j)))^2)/((gamma-1)*M1(j)^2*(sind(beta(i,j)))^2+2));
       if M2(i,j)<1.001
            beta_M2_1(j)=beta(i,j);
            theta_M2_1(j)=theta_x(i,j);
            break;
       end
    end
end

M1=[1.1 1.2 1.3 1.4 1.6 1.8 2 2.25 2.5 3 4 6 8 15 1000];    % M1=1000 equivale a Infinito
beta_T=0:1:90;
for j=1:length(M1)
    for i=1:length(beta_T)
        theta(i,j)=atand((M1(j)^2*(sind(beta_T(i)))^2-1)*2*cotd(beta_T(i))/(gamma*M1(j)^2+M1(j)^2*cosd(2*beta_T(i))+2));
    end
    
end

for j=1:length(M1)
    plot(theta(:,j),beta_T,'-k');
    hold on
end
hold on
plot(theta_M2_1,beta_M2_1,'r');
hold on
[a,b]=max(theta);
plot(a,beta_T(b),'-b');
axis([0 50 0 90]);
grid on;
xlabel(texlabel('theta [º]'));
ylabel(texlabel('beta [º]'),'Rotation',0);

% Curvas correspondientes a valores constantes de M1.
% Representación gráfica de la curva formada por los puntos de ángulo de
% deflexión máxima, theta máxima, correspondientes a cada valor de M1
% (dibujada en rojo)

text(1.40,71, '1.1');
text(3.85,67, '1.2');
text(6.7,65, '1.3');
text(9.55,63, '1.4');
text(14.5,61.5, '1.6');
text(19.0,60.75, '1.8');
text(22.9,60.5,'2');
text(26.51,60.375, '2.25');
text(29.7531,60.4, '2.5');
text(33.9729,61.5, '3');
text(38.6436,63, '4');
text(42.4,64.5,'6');
text(43.8,66,'8');
text(45.8,70,'Inf');
text(19.5,67,'$\theta=\theta_{max}$','interpreter','latex','Color','b');
text(32.5,85,'Onda fuerte, $M_{2}<1$','interpreter','latex');
text(32.5,38,'Onda debil, $M_{2}>1$','interpreter','latex');
text(30,62,'$M_2=1$','interpreter','latex','Color','red');







