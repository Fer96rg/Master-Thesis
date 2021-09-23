clear;
clear all
% close all;

%% Diagrama polar

% Dado un M1_i determinado, se representa el plano p2/p1 para todos los
% valores posibles de beta incidente (de asind(1/M1_i) hasta 90º) respecto
% theta. Destacar que en las relaciones de salto, sólo se tiene en cuenta
% la componente normal de la velocidad y el número de Mach
% (M1n=M1*sind(beta)) al tratarse de una onda de choque oblicua.

gamma=1.4;

theta=@(M1,beta,gamma) atand((M1.^2.*(sind(beta)).^2-1)*2.*cotd(beta)./(gamma*M1.^2+M1.^2.*cosd(2*beta)+2));
p2_p1=@(M1,beta,gamma) (2*gamma*M1.^2.*sind(beta).^2-gamma+1)/(gamma+1); % Sólo se tiene e cuenta la componente normal de M1
M2=@(M1,beta,gamma) sqrt((2+(gamma-1)*M1^2)./(2*gamma*M1^2*sind(beta).^2-gamma+1)+(2*M1^2*(cosd(beta)).^2)./((gamma-1)*M1^2*(sind(beta)).^2+2));

% Onda de choque a M1_i y beta_i_1 determinados:
% INPUTS:
M1_i=7;
% m -> beta_i_1(m);
beta_i_1=asind(1/M1_i):1:90;
beta_i_1(length(beta_i_1))=90;
theta_1=theta(M1_i,beta_i_1,gamma);
p2_p1_1=p2_p1(M1_i,beta_i_1,gamma);
figure()
plot(theta_1,p2_p1_1,'k');
hold on;
plot(-theta_1,p2_p1_1,'k');

% plot(1,1);
% axis([-45 45 1 max(p2_p1_1)]);
% hold on;
% 
% % Eje y (deflexión nula)
% plot([0,0],[0,20],'--k');
% hold on;
% 
% 
% %Ahora para la onda reflejada
% M2_i=M2(M1_i,beta_i_1,gamma); % Es un vector que tiene un valor asociado a cada beta_i_1 
% M1_r=M2_i;
% M1_r(M1_r<1)=NaN;
% M=length(M1_r); % beta_i_1 y M1_r tienen el mismo tamaño
% for j=1:M
%     if isnan(M1_r(j))
%         index=j-1;
%         break
%     end
% end
%     
%     
% % Se calcula la matriz beta_i_2 (filas -> valores posibles de beta_i;
% % columnas -> valores de Mi_r) y, a partir de ella, las matrices theta_2 y
% % p2_p1_2.
% 
% % i -> beta_i
% % j -> M1_r
% 
% beta_i_2=NaN(M);
% for j=1:M
%     rango=length(asind(1/M1_r(j)):1:90);
%     beta_i_2(1:rango,j)=asind(1/M1_r(j)):1:90;
%     if j<=index
%         if beta_i_2(rango,j)<90
%             beta_i_2(rango,j)=90;
%         end
%     end
% end
% theta_2=theta(M1_r,beta_i_2,gamma);
% p2_p1_2=p2_p1(M1_r,beta_i_2,gamma);
% 
% %%
% % Resultados: SE PONE EN COMENTARIO 1. ó 2.  
% 
% % 1. Te resuelve el diagrama polar para un beta_i_1 definido según m
% % (índice del vector).
% 
% % 2. Te resuelve el diagrama polar para el que existe onda reflejada dado
% % un M1_i, es decir, corta con la vertical theta=0.
% 
% % 1.Diagrama polar de la onda reflejada para cierto beta_i_1(m):
% m=12; 
% plot(theta_1(m)+theta_2(:,m),p2_p1_1(m)+p2_p1_2(:,m)-1,'g',theta_1(m)-theta_2(:,m),p2_p1_1(m)+p2_p1_2(:,m)-1,'g');
% % axis([-60 60 0 140]); 
% xlabel(texlabel('theta [º]'));
% ylabel('$\frac{p_2}{p_1}$','interpreter','latex','Rotation',0);
% % text(15,1.5,sprintf('M_1 = %.f', M1_i)); 
% % text(15,3,sprintf('M_2 = %.2f',M1_r(m)));
% % text(15,5,'i'); text(15,7,'r');


% plot(theta_2(:,m),p2_p1_2(:,m),'g',-theta_2(:,m),+p2_p1_2(:,m),'g');

% EJEMPLO DE ONDA DE CHOQUE NOMAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% text(15,5, sprintf('\\beta_{i} = %.2f', beta_i_1(m)));
% axis([-45 45 1 15]);
% xlabel(texlabel('theta [º]'));
% ylabel('$p_2/p_1$','interpreter','latex','Rotation',0);
% text(15,1.5,'$M=3$','interpreter','latex');

% % 2.Diagrama polar de la onda reflejada que corta tangencialmente con la
% % vertical theta=0. Condiciones límites beta_i_1_max:
% %m = beta i
% m=1;
% while m<length(beta_i_1)
%     if theta_1(m)-max(theta_2(:,m))>0
%         m=m-1;
%         break;
%     end
%     m=m+1;
% end
% beta_i_1_max=beta_i_1(m);
% plot(theta_1(m)+theta_2(:,m),p2_p1_1(m)+p2_p1_2(:,m)-1,'g',theta_1(m)-theta_2(:,m),p2_p1_1(m)+p2_p1_2(:,m)-1,'g');
% plot(theta_2(:,m),p2_p1_2(:,m),'g',-theta_2(:,m),p2_p1_2(:,m),'g');
% text(15,9, sprintf('\\beta_{i}^{max} = %.2f', beta_i_1_max));
% axis([-60 60 1 15]);
% xlabel(texlabel('theta [º]'));
% ylabel('$\frac{p_2}{p_1}$','interpreter','latex','Rotation',0);
% text(15,1.5,sprintf('M_1 = %.f', M1_i));
% text(15,3,sprintf('M_2 = %.2f', M1_r(m)));
% text(15,5,'i');
% text(15,7,'r');








