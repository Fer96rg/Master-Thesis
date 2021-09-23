clc
clear all
% close all

%% Resolución de reflexión en pared de onda de choque:

gamma=1.4;

% Se define beta reflejada y las ecuaciones a utilizar:

syms beta_r;
fun1=@(theta,beta,M) tand(theta)-(M^2*(sind(beta))^2-1)*2*cotd(beta)/(gamma*M^2+M^2*cosd(2*beta)+2);
% fun2=@(M2,M1,beta) M2-sqrt((2+(gamma-1)*M1.^2)./(2*gamma*M1.^2.*sind(2*beta)-gamma+1)+(2*M1.^2.*(cosd(beta)).^2)./((gamma-1)*M1.^2.*(sind(beta)).^2+2));

% i -> beta_i
% j -> M1

M1_i=[7 15];
% M1_i=[1.2 5]
N=length(M1_i); 
j=1;

while j<=N
    beta_i_min = asind(1/M1_i(j));
    i=1; %Paso inicial con guess=beta_i para el solver de beta_r
    beta_i(i,j) = beta_i_min + 0.0001;
    while beta_i(i,j)<90
        theta(i,j)=atand((M1_i(j)^2*(sind(beta_i(i,j)))^2-1)*2*cotd(beta_i(i,j))/(gamma*M1_i(j)^2+M1_i(j)^2*cosd(2*beta_i(i,j))+2));
        M2_i(i,j)=sqrt((2+(gamma-1)*M1_i(j)^2)/(2*gamma*M1_i(j)^2*sind(beta_i(i,j))^2-gamma+1)+(2*M1_i(j)^2*(cosd(beta_i(i,j)))^2)/((gamma-1)*M1_i(j)^2*(sind(beta_i(i,j)))^2+2));
        M1_r(i,j) = M2_i(i,j);
        max_theta_M1_r(i,j)=fun_max_theta(M1_r(i,j),gamma);
        
        % IMPORTANTE. Utilizar el solver sólo cuando theta sea menor que
        % theta máxima para un M1_r determinado. Si no, me da error (por
        % eso pongo la condición antes del solver).
        
        % Se sale del bucle while (con break) cuando theta es mayor que el
        % máximo valor de theta para M1_r (graficamente no corta)
        
        if theta(i,j) > max_theta_M1_r(i,j)
            max_theta(j) = max_theta_M1_r(i-1,j);
            %Los valores en esa posición (ya calculados) no se tienen en
            %cuenta: se le da un valor de 0
            beta_i(i,j) = 0;
            theta(i,j) = 0;
            M2_i(i,j) = 0;
            M1_r(i,j) = 0;
            max_theta_M1_r(i,j) = 0;
            break;
        end
        %Se elige para la primera iteración como guess inicial beta_i,
        %y para las sucesivas iteraciones los beta_r previamente
        %calculados.
        if i==1
            beta_r_sol(i,j) = vpasolve(fun1(theta(i,j),beta_r,M1_r(i,j)),beta_r,beta_i(i,j));
        end
        if i>1 
            beta_r_sol(i,j) = vpasolve(fun1(theta(i,j),beta_r,M1_r(i,j)),beta_r,beta_r_sol(i-1,j));
        end
        i=i+1;
        beta_i(i,j)=beta_i(i-1,j)+1;
    end   
    j=j+1;
end


beta_i_max=max(beta_i);
beta_i_max_teorica=beta_i_max(length(beta_i_max));

%Para gamma=1.4
% beta_i(21,4)=NaN;
% %Para gamma=1.2857
% beta_i(,)=NaN;
% %Para gamma=1.66
% beta_i(,)=NaN;


% Acondicionamiento de las matrices resultado: mismo tamaño y teniendo en
% cuenta solo valores que cumplen con las ecuaciones (quitar los últimos
% valores de las columnas de las matrices beta_i, theta, M2_i, M1_r y
% max_theta_M1_r)

[fil,col]=size(beta_r_sol);
beta_i=beta_i(1:fil,:);
theta=real(theta(1:fil,:));
M2_i=real(M2_i(1:fil,:));
M1_r=real(M1_r(1:fil,:));
max_theta_M1_r=double(max_theta_M1_r(1:fil,:));
beta_r_sol=double(beta_r_sol); % Si no se pone double sale en formato symbol

% Acondicionamiento de las matrices resultado: se sustituyen los O´s por
% NaN para así poder representar las columnas gráficamente

beta_i(beta_i==0)=NaN;
theta(theta==0)=NaN;
M2_i(M2_i==0)=NaN;
M1_r(M1_r==0)=NaN;
max_theta_M1_r(max_theta_M1_r==0)=NaN;
beta_r_sol(beta_r_sol==0)=NaN;

% Representación gráfica de los reslutados:
gamma;
beta_i;
theta;
beta_r_sol;
beta_i_max_teorica;

figure()
plot(beta_i(:,:),theta(:,:),'g',beta_i(:,:),beta_r_sol(:,:),'b');
hold on
plot([beta_i_max_teorica,beta_i_max_teorica],[0,80],'--k');
xlabel(texlabel('beta_i [º]'));
text(10,30,'[º]','Color','g','FontSize',12);
text(10,25,'[º]','Color','b','FontSize',12);
text(10,50,'$\theta$','interpreter','latex','Color','g','FontSize',12);
text(10,40,'$\beta_r$','interpreter','latex','Color','b','FontSize',12);
axis([0 60 0 70]);
text(5,65, sprintf('\\gamma = %.1f', gamma));
text(5,60, sprintf('\\beta_{i}^{max} = %.2f', beta_i_max_teorica));
text(50,18,'$M_{1}\uparrow$','interpreter','latex','Color','g','FontSize',12);
text(45,50,'$M_{1}\uparrow$','interpreter','latex','Color','b','FontSize',12);
text(50,65,'$M_{1}=1.2$','interpreter','latex','Color','g','FontSize',12);
text(50,61,'$M_{1}=1.5$','interpreter','latex','Color','g','FontSize',12);
text(50,58,'$M_{1}=2$','interpreter','latex','Color','g','FontSize',12);
text(50,55,'$M_{1}=3$','interpreter','latex','Color','g','FontSize',12);
text(50,53,'$M_{1}=4$','interpreter','latex','Color','g','FontSize',12);
text(50,50,'$M_{1}=5$','interpreter','latex','Color','g','FontSize',12);
text(50,46,'$M_{1}=6$','interpreter','latex','Color','g','FontSize',12);
text(50,41,'$M_{1}=3$','interpreter','latex','Color','b','FontSize',12);
text(50,37,'$M_{1}=4$','interpreter','latex','Color','b','FontSize',12);
text(50,32,'$M_{1}=5$','interpreter','latex','Color','b','FontSize',12);
text(50,27,'$M_{1}=6$','interpreter','latex','Color','b','FontSize',12);


figure()
plot(beta_i_max,M1_i);
hold on
plot([beta_i_max_teorica,beta_i_max_teorica],[1,20],'--k');
xlabel(texlabel('beta_i^{max} [º]'));
ylabel('${M_1}$','interpreter','latex','Rotation',0);
axis([0 90 1 20]);
text(7,5.75, sprintf('\\gamma = %.1f', gamma));
text(7,5.35, sprintf('\\beta_{i}^{max} = %.2f', beta_i_max_teorica));





