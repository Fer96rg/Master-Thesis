clc
clear all 
% close all

load('curva_RH_def_100'); %Cargo los vectores T, R, alpha, P, M1n y M2n

%Diagrama polar de M1_determ
M1_determ=3;
% Para beta determinado:
beta_i=asind(M1_determ^-1):1:90;
for i=1:length(beta_i)
    M1(i,:)=M1n_def./sind(beta_i(i));
    theta(i,:)=atand((R_def-1)./(tand(beta_i(i))+R_def./tand(beta_i(i))));
    M2(i,:)=M2n_def./sind(beta_i(i)-theta(i,:));
end
%Calculo el indice de la matriz M1 que da M1_determ para cada beta (para cada fila)
for i=1:length(beta_i)
    for j=1:length(M1n_def)
        if M1(i,j)>M1_determ && M1(i,j)<(M1_determ+0.1)
            index_M1_determ(i)=j;   
            break;           
        end
    end
    
end
%Utilizo esos indices para calcular theta: calculo el theta de cada beta
%que da M1_determ
for i=1:length(index_M1_determ)
    theta_M1_determ(i)=theta(i,index_M1_determ(i));
end
%Utilizo esos indices para calcular P
for i=1:length(index_M1_determ)
    P_M1_determ(i)=P_def(index_M1_determ(i));
end
%Diagrama polar de M1_determ
figure()
plot(theta_M1_determ,P_M1_determ,'r');
hold on
plot(-theta_M1_determ,P_M1_determ,'r');
% hold on
%M=3
% plot([0,0],[0,30],'--k');
% axis([-60 60 1 15]); 
%M=10
% plot([0,0],[0,150],'--k');
axis([-60 60 0 70]); 
xlabel(texlabel('theta [º]'));
ylabel('$\mathcal{P}$','interpreter','latex','Rotation',0);
text(15,1.5,sprintf('$\\mathcal{M}=%.f$',M1_determ),'interpreter','latex');

%Para y M1=3 y beta_deter=30
beta_deter=30;
M1_inc_reb=M1(11,304);
theta_inc_reb=atand((R_def-1).*tand(beta_i(i))./(R_def+(tand(beta_i(i))).^2));
