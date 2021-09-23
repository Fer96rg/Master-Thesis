function [max_theta] = fun_max_theta(M1,gamma)
beta=0:1:90;
N=length(beta);
for i=1:N
    theta(i)=atand((M1^2*(sind(beta(i)))^2-1)*2*cotd(beta(i))/(gamma*M1^2+M1^2*cosd(2*beta(i))+2));
end
max_theta=max(theta);
