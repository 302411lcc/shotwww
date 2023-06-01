function [y]=myode1(theta,T)
a=theta(1);
b=theta(2);
h=0.005;
N=T/h;
y0 = [2.5,2.5];
% odefun = @(t,y)[a*y(1)-y(1)*y(2);
%    b*y(1)*y(2)-y(2)];
% [~,y1] = ode23s(odefun,tspan,y0);
%     yvec = y1(end,:);

f = @(y)[a*y(1)-y(1)*y(2);
   b*y(1)*y(2)-y(2) ];

yy(:,1)=y0;
for m=2:N+1
    yi= yy(:,m-1);
    K1=f(yi);
    K2=f(yi+h*K1/2);
    K3=f(yi+h*K2/2);
    K4=f(yi+h*K3);
    yy(:,m)=yi+(h/6)*(K1+2*K2+2*K3+K4);
end
y=yy(:,N+1);
end