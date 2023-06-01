%%无噪音下
clear
clc
global csimga ny
np=1;T=1;csimga=0;ny=10;
cs_theta=unifrnd(-5,5,[ny,np]);%y
fid = fopen('theta_data','w');
for i = 1:size(cs_theta,1)
    fprintf(fid,'%g\t',cs_theta(i,:));
    fprintf(fid,'\n');
end
cs_y=[];
for k=1:ny
    theta=[1,cs_theta(k)];
cs_y(k,:)=myode1(theta,T)+normrnd(0,csimga);%x
end
fid = fopen('z_data','w');
for i = 1:size(cs_y,1)
    fprintf(fid,'%g\t',cs_y(i,:));
    fprintf(fid,'\n');
end
%%

x_star=observe(T);
%对待theta求采样和求分布
%由cs数据求核函数

%[x,y]=fmincon(@fun_cs,[2,1])
simga=1;%x(1);
L=sqrt(5);
%simga=x(1);
%L=x(2);
%%
K=zeros(ny,ny);
for k=1:ny
    for j=k:ny
        tt=norm((cs_y(k,:)-cs_y(j,:)),2);
        K(k,j)=simga^2*exp(-tt/2/(L^2));
        K(j,k)=K(k,j);
    end
end

K=K+diag(csimga^2*ones(ny,1));
%%

K_star=zeros(1,ny);
for k=1:ny
    tt=norm((x_star-cs_y(k,:)),2);
    K_star(k)=simga^2*exp(-tt/2/(L^2));
end
%%
 K_star2=simga^2+csimga^2;
 
 ansy=K_star*inv(K)*cs_theta
 var=K_star2-K_star*inv(K)*K_star'
 %%
 figure
 hdata=normrnd(ansy,sqrt(var),100,1);
 h=histogram(hdata,'BinWidth',0.1,'Normalization','pdf');
 t=find(h.Values==max(h.Values));
    y=h.BinEdges(t)
 title('无噪音下b的预测值')