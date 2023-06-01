%生成观测值
clear
clc
h=0.05;
T=[10;20;30;40;50;60;70;80].*h;
observe_y=observe(T);
%%
%由observe值获得信息
meanD =observe_y ;delta=0.5;
A=diag((1/delta)*ones(2,1));


np=2;
theta0=unifrnd(-5,5,[np,1]);
theta(1,:)=theta0;
for k =1: length(T)
    yy(k,:) =myode1(theta0,T(k)); % evaluate model at initial guess
end
chi1 = (meanD-yy);
%%
accept_prob=0;
MClength=1e4;
simga=1;
for n = 1:MClength
    
    theta2(1)=normrnd(theta(n,1),simga);
    theta2(2)=normrnd(theta(n,2),simga);
    
    
    while theta2(1) > 5 || theta2(1) < -5||theta2(2) > 5 || theta2(2) < -5
        theta2(1)=normrnd(theta(n,1),simga);
    theta2(2)=normrnd(theta(n,2),simga);
    end
    
    
    for k =1: length(T)
        yy(k,:) =myode1(theta2,T(k)); % evaluate model at initial guess
    end
    chi2 = (meanD-yy); % candidate sum-of-squares
    e_ratio=exp((diag(-chi2*A*chi2')+diag(chi1*A*chi1'))./2);
    %cumprod(e_ratio)
    ratio = min(1,e_ratio);
    
    if rand(length(T),1) < ratio
        theta(n+1,:) = theta2;
        chi1 = chi2;
        accept_prob = accept_prob + 1/MClength;
    else
        theta(n+1,:) = theta(n,:);
        
    end
    if mod(n,1e3)==0
        disp(['Iteration ', num2str(n)])
        disp(['Acceptance probability = ',num2str(accept_prob)]);
    end
end

total = size(theta,1);
%画图
figure
burn = round(0.2*total); % burn-in, number of iterates to throw away
% Histograms
names = {'a','b'};
for i = 1:2
    subplot(2,1,i)
    h=histogram(theta(burn:end,i),'BinWidth',0.05,'Normalization','pdf');
    t=find(h.Values==max(h.Values));
    y=h.BinEdges(t)
    xlim([-5 5]);
    yl = ylim;
    title(names{i})
    set(gca,'fontsize',14)
end

% a1=tabulate(theta(burn:end,1));
% t=find(a1(:,2)==max(a1(:,2)));
% a=a1(t,1)
% pro_a=a1(t,3)
% 
% a2=tabulate(theta(burn:end,2));
% t2=find(a2(:,3)==max(a2(:,3)));
% pro_b=a2(t,3)
% b=a2(t2,1)

