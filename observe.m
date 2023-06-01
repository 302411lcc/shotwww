function observe_y=observe(T)

theta=[1,1];
y=[];
for k =1: length(T)
    y(k,:) =myode1(theta,T(k)); % evaluate model at initial guess
end
observe_y=y+normrnd(0,0.5,length(T),2);

