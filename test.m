N = 5000;
M = 20;
lam = 1;
z = create_disturbance(N,gain);
x0 = data(:,1,1);
x = z + x0; %Add to one ECG trace
h = lmfir(@sincos,2,M,M,lam.^(2*M:-1:0)); %create filter
zhat = filter(h,1,x); % Filter out disturbance
xhat = x-zhat;
figure(1)
plot(xhat) ; %the cleaned signal
figure(2)
plot(xhat-x0) % the remaining error 
norm(xhat(1000:end)-x0(1000:end))/sqrt((N-1000))