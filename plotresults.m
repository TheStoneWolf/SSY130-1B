% This script plots the HR estimation results for all athletes

load('Data/ath.mat');
load('Data/VisualHr.mat');
close all

fs = 500;

% Determine sizes
[N,Nleads,Nath] = size(data);

% Vector of visually determined HR
%VisualHR = 58*ones(Nath,1); % (This is not the correct value!)

results = zeros(Nath,Nleads);
for ath = 1:Nath
    for lead = 1:Nleads
        %figure(2)
        results(ath,lead) = fbpm(data(:,lead,ath),fs);
    end
end

figure(1)
% Plot esimate for all leads and all athletes
plot(kron(ones(Nleads,1),(1:Nath)'),results(:),'+','MarkerSize',14)
hold on
% Plot the Visual heart rate
plot(1:Nath,VisualHr,'or','MarkerSize',14);
plot(1:Nath,mean(results'),'d','MarkerSize',14);
plot(1:Nath,median(results'),'s','MarkerSize',14);
hold off
%RMS error for the median estimates over  all athletes
norm(VisualHr.'-median(results'))/sqrt(Nath)




