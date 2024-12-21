% Analyze Enriquez et al. (2002) PFD data
% 1b: Data 15 cmab
% 1c: 10 cmab
% 1d: 5 cmab
% 1e: 0 cmab

% clear all; close all; clc
load('enriquez_2002_data/enriquez_fig1_data_cleaned.mat')

z = [15, 10, 5, 0];     % measurement height, cmab

%% Plot time series and extract mean, median, variance
figure(1); clf
for i=1:4
    subplot(2,2,i)
    plot(data(i).t,data(i).PFD,'k.')
    xlabel('time (sec)')
    ylabel('PFD')
    
    PFD_mean(i) = mean(data(i).PFD);
    PFD_med(i) = median(data(i).PFD);
    PFD_var(i) = var(data(i).PFD);
end

%% Plot quantities over depth
figure(2); clf
% subplot(1,2,1)
errorbar(PFD_mean,z,[],[],-sqrt(PFD_var),sqrt(PFD_var),'k.','linewidth',1.5)
hold on

plot(PFD_mean,z,'bo--','markerfacecolor','b','linewidth',1.5)
hold on
plot(PFD_med,z,'rs:','markerfacecolor','r','linewidth',1.5)

set(gca,'fontsize',14)
xlabel('PFD')
ylabel('z (cmab)')

ylim([0 20])

% 
% subplot(1,2,2)
% plot(PFD_var,z,'bo-')

%% What is the PFD at the canopy top (z = 20 cm?)
% p. 893: "80% of the incident light was attenuated in the first 5 cm of
% the canopy"
%       --> PFD_mean(1) = 0.8 * PFD_top
% "only 1% of the irradiance reached the water-sediment interface"
%       --> PFD_mean(4) = 0.01 * PFD_top

PFD_top1 = PFD_mean(1)/.8
PFD_top2 = PFD_mean(4)/.01
