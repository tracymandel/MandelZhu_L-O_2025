% Load and clean digitized data from Enriquez et al. (2002) Fig 1(b-e)
% used https://apps.automeris.io/wpd/ for automated digitization
% 1b: Data 15 cmab
% 1c: 10 cmab
% 1d: 5 cmab
% 1e: 0 cmab

clear all; close all; clc

sub = 'bcde';

data(4) = struct();

for i=1:4
    temp = readtable('enrique_fig1_data.xlsx','Sheet',['Fig 1',sub(i)]);
    data(i).t = table2array(temp(:,1));
    data(i).PFD = table2array(temp(:,2));
end


%%
figure(1); clf
for i=1:4
    subplot(2,2,i)
    plot(data(i).t,data(i).PFD,'k.')
end

% Then need to remove artefacts from plot digitization 
%% Removing spurious points from Fig 1b (top left)

% Ticks at top of figure along horizontal axis
data(1).t(594) = [];
data(1).PFD(594)= [];

data(1).t(483) = [];
data(1).PFD(483)= [];

data(1).t(340) = [];
data(1).PFD(340)= [];

data(1).t(232) = [];
data(1).PFD(232)= [];

data(1).t(120) = [];
data(1).PFD(120)= [];

% from figure label "(b)"
inds = (data(1).t > 54 & data(1).t < 59) & (data(1).PFD > 2150 & data(1).PFD < 2400);
data(1).t(inds) = [];
data(1).PFD(inds) = [];

% remove first and last few second of data to eliminate vertical-axis ticks
inds = data(1).t > 1 & data(1).t < 58;

data(1).t = data(1).t(inds);
data(1).PFD = data(1).PFD(inds);

subplot(2,2,1)
hold on
plot(data(1).t,data(1).PFD,'r.')

%% Removing spurious pts from Fig 1c (top right)
inds = data(2).PFD < 123 & data(2).PFD > 112;
data(2).t(inds) = [];
data(2).PFD(inds) = [];

% remove first and last few second of data to eliminate vertical-axis ticks
inds = data(2).t > 1.5 & data(2).t < 58;

data(2).t = data(2).t(inds);
data(2).PFD = data(2).PFD(inds);

data(2).t(774) = [];
data(2).PFD(774) = [];

subplot(2,2,2)
hold on
plot(data(2).t,data(2).PFD,'r.')


%% Removing spurious pts from Fig 1d (bottom left)
inds = data(3).PFD < 124 & data(3).PFD > 120;
data(3).t(inds) = [];
data(3).PFD(inds) = [];

% from figure label "(d)"
inds = (data(3).t > 9 & data(3).t < 60) & (data(3).PFD > 520 & data(3).PFD < 600);
data(3).t(inds) = [];
data(3).PFD(inds) = [];

% remove first and last few second of data to eliminate vertical-axis ticks
inds = data(3).t > 1.5 & data(3).t < 58;

data(3).t = data(3).t(inds);
data(3).PFD = data(3).PFD(inds);


subplot(2,2,3)
hold on
plot(data(3).t,data(3).PFD,'r.')


%% Removing spurious pts from Fig 1e (bottom right)
inds = data(4).PFD < 122 & data(4).PFD > 119;
data(4).t(inds) = [];
data(4).PFD(inds) = [];

% remove first and last few second of data to eliminate vertical-axis ticks
inds = data(4).t > 1.5 & data(4).t < 58;

data(4).t = data(4).t(inds);
data(4).PFD = data(4).PFD(inds);

subplot(2,2,4)
hold on
plot(data(4).t,data(4).PFD,'r.')

%% Sort data in time order
for i=1:4
    [~,I] = sort(data(i).t);
    data(i).t = data(i).t(I);
    data(i).PFD = data(i).PFD(I);
end

%% Save modified data
save('enriquez_fig1_data_cleaned.mat','data')