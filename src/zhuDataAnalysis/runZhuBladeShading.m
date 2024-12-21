% runZhuBladeShading.m
% Uses Zhu et al. (2020) plant motion dataset and simulates light
% exposure/shading for a range of plant spacing
%
% Mandel & Zhu (2025) L&O
% (c) Tracy Mandel | tracy.mandel@unh.edu
% Last updated 2024/12/21

clear all; close all; clc

load('../../data/zhuData/bladeMotionZhu.mat')

%% read data
h = 1; % [m] water depth 1 m
d = 0.25e-3; % [m] blade thickness 0.25 mm
b = 0.01; % [m] blade width 1 cm

% parameters are saved in dataRMPaper.mat
load('../../data/zhuData/dataRMPaper.mat', 'Ca','EE','H','KC','l','Re','T');
Umax = Re*1e-6/b;
KC_old = KC;
KC = Umax.*T/b;
% Ca is Ca number
% EE is young's modulus
% H is wave height
% KC is KC number
% l is blade length
% Re is Reynolds number
% T is wave period

%% Zhu Buoyancy parameter B
rho = 1000;     % fluid density
rhov = 900;     % blade density
g = 9.81;

d = .25e-3;     % blade thickness
b = 1e-2;       % blade with 
I = b*d^3/12;
B = (rho-rhov) * g * b * d * l.^3./(EE*I);

%% Run blade shading functions

% # of motion cases
nCases = length(H);

% Range of plant spacing
dS_vec = 0.05:.15:0.8;

ctr = 1;

% # of total cases (plant motion x # of plant spacing cases)
M = nCases*length(dS_vec);


for i=1:nCases
    
    fprintf('Plant motion case %i of %i (...%.2f%%) \n',i,nCases,i/nCases*100)
    % 
    blade.x = zhuMotion(i).x/l(i);  % dimensionless
    blade.z = zhuMotion(i).z/l(i);
    blade.t = zhuMotion(i).t;
    nt = length(blade.t);
    ns = length(blade.z(:,1));
    
    for j = 1:length(dS_vec)
        
        dS = dS_vec(j);

        fprintf('    (Overall case %i of %i)\n',ctr,M)

        % Initialize / set to zero before beginning each case
        nhbrShadedPts = zeros(ns,nt);
        selfShadedPts = zeros(ns,nt);
        shadedPts = zeros(ns,nt);
        fracShaded = zeros(1,nt);

        for tc = 1:nt

            x = blade.x(:,tc);
            y = blade.z(:,tc);

            [indSelf,indNhbr,indTotl] = findTotalShading(x,y,dS);

            nhbrShadedPts(:,tc) = indNhbr;
            selfShadedPts(:,tc) = indSelf;
            shadedPts(:,tc) = indTotl;
            fracShaded(tc) = sum(indTotl)/ns;

        end

        % Save for each case

        % Fully time and space resolved shaded points
        zhuShading(ctr).shadedPts = shadedPts;
        zhuShading(ctr).nhbrShadedPts = nhbrShadedPts;
        zhuShading(ctr).selfShadedPts = selfShadedPts;
        zhuShading(ctr).fracShaded = fracShaded;

        % Time and length averaged (self, neighbor, total)

        zhuShading(ctr).avgShaded = mean(mean(shadedPts,2));
        zhuShading(ctr).avgUnshaded = 1-mean(mean(shadedPts,2));
        zhuShading(ctr).avgSelfShaded = mean(mean(selfShadedPts,2));
        zhuShading(ctr).avgNhbrShaded = mean(mean(nhbrShadedPts,2));

        % Phase-averaged, resolved w/depth (same as raw data)
        zhuShading(ctr).shadePhaseAvg = shadedPts; % ns x T array
        
        % Save associated non-dimensional numbers
        zhuShading(ctr).KC = KC(i);
        zhuShading(ctr).Ca = Ca(i);
        zhuShading(ctr).L = l(i)/b/KC(i)*2*pi;
        zhuShading(ctr).R = dS;     % dS/l
        zhuShading(ctr).B = B(i);
        
        % Other data...
        zhuShading(ctr).motionIndex = i;  % associated case index in blade motion structure array
        zhuShading(ctr).nt = nt;      % length of time vector
        zhuShading(ctr).ns = ns;
        zhuShading(ctr).dt = dt;
        
        ctr = ctr+1;
        
%         fprintf('%.2f sec to process overall case %i of %i) \n',toc,ctr,nCases*length(R_vec))
        
    end
end

%% Save data
save('../../data/bladeShadingZhu.mat','zhuShading','-v7.3')
