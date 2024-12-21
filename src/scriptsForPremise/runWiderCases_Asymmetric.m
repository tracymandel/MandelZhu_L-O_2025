% runWiderCases_Asymmetric.m
% Run blade motion simulation for some cases that yield greater asymmetry

clear all; close all; clc
addpath('./originalCode_2DbladeMotion_JGROceans/')
addpath('./shadingFunctions/')

fprintf('Starting code... \n')

load('./Data/bladeMotionZhu.mat')

fprintf('Loaded Zhu motion dataset. \n')

%% Select 40-50 other cases that should yield high asymmetry...
% low E
% large H
% large l
% deep water

thresh = 0.1;   % mean deflection in x-direction > 10% of blade length
mean_asym_inds = [];

for i=1:length(zhuMotion)

    x_mean = mean(zhuMotion(i).x,2);
    if max(abs(x_mean)) > thresh
        mean_asym_inds = [mean_asym_inds, i];
    end

end

nAsymCases = length(mean_asym_inds);
nCases = 50;
indsRandom = randperm(nAsymCases,nCases);

%% Set up cases

h = 1;          % [m] water depth 1 m
d = 0.25e-3;    % [m] blade thickness 0.25 mm
b = 0.01;       % [m] blade width 1 cm
rho_w = 1025;   % fluid density
rho_v = 900;    % blade density
g = 9.81;       % gravity

% Using exact same parameters as Zhu asymmetric cases
for i=1:nCases

    asymCases(i).h = h;
    asymCases(i).l = zhuMotion(indsRandom(i)).l; % use same l, H, T, E as Zhu dataset cases
    asymCases(i).H = zhuMotion(indsRandom(i)).H;
    asymCases(i).T = zhuMotion(indsRandom(i)).T;
    asymCases(i).E = zhuMotion(indsRandom(i)).E;
    asymCases(i).b = b;

end

fprintf('Set up cases. Now running blade motion simulation... \n')

%% Run blade motion simulation

% call function with Longhuan's code
thick = d;
width = b;
rho_structure = rho_v;
rho_water = rho_w;

Delta = 1/sqrt(1700);   % Use just one spacing?

resultsFolder = './Data';
validationCase = 'caseStudyAsym';

num_T_period = 60;
ns = 50;
Dt = 0.01;


% % Run cases that vary about base case
% parfor ii = 1:nCases
% 
%     fprintf('Simulating motion for case %i of %i... \n',ii,nCases)
% 
%     dataFileName = ['motionData_',num2str(ii)];
% 
%     mainFunctionLinearWaveZhu(asymCases(ii).l,width,thick,rho_structure,...
%         asymCases(ii).E,rho_water,h,asymCases(ii).H,asymCases(ii).T,...
%         resultsFolder,validationCase,dataFileName,num_T_period,ns,Dt)
% 
% end



%% Compute nondimensional numbers and save x/l and z/l

% Using index ii because other variables may be part of the saved
% workspace of the .mat files loaded...
for ii=1:nCases

    fprintf('Computing ND #''s for case %i of %i... \n',ii,nCases)

    dataFileName = ['motionData_',num2str(ii)];
    load([resultsFolder,'/',validationCase,'/',dataFileName,'.mat'])

    % Velocity and second moment of inertia
    Uw = mean(Umax(1,:));       % avg velocity at base
    I = asymCases(ii).b*d^3/12;

    asymCases(ii).KC = Uw*asymCases(ii).T/asymCases(ii).b;
    asymCases(ii).Ca = rho_w * asymCases(ii).b * Uw^2 * ...
        asymCases(ii).l^3 / (asymCases(ii).E * I);
    asymCases(ii).L = asymCases(ii).l*2*pi/(asymCases(ii).T*Uw);
    asymCases(ii).B = (rho_w-rho_v)*g*asymCases(ii).b*d*...
        asymCases(ii).l^3 / (asymCases(ii).E * I);

    % Update these variables!! (whatever is saved to .mat file)
    asymCases(ii).x = xPos/asymCases(ii).l;
    asymCases(ii).z = zPos/asymCases(ii).l;

    asymCases(ii).x_dimensional = xPos;
    asymCases(ii).z_dimensional = zPos;

    asymCases(ii).t = tPos;

    asymCases(ii).nPeriods = num_T_period;
    asymCases(ii).ns = 50;
    asymCases(ii).dt = Dt;

    asymCases(ii).Delta = Delta;

    % keyboard

end

%% Compute angles
for i=1:nCases
    x = asymCases(i).x;
    z = asymCases(i).z;

    [ns,nt] = size(x);

    for j=1:nt

        xt = x(:,j);    % x at a given timestep
        zt = z(:,j);    % z at a given timestep
        
        theta(1) = 0;   % assume clamped at base?
        % Estimate theta_i with backward finite difference based on x,z
        theta(2:ns) = atand( (xt(2:end)-xt(1:end-1))./(zt(2:end)-zt(1:end-1)) );
        asymCases(i).theta(:,j) = theta; % in degrees
        
    end
end


%% Run shading

for ii=1:nCases

    fprintf('Computing shading for case %i of %i... \n',ii,nCases)

    dataFileName = ['motionData_',num2str(ii)];
    load([resultsFolder,'/',validationCase,'/',dataFileName,'.mat'])
    clear length

    % Blade position normalized by length
    X = asymCases(ii).x;
    Y = asymCases(ii).z;
    [ns,nt] = size(X);

    % Spacing normalized by length (R = dS/l)
    dS = asymCases(ii).Delta/asymCases(ii).l;

    % Initialize / set to zero before beginning each case
    nhbrShadedPts = zeros(ns,nt);
    selfShadedPts = zeros(ns,nt);
    shadedPts = zeros(ns,nt);
    fracShaded = zeros(1,nt);

    for tc = 1:nt

        x = X(:,tc);
        y = Y(:,tc);

        % figure(1); clf
        % plot(x,y,'k-')
        % axis([-1 1 0 1])
        % drawnow

        [indSelf,indNhbr,indTotl] = findTotalShading(x,y,dS);

        nhbrShadedPts(:,tc) = indNhbr;
        selfShadedPts(:,tc) = indSelf;
        shadedPts(:,tc) = indTotl;
        fracShaded(tc) = sum(indTotl)/ns;

    end

    % Save for each case

    % Fully time and space resolved shaded points
    asymCases(ii).shadedPts = shadedPts;
    asymCases(ii).nhbrShadedPts = nhbrShadedPts;
    asymCases(ii).selfShadedPts = selfShadedPts;
    asymCases(ii).fracShaded = fracShaded;

    % Time and length averaged (self, neighbor, total)

    asymCases(ii).avgShaded = mean(mean(shadedPts,2));
    asymCases(ii).avgUnshaded = 1-mean(mean(shadedPts,2));
    asymCases(ii).avgSelfShaded = mean(mean(selfShadedPts,2));
    asymCases(ii).avgNhbrShaded = mean(mean(nhbrShadedPts,2));

    % Phase-averaged, resolved w/depth (same as raw data)
    asymCases(ii).shadePhaseAvg = shadedPts; % ns x T array

    % Save associated non-dimensional numbers
    % asymCases(ii).KC = asymCases(ii).KC;
    % asymCases(ii).Ca = asymCases(ii).Ca;
    % asymCases(ii).L = asymCases(ii).L;
    asymCases(ii).R = dS;     % dS/l
    % asymCases(ii).B = asymCases(ii).B;

    % Other data...
    asymCases(ii).motionIndex = ii;  % associated case index in blade motion structure array
    asymCases(ii).nt = nt;      % length of time vector
    asymCases(ii).ns = ns;

    % keyboard

end

%% Save data
save([resultsFolder,'/',validationCase,'/caseStudyAsym_MotionAndShading.mat'],...
    'asymCases','-v7.3')

fprintf('Sucessfully saved! \n')