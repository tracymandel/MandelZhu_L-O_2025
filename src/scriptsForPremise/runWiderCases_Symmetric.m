% runValidationCases_widerRange.m
% Run blade motion simulation and shading for the Enriquez "base case" plus
% some variation about that base case.

clear all; close all; clc
addpath('./originalCode_2DbladeMotion_JGROceans/')
addpath('./shadingFunctions/')

fprintf('Starting code... \n')

%% Base case
baseCase = struct();

% baseCase.H = 0.4*2/3;
baseCase.H = 0.4;
baseCase.T = 7;
baseCase.l = 0.20;
baseCase.h = 5;
baseCase.b = 0.01;
baseCase.E = 2400e6;
baseCase.Delta = 1/sqrt(1700);

nCases = 25;

for i=1:nCases
    symCases(i) = baseCase;
end

% Index 1: Base case

% Index 2&5: primary varying H cases for figs
symCases(2).H = 0.4 - 0.2;
symCases(3).H = 0.4 - 0.1;
symCases(4).H = 0.4 + 0.1;
symCases(5).H = 0.4 + 0.2;

% Index 6&9: varying T
symCases(6).T = 7 - 3;
symCases(7).T = 7 - 1.5;
symCases(8).T = 7 + 1.5;
symCases(9).T = 7 + 3;

% Index 10&13: varying l
symCases(10).l = 0.20 - 0.1;
symCases(11).l = 0.20 - 0.05;
symCases(12).l = 0.20 + 0.05;
symCases(13).l = 0.20 + 0.1;

% Index 14&17: varying E
symCases(14).E = 2400e6 - 2000e6;
symCases(15).E = 2400e6 - 1000e6;
symCases(16).E = 2400e6 + 1000e6;
symCases(17).E = 2400e6 + 2000e6;

% Index 18&21: varying spacing Delta
symCases(18).Delta = 1/sqrt(1700) - 0.02;
symCases(19).Delta = 1/sqrt(1700) - 0.01;
symCases(20).Delta = 1/sqrt(1700) + 0.01;
symCases(21).Delta = 1/sqrt(1700) + 0.02;

% Index 25&25: varying water depth h
symCases(22).h = 5 - 3;
symCases(23).h = 5 - 1.5;
symCases(24).h = 5 + 1.5;
symCases(25).h = 5 + 3;


fprintf('Set up cases. Now running blade motion simulation... \n')

nCases = length(symCases);

%% Run blade motion simulation
% keep fixed...
g = 9.81;                   % accel. due to gravity (m/s^2)
d = 0.0003;                 % blade thickness (m)
rho_v = 940;                % density of T. testudinum (kg/m^3)
rho_w = 1025;               % density of (salt) water (kg/m^3)

% call function with Longhuan's code
thick = d;
width = 0.01;
rho_structure = rho_v;
rho_water = rho_w;

resultsFolder = './Data';
validationCase = 'caseStudySym';

num_T_period = 60;
ns = 50;
Dt = 0.01;


% Run cases that vary about base case
parfor ii = 1:nCases

    fprintf('Simulating motion case %i of %i... \n',ii,nCases)

    dataFileName = ['motionData_',num2str(ii)];

    mainFunctionLinearWaveZhu(symCases(ii).l,width,thick,rho_structure,...
        symCases(ii).E,rho_water,symCases(ii).h,symCases(ii).H,symCases(ii).T,...
        resultsFolder,validationCase,dataFileName,num_T_period,ns,Dt)

end



%% Compute nondimensional numbers and save x/l and z/l

% Using index ii because other variables may be part of the saved
% workspace of the .mat files loaded...
for ii=1:nCases

    fprintf('Computing ND #''s for case %i of %i... \n',ii,nCases)

    dataFileName = ['motionData_',num2str(ii)];
    load([resultsFolder,'/',validationCase,'/',dataFileName,'.mat'])

    % Old definition of wave velocity
    % Uw_old = mean(mean(Umax));
    Uw = mean(Umax(1,:));       % avg velocity at base
    I = symCases(ii).b*d^3/12;

    symCases(ii).KC = Uw*symCases(ii).T/symCases(ii).b;
    symCases(ii).Ca = rho_w * symCases(ii).b * Uw^2 * ...
        symCases(ii).l^3 / (symCases(ii).E * I);
    symCases(ii).L = symCases(ii).l*2*pi/(symCases(ii).T*Uw);
    symCases(ii).B = (rho_w-rho_v)*g*symCases(ii).b*d*...
        symCases(ii).l^3 / (symCases(ii).E * I);

    % Update these variables!! (whatever is saved to .mat file)
    symCases(ii).x = xPos/symCases(ii).l;
    symCases(ii).z = zPos/symCases(ii).l;
    symCases(ii).t = tPos;

    symCases(ii).nPeriods = num_T_period;
    symCases(ii).ns = 50;
    symCases(ii).dt = Dt;

    % keyboard

end

%% Compute angles
for i=1:nCases
    x = symCases(i).x;
    z = symCases(i).z;

    [ns,nt] = size(x);

    for j=1:nt

        xt = x(:,j);    % x at a given timestep
        zt = z(:,j);    % z at a given timestep
        
        theta(1) = 0;   % assume clamped at base?
        % Estimate theta_i with backward finite difference based on x,z
        theta(2:ns) = atand( (xt(2:end)-xt(1:end-1))./(zt(2:end)-zt(1:end-1)) );
        symCases(i).theta(:,j) = theta; % in degrees
        
    end
end


%% Run shading

for ii=1:nCases

    fprintf('Computing shading for case %i of %i... \n',ii,nCases)

    dataFileName = ['motionData_',num2str(ii)];
    load([resultsFolder,'/',validationCase,'/',dataFileName,'.mat'])
    clear length

    % Blade position normalized by length
    X = symCases(ii).x;
    Y = symCases(ii).z;
    [ns,nt] = size(X);

    % Spacing normalized by length (R = dS/l)
    dS = symCases(ii).Delta/symCases(ii).l;

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
    symCases(ii).shadedPts = shadedPts;
    symCases(ii).nhbrShadedPts = nhbrShadedPts;
    symCases(ii).selfShadedPts = selfShadedPts;
    symCases(ii).fracShaded = fracShaded;

    % Time and length averaged (self, neighbor, total)

    symCases(ii).avgShaded = mean(mean(shadedPts,2));
    symCases(ii).avgUnshaded = 1-mean(mean(shadedPts,2));
    symCases(ii).avgSelfShaded = mean(mean(selfShadedPts,2));
    symCases(ii).avgNhbrShaded = mean(mean(nhbrShadedPts,2));

    % Phase-averaged, resolved w/depth (same as raw data)
    symCases(ii).shadePhaseAvg = shadedPts; % ns x T array

    % Save associated non-dimensional numbers
    % varyingCases(ii).KC = varyingCases(ii).KC;
    % varyingCases(ii).Ca = varyingCases(ii).Ca;
    % varyingCases(ii).L = varyingCases(ii).L;
    symCases(ii).R = dS;     % dS/l
    % varyingCases(ii).B = varyingCases(ii).B;

    % Other data...
    symCases(ii).motionIndex = ii;  % associated case index in blade motion structure array
    symCases(ii).nt = nt;      % length of time vector
    symCases(ii).ns = ns;

    % keyboard

end


%% Save data
save([resultsFolder,'/',validationCase,'/caseStudySym_MotionAndShading.mat'],...
    'symCases','-v7.3')

fprintf('Sucessfully saved! \n')