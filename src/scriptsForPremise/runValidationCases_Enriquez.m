% runValidationCases_Enriquez.m
% Enriquez et al. (2002) cases
% See "Validation field data ranges.xlsm"
clear all; close all; clc

addpath('./originalCode_2DbladeMotion_JGROceans/')
addpath('./shadingFunctions/')

fprintf('Starting code... \n') 

%%
% Keep vegetation parameters fixed
ll = 0.20;      % leaf length, m
b = 0.01;       % blade width, m
EE = 2.4e9;     % elastic modulus, Pa

% Vary only hydrodynamic parameters (more uncertain)
H = [0.1214, 0.4, 1.1066];  % upper and lower 1/10, plus most frequent wave height
TT = [5, 7, 9];              % wave period, sec
h = [4.85, 5, 5.15];        % water depth, m

% keep fixed...
g = 9.81;                   % accel. due to gravity (m/s^2)
d = 0.0003;                 % blade thickness (m)
rho_v = 940;                % density of T. testudinum (kg/m^3)
rho_w = 1025;               % density of (salt) water (kg/m^3)

nVariables = 3;             % # parameters varying
nCases = 3^nVariables;      % # validation cases

% Plant density later for shading...
Nv = 1700;
Delta = sqrt(1./Nv);       % if unstaggered

%% Combine into 1D array of all combinations...

enriquezCases(nCases) = struct();

ctr = 1;

for i=1:3
    for j=1:3
        for k=1:3
            enriquezCases(ctr).h = h(i);
            enriquezCases(ctr).H = H(j);
            enriquezCases(ctr).T = TT(k);
            enriquezCases(ctr).l = ll;
            enriquezCases(ctr).b = b;
            enriquezCases(ctr).E = EE;
            enriquezCases(ctr).KC = NaN;
            enriquezCases(ctr).L = NaN;
            enriquezCases(ctr).Ca = NaN;
            enriquezCases(ctr).B = NaN;
            ctr = ctr+1;
        end
    end
end

fprintf('Set up cases. Now running blade motion simulation... \n')

%% Call function that contains Longhuan's code

num_T_period = 8;
ns = 50;
Dt = 0.01;

resultsFolder = './Data';
validationCase = 'enriquezValidation';

% parfor ii=1:nCases
% 
%     fprintf('Computing motion for case %i of %i... \n',ii,nCases)
% 
%     Depth_water = enriquezCases(ii).h;
%     Height_wave = enriquezCases(ii).H;
%     T_wave = enriquezCases(ii).T;
% 
%     dataFileName = ['motionData_',num2str(ii)];
% 
%     mainFunctionLinearWaveZhu(ll,b,d,rho_v,EE,rho_w,Depth_water,Height_wave,T_wave,...
%         resultsFolder,validationCase,dataFileName,...
%         num_T_period,ns,Dt)
% 
% end

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
    I = b*d^3/12;

    enriquezCases(ii).KC = Uw*enriquezCases(ii).T/b;
    enriquezCases(ii).Ca = rho_w * b * Uw^2 * ll^3 / (EE * I);
    enriquezCases(ii).L = ll*2*pi/(enriquezCases(ii).T*Uw);
    enriquezCases(ii).B = (rho_w-rho_v)*g*b*d*ll^3 / (EE * I);

    enriquezCases(ii).x = xPos/ll;
    enriquezCases(ii).z = zPos/ll;
    enriquezCases(ii).t = tPos/enriquezCases(ii).T;

    enriquezCases(ii).nPeriods = num_T_period;
    enriquezCases(ii).ns = ns;
    enriquezCases(ii).dt = Dt;

end


%% Adding angle theta...

for i=1:length(enriquezCases)
    x = enriquezCases(i).x;
    z = enriquezCases(i).z;

    [ns,nt] = size(x);

    for j=1:nt

        xt = x(:,j);    % x at a given timestep
        zt = z(:,j);    % z at a given timestep
        
        theta(1) = 0;   % assume clamped at base?
        % Estimate theta_i with backward finite difference based on x,z
        theta(2:ns) = atand( (xt(2:end)-xt(1:end-1))./(zt(2:end)-zt(1:end-1)) );
        enriquezCases(i).theta(:,j) = theta; % in degrees
        
    end
end


%% Run shading on Enriquez cases

for ii=1:nCases

    fprintf('Computing shading for case %i of %i... \n',ii,nCases)

    dataFileName = ['motionData_',num2str(ii)];
    load([resultsFolder,'/',validationCase,'/',dataFileName,'.mat'])
    clear length

    % Blade position normalized by length
    X = enriquezCases(ii).x;
    Y = enriquezCases(ii).z;
    [ns,nt] = size(X);

    % Spacing normalized by length (R = dS/l)
    dS = Delta/enriquezCases(ii).l;

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
    enriquezCases(ii).shadedPts = shadedPts;
    enriquezCases(ii).nhbrShadedPts = nhbrShadedPts;
    enriquezCases(ii).selfShadedPts = selfShadedPts;
    enriquezCases(ii).fracShaded = fracShaded;

    % Time and length averaged (self, neighbor, total)

    enriquezCases(ii).avgShaded = mean(mean(shadedPts,2));
    enriquezCases(ii).avgUnshaded = 1-mean(mean(shadedPts,2));
    enriquezCases(ii).avgSelfShaded = mean(mean(selfShadedPts,2));
    enriquezCases(ii).avgNhbrShaded = mean(mean(nhbrShadedPts,2));

    % Phase-averaged, resolved w/depth (same as raw data)
    enriquezCases(ii).shadePhaseAvg = shadedPts; % ns x T array

    % Save associated non-dimensional numbers
    enriquezCases(ii).R = dS;     % dS/l

    % Other data...
    enriquezCases(ii).motionIndex = ii;  % associated case index in blade motion structure array
    enriquezCases(ii).nt = nt;      % length of time vector
    enriquezCases(ii).ns = ns;

    %         fprintf('%.2f sec to process overall case %i of %i) \n',toc,ctr,nCases*length(R_vec))

end


%% Save results
save([resultsFolder,'/',validationCase,'/','enriquezValidation_MotionAndShading.mat'],...
    'enriquezCases','-v7.3')

fprintf('Sucessfully saved! \n')