% figure_angleTheory.m
% Compute blade bending angles and test simple theory based on bending angle
% 
% Code to generate Figure 8 in paper
% Mandel & Zhu (2025) L&O
% (c) Tracy Mandel | tracy.mandel@unh.edu
% Last updated 2024/12/21


% Case study
load('../../data/caseStudyAsym/caseStudyAsym_MotionAndShading.mat')
load('../../data/caseStudySym/caseStudySym_MotionAndShading.mat')
fprintf('Loaded case study... \n')

% Validation
load('../../data/enriquezValidation/enriquezValidation_MotionAndShading.mat')
fprintf('Loaded validation... \n')

% Zhu dataset
load('../../data/zhuData/bladeMotionZhu.mat')
load('../../data/zhuData/bladeShadingZhu.mat')
fprintf('Loaded Zhu dataset... \n')


%% Decompose bending angle into mean and wave-varying component

for i=1:length(zhuMotion)
    
    theta_st = zhuMotion(i).theta;  % theta as a function of s and t

    ns = zhuMotion(i).ns;
    
    % Time-average
    theta_bar = mean(theta_st,2);
    theta_prime = theta_st - theta_bar;
    
    zhuMotion(i).theta_avg = theta_bar; % time-mean
    zhuMotion(i).theta_osc = theta_prime; % wave-varying
    
    zhuMotion(i).theta_osc_rms = sqrt(mean(theta_prime.^2,2));    % rms bending angle
    zhuMotion(i).theta_osc_max = max(abs(theta_prime),[],2);      % max amplitude of bending angle
    
    % At tip
    zhuMotion(i).theta_avg_tip = theta_bar(end);
    zhuMotion(i).theta_osc_tip = zhuMotion(i).theta_osc_max(end);
    
    % Averaged over the entire blade
    zhuMotion(i).theta_avg_blade = mean(theta_bar);
    zhuMotion(i).theta_osc_blade = mean(zhuMotion(i).theta_osc_max);
    
    % Average 1/theta at the blade tip
    I = mean(1./theta_st(ns,:),2);
    zhuMotion(i).Omega = I; % average of 1/theta
    

    % Average 1/theta along the entire blade??
    I = mean(mean(1./theta_st(2:end,:)));
    zhuMotion(i).Omega_blade = I; % average of 1/theta


    % Max asymmetry relative to blade length
    % --> x-displacement of mean tip position
    xmean = mean(zhuMotion(i).x,2);         % mean x-position
    zhuMotion(i).xmean_tip = xmean(end);
    zhuMotion(i).xmean_max = max(xmean);

    % Blade posture asymmetry from Longhuan's paper (see Zhu 2020 fig 8)
    %   beta_xt = mean(x_tip) / max(abs(x_tip))
    zhuMotion(i).beta = zhuMotion(i).xmean_tip / max(abs(zhuMotion(i).x(end,:)));
    %   beta based on angles...
    %   beta_theta_t = mean(theta) / max(abs(theta))
    zhuMotion(i).beta_tht = abs(theta_bar(end)) / max(abs(theta_st(end,:)));
    
end

%% Same for Enriquez validation
ns = 51;

for i=1:length(enriquezCases)

    theta_st = enriquezCases(i).theta;  % theta as a function of s and t
    % ns = enriquezCases(i).ns;
    
    % Time-average
    theta_bar = mean(theta_st,2);
    theta_prime = theta_st - theta_bar;
    
    enriquezCases(i).theta_avg = theta_bar; % time-mean
    enriquezCases(i).theta_osc = theta_prime; % wave-varying
    
    enriquezCases(i).theta_osc_rms = sqrt(mean(theta_prime.^2,2));    % rms bending angle
    enriquezCases(i).theta_osc_max = max(abs(theta_prime),[],2);      % max amplitude of bending angle
    
    % At tip
    enriquezCases(i).theta_avg_tip = theta_bar(end);
    enriquezCases(i).theta_osc_tip = enriquezCases(i).theta_osc_max(end);
    
    % Averaged over the entire blade
    enriquezCases(i).theta_avg_blade = mean(theta_bar);
    enriquezCases(i).theta_osc_blade = mean(enriquezCases(i).theta_osc_max);
    
    % Average 1/theta at the blade tip
    I = mean(1./theta_st(ns,:),2);
    enriquezCases(i).Omega = I; % average of 1/theta
    

    % Average 1/theta along the entire blade??
    I = mean(mean(1./theta_st(2:end,:)));
    enriquezCases(i).Omega_blade = I; % average of 1/theta



    % Max (mean) asymmetry relative to blade length
    % --> x-displacement of mean tip position
    xmean = mean(enriquezCases(i).x,2);
    enriquezCases(i).xmean_tip = xmean(end);
    enriquezCases(i).xmean_max = max(xmean);

    
    % Blade posture asymmetry from Longhuan's paper (see Zhu 2020 fig 8)
    %   beta_xt = mean(x_tip) / max(abs(x_tip))
    enriquezCases(i).beta = enriquezCases(i).xmean_tip / max(abs(enriquezCases(i).x(end,:)));

    %   beta based on angles...
    %   beta_theta_t = mean(theta) / max(abs(theta))
    enriquezCases(i).beta_tht = abs(theta_bar(end)) / max(abs(theta_st(end,:)));
    
end

%% Same for asymmetrical validation/case study cases
for i=1:length(asymCases)

    theta_st = asymCases(i).theta;  % theta as a function of s and t
    ns = asymCases(i).ns;
    
    % Time-average
    theta_bar = mean(theta_st,2);
    theta_prime = theta_st - theta_bar;
    
    asymCases(i).theta_avg = theta_bar; % time-mean
    asymCases(i).theta_osc = theta_prime; % wave-varying
    
    asymCases(i).theta_osc_rms = sqrt(mean(theta_prime.^2,2));    % rms bending angle
    asymCases(i).theta_osc_max = max(abs(theta_prime),[],2);      % max amplitude of bending angle
    
    % At tip
    asymCases(i).theta_avg_tip = theta_bar(end);
    asymCases(i).theta_osc_tip = asymCases(i).theta_osc_max(end);
    
    % Averaged over the entire blade
    asymCases(i).theta_avg_blade = mean(theta_bar);
    asymCases(i).theta_osc_blade = mean(asymCases(i).theta_osc_max);
    
    % Average 1/theta at the blade tip
    I = mean(1./theta_st(ns,:),2);
    asymCases(i).Omega = I; % average of 1/theta
    

    % Average 1/theta along the entire blade??
    I = mean(mean(1./theta_st(2:end,2:end))); % neglect first column
    asymCases(i).Omega_blade = I; % average of 1/theta


    % Max asymmetry relative to blade length
    % --> x-displacement of mean tip position
    xmean = mean(asymCases(i).x,2);
    asymCases(i).xmean_tip = xmean(end);
    asymCases(i).xmean_max = max(xmean);


    % Blade posture asymmetry from Longhuan's paper (see Zhu 2020 fig 8)
    %   beta_xt = mean(x_tip) / max(abs(x_tip))
    asymCases(i).beta = asymCases(i).xmean_tip / max(abs(asymCases(i).x(end,:)));

        %   beta based on angles...
    %   beta_theta_t = mean(theta) / max(abs(theta))
    asymCases(i).beta_tht = abs(theta_bar(end)) / max(abs(theta_st(end,:)));
    
end

%% Same for main case study cases
ns = 51;

for i=1:length(symCases)

    theta_st = symCases(i).theta;  % theta as a function of s and t
    % ns = enriquezCases(i).ns;
    
    % Time-average
    theta_bar = mean(theta_st,2);
    theta_prime = theta_st - theta_bar;
    
    symCases(i).theta_avg = theta_bar; % time-mean
    symCases(i).theta_osc = theta_prime; % wave-varying
    
    symCases(i).theta_osc_rms = sqrt(mean(theta_prime.^2,2));    % rms bending angle
    symCases(i).theta_osc_max = max(abs(theta_prime),[],2);      % max amplitude of bending angle
    
    % At tip
    symCases(i).theta_avg_tip = theta_bar(end);
    symCases(i).theta_osc_tip = symCases(i).theta_osc_max(end);
    
    % Averaged over the entire blade
    symCases(i).theta_avg_blade = mean(theta_bar);
    symCases(i).theta_osc_blade = mean(symCases(i).theta_osc_max);
    
    % Average 1/theta at the blade tip
    I = mean(1./theta_st(ns,:),2);
    symCases(i).Omega = I; % average of 1/theta
    

    % Average 1/theta along the entire blade??
    I = mean(mean(1./theta_st(2:end,2:end)));
    symCases(i).Omega_blade = I; % average of 1/theta
    

    % Max asymmetry relative to blade length
    % --> x-displacement of mean tip position
    xmean = mean(symCases(i).x,2);
    symCases(i).xmean_tip = xmean(end);
    symCases(i).xmean_max = max(xmean);


    % Blade posture asymmetry from Longhuan's paper (see Zhu 2020 fig 8)
    %   beta_xt = mean(x_tip) / max(abs(x_tip))
    symCases(i).beta = symCases(i).xmean_tip / max(abs(symCases(i).x(end,:)));
    %   beta based on angles...
    %   beta_theta_t = mean(theta) / max(abs(theta))
    symCases(i).beta_tht = abs(theta_bar(end)) / max(abs(theta_st(end,:)));

end

%% Get angle ratios for broader shading dataset as well (should get same angle ratios for same parameters, varying R)

for i=1:length(zhuShading)%nCases
    
    motionInd = zhuShading(i).motionIndex;
    theta_st = zhuMotion(motionInd).theta;  % theta as a function of s and t

    ns = zhuMotion(motionInd).ns;
    
    % Time-average
    theta_bar = mean(theta_st,2);
    theta_prime = theta_st - theta_bar;

    zhuShading(i).theta_osc_rms = sqrt(mean(theta_prime.^2,2));    % rms bending angle
    zhuShading(i).theta_osc_max = max(abs(theta_prime),[],2);      % max amplitude of bending angle

    % At tip
    zhuShading(i).theta_avg_tip = theta_bar(end);
    zhuShading(i).theta_osc_tip = zhuShading(i).theta_osc_max(end);
    
    % Averaged over the entire blade
    zhuShading(i).theta_avg_blade = mean(theta_bar);
    zhuShading(i).theta_osc_blade = mean(zhuShading(i).theta_osc_max);
    
    % Average 1/theta at the blade tip
    I = mean(1./theta_st(ns,:),2);
    zhuShading(i).Omega = I; % average of 1/theta
    

    % Average 1/theta along the entire blade??
    I = mean(mean(1./theta_st(2:end,:)));
    zhuShading(i).Omega_blade = I; % average of 1/theta

    % Max asymmetry relative to blade length
    zhuShading(i).xtip_max = max(abs(zhuMotion(motionInd).x(end,:)));

    % Max asymmetry relative to blade length
    % --> x-displacement of mean tip position
    xmean = mean(zhuMotion(motionInd).x,2);
    zhuShading(i).xmean_tip = xmean(end);
    zhuShading(i).xmean_max = max(xmean);


    % Blade posture asymmetry from Longhuan's paper (see Zhu 2020 fig 8)
    %   beta_xt = mean(x_tip) / max(abs(x_tip))
    zhuShading(i).beta = zhuShading(i).xmean_tip / max(abs(zhuMotion(motionInd).x(end,:)));
    %   beta based on angles...
    %   beta_theta_t = mean(theta) / max(abs(theta))
    zhuShading(i).beta_tht = abs(theta_bar(end)) / max(abs(theta_st(end,:)));
    
end

%% How does Omega scale with all parameters? --> FIGURE USED IN PAPER
KC = [enriquezCases.KC, asymCases.KC, symCases.KC, zhuMotion.KC];
Ca = [enriquezCases.Ca, asymCases.Ca, symCases.Ca, zhuMotion.Ca]; 
L =  [enriquezCases.L,  asymCases.L,  symCases.L,  zhuMotion.L]; 
B =  [enriquezCases.B,  asymCases.B,  symCases.B,   zhuMotion.B];
Omega = [enriquezCases.Omega, asymCases.Omega, symCases.Omega, zhuMotion.Omega];

X = [log10(KC)',log10(Ca)',log10(L)', log10(B)'];
Y = log10(abs(Omega))';

mdl = fitlm(X,Y)
% figure(2); clf
% plot(mdl)

c = mdl.Coefficients.Estimate(1);
b_KC = mdl.Coefficients.Estimate(2);
b_Ca = mdl.Coefficients.Estimate(3);
b_L = mdl.Coefficients.Estimate(4);
b_B = mdl.Coefficients.Estimate(5);

% FIGURE FOR PAPER RELATING TO 1/THETA, P/R, etc.
figure(3); clf

% ==== (b) 1/theta versus nondimensional numbers ====
subplot(1,2,2)

% Zhu data (motion only)
h3 = scatter( 10^c * [zhuMotion.KC].^b_KC .* [zhuMotion.L].^b_L .* [zhuMotion.Ca].^b_Ca .* [zhuMotion.B].^b_B, abs([zhuMotion.Omega]), ...
    15, [zhuMotion.beta_tht], 'filled');
% 15, ([zhuMotion.theta_avg_tip]./[zhuMotion.theta_osc_tip]), 'filled');
hold on

% Validation data
h1 = scatter( 10^c * [enriquezCases.KC].^b_KC .* [enriquezCases.L].^b_L .* [enriquezCases.Ca].^b_Ca .* [enriquezCases.B].^b_B, abs([enriquezCases.Omega]), ...
    20,[enriquezCases.beta_tht], '^','filled','MarkerEdgeColor',[0.8 0.8 0.8],'linewidth',1);
% 20,[enriquezCases.theta_avg_tip]./[enriquezCases.theta_osc_tip], '^','filled','MarkerEdgeColor',[0.8 0.8 0.8],'linewidth',1);

% Case study data
scatter( 10^c * [symCases.KC].^b_KC .* [symCases.L].^b_L .* [symCases.Ca].^b_Ca .* [symCases.B].^b_B, abs([symCases.Omega]), ...
    20, [symCases.beta_tht], 'v','filled','MarkerEdgeColor',[0.8 0.8 0.8],'linewidth',1);
% 20, [symCases.theta_avg_tip]./[symCases.theta_osc_tip], 'v','filled','MarkerEdgeColor',[0.8 0.8 0.8],'linewidth',1);
h2 = scatter( 10^c * [asymCases.KC].^b_KC .* [asymCases.L].^b_L .* [asymCases.Ca].^b_Ca .* [asymCases.B].^b_B, abs([asymCases.Omega]), ...
    20, [asymCases.beta_tht],'v','filled','MarkerEdgeColor',[0.8 0.8 0.8],'linewidth',1);
% 20, [asymCases.theta_avg_tip]./[asymCases.theta_osc_tip],'v','filled','MarkerEdgeColor',[0.8 0.8 0.8],'linewidth',1);

% 1:1 line
h = loglog(logspace(-3,5,10),logspace(-3,5,10),'-.','linewidth',2,'color',[.5 .5 .5]);

hcb = colorbar;
% set(get(hcb,'title'),'string','$\overline{\theta_t}/\widetilde{\theta_t}$','interpreter','latex','fontsize',14)
set(get(hcb,'title'),'string','$\beta_{\theta,t}$','interpreter','latex','fontsize',20)


cmap = colormap(hot(256));
cmap = cmap(1:end-30,:);
colormap(cmap)
% clim([0 15])
clim([0 1])

box on

xlim([10^-3 10^2])
ylim([10^-6 10^4.5])

set(gca,'YScale','log','Xscale','log')
set(gca,'xtick',[10^-2 10^0])
set(gca,'fontsize',14)
% ylabel('$|\overline{1/\theta}|$','interpreter','latex','fontsize',30)
xlabel(sprintf('$%.2f KC^{%.3f} Ca^{%.2f} L^{%.2f}  B^{%.2f}$',10^c,b_KC,b_Ca,b_L,b_B),'interpreter','latex','fontsize',24)

leg = legend([h1,h2,h3,h],{'Validation','Case study','Zhu dataset',['Best-fit, $r^2 =$ ',num2str(mdl.Rsquared.Ordinary,2)]},...
    'interpreter','latex','location','best','fontsize',10);
leg.Position = [0.696    0.197   0.1506    0.1748];

text(10^-2.8,10^3.6,'(c)','interpreter','latex','fontsize',18)


% ==== (a) Plot 1/theta vs. P/R

subplot(1,2,1)

scatter(([zhuShading.avgUnshaded]./[zhuShading.R]),abs([zhuShading.Omega]),15,...
        [zhuShading.beta_tht],'filled')
    % [zhuShading.theta_avg_tip]./[zhuShading.theta_osc_tip],'filled')

hcb = colorbar;
% set(get(hcb,'title'),'string','$\overline{\theta_t}/\widetilde{\theta_t}$','interpreter','latex','fontsize',14)
set(get(hcb,'title'),'string','$\beta_{\theta,t}$','interpreter','latex','fontsize',20)


cmap = colormap(hot(256));
cmap = cmap(1:end-30,:);
colormap(cmap)
% clim([0 15])
clim([-0 1])

hold on

scatter(([enriquezCases.avgUnshaded]./[enriquezCases.R]),abs([enriquezCases.Omega]),20,...
        [enriquezCases.beta_tht],'^','filled','MarkerEdgeColor',[0.8 0.8 0.8],'linewidth',1);
    % [enriquezCases.theta_avg_tip]./[enriquezCases.theta_osc_tip],'^','filled','MarkerEdgeColor',[0.8 0.8 0.8],'linewidth',1);
scatter(([asymCases.avgUnshaded]./[asymCases.R]),abs([asymCases.Omega]),20,...
        [asymCases.beta_tht],'v','filled','MarkerEdgeColor',[0.8 0.8 0.8],'linewidth',1);
    % [asymCases.theta_avg_tip]./[asymCases.theta_osc_tip],'v','filled','MarkerEdgeColor',[0.8 0.8 0.8],'linewidth',1);
scatter(([symCases.avgUnshaded]./[symCases.R]),abs([symCases.Omega]),20,...
        [symCases.beta_tht],'v','filled','MarkerEdgeColor',[0.8 0.8 0.8],'linewidth',1);
    % [symCases.theta_avg_tip]./[symCases.theta_osc_tip],'v','filled','MarkerEdgeColor',[0.8 0.8 0.8],'linewidth',1);


% Linear fit with slope = 1
m = 1.0;
b = -1.8;
h = plot(linspace(10^-0.3,10^1.5,10),10^(b) * linspace(10^-0.3,10^1.5,10).^1,'-.','color',[0.5 0.5 0.5],'linewidth',2);

set(gca,'YScale','log','Xscale','log')
set(gca,'fontsize',14)
xlim([10^-0.3 10^1.5])
ylim([10^-6 10^4.5])
set(gca,'xtick',[10^0, 10^1])

leg = legend(h,['$|\overline{1/\theta_t}| = ',num2str(10^b,2),' (\overline{\langle P \rangle}/R)^1$'],...
    'interpreter','latex','location','best','fontsize',14);
leg.Position = [0.1750 0.2063 0.2267 0.0730];
% legend boxoff

xlabel('$\overline{\langle P \rangle}/R$','interpreter','latex','fontsize',24)
ylabel('$|\overline{1/\theta_t}|$','interpreter','latex','fontsize',24)

box on

text(10^-0.22,10^3.6,'(b)','interpreter','latex','fontsize',18)


set(gcf,'position',[162   457   795   306])

% print(gcf,'./figures/1overTheta_subplots_beta.eps','-depsc')



%% Statistics for these two figures - all data vs. only large angle ratio

% === Fig 11(a) - 1/theta vs. P/R

alldata_X = [[zhuShading.avgUnshaded]./[zhuShading.R], [enriquezCases.avgUnshaded]./[enriquezCases.R], ...
    [asymCases.avgUnshaded]./[asymCases.R], [symCases.avgUnshaded]./[symCases.R]];
alldata_Y = abs([[zhuShading.Omega], [enriquezCases.Omega], ...
    [asymCases.Omega], [symCases.Omega]]);
% alldata_Theta = [[zhuShading.theta_avg_tip]./[zhuShading.theta_osc_tip], [enriquezCases.theta_avg_tip]./[enriquezCases.theta_osc_tip], ...
%     [asymCases.theta_avg_tip]./[asymCases.theta_osc_tip], [symCases.theta_avg_tip]./[symCases.theta_osc_tip]];
alldata_beta = [[zhuShading.beta_tht], [enriquezCases.beta_tht], [asymCases.beta_tht], [symCases.beta_tht]];

N = length(alldata_X);

% Linear fit with slope = 1
m = 1.0;
b = -1.8;
fun = @(x) m*x+b;

% RMSE for all data
RMSEa1 = sqrt( sum( (alldata_Y-fun(alldata_X)).^2)/N );
RMSEa1_norm = RMSEa1 / mean(alldata_Y);

thresh = 0.9;
inds = alldata_beta > thresh;

% RMSE for data with "large" angle ratio
RMSEa2 = sqrt( sum( (alldata_Y(inds)-fun(alldata_X(inds))).^2)/N );
RMSEa2_norm = RMSEa2 / mean(alldata_Y(inds));

fprintf('---- Figure 11(a): average 1/theta vs. P/R ---- \n')
fprintf('RMSE is %.2f for all data and %.2f for data with angle ratio > %.2f \n',...
    RMSEa1,RMSEa2,thresh)
fprintf('Normalized RMSE is %.2f for all data and %.2f for data with angle ratio > %.2f \n',...
    RMSEa1_norm,RMSEa2_norm,thresh)



% === Fig 11(b) - 1/theta vs. KC, Ca, L, etc.
alldata_X = [10^c * [zhuMotion.KC].^b_KC .* [zhuMotion.L].^b_L .* [zhuMotion.Ca].^b_Ca .* [zhuMotion.B].^b_B, ...
    10^c * [enriquezCases.KC].^b_KC .* [enriquezCases.L].^b_L .* [enriquezCases.Ca].^b_Ca .* [enriquezCases.B].^b_B, ...
    10^c * [asymCases.KC].^b_KC .* [asymCases.L].^b_L .* [asymCases.Ca].^b_Ca .* [asymCases.B].^b_B, ...
    10^c * [symCases.KC].^b_KC .* [symCases.L].^b_L .* [symCases.Ca].^b_Ca .* [symCases.B].^b_B];
alldata_Y = abs([[zhuMotion.Omega], [enriquezCases.Omega], ...
    [asymCases.Omega], [symCases.Omega]]);
% alldata_Theta = [[zhuMotion.theta_avg_tip]./[zhuMotion.theta_osc_tip], [enriquezCases.theta_avg_tip]./[enriquezCases.theta_osc_tip], ...
%     [asymCases.theta_avg_tip]./[asymCases.theta_osc_tip], [symCases.theta_avg_tip]./[symCases.theta_osc_tip]];
alldata_beta = [[zhuMotion.beta_tht], [enriquezCases.beta_tht], [asymCases.beta_tht], [symCases.beta_tht]];


N = length(alldata_X);

% RMSE for all data
RMSEb1 = sqrt( sum( (alldata_Y-alldata_X).^2)/N );
RMSEb1_norm = RMSEb1 / mean(alldata_Y);

thresh = 0.9;
inds = alldata_beta > thresh;

% RMSE for data with "large" angle ratio
RMSEb2 = sqrt( sum( (alldata_Y(inds)-alldata_X(inds)).^2)/N );
RMSEb2_norm = RMSEb2 / mean(alldata_Y(inds));

fprintf('---- Figure 11(b): average 1/theta vs. nondimensional numbers ---- \n')
fprintf('RMSE is %.2f for all data and %.2f for data with angle ratio > %.2f \n',...
    RMSEb1,RMSEb2,thresh)
fprintf('Normalized RMSE is %.2f for all data and %.2f for data with angle ratio > %.2f \n',...
    RMSEb1_norm,RMSEb2_norm,thresh)

% Coefficient of determination for all data (r^2 = 1 - SS_res/SS_tot)
SSres = sum( (alldata_Y - alldata_X).^2 );
SStot = sum( (alldata_Y - mean(alldata_Y)).^2 );
Rsq = 1 - SSres/SStot;

SSres = sum( (alldata_Y(inds) - alldata_X(inds)).^2 );
SStot = sum( (alldata_Y(inds) - mean(alldata_Y(inds))).^2 );
Rsq = 1 - SSres/SStot;
