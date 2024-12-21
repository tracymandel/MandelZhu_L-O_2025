% figure_parameterization_allData.m
% Parameterization based on all data
% 
% Code to generate Figure 9 in paper
% Mandel & Zhu (2025) L&O
% (c) Tracy Mandel | tracy.mandel@unh.edu
% Last updated 2024/12/21

clear all; close all; clc

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



%% Figure: Parameterization
% Fit all parameters...
KC = [enriquezCases.KC, asymCases.KC, symCases.KC, zhuShading.KC]; % moreCases.KC, 
Ca = [enriquezCases.Ca, asymCases.Ca, symCases.Ca, zhuShading.Ca]; %moreCases.Ca, 
L =  [enriquezCases.L,  asymCases.L,  symCases.L,  zhuShading.L]; %moreCases.L,  
B =  [enriquezCases.B,  asymCases.B,  symCases.B,   zhuShading.B]; %moreCases.B, 
R =  [enriquezCases.R,  asymCases.R,  symCases.R,  zhuShading.R]; %moreCases.R,  
Phi = [enriquezCases.avgUnshaded, asymCases.avgUnshaded, symCases.avgUnshaded, zhuShading.avgUnshaded]; %moreCases.avgUnshaded, 

% KC = [zhuShading.KC];
% Ca = [zhuShading.Ca];
% L =  [zhuShading.L];
% B =  [zhuShading.B];
% R =  [zhuShading.R];
% KH = [zhuShading.kh];
% SUB = [zhuShading.submRatio];
% Phi = [zhuShading.avgUnshaded];

X = [log10(KC)', log10(Ca)', log10(L)', log10(B)', log10(R)'];
Y = log10(Phi)';

% Piecewise fit: cut off values above around log10(avg_unshaded) = -0.1 (~80%)
cutoff = 0.16; % 69%
XX = X(Y<=-cutoff,:);
YY = Y(Y<=-cutoff);

% XX = X;
% YY = Y;

mdl = fitlm(XX,YY)%,'RobustOpts','on')
% figure(100); clf
% plot(mdl)
% anova(mdl,'summary')

c = mdl.Coefficients.Estimate(1);
b_KC = mdl.Coefficients.Estimate(2);
b_Ca = mdl.Coefficients.Estimate(3);
b_L = mdl.Coefficients.Estimate(4);
b_B = mdl.Coefficients.Estimate(5);
b_R = mdl.Coefficients.Estimate(6);

% Assess model/fit performance
N = length(Phi);
ZM = powerLawModel(KC,Ca,L,B,R,c,b_KC,b_Ca,b_L,b_B,b_R);
RMSE = sqrt(sum((Phi-ZM).^2)/N);
MAPE = 100 * 1/N * sum(abs((Phi-ZM)./Phi));


% ==== FIGURE ====
figure(1); clf

subplot(1,2,1)

% Zhu dataset
h3 = scatter(10^c * [zhuShading.KC].^b_KC .* [zhuShading.Ca].^b_Ca .* ...
    [zhuShading.L].^b_L .* [zhuShading.B].^b_B .* [zhuShading.R].^b_R, ...
    [zhuShading.avgUnshaded],25,[zhuShading.R],'o','filled','MarkerEdgeColor','none');
hold on


% CASE STUDY
h2 = scatter(10^c * [asymCases.KC].^b_KC .* [asymCases.Ca].^b_Ca .* ...
    [asymCases.L].^b_L .* [asymCases.B].^b_B .* [asymCases.R].^b_R, ...
    [asymCases.avgUnshaded],100,[asymCases.R],'v','filled','MarkerEdgeColor','k');

h2 = scatter(10^c * [symCases.KC].^b_KC .* [symCases.Ca].^b_Ca .* ...
    [symCases.L].^b_L .* [symCases.B].^b_B .* [symCases.R].^b_R, ...
    [symCases.avgUnshaded],100,[symCases.R],'v','filled','MarkerEdgeColor','k');

% h2 = scatter(10^c * [moreCases.KC].^b_KC .* [moreCases.Ca].^b_Ca .* ...
%     [moreCases.L].^b_L .* [moreCases.B].^b_B .* [moreCases.R].^b_R, ...
%     [moreCases.avgUnshaded],100,[moreCases.R],'v','filled','MarkerEdgeColor','k');


% NEW VALIDATION
h1 = scatter(10^c * [enriquezCases.KC].^b_KC .* [enriquezCases.Ca].^b_Ca .* ...
    [enriquezCases.L].^b_L .* [enriquezCases.B].^b_B .* [enriquezCases.R].^b_R, ...
    [enriquezCases.avgUnshaded],100,[enriquezCases.R],'^','filled','MarkerEdgeColor','k');


% Parameterization
h = loglog(-0.05:0.01:10,-0.05:0.01:10,'-','linewidth',2,'color','k');
% loglog(1:0.05:100,ones(1,length(1:0.05:100)),'-','linewidth',1.5,'color','k')

% hcb = colorbar;
% set(get(hcb,'title'),'string','$R = \Delta/l$','interpreter','latex','fontsize',20)
colormap summer

set(gca,'YScale','log','Xscale','log')
box on
% grid on

% axis image

xlim([10^-1.5 10^2])
ylim([10^-1.5 1.3])

% axis([10^-1.5 10^2 10^-1.5 10^2])

set(gca,'fontsize',16)
xlabel(sprintf('$%.1f KC^{%.3f} Ca^{%.2f} L^{%.3f} B^{%.2f} R^{%.2f}$',...
    10^c,b_KC,b_Ca,b_L,b_B,b_R),'interpreter','latex','fontsize',22)
    %A^{%.2f}$'
ylabel('Avg. light exposure $\overline{\langle P \rangle}$','interpreter','latex','fontsize',22)




legend([h1,h2,h3,h],{'Validation','Case study','Zhu dataset',['Best-fit, $r^2 =$ ',num2str(mdl.Rsquared.Ordinary,2)]},'interpreter','latex','location','southeast','fontsize',12)
legend boxoff


%% Only Ca, B, R
X = [log10(Ca)', log10(B)', log10(R)'];
Y = log10(Phi)';

% Piecewise fit: cut off values above around log10(avg_unshaded) = -0.1 (~80%)
cutoff = 0.16; % 69%
XX = X(Y<=-cutoff,:);
YY = Y(Y<=-cutoff);

% XX = X;
% YY = Y;

mdl = fitlm(XX,YY)%,'RobustOpts','on')
figure(100)
plot(mdl)

c = mdl.Coefficients.Estimate(1);
b_Ca = mdl.Coefficients.Estimate(2);
b_B = mdl.Coefficients.Estimate(3);
b_R = mdl.Coefficients.Estimate(4);

% Assess model/fit performance
N = length(Phi);
ZM = powerLawModelSimple(Ca,B,R,c,b_Ca,b_B,b_R);
RMSE = sqrt(sum((Phi-ZM).^2)/N);
MAPE = 100 * 1/N * sum(abs((Phi-ZM)./Phi));


% ==== FIGURE ====
figure(1)

subplot(1,2,2)

% Zhu dataset
h3 = scatter(10^c * [zhuShading.Ca].^b_Ca .* ...
    [zhuShading.B].^b_B .* [zhuShading.R].^b_R, ...
    [zhuShading.avgUnshaded],25,[zhuShading.R],'o','filled','MarkerEdgeColor','none');
hold on


% CASE STUDY
h2 = scatter(10^c * [asymCases.Ca].^b_Ca .* ...
    [asymCases.B].^b_B .* [asymCases.R].^b_R, ...
    [asymCases.avgUnshaded],100,[asymCases.R],'v','filled','MarkerEdgeColor','k');

h2 = scatter(10^c * [symCases.Ca].^b_Ca .* ...
    [symCases.B].^b_B .* [symCases.R].^b_R, ...
    [symCases.avgUnshaded],100,[symCases.R],'v','filled','MarkerEdgeColor','k');

% h2 = scatter(10^c * [moreCases.Ca].^b_Ca .* ...
%     [moreCases.B].^b_B .* [moreCases.R].^b_R, ...
%     [moreCases.avgUnshaded],100,[moreCases.R],'v','filled','MarkerEdgeColor','k');


% NEW VALIDATION
h1 = scatter(10^c * [enriquezCases.Ca].^b_Ca .* ...
    [enriquezCases.B].^b_B .* [enriquezCases.R].^b_R, ...
    [enriquezCases.avgUnshaded],100,[enriquezCases.R],'^','filled','MarkerEdgeColor','k');


% Parameterization
h = loglog(-0.05:0.01:1,-0.05:0.01:1,'-','linewidth',2,'color','k');
loglog(1:0.05:100,ones(1,length(1:0.05:100)),'-','linewidth',2,'color','k')

hcb = colorbar;
set(get(hcb,'title'),'string','$R = \Delta/l$','interpreter','latex','fontsize',16)
colormap summer

set(gca,'YScale','log','Xscale','log')
box on
% grid on
xlim([10^-1.5 10^2])
ylim([10^-1.5 1.3])

set(gca,'fontsize',16)
xlabel(sprintf('$%.1f Ca^{%.2f} B^{%.2f} R^{%.2f}$',...
    10^c,b_Ca,b_B,b_R),'interpreter','latex','fontsize',22)
    %A^{%.2f}$'
% ylabel('Avg. light exposure $\overline{\langle P \rangle}$','interpreter','latex','fontsize',22)



% legend([h1,h2,h3,h],{'Validation','Case study','Zhu dataset',['Best-fit, $r^2 =$ ',num2str(mdl.Rsquared.Ordinary,2)]},'interpreter','latex','location','southeast')
legend(h,sprintf('Parameterization, RMSE = %.3f',RMSE),'interpreter','latex','fontsize',10,'location','southeast')
legend boxoff


% [left bottom width height]
sp1 = subplot(1,2,1);
sp2 = subplot(1,2,2);
sp2.Position(3) = sp1.Position(3);

subplot(1,2,1)
text(10^-1.4,10^0,'(a)','Interpreter','latex','fontsize',16)

subplot(1,2,2)
text(10^-1.4,10^0,'(b)','Interpreter','latex','fontsize',16)


set(gcf,'position',[512   963   750   307])


% print(gcf,'./figures/parameterization_dual_subplot.eps','-depsc')




%% Function for computing RMSE to total model
function ZM = powerLawModel(x1,x2,x3,x4,x5,b0,b1,b2,b3,b4,b5)
x_cutoff = 1;       % separation point between two regimes

ZM = 10^b0 * (x1).^b1 .* (x2).^b2 .* (x3).^b3 .* (x4).^b4 .* (x5).^b5;   % power-law model
ZM(ZM >= x_cutoff) = 1;            % regime where data level off

end

%% Function for computing RMSE to simplified model
function ZM = powerLawModelSimple(x1,x2,x3,b0,b1,b2,b3)
x_cutoff = 1;       % separation point between two regimes

ZM = 10^b0 * (x1).^b1 .* (x2).^b2 .* (x3).^b3; % power-law model
ZM(ZM >= x_cutoff) = 1;            % regime where data level off

end