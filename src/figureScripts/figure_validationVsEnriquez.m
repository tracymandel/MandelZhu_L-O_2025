% figure_validationVsEnriquez.m
% Using Enriquez validation cases computed on Premise (60 wave periods)
% 
% Code to generate Figure 2 in paper
% Mandel & Zhu (2025) L&O
% (c) Tracy Mandel | tracy.mandel@unh.edu
% Last updated 2024/12/21

%% ============ Start here if data already simulated! ==============

% Vary only hydrodynamic parameters (more uncertain)
H = [0.1214, 0.4, 1.1066];  % upper and lower 1/10, plus most frequent wave height
TT = [5, 7, 9];              % wave period, sec
h = [4.85, 5, 5.15];        % water depth, m

load('../../data/enriquezValidation/enriquezValidation_MotionAndShading.mat')

%% Enriquez - compare with validation cases

nVariables = 3;             % # parameters varying
nCases = 3^nVariables;      % # validation cases
M = 1*nCases;

load('../../data/enriquezValidation/digitizedData/enriquez_fig1_data_cleaned.mat')

z = [15, 10, 5, 0];     % measurement height, cmab

% Plot time series and extract mean, median, variance
% figure(1); clf
for i=1:4
    % subplot(2,2,i)
    % plot(data(i).t,data(i).PFD,'k.')
    % xlabel('time (sec)')
    % ylabel('PFD')
    
    PFD_mean(i) = mean(data(i).PFD);
    PFD_med(i) = median(data(i).PFD);
    PFD_var(i) = var(data(i).PFD);
    % 
    % hold on
    % plot(0:60,PFD_mean(i)*ones(1,61),'b-','linewidth',1.5)
    % % plot(0:60,PFD_med(i)*ones(1,61),'r-','linewidth',1.5)
end

%% Make main validation figure

figure(3); clf
hold on

cmap = colormap('summer');
linstyle = '-';
linwidth = 3.5;

for i=1:27
    % subplot(3,3,i)
    if enriquezCases(i).H == H(1)
        subplot(1,3,1)
        hold on
    elseif enriquezCases(i).H == H(2)
        subplot(1,3,2)
        hold on
        % clr = cmap(120,:);
    else
        subplot(1,3,3)
        hold on
        % clr = cmap(240,:);
    end

    % if enriquezCases(i).h == h(1)
    %     linstyle = '--';
    % elseif enriquezCases(i).h == h(2)
    %     linstyle = '-';
    % else
    %     linstyle = ':';
    % end

    if enriquezCases(i).h ~= h(2)
        continue
    end

    if enriquezCases(i).T == TT(1)
        % linstyle = ':';
        % linwidth = 1.5;
        clr = cmap(1,:);
    elseif enriquezCases(i).T == TT(2)
        % linstyle = '--';
        % linwidth = 2.5;
        clr = cmap(120,:);
    else
        % linstyle = '-';
        % linwidth = 3.5;
        clr = cmap(240,:);
    end

    phi_z = 1-mean(enriquezCases(i).shadedPts,2);
    z_bar = mean(enriquezCases(i).z,2);
    plot(phi_z,z_bar,linstyle,'linewidth',linwidth,'color',clr);

    if i==10
        h2 = plot(phi_z,z_bar,linstyle,'linewidth',linwidth,'color',clr);
    elseif i==11
        h3 = plot(phi_z,z_bar,linstyle,'linewidth',linwidth,'color',clr);
    elseif i==12
        h4 = plot(phi_z,z_bar,linstyle,'linewidth',linwidth,'color',clr);
    end


end

% Data
for i=1:3
    subplot(1,3,i)
    h1 = plot(PFD_mean/max(PFD_mean),z/20,'ko','linewidth',1.4,'markersize',10,'markerfacecolor','b');

    set(gca,'fontsize',16)
    box on
    ylim([0 1])
    % grid on
    ax = gca;
    ax.LineWidth = 1.2;
    box on
    xticks(0:0.2:1)

    xlabel('$\overline{P}$','interpreter','latex','fontsize',28)
end

subplot(1,3,1)
ylabel('$z/l$','interpreter','latex','fontsize',28)
title(['(a) $H = ',num2str(H(1),2),'$ m'],'interpreter','latex','fontsize',24)
subplot(1,3,2)
title(['(b) $H = ',num2str(H(2),2),'$ m'],'interpreter','latex','fontsize',24)
subplot(1,3,3)
title(['(c) $H = ',num2str(H(3),2),'$ m'],'interpreter','latex','fontsize',24)

% legend([h1 h2 h3 h4],{'Enriquez obs.','Simulated','Enriquez obs. $\langle \overline{\phi} \rangle$',...
%     'Simulated $\langle \overline{\phi} \rangle$'},'interpreter','latex','fontsize',12)
% legend([h1 h2 h3 h4,h5],{'Enriquez observations','Simulated ($H_{min}$)',...
%     'Simulated ($H_{bc}$)','Simulated ($H_{max}$)','Simulated base case',},...
%     'interpreter','latex','fontsize',14,'location','southeast')
subplot(1,3,1)
legend([h1 h2 h3 h4],{'Enriquez obs.',['$T = $ ',num2str(TT(1),2),' sec'],...
    ['$T = $ ',num2str(TT(2),2),' sec'],['$T = $ ',num2str(TT(3),2),' sec']},...
    'interpreter','latex','fontsize',14,'location','southeast')

set(gcf,'position',[85   912   767   340])


% print(gcf,'../../figures/enriquez_validation_H_most_frequent.eps','-depsc')