% figure_PhaseResolvedShading.m
% =============================================================
% Phase-resolved shading throughout one wave period for
% one representative case (use base case of Enriquez validation)
% =============================================================
% 
% Code to generate Figure 3 in paper
% Mandel & Zhu (2025) L&O
% (c) Tracy Mandel | tracy.mandel@unh.edu
% Last updated 2024/12/21

clear all; close all; clc
load('../../data/caseStudySym/caseStudySym_MotionAndShading.mat')

%% Choose a specific case of shading
caseInd = 1;

[~,nt] = size(symCases(caseInd).x);
dt = floor(nt/8);              % show motion every dt timesteps

tlabels = {'0','1/8','1/4','3/8','1/2','5/8','3/4','7/8'};


% Above plant motion: Plot of wave phase
figure(1); clf

h = subplot(3,4,2:3);
t = 0:0.01:1;
u = cos(2*pi*t);
plot(t,u,'-','linewidth',1.5);
hold on
plot(t,zeros(1,length(t)),':','Color',[0.5 0.5 0.5],'linewidth',1.3)
xlim([0 1])
set(gca,'fontsize',12)

theta = 0:1/4:1;
theta_labels = {'0','1/4','1/2','3/4','1'};

xticks(theta)
xticklabels(theta_labels)
box on
ylabel('$U/U_m$','interpreter','latex','fontsize',24)
xlabel('$t/T$','interpreter','latex','fontsize',18)

pos = get(h,'position'); % [left bottom width height]
set(h,'position',[pos(1) pos(2)+0.1 pos(3)+0.015 pos(4)-0.08])


ctr = 5;

% Loop over timesteps
for t = 1:dt:nt-10
    
    x = symCases(caseInd).x(:,t);
    z = symCases(caseInd).z(:,t);

    delta = symCases(caseInd).R;          % spacing
    shaded = logical(symCases(caseInd).shadedPts(:,t));

    frac_unshaded = 1 - mean(shaded);

    sp_handle = subplot(3,4,ctr);
    hold on

    % Define position of neighboring blades
    nx(1,:) = x - delta;
    nx(2,:) = x + delta;
    nz(1,:) = z;
    nz(2,:) = z;
    for ll=1:2
        plot(nx(ll,:),nz(ll,:),'linewidth',4,'color',[0.4660 0.62 0.1880])
    end
    
    % Plot blade of interest and shaded points
    plot(x,z,'linewidth',4,'color',[0.4660 0.68 0.1880])
    plot(x(shaded),z(shaded),'k.','markersize',10)
    
    axis image
    axis([-1.1 1.1 0 1])
    set(gca,'fontsize',12)
%     ylabel('$z/l$','interpreter','latex','fontsize',20)
    
    phase = t/nt;
    title(sprintf('$ t = %s$ $T$',tlabels{ctr-4}),'interpreter','latex','fontsize',14)

    if ctr==5 || ctr==9
        ylabel('$z/l$','interpreter','latex','fontsize',20)
    else
        set(gca,'ytick',[])
    end

    box on
    
    % % Modify subplot position to be more compact
    pos = get(sp_handle, 'Position'); % gives the position of current sub-plot
    % new_pos = pos;
    new_pos = pos + [0 0.02 0.015 0]; % left bottom width height
    if ctr>8
        new_pos = new_pos + [0 -0.01 0 0];
        xlabel('$x/l$','interpreter','latex','fontsize',20)
    end
    set(sp_handle,'Position',new_pos) % set new position of current sub-plot


    text(0.285,0.12,['$\langle P \rangle =$ ',num2str(frac_unshaded,2)],'interpreter','latex','fontsize',11.5)

    % pause

    ctr = ctr+1;

end

set(gcf,'position',[198        1011         931         361])

% % Uncomment to show nondimensional parameters for this case
% fprintf('KC = %.2f, Ca = %.2f, L = %.2f, B = %.2f, R = %.2f \n Phi = %.2f \n',...
%     symCases(caseInd).KC,symCases(caseInd).Ca,symCases(caseInd).L,...
%     symCases(caseInd).B,symCases(caseInd).R,symCases(caseInd).avgUnshaded)

% print(gcf,'../figures/phase_resolved_basecase.eps','-depsc')