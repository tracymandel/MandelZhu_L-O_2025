% figure_shadingAndMotionDimensionalResults
% Plot wider range around Enriquez base case ("case study")
% 
% Code to generate Figures 4-7 and supplemental Figs S1-S3 in paper
% Mandel & Zhu (2025) L&O
% (c) Tracy Mandel | tracy.mandel@unh.edu
% Last updated 2024/12/21

clear all; close all; clc

load('../../data/caseStudySym/caseStudySym_MotionAndShading.mat')
load('../../data/caseStudyAsym/caseStudyAsym_MotionAndShading.mat')


%% "SYMMETRIC" RESULTS - Choose which parameter to vary in figure
% (Figures 4-6 in paper and supplemental Figures S1-S3)

% Set figure case   (change this to create different figures)
figCase = 'Delta';

if strcmp(figCase,'H')
    varCases = [2,5];
elseif strcmp(figCase,'T')
    varCases = [6,9];
elseif strcmp(figCase,'l')
    varCases = [10,13];
elseif strcmp(figCase,'E')
    varCases = [14,17];
elseif strcmp(figCase,'Delta')
    varCases = [18,21];
elseif strcmp(figCase,'h')
    varCases = [22,25];
end


% For labeling pcolor (left side)
theta = 0:1/4:1;
theta_labels = {'0','1/4','1/2','3/4','1'};

figure(1); clf

cmap = colormap('gray');
cmap = flipud(cmap);
colormap(cmap);

% Plot of wave phase
h = subplot(4,2,1);
t = 0:0.01:1;
u = cos(2*pi*t);
plot(t,u,'-','linewidth',1.5);
hold on
plot(t,zeros(1,length(t)),':','Color',[0.5 0.5 0.5],'linewidth',1.3)
xlim([0 1])
set(gca,'fontsize',12)
xticks(theta)
xticklabels(theta_labels)
box on
ylabel('$U/U_m$','interpreter','latex','fontsize',24)

pos = get(h,'position'); % [left bottom width height]
set(h,'position',[pos(1) pos(2) pos(3) pos(4)-0.08])


ctr = 3;

for i = [varCases(1),1,varCases(2)]

    nt = symCases(i).nt;    % # timesteps
    ns = symCases(i).ns;    % spatial resolution

    [~,T] = size(symCases(i).x);
    dt = floor(T/16);                % show motion every dt timesteps

    % For plant visuals (right side)
    % set transparency based on wave phase
    % solid = start of wave phase, transparent = end of wave phase (or
    % half-phase)
    trange = 1:dt:T;
    trange_fullres = 1:T;
    % transp = fliplr(linspace(0.1,1,length(trange)));
    transp = (linspace(0.1,1,length(trange)));


    % ~~ Left column: Phase-resolved shading ~~
    subplot(4,2,ctr)
    pcolor(linspace(0,1,T),linspace(0,1,ns),symCases(i).shadedPts); shading flat
    
    set(gca,'fontsize',12)
    xticks(theta)
    xticklabels(theta_labels)
    box on
    ylabel('$s/l$','interpreter','latex','fontsize',24)
    
%     grid on
%     set(gca,'layer','top')

    ax = gca;
    ax.LineWidth = 1.2;
    box on


    % ~~ Right column: Motion of blade of interest and neighbors + shading ~~
    subplot(4,2,ctr+1);
    dS = symCases(i).R;
    
    ctr2 = 1;
    
    for t = trange
        
        x = symCases(i).x(:,t);
        y = symCases(i).z(:,t);
        
        % Define position of neighboring blades
        nx(1,:) = x - dS;
        nx(2,:) = x + dS;
        ny(1,:) = y;
        ny(2,:) = y;
        for ll=1:2
            plot(nx(ll,:),ny(ll,:),'linewidth',4,'color',[.8 .8 .8, transp(ctr2)])%[0.4660 0.62 0.1880])
            hold on
        end
        
        ctr2 = ctr2 + 1;
        
    end
    
    ctr2 = 1;
    for t = trange
        
        x = symCases(i).x(:,t);
        y = symCases(i).z(:,t);
        
        plot(x,y,'linewidth',4,'color',[0.4660 0.68 0.1880,transp(ctr2)])
       
        ctr2 = ctr2 + 1;
        
    end
    
    ctr2 = 1;
    for t = trange
        
        x = symCases(i).x(:,t);
        y = symCases(i).z(:,t);
        
        shaded = logical(symCases(i).shadedPts(:,t));
        % ~~ issue here -- even if lines are discontinuous, they connect
        % and matlab doesn't support transparency for markers...
        x_shaded = x(shaded);
        y_shaded = y(shaded);

        [x_split,y_split,~] = splitGaps(x,y,shaded);
        % keyboard

        if iscell(x_split)
            N = length(x_split);
            for c = 1:N
                plot(x_split{c},y_split{c},'linewidth',4,'color',[0 0 0 transp(ctr2)])
            end
        else
            plot(x(shaded),y(shaded),'linewidth',4,'color',[0 0 0 transp(ctr2)])%'k.','markersize',9)
        end
        % plot(x(shaded),y(shaded),'.','markersize',3.5,'markerfacecolor',[0 0 0 transp(ctr2)])%'k.','markersize',9)

        hold on
        
        axis equal
        ylim([0 1])
        xlim([-1.2 1.2])
        
        set(gca,'fontsize',12)
        box on
        
        ylabel('$z/l$','interpreter','latex','fontsize',24)
        
       
        ctr2 = ctr2 + 1;
        
    end
    

    if strcmp(figCase,'H')
        text(1.3,0.5,sprintf('%.2f',symCases(i).H),'interpreter','latex','fontsize',18)
    elseif strcmp(figCase,'T')
        text(1.3,0.5,sprintf('%.0f',symCases(i).T),'interpreter','latex','fontsize',18)
    elseif strcmp(figCase,'l')
        text(1.3,0.5,sprintf('%.2f',symCases(i).l),'interpreter','latex','fontsize',18)
    elseif strcmp(figCase,'E')
        text(1.25,0.5,sprintf('%.0f',symCases(i).E*10^(-6)),'interpreter','latex','fontsize',18)
    elseif strcmp(figCase,'Delta')
        text(1.3,0.5,num2str(symCases(i).Delta,2),'interpreter','latex','fontsize',18)
    elseif strcmp(figCase,'h')
        text(1.3,0.5,num2str(symCases(i).h,1),'interpreter','latex','fontsize',18)
    end
   
    if ctr == 3 && strcmp(figCase,'Delta')
        text(0.35,0.15,['$\overline{\langle P \rangle} =$ ',num2str(symCases(i).avgUnshaded,2)],'interpreter','latex','fontsize',12)
    else
        text(0.42,0.15,['$\overline{\langle P \rangle} =$ ',num2str(symCases(i).avgUnshaded,2)],'interpreter','latex','fontsize',12)
    end

    
    ax = gca;
    ax.LineWidth = 1.1;
    box on
    
    ctr = ctr+2;

end


subplot(4,2,7)
xlabel('$t/T$','interpreter','latex','fontsize',24)

subplot(4,2,8)
xlabel('$x/l$','interpreter','latex','fontsize',24)



subplot(4,2,4)
if strcmp(figCase,'H')
    text(1.15,1.15,'$H$ (m):','interpreter','latex','fontsize',20)
elseif strcmp(figCase,'T')
    text(1.15,1.15,'$T$ (sec):','interpreter','latex','fontsize',20)
elseif strcmp(figCase,'l')
    text(1.15,1.15,'$l$ (m):','interpreter','latex','fontsize',20)
elseif strcmp(figCase,'E')
    text(1.0,1.15,'$E$ (MPa):','interpreter','latex','fontsize',20)
elseif strcmp(figCase,'Delta')
    text(1.15,1.15,'$\Delta$ (m):','interpreter','latex','fontsize',20)
elseif strcmp(figCase,'h')
    text(1.1,1.15,'$h$ (m):','interpreter','latex','fontsize',20)
end


set(gcf,'position',[117   968   653   489])%[117        1026         686         431])


% % Comment or uncomment for saving the figure
% if strcmp(figCase,'h')
%     print(gcf,['../figures/dimensionalResults/shading_and_motion_vs_','depth','_enriquez_wider_cases.eps'],'-depsc')
% else
%     print(gcf,['../figures/dimensionalResults/shading_and_motion_vs_',figCase,'_enriquez_wider_cases.eps'],'-depsc')
% end



%% "ASYMMETRIC" RESULTS - Choose which parameter to vary in figure
% (Figure 7 in paper)

% Set figure case
varCases = [47,15,20];

% For labeling pcolor (left side)
theta = 0:1/4:1;
theta_labels = {'0','1/4','1/2','3/4','1'};

figure(1); clf

cmap = colormap('gray');
cmap = flipud(cmap);
colormap(cmap);

% Plot of wave phase
h = subplot(4,2,1);
t = 0:0.01:1;
u = cos(2*pi*t);
plot(t,u,'-','linewidth',1.5);
hold on
plot(t,zeros(1,length(t)),':','Color',[0.5 0.5 0.5],'linewidth',1.3)
xlim([0 1])
set(gca,'fontsize',12)
xticks(theta)
xticklabels(theta_labels)
box on
ylabel('$U/U_m$','interpreter','latex','fontsize',24)

pos = get(h,'position'); % [left bottom width height]
set(h,'position',[pos(1) pos(2) pos(3) pos(4)-0.08])


ctr = 3;

for i = varCases

    nt = asymCases(i).nt;    % # timesteps
    ns = asymCases(i).ns;    % spatial resolution

    [~,T] = size(asymCases(i).x);
    dt = floor(T/16);                % show motion every dt timesteps

    % For plant visuals (right side)
    % set transparency based on wave phase
    % solid = start of wave phase, transparent = end of wave phase (or
    % half-phase)
    trange = 1:dt:T;
    trange_fullres = 1:T;
    % transp = fliplr(linspace(0.1,1,length(trange)));
    transp = (linspace(0.1,1,length(trange)));


    % ~~ Left column: Phase-resolved shading ~~
    subplot(4,2,ctr)
    pcolor(linspace(0,1,T),linspace(0,1,ns),asymCases(i).shadedPts); shading flat
    
    set(gca,'fontsize',12)
    xticks(theta)
    xticklabels(theta_labels)
    box on
    ylabel('$s/l$','interpreter','latex','fontsize',24)
    
%     grid on
%     set(gca,'layer','top')

    ax = gca;
    ax.LineWidth = 1.2;
    box on


    % ~~ Right column: Motion of blade of interest and neighbors + shading ~~
    subplot(4,2,ctr+1);
    dS = asymCases(i).R;
    
    ctr2 = 1;
    
    for t = trange
        
        x = asymCases(i).x(:,t);
        y = asymCases(i).z(:,t);
        
        % Define position of neighboring blades
        nx(1,:) = x - dS;
        nx(2,:) = x + dS;
        ny(1,:) = y;
        ny(2,:) = y;
        for ll=1:2
            plot(nx(ll,:),ny(ll,:),'linewidth',4,'color',[.8 .8 .8, transp(ctr2)])%[0.4660 0.62 0.1880])
            hold on
        end
        
        ctr2 = ctr2 + 1;
        
    end
    
    ctr2 = 1;
    for t = trange
        
        x = asymCases(i).x(:,t);
        y = asymCases(i).z(:,t);
        
        plot(x,y,'linewidth',4,'color',[0.4660 0.68 0.1880,transp(ctr2)])
       
        ctr2 = ctr2 + 1;
        
    end
    
    ctr2 = 1;
    for t = trange
        
        x = asymCases(i).x(:,t);
        y = asymCases(i).z(:,t);
        
        shaded = logical(asymCases(i).shadedPts(:,t));
        % ~~ issue here -- even if lines are discontinuous, they connect
        % and matlab doesn't support transparency for markers...
        x_shaded = x(shaded);
        y_shaded = y(shaded);

        [x_split,y_split,~] = splitGaps(x,y,shaded);
        % keyboard

        if iscell(x_split)
            N = length(x_split);
            for c = 1:N
                plot(x_split{c},y_split{c},'linewidth',4,'color',[0 0 0 transp(ctr2)])
            end
        else
            plot(x(shaded),y(shaded),'linewidth',4,'color',[0 0 0 transp(ctr2)])%'k.','markersize',9)
        end
        % plot(x(shaded),y(shaded),'.','markersize',3.5,'markerfacecolor',[0 0 0 transp(ctr2)])%'k.','markersize',9)

        hold on
        
        axis equal
        ylim([0 1])
        xlim([-1.2 1.2])
        
        set(gca,'fontsize',12)
        box on
        
        ylabel('$z/l$','interpreter','latex','fontsize',24)
        
       
        ctr2 = ctr2 + 1;
        
    end
    

    % if strcmp(figCase,'H')
    %     text(1.3,0.5,sprintf('%.2f',symCases(i).H),'interpreter','latex','fontsize',18)
    % elseif strcmp(figCase,'T')
    %     text(1.3,0.5,sprintf('%.0f',symCases(i).T),'interpreter','latex','fontsize',18)
    % elseif strcmp(figCase,'l')
    %     text(1.3,0.5,sprintf('%.2f',symCases(i).l),'interpreter','latex','fontsize',18)
    % elseif strcmp(figCase,'E')
    %     text(1.25,0.5,sprintf('%.0f',symCases(i).E*10^(-6)),'interpreter','latex','fontsize',18)
    % elseif strcmp(figCase,'Delta')
    %     text(1.3,0.5,num2str(symCases(i).Delta,2),'interpreter','latex','fontsize',18)
    % elseif strcmp(figCase,'h')
    %     text(1.3,0.5,num2str(symCases(i).h,1),'interpreter','latex','fontsize',18)
    % end
    % 
    % if ctr == 3 && strcmp(figCase,'Delta')
    %     text(0.39,0.15,['$\overline{\langle \phi \rangle} =$ ',num2str(symCases(i).avgUnshaded,2)],'interpreter','latex','fontsize',12)
    % else
        text(0.42,0.15,['$\overline{\langle P \rangle} =$ ',num2str(asymCases(i).avgUnshaded,2)],'interpreter','latex','fontsize',12)
    % end

    
    ax = gca;
    ax.LineWidth = 1.1;
    box on
    
    ctr = ctr+2;

end


subplot(4,2,7)
xlabel('$t/T$','interpreter','latex','fontsize',24)

subplot(4,2,8)
xlabel('$x/l$','interpreter','latex','fontsize',24)


subplot(4,2,4)
text(1.35,0.5,'(i)','interpreter','latex','fontsize',18)
subplot(4,2,6)
text(1.35,0.5,'(ii)','interpreter','latex','fontsize',18)
subplot(4,2,8)
text(1.35,0.5,'(iii)','interpreter','latex','fontsize',18)



set(gcf,'position',[117   968   653   489])%[117        1026         686         431])


% % Comment or uncomment for saving the figure
% print(gcf,'../figures/dimensionalResults/shading_and_motion_asymmetric_examples.eps','-depsc')
