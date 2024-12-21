function shaded_inds = findNeighborShading(pos_x,pos_y,dS)
% shaded_ind = findNeighborShading(pos_x,pos_y,dS)
% FINDNEIGHBORSHADING finds the indices where the blade, whose position is
% defined by pos_x and pos_y, is shaded by neighboring blades on either
% side, set a distance +/- dS away from the primary blade
%
% Last updated 2022/08/26 by TLM
% No sun angle yet...


% % Synthetic blade for a test case: Will have complex shading
% % (sinusoid that decays w/height)
% pos_y = linspace(0,1,256);
% pos_x = sin(10*pos_y).*exp(-pos_y);

% % % test case....
% X = vary_Ca(3).X(:,45);
% pos_x = real(X)';    % x-position of blade
% pos_y = imag(X)';    % y-position of blade
% dS = 0.3;

N = length(pos_x);

% Define position of neighboring blades
nx(1,:) = pos_x - dS;
nx(2,:) = pos_x + dS;
ny(1,:) = pos_y;
ny(2,:) = pos_y;

% % Visualize blade positions
% figure(3); clf
% % Blade of interest
% plot(pos_x,pos_y,'bo-','linewidth',1,'markersize',3);
% hold on
% % Neighboring blade(s)
% for j=1:2
%     plot(nx(j,:),ny(j,:),'o-','linewidth',1,'markersize',3,'color',[0.4660 0.6740 0.1880]);
% end

% Loop over both neighboring blades and find overlap...

shaded_inds = zeros(2,N);

for j=1:2
    
    shaded_ind = zeros(1,N);    % indices where the blade is self-shaded
    for k=N-1:-1:1
       ind = ( (pos_y <= ny(j,k)) & (pos_y <= ny(j,k+1))) & ... % below the top segment...
            ( (pos_x >= nx(j,k) & pos_x <= nx(j,k+1)) | (pos_x <= nx(j,k) & pos_x >= nx(j,k+1)) );
        
       shaded_ind = shaded_ind + ind;
        
%         shaded_ind(shaded_ind > 0) = 1;
%         shaded_ind = logical(shaded_ind);
%         
%         % Plot everything
%         figure(3); clf
%         plot(pos_x,pos_y,'bo-','linewidth',1,'markersize',3);
%         hold on
%         plot(nx(j,:),ny(j,:),'o-','linewidth',1,'markersize',3,'color',[0.4660 0.6740 0.1880]);
%         plot(pos_x(shaded_ind),pos_y(shaded_ind),'ks','markersize',5,'markerfacecolor','k')
        
    end
    
    shaded_ind(shaded_ind > 0) = 1;
    shaded_ind = logical(shaded_ind);
    
%     plot(pos_x(shaded_ind),pos_y(shaded_ind),'ks','markersize',5,'markerfacecolor','k')
%     xlim([-1 1])
%     ylim([0 1])
%     set(gca,'fontsize',14)
%     xlabel('$X$','interpreter','latex','fontsize',20)
%     ylabel('$Y$','interpreter','latex','fontsize',20)
    
    shaded_inds(j,:) = shaded_ind;  % save shading by neighbors separately
    
%     drawnow
    
end
