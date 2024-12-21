function shaded_ind = findSelfShading(pos_x,pos_y)
% shaded_ind = findSelfShading(pos_x,pos_y)
% FINDSELFSHADING finds the indices where the blade, whose position is
% defined by pos_x and pos_y, shades itself
%
% Last updated 2022/08/26 by TLM
% No sun angle yet...

% % Synthetic blade for a test case: Will have complex shading
% % (sinusoid that decays w/height)
% pos_y = linspace(0,1,256);
% pos_x = sin(10*y).*exp(-y);

% # points
N = length(pos_x);

% % Plot the blade
% figure(2); clf
% plot(pos_x,pos_y,'bo-','linewidth',1,'markersize',3,'markerfacecolor','b'); hold on

% Find local minima/maxima and endpoints and sequentially shade
[~,loc] = findpeaks(abs(pos_x));

% keyboard

% Add the free end of the blade as well
try
    pks_x = [pos_x(loc),pos_x(end)];
    pks_y = [pos_y(loc),pos_y(end)];
catch
    pks_x = [pos_x(loc)',pos_x(end)];
    pks_y = [pos_y(loc)',pos_y(end)];
end
% plot(pks_x,pks_y,'r*')

% # of peaks/endpoints
npks = length(pks_x);

% Move from top to bottom of blade to find shaded regions
shaded_ind = zeros(N,1);    % indices where the blade is self-shaded
for i = npks:-1:2
    ind = (pos_y <= pks_y(i) & pos_y <= pks_y(i-1)) & ... % below the top segment...
        ( (pos_x >= pks_x(i) & pos_x <= pks_x(i-1)) | (pos_x <= pks_x(i) & pos_x >= pks_x(i-1)) );
        % ... and x-position falls between the two peaks
    shaded_ind = shaded_ind + ind;
end

shaded_ind(shaded_ind > 0) = 1;
shaded_ind = logical(shaded_ind);

% plot(pos_x(shaded_ind),pos_y(shaded_ind),'ks','markersize',5,'markerfacecolor','k')
% xlim([-1 1])
% ylim([0 1])
% set(gca,'fontsize',14)
% xlabel('$X$','interpreter','latex','fontsize',20)
% ylabel('$Y$','interpreter','latex','fontsize',20)

