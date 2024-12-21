function [indSelf,indNhbr,indTotl] = findTotalShading(pos_x,pos_y,dS)
% Uses the findSelfShading and findNeighborShading sub-functions to compute
% shading of the plant by itself and by neighbors, and the combined total
% plant shading

% Find shading by neighbors
% (gives two vectors, representing shading by blades on the left and on the right)
indNhbr_both = findNeighborShading(pos_x',pos_y',dS)';

% Combine into a single neighbor-shading vector
indNhbr = sum(indNhbr_both,2);
indNhbr(indNhbr > 0) = 1;
indNhbr = logical(indNhbr);     % 0 anywhere not shaded, 1 anywhere shaded by either blade

% Find shading by this blade
indSelf = findSelfShading(pos_x,pos_y);

% Total shading (doesn't double-count)
indTotl = sum([indNhbr,indSelf],2);
indTotl(indTotl > 0) = 1;
indTotl = logical(indTotl);     % 0 anywhere not shaded, 1 anywhere shaded by neighbors or itself

% keyboard


% plot(pos_x,pos_y,'b-','linewidth',2)
% hold on
% plot(pos_x+dS,pos_y,'r-','linewidth',2)
% plot(pos_x-dS,pos_y,'r-','linewidth',2)
% plot(pos_x(indTotl),pos_y(indTotl),'k.','markersize',10)
% 
% axis([-1 1 0 1])
% % axis image
% 
% drawnow