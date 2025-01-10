function [J]=Jacob(F_Sys,X)

% ========================================================================
% copyright Longhuan Zhu 2020
% 
% For model details, see
% Zhu, L., Zou, Q.-P., Huguenard, K., and Fredriksson, D. W. (2020). Mechanisms for the Asymmetric
% Motion of Submerged Aquatic Vegetation in Waves: A Consistent-Mass Cable Model. Journal of
% Geophysical Research: Oceans, 125(2):e2019JC015517.
% 
% Primary code repository for Zhu et al. (2020) located at https://github.com/lzhu7/flexibleBladeDynamicsWavesCurrents
% ========================================================================


%input arguments:
%F_Sys--system of nonlinear functions
%X--computing points
%output:
%J--Jacobian matrix

Nx=12;
Dx=1e-6;
J=zeros(6,12);%6 equations 12 variables
%forward finite differences
for j=1:Nx
    Xj=zeros(Nx,1);
    Xj(j)=Dx;
    J(:,j)=(F_Sys(X+Xj)-F_Sys(X))/Dx;
end