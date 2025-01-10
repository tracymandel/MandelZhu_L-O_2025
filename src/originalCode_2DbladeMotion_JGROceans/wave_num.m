function [k]=wave_num(omega,h)

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


g=9.81;
y=omega*omega*h/g;
d(1)=0.66666666667;
d(2)=0.35555555555;
d(3)=0.1608465608;
d(4)=0.0632098765;
d(5)=0.0217540484;
d(6)=0.0065407983;
dy=0;
for i=1:6
    dy=dy+d(i)*y^i;
end
k=sqrt(y*y+y/(1+dy))/h;