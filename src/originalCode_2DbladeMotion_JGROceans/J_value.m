function J_value=J_value(Jacob,Eq_Govern_Done,Xi,Xj,u_wave_old,u_wave_new,v_wave_old,v_wave_new,n,m,w0,Umax)

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


% Umax=zeros(n,1);
%Boundary conditions
    J_value=zeros(6*(n+1),6*(n+1));
    J_value(1,3)=1;
    J_value(2,4)=1;
    J_value(3,6)=1;
    J_value(4,6*n+1)=1;
    J_value(5,6*n+2)=1;
    J_value(6,6*n+5)=1;
                
for k=2:n+1
    X0=Xi(6*(k-2)+1:6*k);
    X=Xj(6*(k-2)+1:6*k);
    Eq=@(X) Eq_Govern_Done(X0,X,u_wave_old(k-1:k),u_wave_new(k-1:k),v_wave_old(k-1:k),v_wave_new(k-1:k),m(k-1),w0(k-1),Umax(k-1));
    J_value(6*(k-1)+1:6*k,6*(k-2)+1:6*k)=Jacob(Eq,X);
end